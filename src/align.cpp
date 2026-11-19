#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace msakit {
namespace {
using Score = std::int64_t;

bool is_residue(char value) {
    return value == 'A' || value == 'C' || value == 'G' || value == 'T' ||
           value == 'U' || value == 'N';
}

Score substitution(char center, char query) {
    if (center == 'N' || query == 'N') return 0;
    if (center == query || (center == 'T' && query == 'U') ||
        (center == 'U' && query == 'T')) return 1;
    return -1;
}

struct Projection {
    std::string at_center;
    std::vector<std::string> insertions;
};

Projection align_pair(std::string_view center, std::string_view query) {
    detail::check_score_bounds(center.size(), query.size());
    const auto height = detail::checked_add(center.size(), 1);
    const auto width = detail::checked_add(query.size(), 1);
    const auto cells = detail::checked_multiply(height, width);
    (void)detail::checked_multiply(cells, sizeof(Score));
    if (cells > std::vector<Score>().max_size()) {
        throw detail::Error(detail::kCoreError, "pairwise matrix exceeds container size");
    }

    // Deliberately use the full, easy-to-read DP table in this reference core.
    // It is local to this pair and is released before processing the next query.
    std::vector<Score> scores(cells);
    const auto cell = [&](std::size_t i, std::size_t j) -> Score& {
        return scores[i * width + j];
    };
    for (std::size_t i = 1; i < height; ++i) {
        cell(i, 0) = static_cast<Score>(i) * detail::kGapScore;
    }
    for (std::size_t j = 1; j < width; ++j) {
        cell(0, j) = static_cast<Score>(j) * detail::kGapScore;
    }
    for (std::size_t i = 1; i < height; ++i) {
        for (std::size_t j = 1; j < width; ++j) {
            cell(i, j) = std::max({
                cell(i - 1, j - 1) + substitution(center[i - 1], query[j - 1]),
                cell(i - 1, j) + detail::kGapScore,
                cell(i, j - 1) + detail::kGapScore});
        }
    }

    std::string aligned_center;
    std::string aligned_query;
    const auto longest_trace = detail::checked_add(center.size(), query.size());
    if (longest_trace > aligned_center.max_size()) {
        throw detail::Error(detail::kCoreError, "pairwise traceback exceeds string size");
    }
    aligned_center.reserve(longest_trace);
    aligned_query.reserve(longest_trace);
    std::size_t i = center.size();
    std::size_t j = query.size();
    while (i != 0 || j != 0) {
        // The order of these branches defines the documented tie rule:
        // consume both, consume center only, then consume query only.
        if (i != 0 && j != 0 && cell(i, j) ==
                cell(i - 1, j - 1) + substitution(center[i - 1], query[j - 1])) {
            aligned_center.push_back(center[--i]);
            aligned_query.push_back(query[--j]);
        } else if (i != 0 && cell(i, j) == cell(i - 1, j) + detail::kGapScore) {
            aligned_center.push_back(center[--i]);
            aligned_query.push_back('-');
        } else if (j != 0 && cell(i, j) == cell(i, j - 1) + detail::kGapScore) {
            aligned_center.push_back('-');
            aligned_query.push_back(query[--j]);
        } else {
            throw detail::Error(detail::kCoreError, "invalid pairwise traceback state");
        }
    }
    std::reverse(aligned_center.begin(), aligned_center.end());
    std::reverse(aligned_query.begin(), aligned_query.end());

    Projection projection;
    projection.at_center.assign(center.size(), '-');
    projection.insertions.resize(height);
    std::size_t consumed = 0;
    for (std::size_t column = 0; column < aligned_center.size(); ++column) {
        if (aligned_center[column] == '-') {
            // Slot s follows exactly s consumed center residues. This also
            // accounts for prefix (0) and suffix (center length) insertions.
            projection.insertions[consumed].push_back(aligned_query[column]);
        } else {
            projection.at_center[consumed++] = aligned_query[column];
        }
    }
    return projection;
}

constexpr std::string_view kProfileSymbols = "ACGTUN-";
constexpr auto kProfileGap = kProfileSymbols.size() - 1;
constexpr Score kProfileGapScore = -2;
using ColumnProfile = std::array<Score, kProfileSymbols.size()>;
constexpr std::size_t kRefinementRadius = 4;
constexpr std::size_t kRefinementBandWidth = 2 * kRefinementRadius + 1;
constexpr Score kUnreachable = std::numeric_limits<Score>::min();

Score placement_gain(char residue, const ColumnProfile& profile) {
    // Compare placing this residue with leaving a gap in this column.
    Score gain = kProfileGapScore * profile[kProfileGap];
    for (std::size_t symbol = 0; symbol < kProfileGap; ++symbol) {
        gain += (substitution(residue, kProfileSymbols[symbol]) - kProfileGapScore) *
                profile[symbol];
    }
    return gain;
}

struct RefinementStep {
    Score score = kUnreachable;
    std::size_t predecessor = 0;
};

struct RefinementBand {
    std::size_t first;
    std::size_t last;
    std::array<RefinementStep, kRefinementBandWidth> steps{};
};

void refine_row(std::string_view sequence, std::string& row,
                const std::vector<ColumnProfile>& profile) {
    if (sequence.size() > std::vector<RefinementBand>().max_size()) {
        throw detail::Error(detail::kCoreError, "refinement matrix exceeds container size");
    }
    std::vector<RefinementBand> bands;
    bands.reserve(sequence.size());
    Score original_score = 0;
    for (std::size_t column = 0; column < row.size(); ++column) {
        if (row[column] == '-') continue;
        original_score += placement_gain(row[column], profile[column]);
        // Form the band without unsigned underflow or overflow at the ends.
        const auto first = column - std::min(column, kRefinementRadius);
        const auto last = column + std::min(row.size() - 1 - column, kRefinementRadius);
        bands.push_back({first, last});
    }

    for (std::size_t residue = 0; residue < sequence.size(); ++residue) {
        auto& band = bands[residue];
        Score best_prefix = residue == 0 ? 0 : kUnreachable;
        std::size_t best_predecessor = 0;
        std::size_t next_predecessor = 0;
        for (std::size_t column = band.first; column <= band.last; ++column) {
            if (residue != 0) {
                const auto& previous = bands[residue - 1];
                // Scan strictly earlier columns, retaining the first column on ties.
                while (next_predecessor <= previous.last - previous.first &&
                       previous.first + next_predecessor < column) {
                    const auto previous_score = previous.steps[next_predecessor].score;
                    if (previous_score > best_prefix) {
                        best_prefix = previous_score;
                        best_predecessor = next_predecessor;
                    }
                    ++next_predecessor;
                }
            }
            if (best_prefix != kUnreachable) {
                band.steps[column - band.first] = {
                    best_prefix + placement_gain(sequence[residue], profile[column]),
                    best_predecessor};
            }
        }
    }

    const auto& last_steps = bands.back().steps;
    const auto best = std::max_element(last_steps.begin(), last_steps.end(),
        [](const auto& a, const auto& b) { return a.score < b.score; });
    if (best->score <= original_score) return;
    auto offset = static_cast<std::size_t>(best - last_steps.begin());
    std::fill(row.begin(), row.end(), '-');
    for (std::size_t end = sequence.size(); end != 0; --end) {
        const auto index = end - 1;
        const auto& band = bands[index];
        row[band.first + offset] = sequence[index];
        offset = band.steps[offset].predecessor;
    }
}

void refine_profile(const std::vector<SequenceRecord>& input, Alignment& alignment) {
    const auto width = alignment.rows.front().size();
    // Every placement contributes at most three points per other row in
    // absolute value. A candidate places at most width residues.
    constexpr auto limit = static_cast<std::uint64_t>(std::numeric_limits<Score>::max() / 3);
    if (alignment.rows.size() - 1 > limit / width) {
        throw detail::Error(detail::kCoreError, "refinement score exceeds signed range");
    }
    if (width > std::vector<ColumnProfile>().max_size()) {
        throw detail::Error(detail::kCoreError, "refinement profile exceeds container size");
    }
    std::vector<ColumnProfile> profile(width);
    const auto update_profile = [&](const std::string& row, Score count) {
        for (std::size_t column = 0; column < width; ++column) {
            profile[column][kProfileSymbols.find(row[column])] += count;
        }
    };
    for (const auto& row : alignment.rows) update_profile(row, 1);

    for (std::size_t i = 1; i < alignment.rows.size(); ++i) {
        auto& row = alignment.rows[i];
        if (input[i].sequence.size() == width) continue;
        update_profile(row, -1);
        refine_row(input[i].sequence, row, profile);
        update_profile(row, 1);
    }

    const auto row_count = static_cast<Score>(alignment.rows.size());
    for (auto& row : alignment.rows) {
        std::size_t kept = 0;
        for (std::size_t column = 0; column < width; ++column) {
            if (profile[column][kProfileGap] != row_count) row[kept++] = row[column];
        }
        row.resize(kept);
    }
}
}  // namespace

// Reference star alignment — replaceable core
// Replace this implementation to use another algorithm. Keep input unchanged,
// preserve row-to-input correspondence, and return only aligned strings.
Alignment align(const std::vector<SequenceRecord>& input) {
    if (input.empty()) {
        throw detail::Error(detail::kCoreError, "the alignment core requires at least one record");
    }
    for (const auto& record : input) {
        if (record.sequence.empty() ||
            !std::all_of(record.sequence.begin(), record.sequence.end(), is_residue)) {
            throw detail::Error(detail::kCoreError,
                                "the alignment core requires nonempty cleaned ACGTUN sequences");
        }
    }
    if (input.size() == 1) return Alignment{{input.front().sequence}};

    const auto& center = input.front().sequence;
    const auto slot_count = detail::checked_add(center.size(), 1);
    std::vector<std::size_t> slot_widths(slot_count, 0);
    std::vector<Projection> projections;
    projections.reserve(input.size() - 1);
    for (std::size_t row = 1; row < input.size(); ++row) {
        auto projection = align_pair(center, input[row].sequence);
        for (std::size_t slot = 0; slot < slot_count; ++slot) {
            slot_widths[slot] = std::max(slot_widths[slot], projection.insertions[slot].size());
        }
        projections.push_back(std::move(projection));
    }

    auto output_width = center.size();
    for (const auto count : slot_widths) output_width = detail::checked_add(output_width, count);
    (void)detail::checked_multiply(input.size(), output_width);
    if (output_width > std::string().max_size()) {
        throw detail::Error(detail::kCoreError, "alignment row exceeds string size");
    }

    Alignment result;
    result.rows.resize(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        auto& output = result.rows[row];
        output.reserve(output_width);
        const Projection* projection = row == 0 ? nullptr : &projections[row - 1];
        for (std::size_t slot = 0; slot < slot_count; ++slot) {
            std::size_t used = 0;
            if (projection != nullptr) {
                const auto& insertion = projection->insertions[slot];
                output.append(insertion);
                used = insertion.size();
            }
            // Each insertion string is left-aligned within its slot. A longer
            // query determines the width; shorter queries receive right gaps.
            output.append(slot_widths[slot] - used, '-');
            if (slot != center.size()) {
                output.push_back(projection == nullptr ? center[slot] : projection->at_center[slot]);
            }
        }
    }
    refine_profile(input, result);
    return result;
}
}  // namespace msakit

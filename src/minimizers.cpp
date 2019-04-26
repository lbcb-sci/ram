/*!
 * @file minimizer.cpp
 *
 * @brief Minimizer class source file
 */

#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <algorithm>

#include "minimizers.hpp"

namespace ram {

std::vector<uint8_t> coder = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255,   0, 255,   1, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255,   3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255,   0, 255,   1, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255,   3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
};

void createMinimizers(std::vector<std::pair<std::uint64_t, std::uint64_t>>& dst,
    const char* sequence, std::uint32_t sequence_length, std::uint32_t id,
    std::uint32_t k, std::uint32_t w) {

    if (k > 32) {
        throw std::invalid_argument("[ram::createMinimizers] error: "
            "invalid kmer size!");
    }

    if (sequence_length < k) {
        return;
    }

    uint64_t mask = (1ULL << (k * 2)) - 1;
    uint64_t shift = (k - 1) * 2;
    uint64_t minimizer = 0, reverse_minimizer = 0;

    std::deque<std::pair<std::uint64_t, std::uint64_t>> window;
    auto window_add = [&window](std::uint64_t value, std::uint64_t location) -> void {
        while (!window.empty() && window.back().first > value) {
            window.pop_back();
        }
        window.emplace_back(value, location);
    };
    auto window_update = [&window](std::uint32_t position) -> void {
        while (!window.empty() && (window.front().second << 32 >> 33) < position) {
            window.pop_front();
        }
    };

    std::uint64_t t = static_cast<std::uint64_t>(id) << 32;

    for (std::uint32_t i = 0; i < sequence_length; ++i) {
        std::uint64_t c = coder[sequence[i]];
        if (c == 255) {
            throw std::invalid_argument("[ram::createMinimizers] error: "
                "invalid character!");
        }
        minimizer = ((minimizer << 2) | c) & mask;
        reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
        if (i >= (k - 1) + w) {
            if (dst.empty() ||
                dst.back().second != window.front().second) {
                dst.emplace_back(window.front());
            }
            window_update(i - (k - 1) - (w - 1));
        }
        if (i >= k - 1) {
            if (minimizer < reverse_minimizer) {
                window_add(minimizer, t | ((i - (k - 1)) << 1 | 0));
            } else if (minimizer > reverse_minimizer) {
                window_add(reverse_minimizer, t | ((i - (k - 1)) << 1 | 1));
            }
        }
    }
}

void sortMinimizers(std::vector<std::pair<std::uint64_t, std::uint64_t>>& src,
    std::uint32_t k) {

    if (src.empty()) {
        throw std::invalid_argument("[ram::sortMinimizers] error: "
            "empty minimizer set");
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> dst(src.size());
    std::uint32_t buckets[0x100] = {};
    std::uint32_t max_shift = ((2 * k + 7) / 8) * 8;

    for (std::uint32_t shift = 0; shift < max_shift; shift += 8) {
        std::uint32_t counts[0x100] = {};
        for (const auto& it: src) {
            ++counts[(it.first >> shift) & 0xFF];
        }
        for (std::uint32_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
            buckets[i] = j;
        }
        for (const auto& it: src) {
            dst[buckets[(it.first >> shift) & 0xFF]++] = it;
        }
        src.swap(dst);
    }
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> map(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& query,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& target,
    const std::unordered_map<std::uint64_t, std::pair<std::uint32_t, std::uint32_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence) {

    std::vector<std::pair<std::uint64_t, std::uint64_t>> matches;

    for (std::uint32_t i = 0; i < query.size(); i++) {

        const auto it = target_hash.find(query[i].first);
        if (it != target_hash.end()) {
            const auto& range = it->second;
            if (range.second >= max_occurence) {
                continue;
            }
            for (std::uint32_t j = range.first; j < range.first + range.second; j++) {
                std::uint64_t strand = (query[i].second & 1) == (target[j].second & 1);

                std::uint64_t query_id = query[i].second >> 32;
                std::uint64_t query_pos = query[i].second << 32 >> 33;

                if (id + 1 == query_id) {
                    query_pos += offset;
                    query_id -= 1;
                }

                std::uint64_t target_id = target[j].second >> 32;
                std::uint64_t target_pos = target[j].second << 32 >> 33;

                // TODO: take minimizers from first/last non-singleton minimizer!
                if (query_id <= target_id) {
                    break;
                }

                std::uint64_t diagonal_diff = !strand ? target_pos + query_pos :
                    (1ULL << 31) + target_pos - query_pos;

                matches.emplace_back(
                    (((target_id << 1) | strand) << 32) | diagonal_diff,
                    (target_pos << 32) | query_pos);
            }
        }
    }

    return matches;
}

bool is_contained(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& query,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& target,
    const std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence,
    const std::vector<std::uint32_t>& sequence_lengths) {

    auto matches = ram::map(query, target, target_hash, id, offset, max_occurence);

    if (matches.size() <= 0) {
        return false;
    }

    std::sort(matches.begin(), matches.end());

    auto op_less = std::less<std::uint64_t>();
    auto op_greater = std::greater<std::uint64_t>();

    for (std::uint32_t i = 1, j = 0; i < matches.size(); ++i) {
        if ((matches[i].first >> 32) != (matches[i - 1].first >> 32) ||
            (matches[i].first << 32 >> 32) - (matches[i - 1].first << 32 >> 32) > 500) {

            if (i - j < 4) {
                j = i;
                continue;
            }

            std::sort(matches.begin() + j, matches.begin() + i,
                [] (const std::pair<std::uint64_t, std::uint64_t>& lhs,
                    const std::pair<std::uint64_t, std::uint64_t>& rhs) {
                    return lhs.second < rhs.second;
                }
            );

            std::vector<std::uint32_t> indices;
            bool strand = (matches[i - 1].first >> 32) & 1;
            if (strand) {
                indices = ram::longestSubsequence(matches.begin() + j,
                    matches.begin() + i, op_less);
            } else {
                indices = ram::longestSubsequence(matches.begin() + j,
                    matches.begin() + i, op_greater);
            }

            // TODO: only check for groups in diag_eps of 500
            // TODO: check lis if there are multiple target_pos
            std::uint64_t target_begin = matches[j + indices.front()].second >> 32;
            std::uint64_t target_end = 15 + (matches[j + indices.back()].second >> 32);
            std::uint64_t query_begin = strand ?
                matches[j + indices.front()].second << 32 >> 32 :
                matches[j + indices.back()].second << 32 >> 32;
            std::uint64_t query_end = 15 + (strand ?
                matches[j + indices.back()].second << 32 >> 32 :
                matches[j + indices.front()].second << 32 >> 32);

            std::uint32_t target_id = matches[i - 1].first >> 33;
            std::uint32_t overhang = std::min(target_begin, query_begin) + std::min(
                sequence_lengths[id] - query_end,
                sequence_lengths[target_id] - target_end);

            if (target_end - target_begin > (target_end - target_begin + overhang) * 0.875 &&
                query_end - query_begin > (query_end - query_begin + overhang) * 0.875 &&
                sequence_lengths[target_id] - target_end >= sequence_lengths[id] - query_end &&
                target_begin >= query_begin) {

                return true;
            }

            j = i;
        }
    }

    return false;
}

}

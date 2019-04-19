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

void create_clusters(
    std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> &clusters,
    std::vector<std::uint32_t> &indices,
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator begin) {

    if (indices.size() <= 0) {
        return;
    }

    std::uint32_t start_ref = (begin+indices[0])->second >> 32;
    std::uint32_t len_ref = 0;
    std::uint32_t start_read = ((begin+indices[0])->second << 32) >> 32;
    std::uint32_t len_read = 0;

    auto previous_it = indices.begin();
    auto current_it = previous_it + 1;

    while(current_it != indices.end()) {
        auto match = (begin + (*current_it));
        auto previous_match = (begin + (*previous_it));

        std::uint32_t strand = ((previous_match->first >> 32) & 1);

        std::uint64_t ref_len = (match->second >> 32) - (previous_match->second >> 32);
        std::uint64_t read_len = strand ? ((match->second << 32) >> 32) - ((previous_match->second << 32) >> 32):
            ((previous_match->second << 32) >> 32) - ((match->second << 32) >> 32);

        if (ref_len > 1000) {
            auto adjusted_read_start = !strand ? start_read : ((match->second << 32) >> 32);
            clusters.emplace_back(start_ref, len_ref+15, adjusted_read_start, len_read+15);

            start_ref = match->second >> 32;
            start_read = (match->second << 32) >> 32;
            len_ref = 0;
            len_read = 0;
        } else {
            len_ref += ref_len;
            len_read += read_len;
        }

        current_it += 1;
        previous_it += 1;
    }

    if (len_ref > 0) {
        clusters.emplace_back(start_ref, len_ref+15, start_read, len_read+15);
    }

    return;
}

bool are_clusters_contained(std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> &clusters) {

    auto start_it = clusters.begin();

    std::uint32_t aligned_read_len = 0;

    std::uint32_t min_query_pos = (0-1);
    std::uint32_t min_target_pos = (0-1);

    std::uint32_t max_query_pos = 0;
    std::uint32_t max_target_pos = 0;

    while (start_it != clusters.end()) {
        min_query_pos = std::min(min_query_pos, std::get<2>(*start_it));
        min_target_pos = std::min(min_target_pos, std::get<0>(*start_it));

        max_query_pos = std::max(max_query_pos, std::get<2>(*start_it) + std::get<3>(*start_it));
        max_target_pos = std::max(max_target_pos, std::get<0>(*start_it) + std::get<1>(*start_it));

        aligned_read_len += std::get<3>(*start_it);

        start_it += 1;
    }

    return min_query_pos < min_target_pos && max_query_pos < max_target_pos && aligned_read_len > 1300;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> map(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& query,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& target,
    std::unordered_map<std::uint64_t, std::pair<std::uint32_t, std::uint32_t>>& target_hash,
    std::uint32_t second_sequence_offset, std::uint32_t id,
    std::uint32_t max_occurence) {

    std::vector<std::pair<std::uint64_t, std::uint64_t>> matches;

    for (std::uint32_t i = 0; i < query.size(); i++) {

        const auto it = target_hash.find(query[i].first);
        if (it != target_hash.end()) {
            const auto& range = it->second;
            if (false && range.second >= max_occurence) {
                continue;
            }
            for (std::uint32_t j = range.first; j < range.first + range.second; j++) {
                std::uint64_t strand = (query[i].second & 1) == (target[j].second & 1);

                std::uint64_t query_id = query[i].second >> 32;
                std::uint64_t query_pos = query[i].second << 32 >> 33;

                if (id + 1 == query_id) {
                    query_pos += second_sequence_offset;
                    query_id -= 1;
                }

                std::uint64_t target_id = target[j].second >> 32;
                std::uint64_t target_pos = target[j].second << 32 >> 33;

                // fix to check lenghts or do nothing, idk
                //if (query_id >= target_id) {
                //    continue;
                //}

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

bool is_read_contained(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& query,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& target,
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>>& target_hash,
    std::uint32_t second_sequence_offset, std::uint32_t id,
    std::uint32_t max_occurence) {

    auto matches = ram::map(query, target, target_hash, second_sequence_offset,
        id, max_occurence);

    if (matches.size() <= 0) {
        return false;
    }

    std::sort(matches.begin(), matches.end());

    auto op_less = std::less<std::uint64_t>();
    auto op_greater = std::greater<std::uint64_t>();

    std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> clusters;

    for (unsigned i = 1, j = 0; i < matches.size(); ++i) {
        if ((matches[i].first >> 32) != matches[i - 1].first >> 32 ||
            (matches[i].first << 32 >> 32) - (matches[i - 1].first << 32 >> 32) > 500 ||
            i == matches.size() - 1) {

            std::sort(matches.begin() + j, matches.begin() + i,
                [](const std::pair<std::uint64_t, std::uint64_t>& lhs,
                   const std::pair<std::uint64_t, std::uint64_t>& rhs) {
                    return lhs.second < rhs.second;
                }
            );

            std::vector<std::uint32_t> indices;
            if ((matches[i - 1].first >> 32) & 1) {
                indices = ram::longestSubsequence(matches.begin() + j,
                    matches.begin() + i, op_less);
            } else {
                indices = ram::longestSubsequence(matches.begin() + j,
                    matches.begin() + i, op_greater);
            }

            create_clusters(clusters, indices, matches.begin() + j);

            if ((matches[i].first >> 32) != matches[i - 1].first >> 32) {
                if (are_clusters_contained(clusters)) {
                    return true;
                }
                clusters.clear();
            }

            j = i;
        }
    }

    return false;
}

}

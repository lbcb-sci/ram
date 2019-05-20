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

#include "thread_pool/thread_pool.hpp"

#include "minimizers.hpp"

namespace ram {

std::vector<std::uint8_t> coder = {
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

void createMinimizers(std::vector<uint128_t>& dst, const char* sequence,
    std::uint32_t sequence_length, std::uint32_t sequence_id,
    std::uint32_t k, std::uint32_t w) {

    if (k > 32) {
        throw std::invalid_argument("[ram::createMinimizers] error: "
            "invalid kmer size!");
    }

    if (sequence_length < k) {
        return;
    }

    std::uint64_t mask = (1ULL << (k * 2)) - 1;
    std::uint64_t shift = (k - 1) * 2;
    std::uint64_t minimizer = 0, reverse_minimizer = 0;

    std::deque<uint128_t> window;
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

    std::uint64_t id = static_cast<std::uint64_t>(sequence_id) << 32;

    for (std::uint32_t i = 0; i < sequence_length; ++i) {
        std::uint64_t c = coder[sequence[i]];
        if (c == 255) {
            throw std::invalid_argument("[ram::createMinimizers] error: "
                "invalid character!");
        }
        minimizer = ((minimizer << 2) | c) & mask;
        reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
        if (i >= k - 1) {
            if (minimizer < reverse_minimizer) {
                window_add(minimizer, id | ((i - (k - 1)) << 1 | 0));
            } else if (minimizer > reverse_minimizer) {
                window_add(reverse_minimizer, id | ((i - (k - 1)) << 1 | 1));
            }
        }
        if (i >= (k - 1) + (w - 1)) {
            if (dst.empty() ||
                dst.back().second != window.front().second) {
                dst.emplace_back(window.front());
            }
            window_update(i - (k - 1) - (w - 1) + 1);
        }
    }
}

void sortMinimizers(std::vector<uint128_t>& src, std::uint32_t k) noexcept {

    if (src.empty()) {
        return;
    }

    std::vector<uint128_t> dst(src.size());
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

void transformMinimizers(std::vector<std::vector<uint128_t>>& hash,
    std::vector<std::unordered_map<std::uint64_t, uint128_t>>& index,
    std::vector<std::vector<uint128_t>>& minimizers, std::uint32_t k,
    const std::unique_ptr<thread_pool::ThreadPool>& thread_pool) {

    hash.clear();
    hash.resize(0x4000);

    index.clear();
    index.resize(hash.size());

    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
        for (const auto& it: minimizers[i]) {
            hash[it.first & 0x3FFF].emplace_back(it);
        }
        std::vector<uint128_t>().swap(minimizers[i]);
    }

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < hash.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&hash, &index, k] (std::uint32_t i) -> void {
                if (hash[i].empty()) {
                    return;
                }
                ram::sortMinimizers(hash[i], k);
                for (std::uint32_t j = 0, count = 0; j < hash[i].size(); ++j) {
                    if (j > 0 && hash[i][j - 1].first != hash[i][j].first) {
                        index[i][hash[i][j - 1].first] = std::make_pair(j - count, count);
                        count = 0;
                    }
                    if (j == hash[i].size() - 1) {
                        index[i][hash[i][j].first] = std::make_pair(j - count, count + 1);
                    }
                    ++count;
                }
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();
}

std::vector<uint128_t> map(const std::vector<uint128_t>& query,
    const std::vector<std::vector<uint128_t>>& target,
    const std::vector<std::unordered_map<std::uint64_t, uint128_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence,
    bool diagonal, bool triangle) {

    std::vector<uint128_t> matches;

    for (std::uint32_t i = 0; i < query.size(); i++) {
        std::uint32_t bin = query[i].first & 0x3FFF;
        const auto it = target_hash[bin].find(query[i].first);
        if (it != target_hash[bin].end()) {
            const auto& range = it->second;
            if (range.second >= max_occurence) {
                continue;
            }
            for (std::uint32_t j = range.first; j < range.first + range.second; j++) {
                std::uint64_t strand = (query[i].second & 1) == (target[bin][j].second & 1);

                std::uint64_t query_id = query[i].second >> 32;
                std::uint64_t query_pos = query[i].second << 32 >> 33;

                if (id + 1 == query_id) {
                    query_pos += offset;
                    query_id -= 1;
                }

                std::uint64_t target_id = target[bin][j].second >> 32;
                std::uint64_t target_pos = target[bin][j].second << 32 >> 33;

                if (diagonal && query_id == target_id) {
                    continue;
                }
                if (triangle && query_id < target_id) {
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

bool is_contained(const std::vector<uint128_t>& query,
    const std::vector<std::vector<uint128_t>>& target,
    const std::vector<std::unordered_map<std::uint64_t, uint128_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence,
    std::uint32_t query_length, const std::vector<std::uint32_t>& sequence_lengths) {

    auto matches = ram::map(query, target, target_hash, id, offset,
        max_occurence, true, true);

    if (matches.empty()) {
        return false;
    }

    std::sort(matches.begin(), matches.end());

    auto op_less = std::less<std::uint64_t>();
    auto op_greater = std::greater<std::uint64_t>();

    matches.emplace_back(-1, -1); // stop dummy
    for (std::uint32_t i = 1, j = 0; i < matches.size(); ++i) {
        if ((matches[i].first >> 32) != (matches[i - 1].first >> 32) ||
            (matches[i].first << 32 >> 32) - (matches[i - 1].first << 32 >> 32) > 500) {

            if (i - j < 4) {
                j = i;
                continue;
            }

            std::sort(matches.begin() + j, matches.begin() + i,
                [] (const uint128_t& lhs, const uint128_t& rhs) {
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
                query_length - query_end, sequence_lengths[target_id] - target_end);

            if (target_end - target_begin > (target_end - target_begin + overhang) * 0.875 &&
                query_end - query_begin > (query_end - query_begin + overhang) * 0.875 &&
                sequence_lengths[target_id] - target_end >= query_length - query_end &&
                target_begin >= query_begin) {

                return true;
            }

            j = i;
        }
    }

    return false;
}

}

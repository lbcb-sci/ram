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

    uint64_t mask = (1 << (k * 2)) - 1;
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

    for (std::uint32_t i = 0; i < sequence_length - k + 1; ++i) {
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
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& lhs,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& rhs,
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>>& hash) {

    std::vector<std::pair<std::uint64_t, std::uint64_t>> matches;

    for (std::uint32_t i = 0; i < lhs.size(); i++) {

        auto value = hash.find(lhs[i].first);
        if (value != hash.end()) {
            auto range = value->second;
            for (std::uint32_t j = range.first; j < range.first + range.second; j++) {
                std::uint32_t strand = (lhs[i].second & 1) == (rhs[j].second & 1);
                std::uint64_t query_id = lhs[i].second >> 32;
                std::uint64_t target_id = rhs[j].second >> 32;

                if (query_id >= target_id) {
                    continue;
                }

                std::uint32_t query_pos = (static_cast<std::uint32_t>(lhs[i].second) >> 1);
                std::uint32_t target_pos = (static_cast<std::uint32_t>(rhs[j].second) >> 1);

                // std::uint32_t bounds = (1<<31);
                // std::cout <<  bounds << std::endl;
                std::uint32_t diagonal_diff = strand ? (100000 + target_pos - query_pos) : target_pos+query_pos;

                // std::cout << "strand " << strand << std::endl;
                // std::cout << "(target_id << 1) " << (target_id << 1) << std::endl;
                std::uint64_t match_first = (((target_id << 1) | strand) << 32) | diagonal_diff;
                std::uint64_t match_second = (static_cast<std::uint64_t>(target_pos) << 32) | static_cast<std::uint64_t>(query_pos);

                // std::cout << match_first << " " << match_second << std::endl;
                matches.emplace_back(match_first, match_second);
            }
        }
    }
    return matches;
}

}

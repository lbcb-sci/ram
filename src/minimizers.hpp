/*!
 * @file minimizer.hpp
 *
 * @brief Minimizer class header file
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <unordered_map>

namespace thread_pool {
    class ThreadPool;
}

namespace ram {

using uint128_t = std::pair<std::uint64_t, std::uint64_t>;

// [127: 64] := minimizer
// [ 63: 32] := identifier
// [ 31:  1] := position
// [  0:  0] := strand
void createMinimizers(std::vector<uint128_t>& dst, const char* sequence,
    std::uint32_t sequence_length, std::uint32_t sequence_id,
    std::uint32_t k, std::uint32_t w);

void sortMinimizers(std::vector<uint128_t>& src, std::uint32_t k) noexcept;

void transformMinimizers(std::vector<std::vector<uint128_t>>& hash,
    std::vector<std::unordered_map<std::uint64_t, uint128_t>>& index,
    std::vector<std::vector<uint128_t>>& minimizers, std::uint32_t k,
    const std::unique_ptr<thread_pool::ThreadPool>& thread_pool);

template<typename Operator>
inline std::vector<std::uint32_t> longestSubsequence(
    std::vector<uint128_t>::const_iterator begin,
    std::vector<uint128_t>::const_iterator end,
    const Operator& op) {

    if (begin >= end) {
        throw std::invalid_argument("[ram::longestSubsequence] error: "
            "empty match set");
    }

    std::vector<std::uint32_t> smallest(end - begin + 1, 0);
    std::vector<std::uint32_t> predecessor(end - begin, 0);

    // TODO: not sure about double lis
    auto op_less = std::less<std::uint64_t>();

    std::uint32_t length = 0;
    for (auto it = begin; it != end; ++it) {
        std::uint32_t l = 1, h = length;
        while (l <= h) {
            std::uint32_t m = (l + h) >> 1;
            if (op_less((begin + smallest[m])->second >> 32, it->second >> 32) &&
                op((begin + smallest[m])->second << 32 >> 32, it->second << 32 >> 32)) {
                l = m + 1;
            } else {
                h = m - 1;
            }
        }

        predecessor[it - begin] = smallest[l - 1];
        smallest[l] = it - begin;
        length = std::max(length, l);
    }

    std::vector<std::uint32_t> dst;
    dst.reserve(length);
    for (std::uint32_t i = 0, j = smallest[length]; i < length; ++i, j = predecessor[j]) {
        dst.emplace_back(j);
    }
    std::reverse(dst.begin(), dst.end());

    return dst;
}

// [127: 97] := identifier
// [ 96: 96] := relative strand
// [ 95: 64] := diagonal identifier
// [ 63: 32] := target position
// [ 31:  0] := query position
std::vector<uint128_t> map(const std::vector<uint128_t>& query,
    const std::vector<std::vector<uint128_t>>& target,
    const std::vector<std::unordered_map<std::uint64_t, uint128_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence,
    bool diagonal, bool triangle);

bool is_contained(const std::vector<uint128_t>& query,
    const std::vector<std::vector<uint128_t>>& target,
    const std::vector<std::unordered_map<std::uint64_t, uint128_t>>& target_hash,
    std::uint32_t id, std::uint32_t offset, std::uint32_t max_occurence,
    std::uint32_t query_length, const std::vector<std::uint32_t>& sequence_lengths);

}

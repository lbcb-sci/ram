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

namespace ram {

void createMinimizers(std::vector<std::pair<std::uint64_t, std::uint64_t>>& dst,
    const char* sequence, std::uint32_t sequence_length, std::uint32_t id,
    std::uint32_t k, std::uint32_t w);

void sortMinimizers(std::vector<std::pair<std::uint64_t, std::uint64_t>>& src,
    std::uint32_t k);

template<typename Operator>
inline std::vector<std::uint32_t> longestIncreasingSubsequence(
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator begin,
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator end,
    Operator op) {

    if (begin == end) {
        throw std::invalid_argument("[ram::longestIncreasingSubsequence] error: "
            "empty match set");
    }

    std::vector<std::uint32_t> smallest(end - begin + 1, 0);
    std::vector<std::uint32_t> predecessor(end - begin, 0);

    std::uint32_t length = 0;
    for (auto it = begin; it != end; ++it) {
        std::uint32_t l = 1, h = length;
        while (l <= h) {
            std::uint32_t m = (l + h) >> 1;
            if (op((begin + smallest[m])->second << 32 >> 32, it->second << 32 >> 32)) {
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

std::vector<std::pair<std::uint64_t, std::uint64_t>> map(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& lhs,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& rhs,
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>>& hash);

}

struct Minimizer {
    std::uint64_t value;
    std::uint32_t location;
    bool strand;

    Minimizer(std::uint64_t value, std::uint64_t location, bool strand) {
        this->value = value;
        this->location = location;
        this->strand = strand;
    }

    bool operator<(Minimizer other) const {
        return value < other.value;
    }

    std::string to_string() const {
        return "Location: " + std::to_string(location >> 1) + ", value: " +
            std::to_string(value) + ", strand: " + std::to_string(location & 1);
    }
};

struct Match {
    std::uint64_t ref_start;
    std::uint64_t read_start;
    bool strand;

    Match(Minimizer read_minimizer, Minimizer ref_minimizer) {
        this->ref_start = ref_minimizer.location;
        this->read_start = read_minimizer.location;
        this->strand = read_minimizer.strand != ref_minimizer.strand;
    }

    bool operator<(Match other) const {
        return ref_start < other.ref_start;
    }

    std::string to_string() const {
        return "ReadStart: " + std::to_string(read_start) + ", refStart: " +
            std::to_string(ref_start) + ", strand: " + std::to_string(strand);
    }
};

struct Cluster {
    std::uint64_t ref_start;
    std::uint64_t read_start;
    std::uint64_t ref_length;
    std::uint64_t read_length;
    bool strand;

    Cluster(std::uint64_t ref_start, std::uint64_t ref_length,
        std::uint64_t read_start, std::uint64_t read_length, bool strand) {

        this->ref_start = ref_start;
        this->read_start = read_start;
        this->ref_length = ref_length;
        this->read_length = read_length;
        this->strand = strand;
    }

    std::string to_string() const {
        return "RefStart: " + std::to_string(read_start) + " refLen: " +
            std::to_string(read_length) + ", readStart: " + std::to_string(ref_start) +
            ", readLen: " + std::to_string(ref_length) + ", strand: " +
            std::to_string(strand);
    }
};

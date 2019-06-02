/*!
 * @file minimizer_engine.hpp
 *
 * @brief MinimizerEngine class header file
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <unordered_map>

namespace thread_pool {
    class ThreadPool;
}

namespace ram {

struct Sequence {
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length);
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length,
        const char*, std::uint32_t);
    ~Sequence() = default;

    static std::uint32_t sequence_id;

    std::uint32_t id;
    std::string name;
    std::string data;
};

struct Overlap {
    Overlap(std::uint32_t t_id, std::uint32_t t_begin, std::uint32_t t_end,
        std::uint32_t q_begin, std::uint32_t q_end, std::uint32_t strand);
    ~Overlap() = default;

    std::uint32_t t_id;
    std::uint32_t t_begin;
    std::uint32_t t_end;
    std::uint32_t q_begin;
    std::uint32_t q_end;
    std::uint32_t strand;
};

using uint128_t = std::pair<std::uint64_t, std::uint64_t>;

class MinimizerEngine {
public:
    MinimizerEngine(std::uint8_t k, std::uint8_t w,
        std::uint32_t num_threads = 1);
    ~MinimizerEngine() = default;

    void minimize(const std::vector<std::unique_ptr<Sequence>>& src,
        double f = 0.001);
    void minimize(
        std::vector<std::unique_ptr<Sequence>>::const_iterator begin,
        std::vector<std::unique_ptr<Sequence>>::const_iterator end,
        double f = 0.001);

    std::vector<Overlap> map(const std::unique_ptr<Sequence>& src,
        bool diagonal, bool triangle, std::uint32_t e = -1) const;
private:
    std::vector<uint128_t> minimize(const std::unique_ptr<Sequence>& src,
        std::uint32_t e = -1) const;

    void radix_sort(
        std::vector<uint128_t>::iterator begin,
        std::vector<uint128_t>::iterator end,
        bool first = true) const;

    std::vector<uint32_t> longest_subsequence(
        std::vector<uint128_t>::const_iterator begin,
        std::vector<uint128_t>::const_iterator end,
        bool increasing = true) const;

    struct MinimizerHash {
        void initialize();

        inline std::pair<std::vector<uint128_t>::const_iterator,
            std::vector<uint128_t>::const_iterator> operator[](std::uint64_t m) const;

        std::vector<std::vector<uint128_t>> minimizers;
        std::vector<std::unordered_map<std::uint64_t, uint128_t>> index;
    };

    std::uint8_t k_;
    std::uint8_t w_;
    std::uint32_t s_;
    MinimizerHash hash_;
    std::unique_ptr<thread_pool::ThreadPool> thread_pool_;
};

}

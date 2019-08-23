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
    class Threadpool;
}

namespace ram {

class Sequence;
class Overlap;
class MinimizerEngine;

using uint128_t = std::pair<std::uint64_t, std::uint64_t>;

std::unique_ptr<MinimizerEngine> createMinimizerEngine(std::uint8_t k,
    std::uint8_t w, std::shared_ptr<thread_pool::ThreadPool> thread_pool);

class MinimizerEngine {
public:
    ~MinimizerEngine() = default;

    /*!
     * @brief Transforms a set of sequences to a hash of minimizers without
     * the most frequent f.
     */
    void minimize(const std::vector<std::unique_ptr<Sequence>>& sequences,
        double f);
    void minimize(
        typename std::vector<std::unique_ptr<Sequence>>::const_iterator begin,
        typename std::vector<std::unique_ptr<Sequence>>::const_iterator end,
        double f);

    /*!
     * @brief Maps a sequence to a preconstructed minimizer hash while ignoring
     * self overlaps if d, and dual overlaps if t (both based on sequence id)
     */
    std::vector<Overlap> map(const std::unique_ptr<Sequence>& sequence, bool d,
        bool t, std::uint32_t e = -1) const;

    friend std::unique_ptr<MinimizerEngine> createMinimizerEngine(std::uint8_t k,
        std::uint8_t w, std::shared_ptr<thread_pool::ThreadPool> thread_pool);
private:
    MinimizerEngine(std::uint8_t k, std::uint8_t w,
        std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    MinimizerEngine(const MinimizerEngine&) = delete;
    const MinimizerEngine& operator=(const MinimizerEngine&) = delete;

    std::vector<uint128_t> minimize(const std::unique_ptr<Sequence>& sequence,
        std::uint32_t e = -1) const;

    std::uint8_t k_;
    std::uint8_t w_;
    std::uint64_t o_;
    std::vector<std::vector<uint128_t>> minimizers_;
    std::vector<std::unordered_map<std::uint64_t, uint128_t>> index_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}

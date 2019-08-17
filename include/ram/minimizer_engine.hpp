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

#include "thread_pool/thread_pool.hpp"

namespace ram {

struct Overlap {
    Overlap(std::uint32_t q_id, std::uint32_t q_begin, std::uint32_t q_end,
        std::uint32_t t_id, std::uint32_t t_begin, std::uint32_t t_end,
        std::uint32_t strand, std::uint32_t minimizers)
            : q_id(q_id), q_begin(q_begin), q_end(q_end), t_id(t_id),
            t_begin(t_begin), t_end(t_end), strand(strand), minimizers(minimizers) {
    }

    std::uint32_t q_id;
    std::uint32_t q_begin;
    std::uint32_t q_end;
    std::uint32_t t_id;
    std::uint32_t t_begin;
    std::uint32_t t_end;
    std::uint16_t strand;
    std::uint16_t minimizers;
}; // uint224_t

using uint128_t = std::pair<std::uint64_t, std::uint64_t>;

class MinimizerEngine;
std::unique_ptr<MinimizerEngine> createMinimizerEngine(std::uint8_t k,
    std::uint8_t w, std::shared_ptr<thread_pool::ThreadPool> thread_pool);

class MinimizerEngine {
public:
    ~MinimizerEngine() = default;

    /*!
     * @brief Transforms a set of sequences to a hash of minimizers without
     * the most frequent f.
     */
    template<typename T>
    void minimize(const std::vector<T>& sequences, double f) {
        minimize<T>(sequences.begin(), sequences.end(), f);
    }
    template<typename T>
    void minimize(
        typename std::vector<T>::const_iterator begin,
        typename std::vector<T>::const_iterator end,
        double f) {

        hash_.clear();

        if (begin >= end) {
            return;
        }

        std::vector<std::vector<uint128_t>> minimizers(end - begin);
        std::vector<std::future<void>> thread_futures;
        for (auto it = begin; it != end; ++it) {
            thread_futures.emplace_back(thread_pool_->submit(
                [&] (typename std::vector<T>::const_iterator it) -> void {
                    minimizers[it - begin] = minimize((*it)->id(), (*it)->data());
                }
            , it));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }
        thread_futures.clear();

        return transform(minimizers, f);
    }

    /*!
     * @brief Maps a sequence to a preconstructed minimizer hash while ignoring
     * self overlaps if d, and dual overlaps if t
     */
    template<typename T>
    std::vector<Overlap> map(const T& sequence, bool d, bool t, std::uint32_t e = -1) const {
        auto minimizers = minimize(sequence->id(), sequence->data(), e);
        return match(minimizers, d, t);
    }

    friend std::unique_ptr<MinimizerEngine> createMinimizerEngine(std::uint8_t k,
        std::uint8_t w, std::shared_ptr<thread_pool::ThreadPool> thread_pool);
private:
    MinimizerEngine(std::uint8_t k, std::uint8_t w,
        std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    MinimizerEngine(const MinimizerEngine&) = delete;
    const MinimizerEngine& operator=(const MinimizerEngine&) = delete;

    std::vector<uint128_t> minimize(std::uint32_t id, const std::string& data,
        std::uint32_t e = -1) const;

    void transform(std::vector<std::vector<uint128_t>>& minimizers, double f);

    std::vector<Overlap> match(std::vector<uint128_t>& minimizers, bool d, bool t) const;

    struct MinimizerHash {
        MinimizerHash(std::uint16_t num_bins);
        ~MinimizerHash() = default;

        void clear();

        inline std::pair<std::vector<uint128_t>::const_iterator,
            std::vector<uint128_t>::const_iterator> operator[](std::uint64_t minimizer) const;

        std::uint64_t m_;
        std::vector<std::vector<uint128_t>> minimizers;
        std::vector<std::unordered_map<std::uint64_t, uint128_t>> index;
    };

    std::uint8_t k_;
    std::uint8_t w_;
    std::uint64_t o_;
    MinimizerHash hash_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}

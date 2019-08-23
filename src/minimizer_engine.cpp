/*!
 * @file minimizer_engine.cpp
 *
 * @brief MinimizerEngine class source file
 */

#include <cstdlib>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <queue>

#include "thread_pool/thread_pool.hpp"

#include "ram/minimizer_engine.hpp"
#include "ram/overlap.hpp"
#include "ram/sequence.hpp"

namespace ram {

std::vector<std::uint8_t> kCoder = {
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

std::uint32_t Sequence::num_objects = 0;

inline std::uint64_t uint128_t_first(const uint128_t& src) {
    return src.first;
}
inline std::uint64_t uint128_t_second(const uint128_t& src) {
    return src.second;
}

template<typename F>
void radix_sort(std::vector<uint128_t>::iterator begin,
    std::vector<uint128_t>::iterator end, std::uint8_t num_bits, F func) {

    if (begin >= end) {
        return;
    }

    std::vector<uint128_t> dst(end - begin);
    auto b = dst.begin(), e = dst.end();
    std::uint64_t buckets[0x100] = {};

    std::uint8_t shift = 0;
    for (; shift < num_bits; shift += 8) {
        std::uint64_t counts[0x100] = {};
        for (auto it = begin; it != end; ++it) {
            ++counts[func(*it) >> shift & 0xFF];
        }
        for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
            buckets[i] = j;
        }
        for (auto it = begin; it != end; ++it) {
            *(b + buckets[func(*it) >> shift & 0xFF]++) = *it;
        }
        std::swap(b, begin);
        std::swap(e, end);
    }
    if (shift / 8 & 1) {
        for (; begin != end; ++begin, ++b) {
            *b = *begin;
        }
    }
}

template<typename F>
std::vector<uint64_t> longest_subsequence(std::vector<uint128_t>::const_iterator begin,
    std::vector<uint128_t>::const_iterator end, F func) {

    std::vector<std::uint64_t> dst;
    if (begin >= end) {
        return dst;
    }

    std::vector<std::uint64_t> smallest(end - begin + 1, 0);
    std::vector<std::uint64_t> predecessor(end - begin, 0);

    // TODO: not sure about double lis
    std::uint64_t length = 0;
    for (auto it = begin; it != end; ++it) {
        std::uint64_t l = 1, h = length;
        while (l <= h) {
            std::uint64_t m = (l + h) >> 1;
            if ((begin + smallest[m])->second >> 32 < it->second >> 32 &&
                func((begin + smallest[m])->second << 32 >> 32, it->second << 32 >> 32)) {
                l = m + 1;
            } else {
                h = m - 1;
            }
        }

        predecessor[it - begin] = smallest[l - 1];
        smallest[l] = it - begin;
        length = std::max(length, l);
    }

    dst.reserve(length);
    for (std::uint64_t i = 0, j = smallest[length]; i < length; ++i, j = predecessor[j]) {
        dst.emplace_back(j);
    }
    std::reverse(dst.begin(), dst.end());

    return dst;
}

std::unique_ptr<MinimizerEngine> createMinimizerEngine(std::uint8_t k,
    std::uint8_t w, std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    if (k == 0 || k > 32) {
        throw std::invalid_argument("[ram::createMinimizerEngine] error: "
            "invalid k-mer length!");
    }
    if (w == 0) {
        throw std::invalid_argument("[ram::createMinimizerEngine] error: "
            "invalid window length!");
    }
    if (thread_pool == nullptr) {
        throw std::invalid_argument("[ram::createMinimizerEngine] error: "
            "thread_pool is nullptr!");
    }

    return std::unique_ptr<MinimizerEngine>(new MinimizerEngine(k, w, thread_pool));
}

MinimizerEngine::MinimizerEngine(std::uint8_t k, std::uint8_t w,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
        : k_(k), w_(w), o_(-1), minimizers_(1U << std::min(14, 2 * k)),
        index_(minimizers_.size()), thread_pool_(thread_pool) {
}

void MinimizerEngine::minimize(
    const std::vector<std::unique_ptr<Sequence>>& sequences,
    double f) {

    minimize(sequences.begin(), sequences.end(), f);
}

void MinimizerEngine::minimize(
    typename std::vector<std::unique_ptr<Sequence>>::const_iterator begin,
    typename std::vector<std::unique_ptr<Sequence>>::const_iterator end,
    double f) {

    for (auto& it: minimizers_) {
        it.clear();
    }
    for (auto& it: index_) {
        it.clear();
    }

    if (begin >= end) {
        return;
    }

    std::vector<std::vector<uint128_t>> storage(end - begin);
    std::vector<std::future<void>> thread_futures;
    for (auto it = begin; it != end; ++it) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::vector<std::unique_ptr<Sequence>>::const_iterator it) -> void {
                storage[it - begin] = minimize((*it));
            }
        , it));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    std::uint64_t mask = minimizers_.size() - 1;

    for (std::uint64_t i = 0; i < storage.size(); ++i) {
        for (const auto& it: storage[i]) {
            minimizers_[it.first & mask].emplace_back(it);
        }
        std::vector<uint128_t>().swap(storage[i]);
    }
    std::vector<std::vector<uint128_t>>().swap(storage);

    for (std::uint32_t i = 0; i < minimizers_.size(); ++i) {
        if (minimizers_[i].empty()) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t i) -> void {
                auto& minimizers = minimizers_[i];
                auto& index = index_[i];

                radix_sort(minimizers.begin(), minimizers.end(), k_ * 2,
                    uint128_t_first);

                for (std::uint64_t j = 0, count = 0; j < minimizers.size(); ++j) {
                    if (j > 0 && minimizers[j - 1].first != minimizers[j].first) {
                        index[minimizers[j - 1].first] = std::make_pair(j - count, count);
                        count = 0;
                    }
                    if (j == minimizers.size() - 1) {
                        index[minimizers[j].first] = std::make_pair(j - count, count + 1);
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

    if (f == 0) {
        o_ = -1;
        return;
    }

    std::vector<std::uint64_t> occurrences;
    for (std::uint64_t i = 0; i < index_.size(); ++i) {
        for (const auto& it: index_[i]) {
            occurrences.emplace_back(it.second.second);
        }
    }
    std::nth_element(occurrences.begin(), occurrences.begin() +
        (1 - f) * occurrences.size(), occurrences.end());
    o_ = occurrences[(1 - f) * occurrences.size()] + 1;
}

std::vector<Overlap> MinimizerEngine::map(const std::unique_ptr<Sequence>& sequence,
    bool d, bool t, std::uint32_t e) const {

    auto minimizers = minimize(sequence, e);

    std::vector<Overlap> dst;
    if (minimizers.empty()) {
        return dst;
    }

    std::uint64_t mask = minimizers_.size() - 1;
    std::vector<uint128_t> matches;

    for (const auto& it: minimizers) {

        std::uint32_t bin = it.first & mask;
        auto info = index_[bin].find(it.first);
        if (info == index_[bin].end()) {
            continue;
        }
        if (info->second.second > o_) {
            continue;
        }

        auto jt = minimizers_[bin].begin() + info->second.first;
        auto end = jt + info->second.second;

        for (; jt != end; ++jt) {

            std::uint64_t t_id = jt->second >> 32;

            if (d && sequence->id == t_id) {
                continue;
            }
            if (t && sequence->id > t_id) {
                continue;
            }

            std::uint64_t strand = (it.second & 1) == (jt->second & 1);
            std::uint64_t q_pos = it.second << 32 >> 33;
            std::uint64_t t_pos = jt->second << 32 >> 33;

            std::uint64_t diff = !strand ? t_pos + q_pos :
                (1ULL << 31) + t_pos - q_pos;

            matches.emplace_back((((t_id << 1) | strand) << 32) | diff,
                (q_pos << 32) | t_pos);
        }
    }

    radix_sort(matches.begin(), matches.end(), 64, uint128_t_first);
    matches.emplace_back(-1, -1); // stop dummy

    for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {
        if (matches[i].first - matches[i - 1].first > 500) {

            if (i - j < 4) {
                j = i;
                continue;
            }

            radix_sort(matches.begin() + j, matches.begin() + i, 64,
                uint128_t_second);

            std::vector<std::uint64_t> indices;
            if (matches[i - 1].first >> 32 & 1) {
                indices = longest_subsequence(matches.begin() + j,
                    matches.begin() + i, std::less<std::uint64_t>());
            } else {
                indices = longest_subsequence(matches.begin() + j,
                    matches.begin() + i, std::greater<std::uint64_t>());
            }

            if (indices.size() < 4) {
                j = i;
                continue;
            }

            bool strand = matches[i - 1].first >> 32 & 1;
            std::uint32_t q_matches = 0, q_begin = 0, q_end = 0;
            std::uint32_t t_matches = 0, t_begin = 0, t_end = 0;
            for (std::uint32_t k = 0; k < indices.size(); ++k) {
                std::uint32_t q_pos = matches[j + indices[k]].second >> 32;
                if (q_pos > q_end) {
                    q_matches += q_end - q_begin;
                    q_begin = q_pos;
                }
                q_end = q_pos + k_;

                std::uint32_t t_pos = matches[j + indices[k]].second << 32 >> 32;
                t_pos = strand ? t_pos : (1U << 31) - (t_pos + k_ - 1);
                if (t_pos > t_end) {
                    t_matches += t_end - t_begin;
                    t_begin = t_pos;
                }
                t_end = t_pos + k_;
            }
            q_matches += q_end - q_begin;
            t_matches += t_end - t_begin;
            if (q_matches < 100 || t_matches < 100) {
                j = i;
                continue;
            }

            dst.emplace_back(sequence->id,
                matches[j + indices.front()].second >> 32,
                k_ + (matches[j + indices.back()].second >> 32),
                matches[i - 1].first >> 33,
                strand ? matches[j + indices.front()].second << 32 >> 32 :
                    matches[j + indices.back()].second << 32 >> 32,
                k_ + (strand ? matches[j + indices.back()].second << 32 >> 32 :
                    matches[j + indices.front()].second << 32 >> 32),
                strand, std::min(q_matches, t_matches));

            j = i;
        }
    }

    return dst;
}

std::vector<uint128_t> MinimizerEngine::minimize(
    const std::unique_ptr<Sequence>& sequence, std::uint32_t e) const {

    std::vector<uint128_t> dst;
    if (sequence->data.size() < k_) {
        return dst;
    }

    e = std::min(e, static_cast<std::uint32_t>(sequence->data.size()));
    std::uint32_t b = sequence->data.size() - e;
    std::swap(b, e);

    std::uint64_t mask = (1ULL << (k_ * 2)) - 1;
    std::uint64_t shift = (k_ - 1) * 2;
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

    std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
    std::uint64_t is_stored = 1ULL << 63;

    for (std::uint32_t i = 0; i < sequence->data.size(); ++i) {
        std::uint64_t c = kCoder[sequence->data[i]];
        if (c == 255) {
            throw std::invalid_argument("[ram::MinimizerEngine::minimize] error: "
                "invalid character!");
        }
        minimizer = ((minimizer << 2) | c) & mask;
        reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
        if (i >= k_ - 1U) {
            if (minimizer < reverse_minimizer) {
                window_add(std::hash<std::uint64_t>{}(minimizer), (i - (k_ - 1U)) << 1 | 0);
            } else if (minimizer > reverse_minimizer) {
                window_add(std::hash<std::uint64_t>{}(reverse_minimizer), (i - (k_ - 1U)) << 1 | 1);
            }
        }
        if (i >= (k_ - 1U) + (w_ - 1U)) {
            if (i < b || i - (k_ - 1U) >= e) {
                for (auto it = window.begin(); it != window.end(); ++it) {
                    if (it->first != window.front().first) {
                        break;
                    }
                    if (it->second & is_stored) {
                        continue;
                    }
                    dst.emplace_back(it->first, id | it->second);
                    it->second |= is_stored;
                }
            }
            window_update(i - (k_ - 1U) - (w_ - 1U) + 1);
        }
    }

    return dst;
}

}
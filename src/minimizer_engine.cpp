/*!
 * @file minimizer_engine.cpp
 *
 * @brief MinimizerEngine class source file
 */

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <queue>

#include "thread_pool/thread_pool.hpp"

#include "ram/minimizer_engine.hpp"

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

Sequence::Sequence(const char* name, std::uint32_t name_length,
    const char* data, std::uint32_t data_length)
        : id(sequence_id++), name(name, name_length), data(data, data_length) {
}

Sequence::Sequence(const char* name, std::uint32_t name_length,
    const char* data, std::uint32_t data_length,
    const char*, std::uint32_t)
        : Sequence(name, name_length, data, data_length) {
}

Overlap::Overlap(std::uint64_t t_id, std::uint32_t t_begin, std::uint32_t t_end,
    std::uint32_t q_begin, std::uint32_t q_end, std::uint32_t strand)
        : t_id(t_id), t_begin(t_begin), t_end(t_end), q_begin(q_begin),
        q_end(q_end), strand(strand) {
}

std::uint64_t Sequence::sequence_id = 0;

MinimizerEngine::MinimizerHash::MinimizerHash(std::uint16_t num_bins)
        : m_(num_bins - 1), minimizers(num_bins), index(num_bins) {
}

void MinimizerEngine::MinimizerHash::clear() {
    for (auto& it: minimizers) {
        it.clear();
    }
    for (auto& it: index) {
        it.clear();
    }
}

inline std::pair<std::vector<uint128_t>::const_iterator, std::vector<uint128_t>::const_iterator>
    MinimizerEngine::MinimizerHash::operator[](std::uint64_t minimizer) const {

    std::uint32_t bin = minimizer & m_;
    auto it = index[bin].find(minimizer);
    if (it == index[bin].end()) {
        return std::make_pair(minimizers.front().begin(), minimizers.front().begin());
    }
    return std::make_pair(minimizers[bin].begin() + it->second.first,
        minimizers[bin].begin() + it->second.second);
}

MinimizerEngine::MinimizerEngine(std::uint8_t k, std::uint8_t w,
    std::uint32_t num_threads)
        : k_(k), w_(w), s_(-1), hash_(1U << std::min(14, 2 * k)),
        thread_pool_(thread_pool::createThreadPool(num_threads)) {
}

void MinimizerEngine::minimize(const std::vector<std::unique_ptr<Sequence>>& src) {
    return minimize(src.begin(), src.end());
}

void MinimizerEngine::minimize(
    std::vector<std::unique_ptr<Sequence>>::const_iterator begin,
    std::vector<std::unique_ptr<Sequence>>::const_iterator end) {

    hash_.clear();

    if (begin >= end) {
        return;
    }

    std::vector<std::vector<uint128_t>> minimizers(end - begin);
    std::vector<std::future<void>> thread_futures;
    for (auto it = begin; it != end; ++it) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint64_t i) -> void {
                minimizers[i] = minimize(*(begin + i));
            }
        , it - begin));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    for (std::uint64_t i = 0; i < minimizers.size(); ++i) {
        for (const auto& it: minimizers[i]) {
            hash_.minimizers[it.first & hash_.m_].emplace_back(it);
        }
        std::vector<uint128_t>().swap(minimizers[i]);
    }

    for (std::uint32_t i = 0; i < hash_.minimizers.size(); ++i) {
        if (hash_.minimizers[i].empty()) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t i) -> void {
                auto& minimizers = hash_.minimizers[i];
                auto& index = hash_.index[i];

                radix_sort(minimizers.begin(), minimizers.end(), k_ * 2,
                    uint128_t_first);

                for (std::uint64_t j = 0, count = 0; j < minimizers.size(); ++j) {
                    if (j > 0 && minimizers[j - 1].first != minimizers[j].first) {
                        index[minimizers[j - 1].first] = std::make_pair(j - count, j);
                        count = 0;
                    }
                    if (j == minimizers.size() - 1) {
                        index[minimizers[j].first] = std::make_pair(j - count, j + 1);
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

void MinimizerEngine::filter(double f) {

    if (f == 0) {
        return;
    }

    std::vector<std::uint64_t> frequencies;
    for (std::uint64_t i = 0; i < hash_.index.size(); ++i) {
        for (const auto& it: hash_.index[i]) {
            frequencies.emplace_back(it.second.second - it.second.first);
        }
    }
    std::nth_element(frequencies.begin(), frequencies.begin() +
        (1 - f) * frequencies.size(), frequencies.end());
    s_ = frequencies[(1 - f) * frequencies.size()];
}

std::vector<Overlap> MinimizerEngine::map(const std::unique_ptr<Sequence>& src,
    bool diagonal, bool triangle, std::uint32_t e) const {

    auto minimizers = minimize(src, e);
    radix_sort(minimizers.begin(), minimizers.end(), k_ * 2, uint128_t_first);

    std::vector<uint128_t> matches;
    for (const auto& it: minimizers) {
        auto m = hash_[it.first];
        if (m.second - m.first >= s_) {
            continue;
        }
        for (auto jt = m.first; jt != m.second; ++jt) {

            std::uint64_t query_id = it.second >> 32;
            std::uint64_t target_id = jt->second >> 32;

            if (diagonal && query_id == target_id) {
                continue;
            }
            if (triangle && query_id < target_id) {
                continue;
            }

            std::uint64_t strand = (it.second & 1) == (jt->second & 1);
            std::uint64_t query_pos = it.second << 32 >> 33;
            std::uint64_t target_pos = jt->second << 32 >> 33;

            std::uint64_t diagonal_diff = !strand ? target_pos + query_pos :
                (1ULL << 31) + target_pos - query_pos;

            matches.emplace_back(
                (((target_id << 1) | strand) << 32) | diagonal_diff,
                (target_pos << 32) | query_pos);
        }
    }

    radix_sort(matches.begin(), matches.end(), 64, uint128_t_first);
    matches.emplace_back(-1, -1); // stop dummy

    std::vector<Overlap> dst;

    for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {
        if ((matches[i].first >> 32) != (matches[i - 1].first >> 32) ||
            (matches[i].first << 32 >> 32) - (matches[i - 1].first << 32 >> 32) > 500) {

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
            dst.emplace_back(matches[i - 1].first >> 33,
                matches[j + indices.front()].second >> 32,
                k_ + (matches[j + indices.back()].second >> 32),
                strand ? matches[j + indices.front()].second << 32 >> 32 :
                    matches[j + indices.back()].second << 32 >> 32,
                k_ + (strand ? matches[j + indices.back()].second << 32 >> 32 :
                    matches[j + indices.front()].second << 32 >> 32),
                strand);

            j = i;
        }
    }

    return dst;
}

std::vector<uint128_t> MinimizerEngine::minimize(const std::unique_ptr<Sequence>& src,
    std::uint32_t e) const {

    std::vector<uint128_t> dst;
    if (src->data.size() < k_) {
        return dst;
    }

    e = std::min(e, static_cast<std::uint32_t>(src->data.size()));
    std::uint32_t b = src->data.size() - e;
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

    std::uint64_t id = static_cast<std::uint64_t>(src->id) << 32;

    for (std::uint32_t i = 0; i < src->data.size(); ++i) {
        std::uint64_t c = coder[src->data[i]];
        if (c == 255) {
            throw std::invalid_argument("[ram::MinimizerEngine::minimize] error: "
                "invalid character!");
        }
        minimizer = ((minimizer << 2) | c) & mask;
        reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
        if (i >= k_ - 1U) {
            if (minimizer < reverse_minimizer) {
                window_add(minimizer, id | ((i - (k_ - 1U)) << 1 | 0));
            } else if (minimizer > reverse_minimizer) {
                window_add(reverse_minimizer, id | ((i - (k_ - 1U)) << 1 | 1));
            }
        }
        if (i >= (k_ - 1U) + (w_ - 1U)) {
            if ((dst.empty() || dst.back().second != window.front().second) &&
                (i < b || i - (k_ - 1U) >= e)) {
                dst.emplace_back(window.front());
            }
            window_update(i - (k_ - 1U) - (w_ - 1U) + 1);
        }
    }

    return dst;
}

}

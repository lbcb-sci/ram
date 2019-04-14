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
        while (!window.empty() && (static_cast<std::uint32_t>(window.front().second) >> 1) < position) {
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


std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t> map(
    const std::vector<std::pair<std::uint64_t, std::uint32_t>>& lhs,
    const std::vector<std::pair<std::uint64_t, std::uint32_t>>& rhs) {

    std::vector<std::pair<std::uint32_t, std::uint32_t>> matches;

    for (std::uint32_t i = 0, j = 0; i < lhs.size() && j < rhs.size();) {
        if (lhs[i].first < rhs[j].first) {
            ++i;
        } else if (lhs[i].first == rhs[j].first) {
            for (std::uint32_t k = j; k < rhs.size(); ++k) {
                matches.emplace_back(lhs[i].second, rhs[j].second);
                if (lhs[i].first != rhs[k].first) {
                    break;
                }
            }
            ++i;
        } else {
            ++j;
        }
    }

    //for (const auto& it: matches) {
        //std::cerr << (it.first >> 1) << " " << (it.second >> 1) << std::endl;
    //}

    return std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>();
}


}

std::uint64_t base_to_code(char base) {
    switch (base) {
        case 'C': return 0;
        case 'A': return 1;
        case 'T': return 2;
        case 'G': return 3;
        default: return 4;
    }
}

std:: uint64_t base_to_code_complement(char base) {
    switch (base) {
        case 'C': return 3;
        case 'A': return 2;
        case 'T': return 1;
        case 'G': return 0;
        default: return 4;
    }
}

std::uint64_t get_next_minimzer(const char *sequence, unsigned int position,
    unsigned int kmer_length, unsigned int current_minimizer) {

    unsigned int value = 0;
    unsigned int mask = (1 << (2 * kmer_length)) - 1;
    value = ((current_minimizer << 2) & mask) | base_to_code(sequence[position + kmer_length - 1]);
    return value;
}

std::uint64_t get_next_minimzer_reversed(const char *sequence, unsigned int position,
    unsigned int kmer_length, unsigned int current_minimizer) {

    unsigned int value = 0;
    value = ((base_to_code_complement(sequence[position + kmer_length - 1]) << 2 * kmer_length) | current_minimizer) >> 2;
    return value;
}


void get_minimizers(const char *sequence, std::uint64_t sequence_length,
    std::uint64_t kmer_length, std::uint64_t window_length,
    std::set<Minimizer> &minimizers_set) {

    std::deque<Minimizer> queue;

    std::uint64_t forward_minimizer_value = 0;
    std::uint64_t reversed_minimizer_value = 0;

    for (std::uint64_t i = 0; i < kmer_length; i++) {
        forward_minimizer_value = forward_minimizer_value << 2;
        forward_minimizer_value = forward_minimizer_value | base_to_code(sequence[i]);
    }

    for (std::uint64_t i = kmer_length; i > 0; i--) {
        reversed_minimizer_value = reversed_minimizer_value << 2;
        reversed_minimizer_value = reversed_minimizer_value | base_to_code_complement(sequence[i-1]);
    }

    std::uint64_t current_minimizer_value = forward_minimizer_value < reversed_minimizer_value ? forward_minimizer_value : reversed_minimizer_value;
    int current_minimizer_strand = forward_minimizer_value < reversed_minimizer_value ? 0 : 1;
    queue.emplace_back(current_minimizer_value, 0, current_minimizer_strand);

    minimizers_set.insert(queue.front());

    for (std::uint64_t i = 1; i < window_length - 1; i++) {
        forward_minimizer_value = get_next_minimzer(sequence, i, kmer_length, forward_minimizer_value);
        reversed_minimizer_value = get_next_minimzer_reversed(sequence, i, kmer_length, reversed_minimizer_value);

        std::uint64_t current_minimizer_value = forward_minimizer_value < reversed_minimizer_value ? forward_minimizer_value : reversed_minimizer_value;
        int current_minimizer_strand = forward_minimizer_value < reversed_minimizer_value ? 0 : 1;
        queue.emplace_back(current_minimizer_value, 0, current_minimizer_strand);
    }

    std::uint64_t upper_bound = sequence_length - (kmer_length-1);

    for (std::uint64_t beginning_position = window_length - 1; beginning_position <= upper_bound; beginning_position++) {

        forward_minimizer_value = get_next_minimzer(sequence, beginning_position, kmer_length, forward_minimizer_value);
        reversed_minimizer_value = get_next_minimzer_reversed(sequence, beginning_position, kmer_length, reversed_minimizer_value);
        uint64_t current_minimizer_value = forward_minimizer_value < reversed_minimizer_value ? forward_minimizer_value : reversed_minimizer_value;
        int current_minimizer_strand = forward_minimizer_value < reversed_minimizer_value ? 0 : 1;

        while (!queue.empty()) {
            Minimizer current_minimizer = queue.back();
            if (current_minimizer.value >= current_minimizer_value) {
                queue.pop_back();
            } else {
                break;
            }
        }

        queue.emplace_back(current_minimizer_value, beginning_position, current_minimizer_strand);
        Minimizer front_element = queue.front();
        minimizers_set.insert(front_element);

        if (front_element.location + (window_length - 1) <= beginning_position) {
            queue.pop_front();
        }
    }

    if(queue.size() > 0) {
        minimizers_set.insert(queue.back());
    }
}

std::vector<Minimizer> createMinimizers(const char *sequence, std::uint64_t sequence_length,
    std::uint64_t kmer_length, std::uint64_t window_length) {

    std::vector<Minimizer> minimizers_vector;
    std::set<Minimizer> minimizers_set;

    get_minimizers(sequence, sequence_length, kmer_length, window_length, minimizers_set);

    for (auto minimizer: minimizers_set) {
        minimizers_vector.emplace_back(minimizer);
    }

    return minimizers_vector;
}

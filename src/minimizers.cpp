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

inline std::uint8_t coder(char c) {
    switch (c) {
        case 'A':
        case 'a': return 1;
        case 'C':
        case 'c': return 0;
        case 'G':
        case 'g': return 3;
        case 'T':
        case 't': return 2;
        default: throw std::invalid_argument("[ram::coder] error: "
            "invalid character!");
    }
}

inline std::uint8_t reverseCoder(char c) {
    switch (c) {
        case 'A':
        case 'a': return 2;
        case 'C':
        case 'c': return 3;
        case 'G':
        case 'g': return 0;
        case 'T':
        case 't': return 1;
        default: throw std::invalid_argument("[ram::reverseCoder] error: "
            "invalid character!");
    }
}

std::vector<std::pair<std::uint64_t, std::uint32_t>> createMinimizers(
    const char* sequence, std::uint32_t sequence_length, std::uint32_t k,
    std::uint32_t w) {

    if (k > 31) {
        throw std::invalid_argument("[ram::createMinimizers] error: "
            "invalid kmer size!");
    }

    std::vector<std::pair<std::uint64_t, std::uint32_t>> minimizers;
    if (sequence_length < k) {
        return minimizers;
    }

    uint64_t mask = (1 << (k * 2)) - 1;
    uint64_t minimizer = 0, reverse_minimizer = 0;

    std::deque<std::pair<std::uint64_t, std::uint32_t>> window;
    auto window_add = [&window](std::uint64_t value, std::uint32_t location) -> void {
        while (!window.empty() && window.back().first > value) {
            window.pop_back();
        }
        window.emplace_back(value, location);
    };
    auto window_update = [&window](std::uint32_t position) -> void {
        while (!window.empty() && (window.front().second >> 1) < position) {
            window.pop_front();
        }
    };

    for (std::uint32_t i = 0; i < sequence_length - k + 1; ++i) {
        minimizer = (minimizer << 2 | coder(sequence[i])) & mask;
        reverse_minimizer = (reverse_minimizer << 2 |
            reverseCoder(sequence[sequence_length - i - 1])) & mask;
        if (i >= (k - 1) + w) {
            if (minimizers.empty() ||
                minimizers.back().second != window.front().second) {
                minimizers.emplace_back(window.front());
            }
            window_update(i - (k - 1) - (w - 1));
        }
        if (i >= k - 1) {
            if (minimizer < reverse_minimizer) {
                window_add(minimizer, (i - (k - 1)) << 1 | 0);
            } else if (minimizer > reverse_minimizer) {
                window_add(reverse_minimizer, (i - (k - 1)) << 1 | 1);
            }
        }
    }

    std::sort(minimizers.begin(), minimizers.end(), [](
        const std::pair<std::uint64_t, std::uint32_t>& lhs,
        const std::pair<std::uint64_t, std::uint32_t>& rhs) {
            if (lhs.first < rhs.first) return true;
            if (lhs.first == rhs.second) return lhs.second < rhs.second;
            return false;
        });

    return minimizers;
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

    for (const auto& it: matches) {
        std::cerr << (it.first >> 1) << " " << (it.second >> 1) << std::endl;
    }

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

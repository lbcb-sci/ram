#include <queue>
#include <set>

#include "minimizers.hpp"

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

std::vector<Minimizer> create_minimizers(const char *sequence, std::uint64_t sequence_length,
    std::uint64_t kmer_length, std::uint64_t window_length) {

    std::vector<Minimizer> minimizers_vector;
    std::set<Minimizer> minimizers_set;

    get_minimizers(sequence, sequence_length, kmer_length, window_length, minimizers_set);

    for (auto minimizer: minimizers_set) {
      minimizers_vector.emplace_back(minimizer);
    }

    return minimizers_vector;

}

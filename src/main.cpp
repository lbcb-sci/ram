#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <unordered_map>

#include "minimizers.hpp"
#include "aligner.hpp"

#include "bioparser/bioparser.hpp"
#include "logger/logger.hpp"

static const std::string version = "v0.0.1";

static struct option options[] = {
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

void help();

struct Sequence {
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length)
            : name(name, name_length), data(data, data_length) {
    }

    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length,
        const char*, std::uint32_t)
            : Sequence(name, name_length, data, data_length) {
    }

    ~Sequence() {
    }

    std::string name;
    std::string data;
};

uint64_t get_complement_letter_value(char letter) {
    switch (letter) {
        case 'C': return 'G';
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        default: return  0;
    }
}

std::uint64_t value_for_minimizer_location(std::uint64_t minimizer_location) {
    return static_cast<std::uint64_t>((static_cast<std::uint32_t>(minimizer_location) >> 1));
}

std::uint32_t sequence_for_minimizer_location(std::uint64_t minimizer_location) {
    return minimizer_location >> 32;
}

int get_ceil_index(std::vector<std::pair<std::uint64_t, std::uint64_t>> matches,
    std::vector<std::uint64_t> &tail_indices,
    std::uint64_t left, std::uint64_t right, std::uint64_t key) {
    while (right - left > 1) {
        int middle = left + (right - left)/2;
        if (value_for_minimizer_location(matches[tail_indices[middle]].second) >= key) {
            right = middle;
        } else {
            left = middle;
        }
    }
    return right;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> find_lis_matches_reversed(std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches) {

    std::vector<std::uint64_t> tail_indices(matches.size(), 0);
    std::vector<std::uint64_t> prev_indices(matches.size(), -1);

    std::cout << "in lis" << std::endl;

    std::uint64_t len = 1;
    for (std::uint64_t i = 1; i < matches.size(); i++) {
        if (value_for_minimizer_location(matches[i].second) < value_for_minimizer_location(matches[tail_indices[0]].second)) {
            tail_indices[0] = i;
        } else if (value_for_minimizer_location(matches[i].second) > value_for_minimizer_location(matches[tail_indices[len-1]].second)) {
            prev_indices[i] = tail_indices[len-1];
            tail_indices[len++] = i;
        } else {
            std::uint64_t pos = get_ceil_index(matches, tail_indices, -1, len-1, value_for_minimizer_location(matches[i].second));
            prev_indices[i] = tail_indices[pos-1];
            tail_indices[pos] = i;
        }
    }

    std::cout << "after for" << std::endl;

    std::vector<std::pair<std::uint64_t, std::uint64_t>> best_matches;

    for (int i = tail_indices[len-1]; i >= 0; i = prev_indices[i]) {
        best_matches.push_back(matches[i]);
    }

    std::cout << "before return" << std::endl;

    return best_matches;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> find_lis_matches(std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches) {

    std::vector<std::uint64_t> tail_indices(matches.size(), 0);
    std::vector<std::uint64_t> prev_indices(matches.size(), -1);

    std::uint64_t len = 1;
    for (std::uint64_t i = 1; i < matches.size(); i++) {
        if (value_for_minimizer_location(matches[i].second) > value_for_minimizer_location(matches[tail_indices[0]].second)) {
            tail_indices[0] = i;
        } else if (value_for_minimizer_location(matches[i].second) < value_for_minimizer_location(matches[tail_indices[len-1]].second)) {
            prev_indices[i] = tail_indices[len-1];
            tail_indices[len++] = i;
        } else {
            std::uint64_t pos = get_ceil_index(matches, tail_indices, -1, len-1, value_for_minimizer_location(matches[i].second));
            prev_indices[i] = tail_indices[pos-1];
            tail_indices[pos] = i;
        }
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> best_matches;

    for (int i = tail_indices[len-1]; i >= 0; i = prev_indices[i]) {
        best_matches.push_back(matches[i]);
    }
    return best_matches;
}

bool sortbysequence(const std::pair<std::uint64_t, std::uint64_t> &a,
              const std::pair<std::uint64_t, std::uint64_t> &b) {
                  return a.second > b.second;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> find_best_matches(std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches) {
    std::sort(matches.begin(), matches.end(), sortbysequence);

    std::vector<std::pair<std::uint64_t, std::uint64_t>> sequence_matches;

    if (matches.size() <= 0) {
        return matches;
    }

    std::uint32_t best_sequence_index = 0;
    std::uint32_t best_sequence_size = 1;
    std::uint64_t previous_sequence = matches[0].second >> 32;

    std::uint32_t current_sequence_size = 1;
    std::uint32_t current_sequence_index = 0;

    for (std::uint32_t i = 1; i < matches.size(); i++) {
        auto match = matches[i];
        std::uint64_t current_sequence = match.second >> 32;

        // std::cout << current_sequence << " tt " << (match.first >> 32) << std::endl;
        if (current_sequence != previous_sequence) {
            if (current_sequence_size > best_sequence_size && previous_sequence != match.first >> 32) {
                best_sequence_size = current_sequence_size;
                best_sequence_index = current_sequence_index;
            }
            current_sequence_size = 1;
            current_sequence_index = i;
            previous_sequence = current_sequence;
        } else {
            current_sequence_size += 1;
        }
    }

    if (current_sequence_size > best_sequence_size && previous_sequence != matches[matches.size()-1].first >> 32) {
        best_sequence_size = current_sequence_size;
        best_sequence_index = current_sequence_index;
    }

    // std::cout << "best_sequence_index " << best_sequence_index << std::endl;
    // std::cout << "best_sequence_size " << best_sequence_size << std::endl;

    for (std::uint32_t i = best_sequence_index; i < (best_sequence_index + best_sequence_size); i++) {
        sequence_matches.emplace_back(matches[i]);
    }

    std::cout << "sequence_matches " << sequence_matches.size() << std::endl;

    std::vector<std::pair<std::uint64_t, std::uint64_t>> forward_matches;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> reversed_matches;

    for (auto const &match: sequence_matches) {
        if ((match.first & 1) == (match.second & 1)) {
            forward_matches.push_back(match);
        } else {
            reversed_matches.push_back(match);
        }
    }

    std::cout << "forward_matches " << forward_matches.size() << std::endl;
    std::cout << "reversed_matches " << reversed_matches.size() << std::endl;

    std::reverse(reversed_matches.begin(), reversed_matches.end());
    std::vector<std::pair<std::uint64_t, std::uint64_t>> best_matches_reversed;
    if (reversed_matches.size() > 0) {
        best_matches_reversed = find_lis_matches_reversed(reversed_matches);
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> best_matches_forward;
    if (forward_matches.size() > 0) {
        best_matches_forward = find_lis_matches(forward_matches);
    }

    std::cout << "Best forward matches " << best_matches_forward.size() << std::endl;
    std::cout << "Best reversed matches " << best_matches_reversed.size() << std::endl;

    if (best_matches_forward.size() > best_matches_reversed.size()) {
        return best_matches_forward;
    } else {
        return best_matches_reversed;
    }
}

std::vector<Cluster> create_clusters(std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches, std::uint64_t kmer_length) {

    std::vector<Cluster> clusters;

    if (matches.size() <= 0) {
        return clusters;
    }

    std::uint64_t start_ref = value_for_minimizer_location(matches[0].first);
    std::uint64_t len_ref = 0;
    std::uint64_t start_read = value_for_minimizer_location(matches[0].second);
    std::uint64_t len_read = 0;

    for (std::uint64_t i = 1; i < matches.size(); i++) {
        auto match = matches[i];
        auto previous_match = matches[i-1];

        // std::cout << "match location " << value_for_minimizer_location(match.first) << " " << value_for_minimizer_location(match.second) << std::endl;
        // std::cout << "match sequence " << sequence_for_minimizer_location(match.first) << " " << sequence_for_minimizer_location(match.second) << std::endl;

        bool strand = (match.first & 1) == (match.second & 1);

        std::uint64_t ref_len = value_for_minimizer_location(match.first) - value_for_minimizer_location(previous_match.first);
        std::uint64_t read_len = strand ?
            value_for_minimizer_location(match.second) - value_for_minimizer_location(previous_match.second):
            value_for_minimizer_location(previous_match.second) - value_for_minimizer_location(match.second);

        if (read_len < 15 && ref_len > 100) {
            std::uint64_t adjusted_start_read = strand ? start_read : start_read-len_read;
            clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, strand);

            start_ref = value_for_minimizer_location(match.first);
            start_read = value_for_minimizer_location(match.second);
            len_ref = 0;
            len_read = 0;
        } else {
            len_ref += ref_len;
            len_read += read_len;
        }
    }

    // std::cout << "len_ref " << len_ref << std::endl;
    if (len_ref > 0) {
        bool strand = (matches[0].first & 1) == (matches[0].second & 1);

        std::uint64_t adjusted_start_read = strand ? start_read-len_read : start_read;
        clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, strand);
    }

  return clusters;
}

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool sortbylength(const std::unique_ptr<Sequence> &a,
              const std::unique_ptr<Sequence> &b) {
                  return a->data.size() < b->data.size();
}

int main(int argc, char** argv) {

    std::vector<std::string> input_paths;

    std::uint64_t kmer_length = 15;
    std::uint64_t window_length = 5;

    char argument;
    while ((argument = getopt_long(argc, argv, "h", options, nullptr)) != -1) {
        switch (argument) {
            case 'v': std::cout << version << std::endl; return 0;
            case 'h': help(); return 0;
            default: return 1;
        }
    }

    for (std::int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.empty()) {
        std::cerr << "[ram::] error: missing input file(s)!" << std::endl;
        help();
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr;

    if (isSuffix(input_paths[0], ".fasta") || isSuffix(input_paths[0], ".fa") ||
        isSuffix(input_paths[0], ".fasta.gz") || isSuffix(input_paths[0], ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            input_paths[0]);
    } else if (isSuffix(input_paths[0], ".fastq") || isSuffix(input_paths[0], ".fq") ||
               isSuffix(input_paths[0], ".fastq.gz") || isSuffix(input_paths[0], ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            input_paths[0]);
    } else {
        std::cerr << "[ram::] error: file " << input_paths[0] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
        std::endl;
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> tparser = nullptr;

    if (input_paths.size() > 1) {
        if (isSuffix(input_paths[1], ".fasta") || isSuffix(input_paths[1], ".fa") ||
            isSuffix(input_paths[1], ".fasta.gz") || isSuffix(input_paths[1], ".fa.gz")) {
            tparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
                input_paths[1]);
        } else if (isSuffix(input_paths[1], ".fastq") || isSuffix(input_paths[1], ".fq") ||
                   isSuffix(input_paths[1], ".fastq.gz") || isSuffix(input_paths[1], ".fq.gz")) {
            tparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
                input_paths[1]);
        } else {
            std::cerr << "[ram::] error: file " << input_paths[1] <<
                " has unsupported format extension (valid extensions: .fasta, "
                ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
            std::endl;
            return 1;
        }
    }

    std::vector<std::unique_ptr<Sequence>> containee_reads;
    sparser->parse(containee_reads, -1);

    std::vector<std::unique_ptr<Sequence>> contained_reads;
    tparser->parse(contained_reads, -1);

    auto logger = logger::Logger();
    logger.log();

    std::unordered_map<uint64_t, std::vector<uint64_t>> minimizer_hash;

    sort(contained_reads.begin(), contained_reads.end(), sortbylength);

    uint32_t id = 0;
    std::vector<std::pair<uint64_t, uint64_t>> minimizers;
    for (const auto& it: contained_reads) {
        ram::createMinimizers(minimizers, it->data.c_str(), it->data.size(), id, 15, 5);
        ++id;
    }

    std::cerr << minimizers.size() << std::endl;

    logger.log("[ram::] collected minimizers in");
    logger.log();

    ram::sortMinimizers(minimizers, 15);

    logger.log("[ram::] sorted minimizers in");
    logger.log();

    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
        if (i != 0 && minimizers[i].first == minimizers[i-1].first &&
            minimizers[i].second < minimizers[i - 1].second) {
            // std::cerr << minimizers[i].first << " " << (minimizers[i].second >> 32) << " " << (static_cast<std::uint32_t>(minimizers[i].second) >> 1) << " " << (minimizers[i].second & 1) << std::endl;
            // std::cerr << minimizers[i - 1].first << " " << (minimizers[i - 1].second >> 32) << " " << (static_cast<std::uint32_t>(minimizers[i - 1].second) >> 1) << " " << (minimizers[i - 1].second & 1) << std::endl;
            //throw std::logic_error("minimizers");
        }
    }

    std::uint32_t num_minimizers = 0;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
        if (i != 0 && minimizers[i - 1].first != minimizers[i].first) {
            ++num_minimizers;
        }
    }

    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> hash;
    hash.reserve(num_minimizers);

    std::vector<std::uint32_t> counts;
    counts.reserve(num_minimizers);

    std::uint32_t count = 0;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
        if (i != 0 && minimizers[i - 1].first != minimizers[i].first) {
            hash[minimizers[i - 1].first] = std::make_pair(i - count, count);
            counts.emplace_back(count);
            count = 0;
        }
        ++count;
    }

    logger.log("[ram::] created hash in");
    logger.log();

    std::sort(counts.begin(), counts.end());

    std::cerr << "Hash size: " << hash.size() << std::endl;
    std::cerr << counts[(1 - 0.001) * counts.size()] << std::endl;

    logger.log("[ram::] found occurences in");
    logger.log();

    auto indices = ram::longestIncreasingSubsequence(minimizers.begin(), minimizers.end());

    logger.log("[ram::] found lis in");

    return 0;

    uint32_t read_shorter_than_2000 = 0;

    id = 0;
    for (const auto& it: contained_reads) {

        if (it->data.size() <= 2000) {
            read_shorter_than_2000 += 1;
            id += 1;
            continue;
        }

        std::cout << "----------------" << std::endl;
        std::vector<std::pair<uint64_t, uint64_t>> read_minimizers;

        std::uint32_t end_read_size = 1000;

        ram::createMinimizers(read_minimizers, it->data.c_str(), end_read_size, id, 15, 5);
        ram::createMinimizers(read_minimizers, it->data.c_str() + (it->data.size() - (end_read_size+1)), end_read_size, id, 15, 5);
        ram::sortMinimizers(read_minimizers, 15);

        auto matches = ram::map(read_minimizers, minimizers, hash);
        std::cout << "matches size " << matches.size() << std::endl;

        std::vector<std::pair<std::uint64_t, std::uint64_t>> best_matches_start = find_best_matches(matches);
        std::cout << "best matches size " << best_matches_start.size() << std::endl;

        std::vector<Cluster> clusters = create_clusters(best_matches_start, kmer_length);

        if(clusters.size() > 0) {
            std::cout << "clusters size " << clusters.size() << std::endl;
            for (auto const &cluster: clusters) {
                std::cout << cluster.to_string() << std::endl;
            }

            std::cout << "size of sequence " << contained_reads[best_matches_start[0].second >> 32]->data.size() << std::endl;
            std::cout << "read size " << it->data.size() << std::endl;

            int back_clipping = it->data.size() - (clusters[0].ref_start + clusters[0].ref_length);
            bool isStartContained = clusters[0].read_start > clusters[0].ref_start;
            bool isEndContained = (clusters[0].ref_start + clusters[0].ref_length + back_clipping) < contained_reads[best_matches_start[0].second >> 32]->data.size();

            std::cout << isStartContained << std::endl;
            std::cout << isEndContained << std::endl;

            if ((double)clusters[0].read_length / (double)it->data.size() > 0.5 && isStartContained && isEndContained) {
                std::cout << "Read " << id << " contained in " << (best_matches_start[0].second >> 32) << std::endl;
            }
        }

        id += 1;
    }

    std::cout << "reads shorter than 2000 " << read_shorter_than_2000 << std::endl;

    return 0;
}

void help() {
    std::cout <<
        "usage: ram [options ...] <sequences> [<target sequences>]\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "    <target sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing target sequences\n"
        "\n"
        "    options:\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}

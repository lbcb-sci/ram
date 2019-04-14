#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>

#include "minimizers.hpp"
#include "aligner.hpp"

#include "bioparser/bioparser.hpp"

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

int get_ceil_index(std::vector<Match> matches, std::vector<std::uint64_t> &tail_indices,
    std::uint64_t left, std::uint64_t right, std::uint64_t key) {
    while (right - left > 1) {
        int middle = left + (right - left)/2;
        if (matches[tail_indices[middle]].read_start >= key) {
            right = middle;
        } else {
            left = middle;
        }
    }
    return right;
}

std::vector<Match> find_lis_matches(std::vector<Match> &matches) {

    std::vector<std::uint64_t> tail_indices(matches.size(), 0);
    std::vector<std::uint64_t> prev_indices(matches.size(), -1);

    std::uint64_t len = 1;
    for (std::uint64_t i = 1; i < matches.size(); i++) {
        if (matches[i].read_start < matches[tail_indices[0]].read_start) {
            tail_indices[0] = i;
        } else if (matches[i].read_start > matches[tail_indices[len-1]].read_start) {
            prev_indices[i] = tail_indices[len-1];
            tail_indices[len++] = i;
        } else {
            std::uint64_t pos = get_ceil_index(matches, tail_indices, -1, len-1, matches[i].read_start);
            prev_indices[i] = tail_indices[pos-1];
            tail_indices[pos] = i;
        }
    }

    std::vector<Match> best_matches;

    for (int i = tail_indices[len-1]; i >= 0; i = prev_indices[i]) {
        best_matches.push_back(matches[i]);
    }
    return best_matches;
}

std::vector<Match> find_best_matches(std::vector<Match> &matches) {
    std::sort(matches.begin(), matches.end());

    std::vector<Match> forward_matches;
    std::vector<Match> reversed_matches;

    for (auto const &match: matches) {
        if (match.strand) {
            reversed_matches.push_back(match);
        } else {
            forward_matches.push_back(match);
        }
    }

    std::reverse(reversed_matches.begin(), reversed_matches.end());
    std::vector<Match> best_matches_reversed;
    if (reversed_matches.size() > 0) {
        best_matches_reversed = find_lis_matches(reversed_matches);
    }

    std::vector<Match> best_matches_forward;
    if (forward_matches.size() > 0) {
        best_matches_forward = find_lis_matches(forward_matches);
    }

    if (best_matches_forward.size() > best_matches_reversed.size()) {
        return best_matches_forward;
    } else {
        return best_matches_reversed;
    }
}

std::vector<Match> get_matches(std::vector<Minimizer> &first_sequence, std::vector<Minimizer> &second_sequence) {
    std::vector<Match> matches;
    std::uint64_t first_minimizers_index = 0;
    std::uint64_t second_minimizers_index = 0;

    while (first_minimizers_index < first_sequence.size() && second_minimizers_index < second_sequence.size()) {
        Minimizer first_minimizer = first_sequence[first_minimizers_index];
        Minimizer second_minimizer = second_sequence[second_minimizers_index];

        if (first_minimizer.value == second_minimizer.value) {
            std::uint64_t first_minimizers_inner_index = first_minimizers_index;
            std::uint64_t second_minimizers_inner_index = second_minimizers_index;

            while (first_minimizers_inner_index < first_sequence.size()) {
                Minimizer first_minimizer_inner = first_sequence[first_minimizers_inner_index];

                if (first_minimizer.value != first_minimizer_inner.value) {
                    break;
                }

                second_minimizers_inner_index = second_minimizers_index;
                while (second_minimizers_inner_index < second_sequence.size()) {
                    Minimizer second_minimizer_inner = second_sequence[second_minimizers_inner_index];
                    if (second_minimizer_inner.value == first_minimizer_inner.value) {
                        matches.emplace_back(first_minimizer_inner, second_minimizer_inner);
                        second_minimizers_inner_index += 1;
                    } else {
                        break;
                    }
                }
                first_minimizers_inner_index += 1;
            }

            first_minimizers_index = first_minimizers_inner_index;
            second_minimizers_index = second_minimizers_inner_index;
        } else {
            if (first_minimizer.value < second_minimizer.value) {
                first_minimizers_index += 1;
            } else {
                second_minimizers_index += 1;
            }
        }
    }

    return matches;
}

std::vector<Cluster> create_clusters(std::vector<Match> matches, std::uint64_t kmer_length) {

    std::vector<Cluster> clusters;

    if (matches.size() <= 0) {
        return clusters;
    }

    std::uint64_t start_ref = matches[0].ref_start;
    std::uint64_t len_ref = 0;
    std::uint64_t start_read = matches[0].read_start;
    std::uint64_t len_read = 0;

    for (std::uint64_t i = 1; i < matches.size(); i++) {
        Match match = matches[i];
        Match previous_match = matches[i-1];

        std::uint64_t ref_len = match.ref_start - previous_match.ref_start;
        std::uint64_t read_len = match.strand ?
            previous_match.read_start - match.read_start :
            match.read_start - previous_match.read_start;

        if (read_len < 15 && ref_len > 100) {
            std::uint64_t adjusted_start_read = match.strand ? start_read-len_read : start_read;
            clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, match.strand);

            start_ref = match.ref_start;
            start_read = match.read_start;
            len_ref = 0;
            len_read = 0;
        } else {
            len_ref += ref_len;
            len_read += read_len;
        }
    }
    if (len_ref > 0) {
        std::uint64_t adjusted_start_read = matches[0].strand ? start_read-len_read : start_read;
        clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, matches[0].strand);
    }

  return clusters;
}

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
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
    sparser->parse_objects(containee_reads, -1);

    std::string containee = containee_reads.front()->data;

    std::vector<std::unique_ptr<Sequence>> contained_reads;
    tparser->parse_objects(contained_reads, -1);
    std::string& contained = contained_reads.front()->data;

    auto containee_mini = ram::createMinimizers(containee.c_str(), contained.size(), 15, 5);
    auto contained_mini = ram::createMinimizers(contained.c_str(), contained.size(), 15, 5);

    ram::map(containee_mini, contained_mini);

    return 0;

    std::cout << "Containee size: " << containee.length() << std::endl;
    std::cout << "Contained size: " << contained.length() << std::endl;

    std::vector<Minimizer> minimizers_contained;
    minimizers_contained = createMinimizers(contained.c_str(),
        contained.length(), kmer_length, window_length);

    std::cout << "Minimizers contained: " << minimizers_contained.size() << std::endl;

    std::vector<Minimizer> start_minimizers;
    std::vector<Minimizer> end_minimizers;

    for (auto const &minimizer: minimizers_contained) {
      if(minimizer.location < 1000) {
        start_minimizers.push_back(minimizer);
      }
      if(minimizer.location > (contained.length()-1000)) {
        end_minimizers.push_back(minimizer);
      }
    }

    std::cout << "Start minimizers size: " << start_minimizers.size() << std::endl;
    std::cout << "End minimizers size: " << end_minimizers.size() << std::endl;

    std::vector<Minimizer> minimizers_containee;
    minimizers_containee = createMinimizers(containee.c_str(),
        containee.length(), kmer_length, window_length);

    std::vector<Match> matches_start = get_matches(start_minimizers, minimizers_containee);
    std::vector<Match> result_matches_start = find_best_matches(matches_start);
    std::vector<Cluster> clusters_start = create_clusters(result_matches_start, kmer_length);

    for (auto const &match: result_matches_start) {
      std::cout << match.to_string() << std::endl;
    }

    for (auto const &cluster: clusters_start) {
      std::cout << cluster.to_string() << std::endl;
    }

    std::vector<Match> matches_end = get_matches(end_minimizers, minimizers_containee);
    std::vector<Match> result_matches_end = find_best_matches(matches_end);
    std::vector<Cluster> clusters_end = create_clusters(result_matches_end, kmer_length);

    for (auto const &match: result_matches_end) {
      std::cout << match.to_string() << std::endl;
    }

    for (auto const &cluster: clusters_end) {
      std::cout << cluster.to_string() << std::endl;
    }

    // todo alignment
    std::string cigar;
    std::uint64_t target_begin = 0;

    Cluster cluster = clusters_start[0];

    std::string t_sub = containee.substr(cluster.ref_start, cluster.ref_length);
    std::string q_sub = contained.substr(cluster.read_start, cluster.read_length);

    char q_sub_reverse[q_sub.size()];
    for (std::uint64_t i = 0; i < q_sub.size(); i++) {
      q_sub_reverse[(q_sub.size()-1)-i] = get_complement_letter_value(q_sub[i]);
    }

    const char* t_sub_char = t_sub.c_str();

    pairwise_alignment(&q_sub_reverse[0], q_sub.size(), t_sub_char, t_sub.size(),
        2, -1, -2, cigar, target_begin);

    cigar = std::string(cigar.rbegin(), cigar.rend());

    std::cout << cigar << std::endl;

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

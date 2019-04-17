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

void create_clusters(
    std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> &clusters,
    std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches,
    std::vector<std::uint32_t> &matches_indexes) {

    if (matches.size() <= 0) {
        return;
    }

    std::uint32_t start_ref = matches[0].second >> 32;
    std::uint32_t len_ref = 0;
    std::uint32_t start_read = (matches[0].second << 32) >> 32;
    std::uint32_t len_read = 0;

    auto previous_it = matches_indexes.begin();
    auto current_it = previous_it + 1;

    while(current_it != matches_indexes.end()) {
        auto match = matches[*current_it];
        auto previous_match = matches[*previous_it];

        std::cout << (match.second >> 32) << " --- " << ((match.second << 32) >> 32) << std::endl;

        // bool strand = (match.first & 1) == (match.second & 1);

        std::uint64_t ref_len = (match.second >> 32) - (previous_match.second >> 32);
        std::uint64_t read_len = ((match.second << 32) >> 32) -  ((previous_match.second << 32) >> 32);

        if (ref_len > 1000) {
            clusters.emplace_back(start_ref, len_ref+15, start_read, len_read+15);

            start_ref = match.second >> 32;
            start_read = (match.second << 32) >> 32;
            len_ref = 0;
            len_read = 0;
        } else {
            len_ref += ref_len;
            len_read += read_len;
        }

        // std::cout << start_ref << " " << len_ref << " " << start_read << " " << len_read << std::endl;

        current_it += 1;
        previous_it += 1;
    }

    std::cout << "len_ref " << len_ref << std::endl;
    if (len_ref > 0) {
        clusters.emplace_back(start_ref, len_ref+15, start_read, len_read+15);
    }

    // std::cout << "Number of clusters " << clusters.size() << std::endl;

    return;
}

bool are_clusters_contained(std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> &clusters) {

    auto start_it = clusters.begin();

    std::uint32_t min_query_pos = (0-1);
    std::uint32_t min_target_pos = (0-1);

    std::uint32_t max_query_pos = 0;
    std::uint32_t max_target_pos = 0;

    while (start_it != clusters.end()) {
        min_query_pos = std::min(min_query_pos, std::get<2>(*start_it));
        min_target_pos = std::min(min_target_pos, std::get<1>(*start_it));

        max_query_pos = std::max(max_target_pos, std::get<2>(*start_it) + std::get<3>(*start_it));
        max_target_pos = std::max(max_target_pos, std::get<0>(*start_it) + std::get<1>(*start_it));

        start_it += 1;
    }

    std::cout << min_query_pos << " " << min_target_pos << " " << max_query_pos << " " << max_target_pos << std::endl;

    return min_query_pos < min_target_pos && max_query_pos < max_target_pos;
}

bool find_best_matches(std::vector<std::pair<std::uint64_t, std::uint64_t>> &matches) {

    std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t>> clusters;

    if (matches.size() <= 0) {
        return false;
    }

    std::sort(matches.begin(), matches.end());

    std::cout << "all matches" << std::endl;
    for(auto &match : matches) {
        // std::cout << match.first << " " << match.second << std::endl;
        std::cout << (match.first >> 33) << ":  " << (match.first << 32 >> 32) << " " << (match.second >> 32)  << " " << (match.second << 32 >> 32) << " - " << (match.first >> 32 & 1) << std::endl;
    }

    bool found_contained = false;
    auto start_it = matches.begin();
    auto previous_it = start_it;
    auto current_it = previous_it + 1;

    while(current_it != matches.end() && !found_contained) {
        if ((previous_it->first >> 32) != (current_it->first >> 32)) {
            // std::cout << "----" << std::endl;
            // auto tmp_it = start_it;
            // while(tmp_it != current_it) {
            //     std::cout << tmp_it->first << " " << tmp_it->second << std::endl;
            //     tmp_it += 1;
            // }
            // std::cout << "----" << std::endl;
            auto matches_indexes = ram::longestIncreasingSubsequence(start_it, current_it);

            std::cout << "lis matches" << std::endl;
            for(auto &mi : matches_indexes) {
                std::cout << mi << " - " << (matches[mi].first >> 33) << ":  " << (matches[mi].first << 32 >> 32) << "  " << (matches[mi].second >> 32) << "  " << (matches[mi].second << 32 >> 32) << " " << " - " << (matches[mi].first >> 32 & 1) << std::endl;
            }

            create_clusters(clusters, matches, matches_indexes);

            // std::cout << "----xxxx" << std::endl;

            auto are_contained = are_clusters_contained(clusters);
            std::cout << "number of clusters " << clusters.size() << std::endl;
            if (are_contained) {
                return true;
            }

            clusters.empty();

            start_it = current_it;
        } else if ((((current_it->first << 32) >> 32) - ((previous_it->first << 32) >> 32)) > 200) {
            // std::cout << "xxxx" << std::endl;
            // auto tmp_it = start_it;
            // while(tmp_it != current_it) {
            //     std::cout << tmp_it->first << " " << tmp_it->second << std::endl;
            //     tmp_it += 1;
            // }
            // std::cout << "xxxx" << std::endl;
            auto matches_indexes = ram::longestIncreasingSubsequence(start_it, current_it);
            create_clusters(clusters, matches, matches_indexes);

            start_it = current_it;
        }
        previous_it = current_it;
        current_it += 1;
    }

    // std::cout << "----" << std::endl;
    // auto tmp_it = start_it;
    // while(tmp_it != current_it) {
    //     std::cout << tmp_it->first << " " << tmp_it->second << std::endl;
    //     tmp_it += 1;
    // }
    // std::cout << "----" << std::endl;

    // std::cout << "gotovo" << std::endl;

    auto matches_indexes = ram::longestIncreasingSubsequence(start_it, current_it);
    create_clusters(clusters, matches, matches_indexes);

    std::cout << "number of clusters " << clusters.size() << std::endl;

    auto are_contained = are_clusters_contained(clusters);
    return are_contained;
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

        std::sort(read_minimizers.begin(), read_minimizers.end(), [](
            const std::pair<std::uint64_t, std::uint32_t>& lhs,
            const std::pair<std::uint64_t, std::uint32_t>& rhs) {
                if (lhs.first < rhs.first) return true;
                if (lhs.first == rhs.first) return lhs.second < rhs.second;
                return false;
            });

        auto matches = ram::map(read_minimizers, minimizers, hash);
        std::cout << "matches size " << matches.size() << std::endl;

        auto is_contained = find_best_matches(matches);
        std::cout << "is contained: " << is_contained << std::endl;

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

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

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char** argv) {

    std::uint32_t l = 1000;
    std::uint32_t k = 15;
    std::uint32_t w = 5;
    double f = 0.001;

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

    sort(containee_reads.begin(), containee_reads.end(),
        [](const std::unique_ptr<Sequence>& lhs,
           const std::unique_ptr<Sequence>& rhs) {
            return lhs->data.size() < rhs->data.size();
        }
    );
    sort(contained_reads.begin(), contained_reads.end(),
        [](const std::unique_ptr<Sequence>& lhs,
           const std::unique_ptr<Sequence>& rhs) {
            return lhs->data.size() < rhs->data.size();
        }
    );

    uint32_t id = 0;
    std::vector<std::pair<uint64_t, uint64_t>> minimizers;
    for (const auto& it: containee_reads) {
        ram::createMinimizers(minimizers, it->data.c_str(), it->data.size(),
            id++, k, w);
    }

    std::cerr << "Number of minimizers: " << minimizers.size() << std::endl;

    logger.log("[ram::] collected minimizers in");
    logger.log();

    ram::sortMinimizers(minimizers, k);

    logger.log("[ram::] sorted minimizers in");
    logger.log();

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
    std::uint32_t max_occurence = counts[(1 - f) * counts.size()];

    std::cerr << "Hash size: " << hash.size() << std::endl;
    std::cerr << "Counts size " << counts.size() << std::endl;
    std::cerr << "Max occurence " << max_occurence << std::endl;

    logger.log("[ram::] found occurences in");
    logger.log();

    uint32_t num_short_reads = 0;
    uint32_t num_contained = 0;

    id = 0;
    for (const auto& it: contained_reads) {

        if (it->data.size() <= 2 * l) {
            ++num_short_reads;
            ++id;
            continue;
        }

        std::vector<std::pair<uint64_t, uint64_t>> read_minimizers;

        ram::createMinimizers(read_minimizers, it->data.c_str(), l, id, k, w);
        ram::createMinimizers(read_minimizers, it->data.c_str() +
            (it->data.size() - l - 1), l, id + 1, k, w);

        std::sort(read_minimizers.begin(), read_minimizers.end());

        num_contained += ram::is_read_contained(read_minimizers, minimizers,
            hash, it->data.size() - 2 * l, id, max_occurence);

        ++id;
    }

    std::cout << "Number of reads: " << contained_reads.size() << std::endl;
    std::cout << "Number of short reads: " << num_short_reads << std::endl;
    std::cout << "Number of contained reads: " << num_contained << std::endl;

    logger.log("[ram::] mapped in");
    logger.total("[ram::] total time");

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

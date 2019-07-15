#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"

#include "ram/ram.hpp"

static const std::string version = "v0.0.3";

static struct option options[] = {
    {"kmer-length", required_argument, nullptr, 'k'},
    {"window-length", required_argument, nullptr, 'w'},
    {"filter-threshold", required_argument, nullptr, 'f'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

constexpr std::uint32_t kChunkSize = 1024 * 1024 * 1024; // ~1GB

void help();

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

template<typename T>
void shrinkToFit(std::vector<T>& src, std::uint64_t begin) {

    std::uint64_t i = begin;
    for (std::uint64_t j = begin; i < src.size(); ++i) {
        if (src[i] != nullptr) {
            continue;
        }

        j = std::max(j, i);
        while (j < src.size() && src[j] == nullptr) {
            ++j;
        }

        if (j >= src.size()) {
            break;
        } else if (i != j) {
            std::swap(src[i], src[j]);
        }
    }
    if (i < src.size()) {
        src.resize(i);
    }
}

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path) {

    if (isSuffix(path, ".fasta")    || isSuffix(path, ".fa") ||
        isSuffix(path, ".fasta.gz") || isSuffix(path, ".fa.gz")) {
        return bioparser::createParser<bioparser::FastaParser, ram::Sequence>(path);
    }
    if (isSuffix(path, ".fastq")    || isSuffix(path, ".fq") ||
        isSuffix(path, ".fastq.gz") || isSuffix(path, ".fq.gz")) {
        return bioparser::createParser<bioparser::FastqParser, ram::Sequence>(path);
    }

    std::cerr << "[ram::] error: file " << path
              << " has unsupported format extension (valid extensions: .fasta, "
              << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!"
              << std::endl;
    return nullptr;
}

int main(int argc, char** argv) {

    std::uint32_t k = 15;
    std::uint32_t w = 5;
    double f = 0.001;
    std::uint32_t num_threads = 1;

    std::vector<std::string> input_paths;

    char argument;
    while ((argument = getopt_long(argc, argv, "k:w:f:t:h", options, nullptr)) != -1) {
        switch (argument) {
            case 'k': k = atoi(optarg); break;
            case 'w': w = atoi(optarg); break;
            case 'f': f = atof(optarg); break;
            case 't': num_threads = atoi(optarg); break;
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

    std::unique_ptr<bioparser::Parser<ram::Sequence>> tparser = createParser(input_paths[0]);
    if (tparser == nullptr) {
        return 1;
    }

    auto thread_pool = thread_pool::createThreadPool(num_threads);
    auto logger = logger::Logger();
    ram::MinimizerEngine minimizer_engine(k, w, num_threads);

    logger.log();

    std::vector<std::unique_ptr<ram::Sequence>> sequences;
    tparser->parse(sequences, -1);

    logger.log("[ram::] parsed targets in");
    logger.log();

    minimizer_engine.minimize(sequences.begin(), sequences.end());
    minimizer_engine.filter(f);

    logger.log("[ram::] created minimizers in");

    std::vector<std::unique_ptr<bioparser::Parser<ram::Sequence>>> sparsers;

    if (input_paths.size() > 1) {
        for (std::uint32_t i = 1; i < input_paths.size(); ++i) {
            sparsers.emplace_back(createParser(input_paths[i]));
            if (sparsers.back() == nullptr) {
                return 1;
            }
        }
    } else {
        sparsers.emplace_back(createParser(input_paths[0]));
    }

    for (const auto& it: sparsers) {
        while (true) {
            std::uint32_t l = sequences.size();

            logger.log();

            bool status = it->parse(sequences, kChunkSize);

            logger.log("[ram::] parsed chunk of sequences in");
            logger.log();

            std::vector<std::future<void>> thread_futures;
            for (std::uint32_t i = l; i < sequences.size(); ++i) {
                thread_futures.emplace_back(thread_pool->submit(
                    [&] (std::uint32_t i) -> void {
                        auto overlaps = minimizer_engine.map(sequences[i], true, true);
                    }
                , i));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }

            logger.log("[ram::] mapped chunk of sequences in");
            logger.log();

            sequences.resize(l);

            if (!status) {
                break;
            }
        }
    }

    logger.total("[ram::] total time");

    return 0;
}

void help() {
    std::cout <<
        "usage: ram [options ...] <sequences> [<sequences> ...]\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "\n"
        "    options:\n"
        "        -k, --kmer-length <int>\n"
        "            default: 15\n"
        "            length of minimizers\n"
        "        -w, --window-length <int>\n"
        "            default: 5\n"
        "            window length from which minimizers are found\n"
        "        -f, --filter-threshold <float>\n"
        "            default: 0.001\n"
        "            threshold for ignoring most frequent minimizers\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}

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

static const std::string version = "v0.0.7";

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

    std::unique_ptr<bioparser::Parser<ram::Sequence>> tparser =
        createParser(input_paths[0]);
    if (tparser == nullptr) {
        return 1;
    }

    std::shared_ptr<thread_pool::ThreadPool> thread_pool;
    try {
        thread_pool = thread_pool::createThreadPool(num_threads);
    } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    auto logger = logger::Logger();
    auto minimizer_engine = ram::createMinimizerEngine(k, w, thread_pool);

    logger.log();

    std::vector<std::unique_ptr<ram::Sequence>> sequences;
    try {
        tparser->parse(sequences, -1);
    } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    if (sequences.size() == 0) {
        std::cerr << "[ram::] empty sequence file" << std::endl;
        return 1;
    }

    logger.log("[ram::] parsed targets in");
    logger.log();

    minimizer_engine->minimize(sequences, f);

    logger.log("[ram::] created minimizers in");

    std::vector<std::unique_ptr<bioparser::Parser<ram::Sequence>>> sparsers;
    std::vector<std::uint8_t> is_equal_sparser;

    if (input_paths.size() > 1) {
        for (std::uint32_t i = 1; i < input_paths.size(); ++i) {
            sparsers.emplace_back(createParser(input_paths[i]));
            if (sparsers.back() == nullptr) {
                return 1;
            }
            is_equal_sparser.emplace_back(input_paths[0].compare(input_paths[i]) == 0 ?
                true : false);
        }
    } else {
        sparsers.emplace_back(createParser(input_paths[0]));
        is_equal_sparser.emplace_back(true);
    }

    for (std::uint32_t i = 0; i < sparsers.size(); ++i) {

        ram::Sequence::num_objects = is_equal_sparser[i] ? 0 : sequences.size();

        while (true) {
            std::uint32_t l = sequences.size();

            logger.log();

            bool status;
            try {
                status = sparsers[i]->parse(sequences, kChunkSize);
            } catch (std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
                return 1;
            }

            logger.log("[ram::] parsed chunk of sequences in");
            logger.log();

            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t j = l; j < sequences.size(); ++j) {
                thread_futures.emplace_back(thread_pool->submit(
                    [&] (std::uint32_t j) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[j], is_equal_sparser[i],
                            is_equal_sparser[i]);
                    }
                , j));
            }
            for (std::uint32_t j = 0; j < thread_futures.size(); ++j) {
                thread_futures[j].wait();
                auto overlaps = thread_futures[j].get();
                for (const auto& it: overlaps) {
                    std::cout << sequences[l + j]->name << "\t"
                              << sequences[l + j]->data.size() << "\t"
                              << it.q_begin << "\t"
                              << it.q_end << "\t"
                              << (it.strand ? "+" : "-") << "\t"
                              << sequences[it.t_id]->name << "\t"
                              << sequences[it.t_id]->data.size() << "\t"
                              << it.t_begin << "\t"
                              << it.t_end << "\t"
                              << it.matches << "\t"
                              << std::max(it.q_end - it.q_begin, it.t_end - it.t_begin)<< "\t"
                              << 255
                              << std::endl;
                }
            }

            logger.log("[ram::] mapped chunk of sequences in");

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

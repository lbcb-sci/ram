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

static const std::string version = "v0.0.17";

static struct option options[] = {
    {"kmer-length", required_argument, nullptr, 'k'},
    {"window-length", required_argument, nullptr, 'w'},
    {"filter-threshold", required_argument, nullptr, 'f'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path);

void help();

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

    bool is_equal = false;
    std::unique_ptr<bioparser::Parser<ram::Sequence>> sparser = nullptr;
    if (input_paths.size() > 1) {
        sparser = createParser(input_paths[1]);
        if (sparser == nullptr) {
            return 1;
        }
        is_equal = input_paths[0] == input_paths[1];
    } else {
        sparser = createParser(input_paths[0]);
        is_equal = true;
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

    while (true) {
        std::vector<std::unique_ptr<ram::Sequence>> targets;

        bool status;
        try {
            status = tparser->parse(targets, 1U << 30);
        } catch (std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return 1;
        }

        std::uint32_t num_targets = ram::Sequence::num_objects;
        std::uint32_t t_offset = targets.empty() ? 0 : targets.front()->id;

        ram::Sequence::num_objects = 0;

        logger.log();

        minimizer_engine->minimize(targets);
        minimizer_engine->filter(f);

        logger.log("[ram::] created minimizers of " + std::to_string(targets.size()) + " sequences in");

        while (true) {
            std::vector<std::unique_ptr<ram::Sequence>> sequences;

            bool status_s;
            try {
                status_s = sparser->parse(sequences, 1ULL << 32);
            } catch (std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
                return 1;
            }

            std::uint32_t q_offset = sequences.empty() ? 0 : sequences.front()->id;

            logger.log();

            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t i = 0; i < sequences.size(); ++i) {
                thread_futures.emplace_back(thread_pool->submit(
                    [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[i], is_equal,
                            is_equal);
                    }
                , i));
            }
            for (auto& it: thread_futures) {
                for (const auto& jt: it.get()) {
                    std::cout << sequences[jt.q_id - q_offset]->name << "\t"
                              << sequences[jt.q_id - q_offset]->data.size() << "\t"
                              << jt.q_begin << "\t"
                              << jt.q_end << "\t"
                              << (jt.strand ? "+" : "-") << "\t"
                              << targets[jt.t_id - t_offset]->name << "\t"
                              << targets[jt.t_id - t_offset]->data.size() << "\t"
                              << jt.t_begin << "\t"
                              << jt.t_end << "\t"
                              << jt.matches << "\t"
                              << std::max(jt.q_end - jt.q_begin, jt.t_end - jt.t_begin)<< "\t"
                              << 255
                              << std::endl;
                }
            }

            logger.log("[ram::] mapped " + std::to_string(sequences.size()) + " sequences in");

            if (!status_s) {
                sparser->reset();
                break;
            }
        }

        ram::Sequence::num_objects = num_targets;

        if (!status) {
            break;
        }
    }

    logger.total("[ram::]");

    return 0;
}

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path) {

    auto is_suffix = [] (const std::string& src, const std::string& suffix) {
        return src.size() < suffix.size() ? false :
            src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
        is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
        return bioparser::createParser<bioparser::FastaParser, ram::Sequence>(path);
    }
    if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
        is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
        return bioparser::createParser<bioparser::FastqParser, ram::Sequence>(path);
    }

    std::cerr << "[ram::] error: file " << path
              << " has unsupported format extension (valid extensions: .fasta, "
              << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!"
              << std::endl;
    return nullptr;
}

void help() {
    std::cout <<
        "usage: ram [options ...] <sequences> [<sequences>]\n"
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

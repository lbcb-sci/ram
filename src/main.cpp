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

static const std::string version = "v0.0.2";

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

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

int main(int argc, char** argv) {

    std::uint32_t e = 1000;
    std::uint32_t k = 15;
    std::uint32_t w = 5;
    double f = 0.0001;
    std::uint32_t num_threads = 1;

    std::vector<std::string> input_paths;

    char argument;
    while ((argument = getopt_long(argc, argv, "t:h", options, nullptr)) != -1) {
        switch (argument) {
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

    std::unique_ptr<bioparser::Parser<ram::Sequence>> sparser = nullptr;

    if (isSuffix(input_paths[0], ".fasta") || isSuffix(input_paths[0], ".fa") ||
        isSuffix(input_paths[0], ".fasta.gz") || isSuffix(input_paths[0], ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, ram::Sequence>(
            input_paths[0]);
    } else if (isSuffix(input_paths[0], ".fastq") || isSuffix(input_paths[0], ".fq") ||
               isSuffix(input_paths[0], ".fastq.gz") || isSuffix(input_paths[0], ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, ram::Sequence>(
            input_paths[0]);
    } else {
        std::cerr << "[ram::] error: file " << input_paths[0] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
            std::endl;
        return 1;
    }

    std::vector<std::unique_ptr<bioparser::Parser<ram::Sequence>>> tparsers;

    if (input_paths.size() > 1) {
        for (std::uint32_t i = 1; i < input_paths.size(); ++i) {
            if (isSuffix(input_paths[i], ".fasta") || isSuffix(input_paths[i], ".fa") ||
                isSuffix(input_paths[i], ".fasta.gz") || isSuffix(input_paths[i], ".fa.gz")) {
                tparsers.emplace_back(bioparser::createParser<bioparser::FastaParser, ram::Sequence>(
                    input_paths[i]));
            } else if (isSuffix(input_paths[i], ".fastq") || isSuffix(input_paths[i], ".fq") ||
                       isSuffix(input_paths[i], ".fastq.gz") || isSuffix(input_paths[i], ".fq.gz")) {
                tparsers.emplace_back(bioparser::createParser<bioparser::FastqParser, ram::Sequence>(
                    input_paths[i]));
            } else {
                std::cerr << "[ram::] error: file " << input_paths[i] <<
                    " has unsupported format extension (valid extensions: .fasta, "
                    ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
                    std::endl;
                return 1;
            }
        }
    }

    auto thread_pool = thread_pool::createThreadPool(num_threads);

    auto logger = logger::Logger();

    ram::MinimizerEngine minimizer_engine(k, w, num_threads);

    std::vector<std::unique_ptr<ram::Sequence>> sequences;
    std::uint32_t num_chunks = 0;
    while (true) {
        ++num_chunks;

        logger.log();

        bool status = sparser->parse(sequences, -1);
        std::sort(sequences.begin(), sequences.end(),
            [] (const std::unique_ptr<ram::Sequence>& lhs,
                const std::unique_ptr<ram::Sequence>& rhs) -> bool {
                return lhs->data.size() > rhs->data.size();
            });
        for (std::uint32_t i = 0; i < sequences.size(); ++i) {
            sequences[i]->id = i;
        }

        logger.log("[ram::] parsed sequences in");
        logger.log();

        minimizer_engine.minimize(sequences.begin(),
            sequences.begin() + 0.1 * sequences.size(), f);

        logger.log("[ram::] created minimizers in");
        logger.log();

        std::vector<uint8_t> is_valid(sequences.size(), 1);

        std::vector<std::future<void>> thread_futures;
        for (std::uint32_t i = 0; i < sequences.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit(
                [&] (std::uint32_t i) -> void {
                    auto overlaps = minimizer_engine.map(sequences[i], true, true, e);
                    for (const auto& it: overlaps) {
                        std::uint32_t overhang = std::min(it.t_begin, it.q_begin) +
                            std::min(sequences[i]->data.size() - it.q_end,
                            sequences[it.t_id]->data.size() - it.t_end);

                        if (it.t_end - it.t_begin > (it.t_end - it.t_begin + overhang) * 0.875 &&
                            it.q_end - it.q_begin > (it.q_end - it.q_begin + overhang) * 0.875 &&
                            sequences[it.t_id]->data.size() - it.t_end >= sequences[i]->data.size() - it.q_end &&
                            it.t_begin >= it.q_begin) {

                            is_valid[i] = 0;
                            break;
                        }
                    }
                }
            , i));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }

        logger.log("[ram::] mapped in");

        for (std::uint32_t i = 0; i < sequences.size(); ++i) {
            if (!is_valid[i]) {
                sequences[i].reset();
            }
        }
        shrinkToFit(sequences, 0);

        if (!status) {
            break;
        }
    }

    std::cerr << "[ram::] num uncontained sequences = " << sequences.size() << std::endl;
    logger.total("[ram::] total time");

    /*for (const auto& it: sequences) {
        std::cout << ">" << it->name << std::endl;
        std::cout << it->data << std::endl;
    }*/

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
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}

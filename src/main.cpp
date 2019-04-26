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
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"

static const std::string version = "v0.0.1";

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
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

template<typename T>
void shrinkToFit(std::vector<T>& src, uint64_t begin) {

    uint64_t i = begin;
    for (uint64_t j = begin; i < src.size(); ++i) {
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

void removeContained(std::vector<std::unique_ptr<Sequence>>& sequences,
    std::uint32_t b, std::uint32_t e, std::uint32_t k, std::uint32_t w,
    double f, const std::unique_ptr<thread_pool::ThreadPool>& thread_pool) {

    logger::Logger logger;
    logger.log();

    std::sort(sequences.begin() + b, sequences.end(),
        [] (const std::unique_ptr<Sequence>& lhs,
            const std::unique_ptr<Sequence>& rhs) {
            return lhs->data.size() > rhs->data.size();
        }
    );

    std::vector<std::pair<uint64_t, uint64_t>> minimizers;
    for (std::uint32_t i = b; i < sequences.size(); ++i) {
        ram::createMinimizers(minimizers, sequences[i]->data.c_str(),
            sequences[i]->data.size(), i - b, k, w);
    }

    std::cerr << "Num minimizers: " << minimizers.size() << std::endl;

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
            if (count > 1) {
                hash[minimizers[i - 1].first] = std::make_pair(i - count, count);
            }
            counts.emplace_back(count);
            count = 0;
        }
        ++count;
    }

    logger.log("[ram::] created hash in");
    logger.log();

    std::nth_element(counts.begin(), counts.begin() + (1 - f) * counts.size(),
        counts.end());
    std::uint32_t max_occurence = counts[(1 - f) * counts.size()];

    std::cerr << "Hash size: " << hash.size() << std::endl;
    std::cerr << "Counts size: " << counts.size() << std::endl;
    std::cerr << "Max occurence: " << max_occurence << std::endl;

    logger.log("[ram::] found occurences in");
    logger.log();

    std::vector<std::uint32_t> sequence_lengths(sequences.size() - b);
    for (std::uint32_t i = b; i < sequences.size(); ++i) {
        sequence_lengths[i - b] = sequences[i]->data.size();
    }

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = b; i < sequences.size(); ++i) {

        thread_futures.emplace_back(thread_pool->submit(
            [&](std::uint32_t j) -> void {
                std::vector<std::pair<uint64_t, uint64_t>> sequence_minimizers;
                if (sequences[j]->data.size() <= 2 * e) {
                    ram::createMinimizers(sequence_minimizers, sequences[j]->data.c_str(),
                        sequences[j]->data.size(), j - b, k, w);
                } else {
                    ram::createMinimizers(sequence_minimizers, sequences[j]->data.c_str(),
                        e, j - b, k, w);
                    ram::createMinimizers(sequence_minimizers, sequences[j]->data.c_str() +
                        sequences[j]->data.size() - e, e, j - b + 1, k, w);
                }

                std::sort(sequence_minimizers.begin(), sequence_minimizers.end());

                bool ic = ram::is_contained(sequence_minimizers, minimizers,
                    hash, j - b, sequences[j]->data.size() - e, max_occurence,
                    sequence_lengths);

                if (ic) {
                    sequences[j].reset();
                }
            }, i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }

    shrinkToFit(sequences, b);

    logger.log("[ram::] mapped in");
}

int main(int argc, char** argv) {

    std::uint32_t e = 1000;
    std::uint32_t k = 15;
    std::uint32_t w = 5;
    double f = 0.001;
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

    std::vector<std::unique_ptr<bioparser::Parser<Sequence>>> tparsers;

    if (input_paths.size() > 1) {
        for (std::uint32_t i = 1; i < input_paths.size(); ++i) {
            if (isSuffix(input_paths[i], ".fasta") || isSuffix(input_paths[i], ".fa") ||
                isSuffix(input_paths[i], ".fasta.gz") || isSuffix(input_paths[i], ".fa.gz")) {
                tparsers.emplace_back(bioparser::createParser<bioparser::FastaParser, Sequence>(
                    input_paths[i]));
            } else if (isSuffix(input_paths[i], ".fastq") || isSuffix(input_paths[i], ".fq") ||
                       isSuffix(input_paths[i], ".fastq.gz") || isSuffix(input_paths[i], ".fq.gz")) {
                tparsers.emplace_back(bioparser::createParser<bioparser::FastqParser, Sequence>(
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

    std::vector<std::unique_ptr<Sequence>> sequences;
    std::uint32_t num_chunks = 0;
    while (true) {
        ++num_chunks;

        logger.log();

        std::uint32_t l = sequences.size();
        bool status = sparser->parse(sequences, 256 * 1024 * 1024);

        std::cerr << "Num sequences: " << sequences.size() - l << std::endl;

        logger.log("[ram::] parsed sequences in");

        removeContained(sequences, l, e, k, w, f, thread_pool);

        std::uint64_t remaining_size = 0;
        for (std::uint32_t i = 0; i < sequences.size(); ++i) {
            remaining_size += sequences[i]->data.size();
        }

        if (remaining_size > 256 * 1024 * 1024 || (!status && num_chunks > 1)) {
            removeContained(sequences, 0, e, k, w, f, thread_pool);
        }

        if (!status) {
            break;
        }
    }

    std::cerr << "Num suriving reads: " << sequences.size() << std::endl;
    logger.total("[ram::] total time");

    for (const auto& it: sequences) {
        std::cout << ">" << it->name << std::endl;
        std::cout << it->data << std::endl;
    }

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

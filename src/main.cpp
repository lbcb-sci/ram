#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <unordered_map>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"

#include "minimizers.hpp"
#include "aligner.hpp"

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
            : id(sequence_id++), name(name, name_length), data(data, data_length) {
    }

    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length,
        const char*, std::uint32_t)
            : Sequence(name, name_length, data, data_length) {
    }

    ~Sequence() {
    }

    static std::uint32_t sequence_id;

    std::uint32_t id;
    std::string name;
    std::string data;
};

std::uint32_t Sequence::sequence_id = 0;

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

void removeContained(std::vector<std::unique_ptr<Sequence>>& sequences,
    std::uint32_t b, std::uint32_t lid, std::uint32_t e, std::uint32_t k,
    std::uint32_t w, double f,
    const std::unique_ptr<thread_pool::ThreadPool>& thread_pool) {

    logger::Logger logger;
    logger.log();

    std::sort(sequences.begin() + b, sequences.end(),
        [] (const std::unique_ptr<Sequence>& lhs,
            const std::unique_ptr<Sequence>& rhs) -> bool {
            return lhs->data.size() > rhs->data.size();
        }
    );

    std::vector<std::vector<ram::uint128_t>> minimizers(sequences.size() - b);
    std::vector<std::future<void>> thread_futures;
    std::uint32_t num_minimizers = 0;
    for (std::uint32_t i = b; i < sequences.size(); ++i) {
        if (sequences[i]->data.size() < 10000) {
            break;
        }
        thread_futures.emplace_back(thread_pool->submit(
            [&minimizers, &sequences, b, k, w] (std::uint32_t i) -> void {
                ram::createMinimizers(minimizers[i - b], sequences[i]->data.c_str(),
                    sequences[i]->data.size(), i, k, w);
            }
        , i));
    }
    for (std::uint32_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();
        num_minimizers += minimizers[i].size();
    }
    thread_futures.clear();

    std::cerr << "Num minimizers: " << num_minimizers << std::endl;

    logger.log("[ram::] collected minimizers in");
    logger.log();

    std::vector<std::vector<ram::uint128_t>> hash;
    std::vector<std::unordered_map<std::uint64_t, ram::uint128_t>> index;
    ram::transformMinimizers(hash, index, minimizers, k, thread_pool);

    logger.log("[ram::] created hash in");
    logger.log();

    std::vector<std::uint32_t> counts;
    for (std::uint32_t i = 0; i < index.size(); ++i) {
        for (const auto& it: index[i]) {
            counts.emplace_back(it.second.second);
        }
    }

    std::nth_element(counts.begin(), counts.begin() + (1 - f) * counts.size(),
        counts.end());
    std::uint32_t max_occurence = counts[(1 - f) * counts.size()];

    std::cerr << "Counts size: " << counts.size() << std::endl;
    std::cerr << "Max occurence: " << max_occurence << std::endl;

    logger.log("[ram::] found occurences in");
    logger.log();

    std::vector<std::uint32_t> sequence_lengths;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        sequence_lengths.emplace_back(sequences[i]->data.size());
    }

    std::vector<std::uint32_t> id_map;
    for (std::uint32_t i = 0, j = b; i < b; ++i) {
        while (j < sequences.size() && sequences[i]->data.size() < sequences[j]->data.size()) {
            ++j;
        }
        id_map.emplace_back(j);
    }

    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&sequences, &sequence_lengths, &hash, &index, &id_map, b, e, k, w, max_occurence] (std::uint32_t i) -> void {
                std::vector<ram::uint128_t> sequence_minimizers;
                if (sequences[i]->data.size() <= 2 * e) {
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str(),
                        sequences[i]->data.size(), i < b ? id_map[i] : i, k, w);
                } else {
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str(),
                        e, i < b ? id_map[i] : i, k, w);
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str() +
                        sequences[i]->data.size() - e, e, (i < b ? id_map[i] : i) + 1, k, w);
                }

                std::sort(sequence_minimizers.begin(), sequence_minimizers.end());

                bool ic = ram::is_contained(sequence_minimizers, hash, index,
                    i < b ? id_map[i] : i, sequences[i]->data.size() - e, max_occurence,
                    sequences[i]->data.size(), sequence_lengths);

                if (ic) {
                    sequences[i].reset();
                }
            }, i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    shrinkToFit(sequences, 0);

    logger.log("[ram::] mapped in");
    logger.log();

    // map new to old
    for (b = 0; b < sequences.size() && sequences[b]->id < lid; ++b);
    if (b == 0) {
        return;
    }
    minimizers.resize(b);

    for (std::uint32_t i = 0; i < b; ++i) {
        if (sequences[i]->data.size() < 10000) {
            break;
        }
        thread_futures.emplace_back(thread_pool->submit(
            [&minimizers, &sequences, k, w] (std::uint32_t i) -> void {
                ram::createMinimizers(minimizers[i], sequences[i]->data.c_str(),
                    sequences[i]->data.size(), i, k, w);
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    logger.log("[ram::] collected 2nd part of minimizers in");
    logger.log();

    ram::transformMinimizers(hash, index, minimizers, k, thread_pool);

    logger.log("[ram::] created 2nd part of hash in");
    logger.log();

    counts.clear();
    for (std::uint32_t i = 0; i < index.size(); ++i) {
        for (const auto& it: index[i]) {
            counts.emplace_back(it.second.second);
        }
    }

    std::nth_element(counts.begin(), counts.begin() + (1 - f) * counts.size(),
        counts.end());
    max_occurence = counts[(1 - f) * counts.size()];

    std::cerr << "Counts size: " << counts.size() << std::endl;
    std::cerr << "Max occurence: " << max_occurence << std::endl;

    logger.log("[ram::] found 2nd part of occurences in");
    logger.log();

    sequence_lengths.clear();
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        sequence_lengths.emplace_back(sequences[i]->data.size());
    }

    id_map.clear();
    for (std::uint32_t i = b, j = 0; i < sequences.size(); ++i) {
        while (j < b && sequences[i]->data.size() < sequences[j]->data.size()) {
            ++j;
        }
        id_map.emplace_back(j);
    }

    for (std::uint32_t i = b; i < sequences.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&sequences, &sequence_lengths, &hash, &index, &id_map, b, e, k, w, max_occurence] (std::uint32_t i) -> void {
                std::vector<ram::uint128_t> sequence_minimizers;
                if (sequences[i]->data.size() <= 2 * e) {
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str(),
                        sequences[i]->data.size(), id_map[i - b], k, w);
                } else {
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str(),
                        e, id_map[i - b], k, w);
                    ram::createMinimizers(sequence_minimizers, sequences[i]->data.c_str() +
                        sequences[i]->data.size() - e, e, id_map[i - b] + 1, k, w);
                }

                std::sort(sequence_minimizers.begin(), sequence_minimizers.end());

                bool ic = ram::is_contained(sequence_minimizers, hash, index,
                    id_map[i - b], sequences[i]->data.size() - e, max_occurence,
                    sequences[i]->data.size(), sequence_lengths);

                if (ic) {
                    sequences[i].reset();
                }
            }, i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    shrinkToFit(sequences, 0);

    logger.log("[ram::] 2nd part mapped in");
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
        std::uint32_t lid = Sequence::sequence_id;
        bool status = sparser->parse(sequences, 512 * 1024 * 1024);

        std::cerr << "Num sequences: " << sequences.size() - l << std::endl;

        logger.log("[ram::] parsed sequences in");

        removeContained(sequences, l, lid, e, k, w, f, thread_pool);

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

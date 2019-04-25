#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <algorithm>

#include "bioparser/bioparser.hpp"

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

    if (argc < 2) {
        std::cerr << "[ram::] error: file missing input files!" << std::endl;
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> fparser = nullptr;
    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr;

    if (isSuffix(argv[1], ".fasta") || isSuffix(argv[1], ".fa") ||
        isSuffix(argv[1], ".fasta.gz") || isSuffix(argv[1], ".fa.gz")) {
        fparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            argv[1]);
    } else if (isSuffix(argv[1], ".fastq") || isSuffix(argv[1], ".fq") ||
               isSuffix(argv[1], ".fastq.gz") || isSuffix(argv[1], ".fq.gz")) {
        fparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            argv[1]);
    } else {
        std::cerr << "[ram::] error: file " << argv[1] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
            std::endl;
        return 1;
    }

    if (isSuffix(argv[2], ".fasta") || isSuffix(argv[2], ".fa") ||
        isSuffix(argv[2], ".fasta.gz") || isSuffix(argv[2], ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            argv[2]);
    } else if (isSuffix(argv[2], ".fastq") || isSuffix(argv[2], ".fq") ||
               isSuffix(argv[2], ".fastq.gz") || isSuffix(argv[2], ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            argv[2]);
    } else {
        std::cerr << "[ram::] error: file " << argv[2] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
            std::endl;
        return 1;
    }

    std::vector<std::unique_ptr<Sequence>> f;
    fparser->parse(f, -1);

    std::vector<std::unique_ptr<Sequence>> s;
    sparser->parse(s, -1);

    std::sort(f.begin(), f.end(),
        [] (const std::unique_ptr<Sequence>& lhs,
            const std::unique_ptr<Sequence>& rhs) {
            return lhs->name.compare(rhs->name) < 0;
        });

    std::sort(s.begin(), s.end(),
        [] (const std::unique_ptr<Sequence>& lhs,
            const std::unique_ptr<Sequence>& rhs) {
            return lhs->name.compare(rhs->name) < 0;
        });

    std::cerr << "First size = " << f.size() << std::endl;
    std::cerr << "Second size = " << s.size() << std::endl;

    std::uint32_t num_f_in_s = 0;
    for (std::uint32_t i = 0; i < f.size(); ++i) {
        bool is_present = false;
        for (std::uint32_t j = 0; j < s.size(); ++j) {
            auto c = f[i]->name.compare(s[j]->name);
            if (c == 0) {
                is_present = true;
                ++num_f_in_s;
                break;
            } else if (c < 0) {
                //break;
            }
        }
        if (!is_present) {
            // std::cerr << "Missing read " << f[i]->name << " with size of " <<
            //     f[i]->data.size() << std::endl;
        }
    }

    std::cerr << "Num of first in second = " << num_f_in_s << std::endl;

    return 0;
}

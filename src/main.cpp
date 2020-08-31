// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <cstdlib>
#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"

#include "ram/minimizer_engine.hpp"

std::atomic<std::uint64_t> biosoup::Sequence::num_objects{0};

namespace {

const char* ram_version = RAM_VERSION;

static struct option options[] = {
  {"kmer-length", required_argument, nullptr, 'k'},
  {"window-length", required_argument, nullptr, 'w'},
  {"frequency-threshold", required_argument, nullptr, 'f'},
  {"micromize", no_argument, nullptr, 'm'},
  {"wildcards", required_argument, nullptr, 'W'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"rna", no_argument, nullptr, 'r'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>>
    CreateParser(const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[ram::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: ram [options ...] <sequences> [<sequences>]\n"
      "\n"
      "  # default output is stdout\n"
      "  <sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -k, --kmer-length <int>\n"
      "      default: 15\n"
      "      length of minimizers\n"
      "    -w, --window-length <int>\n"
      "      default: 5\n"
      "      window length from which minimizers are found\n"
      "    -f, --frequency-threshold <float>\n"
      "      default: 0.001\n"
      "      threshold for ignoring most frequent minimizers\n"
      "    -m, --micromize\n"
      "      use only a portion of all minimizers\n"
      "    --wildcards <string>\n"
      "      default: 111111111111111\n"
      "      allow wildcards (0) in minimizers\n"
      "    --rna\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::uint32_t k = 15;
  std::uint32_t w = 5;
  double frequency = 0.001;
  bool micromize = false;
  std::string wildcards;
  std::uint64_t wildcard_mask = -1;
  std::uint32_t num_threads = 1;
  std::uint32_t is_rna = 0;

  std::vector<std::string> input_paths;

  const char* optstr = "k:w:f:mt:h";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, nullptr)) != -1) {
    switch (arg) {
      case 'k': k = std::atoi(optarg); break;
      case 'w': w = std::atoi(optarg); break;
      case 'f': frequency = std::atof(optarg); break;
      case 'm': micromize = true; break;
      case 'W': wildcards = optarg; break;
      case 'r': is_rna = 1; break;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << ram_version << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  for (auto i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.empty()) {
    std::cerr << "[ram::] error: missing sequences file" << std::endl;
    return 1;
  }
  if (!wildcards.empty()) {
    if (wildcards.size() != k ||
        wildcards.find_first_not_of("01") != std::string::npos) {
      std::cerr << "[ram::] error: invalid wildcards string" << std::endl;
      return 1;
    }
    for (std::uint32_t i = 0; i < k; ++i) {
      if (wildcards[i] == '0') {
        wildcard_mask ^= 3ULL << (2 * (k - i - 1));
      }
    }
  }

  auto tparser = CreateParser(input_paths[0]);
  if (tparser == nullptr) {
    return 1;
  }

  bool is_ava = false;
  std::unique_ptr<bioparser::Parser<biosoup::Sequence>> sparser = nullptr;
  if (input_paths.size() > 1) {
    sparser = CreateParser(input_paths[1]);
    if (sparser == nullptr) {
      return 1;
    }
    is_ava = input_paths[0] == input_paths[1];
  } else {
    sparser = CreateParser(input_paths[0]);
    is_ava = true;
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
  ram::MinimizerEngine minimizer_engine{k, w, thread_pool, wildcard_mask};

  biosoup::Timer timer{};

  auto byte_len = is_rna ? -1 : 1U << 30;

  while (true) {
    timer.Start();

    std::vector<std::unique_ptr<biosoup::Sequence>> targets;
    try {
      targets = tparser->Parse(byte_len);
    } catch (std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (targets.empty()) {
      break;
    }

    std::cerr << "[ram::] parsed " << targets.size() << " sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    minimizer_engine.Minimize(targets.begin(), targets.end());
    minimizer_engine.Filter(frequency);

    std::cerr << "[ram::] minimized sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::uint64_t num_targets = biosoup::Sequence::num_objects;
    biosoup::Sequence::num_objects = 0;

    bool should_print = true;

    while (true) {
      timer.Start();

      std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
      try {
        sequences = sparser->Parse(byte_len);
      } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
      }

      if (sequences.empty()) {
        break;
      }

      if(is_rna) {
        std::vector<std::future<std::string>> futures;
        for (const auto& it : sequences) {
          futures.emplace_back(thread_pool->Submit(
              [&] (const std::unique_ptr<biosoup::Sequence>& sequence)
                  -> std::string {
                  return minimizer_engine.MapRNA(sequence, is_ava, is_ava, should_print, targets, micromize);
              },
              std::ref(it)));
        }

        biosoup::ProgressBar bar{
            static_cast<std::uint32_t>(sequences.size()), 16};

        for (auto& it : futures) {
          // for (const auto& jt : it.get()) {
          std::cout << it.get();
          // }

          if (++bar) {
            std::cerr << "[ram::] mapped " << bar.event_counter() << " sequences "
                      << "[" << bar << "] "
                      << std::fixed << timer.Lap() << "s"
                      << "\r";
          }
        }
        std::cerr << std::endl;
        timer.Stop();

        should_print = true;

        if (is_ava && biosoup::Sequence::num_objects == num_targets) {
          break;
        }
      } else {
        std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
        for (const auto& it : sequences) {
          futures.emplace_back(thread_pool->Submit(
              [&] (const std::unique_ptr<biosoup::Sequence>& sequence)
                  -> std::vector<biosoup::Overlap> {
                  return minimizer_engine.Map(sequence, is_ava, is_ava, micromize);
              },
              std::ref(it)));
        }

        biosoup::ProgressBar bar{
            static_cast<std::uint32_t>(sequences.size()), 16};

        for (auto& it : futures) {
          std::uint64_t rhs_offset = targets.front()->id;
          std::uint64_t lhs_offset = sequences.front()->id;
          for (const auto& jt : it.get()) {
            std::cout << sequences[(jt.lhs_id >> 1) - lhs_offset]->name << "\t"
                      << sequences[(jt.lhs_id >> 1) - lhs_offset]->data.size() << "\t"  // NOLINT
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.lhs_id & 1 ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->data.size() << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                            jt.lhs_end - jt.lhs_begin,
                            jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          }

          if (++bar) {
            std::cerr << "[ram::] mapped " << bar.event_counter() << " sequences "
                      << "[" << bar << "] "
                      << std::fixed << timer.Lap() << "s"
                      << "\r";
          }
        }
        std::cerr << std::endl;
        timer.Stop();

        if (is_ava && biosoup::Sequence::num_objects == num_targets) {
          break;
        }
      }
    }

    sparser->Reset();
    biosoup::Sequence::num_objects = num_targets;
  }

  std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;

  return 0;
}

// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <bitset>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"
#include "edlib.h"  // NOLINT

#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"kmer-length", required_argument, nullptr, 'k'},
  {"window-length", required_argument, nullptr, 'w'},
  {"frequency-threshold", required_argument, nullptr, 'f'},
  {"bandwidth", required_argument, nullptr, 'b'},
  {"chain", required_argument, nullptr, 'c'},
  {"matches", required_argument, nullptr, 'm'},
  {"gap", required_argument, nullptr, 'g'},
  {"minhash", no_argument, nullptr, 'M'},
  {"begin", required_argument, nullptr, 'B'},
  {"end", required_argument, nullptr, 'E'},
  {"offset", required_argument, nullptr, 'O'},
  {"ignore", required_argument, nullptr, 'I'},
  {"dummy", required_argument, nullptr, 'D'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>>
    CreateParser(const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz")   ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[ram::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: ram [options ...] <target> [<sequences>]\n"
      "\n"
      "  # default output is stdout\n"
      "  <target>/<sequences> \n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -k, --kmer-length <int>\n"
      "      default: 15\n"
      "      length of minimizers\n"
      "    -w, --window-length <int>\n"
      "      default: 5\n"
      "      length of sliding window from which minimizers are sampled\n"
      "    -f, --frequency-threshold <float>\n"
      "      default: 0.001\n"
      "      threshold for ignoring most frequent minimizers\n"
      "    --bandwidth <int>\n"
      "      default: 500\n"
      "      size of bandwidth in which minimizer hits can be chained\n"
      "    --chain <int>\n"
      "      default: 4\n"
      "      minimal number of chained minimizer hits in overlap\n"
      "    --matches <int>\n"
      "      default: 100\n"
      "      minimal number of matching bases in overlap\n"
      "    --gap <int>\n"
      "      default: 10000\n"
      "      maximal gap between minimizer hits in a chain\n"
      "    --minhash\n"
      "      use only a portion of all minimizers\n"
      "    -B, --begin <int>\n"
      "      default: 12000\n"
      "      region of interest begin position\n"
      "    -E, --end <int>\n"
      "      default: 35355\n"
      "      region of interest end position\n"
      "    -O, --offset <int>\n"
      "      default: 100\n"
      "      number of flanking bases left and right of region of interest\n"
      "    -I, --ignore <int>\n"
      "      default: 10\n"
      "      number of bases to skip before offseting\n"
      "    -D, --dummy <int>\n"
      "      default: 5000\n"
      "      size of dummy sequence used for padding\n"
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
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t matches = 100;
  std::uint32_t gap = 10000;
  double frequency = 0.001;
  bool minhash = false;
  std::uint32_t num_threads = 1;

  std::uint32_t begin = 12000;
  std::uint32_t end = 35355;
  std::uint32_t offset = 100;
  std::uint32_t ignore = 10;
  std::uint32_t dummy = 5000;

  std::vector<std::string> input_paths;

  const char* optstr = "k:w:f:t:B:E:O:I:D:h";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, nullptr)) != -1) {
    switch (arg) {
      case 'k': k = std::atoi(optarg); break;
      case 'w': w = std::atoi(optarg); break;
      case 'b': bandwidth = std::atoi(optarg); break;
      case 'c': chain = std::atoi(optarg); break;
      case 'm': matches = std::atoi(optarg); break;
      case 'g': gap = std::atoi(optarg); break;
      case 'f': frequency = std::atof(optarg); break;
      case 'M': minhash = true; break;
      case 'B': begin = std::atoi(optarg); break;
      case 'E': end = std::atoi(optarg); break;
      case 'O': offset = std::atoi(optarg); break;
      case 'I': ignore = std::atoi(optarg); break;
      case 'D': dummy = std::atoi(optarg); break;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
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
    std::cerr << "[ram::] error: missing target file" << std::endl;
    return 1;
  }

  auto tparser = CreateParser(input_paths[0]);
  if (tparser == nullptr) {
    return 1;
  }

  bool is_ava = false;
  std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> sparser = nullptr;
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
  ram::MinimizerEngine minimizer_engine{
      thread_pool,
      k,
      w,
      bandwidth,
      chain,
      matches,
      gap};

  biosoup::Timer timer{};

  std::ofstream os("trimmed.fa");

  while (true) {
    timer.Start();

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
    try {
      targets = tparser->Parse(1ULL << 32);
    } catch (std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (targets.empty()) {
      break;
    }

    std::cerr << "[ram::] parsed " << targets.size() << " targets "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    minimizer_engine.Minimize(targets.begin(), targets.end(), minhash);
    minimizer_engine.Filter(frequency);

    std::cerr << "[ram::] minimized targets "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
    biosoup::NucleicAcid::num_objects = 0;

    while (true) {
      timer.Start();

      std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
      try {
        sequences = sparser->Parse(1U << 30);
      } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
      }

      if (sequences.empty()) {
        break;
      }

      std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
      for (const auto& it : sequences) {
        if (is_ava && it->id >= num_targets) {
          continue;
        }
        futures.emplace_back(thread_pool->Submit(
            [&] (const std::unique_ptr<biosoup::NucleicAcid>& sequence)
                -> std::vector<biosoup::Overlap> {
              return minimizer_engine.Map(sequence, is_ava, is_ava, minhash);
            },
            std::ref(it)));
      }

      biosoup::ProgressBar bar{
          static_cast<std::uint32_t>(futures.size()), 16};

      // std::uint64_t rhs_offset = targets.front()->id;
      std::uint64_t lhs_offset = sequences.front()->id;
      for (auto& it : futures) {
        bool c = false, l = false, r = false;
        auto overlaps = it.get();
        for (const auto& jt : overlaps) {
          // 11900 - 12000 and 35355 - 35455
          if (!(jt.rhs_end < begin || jt.rhs_begin > end)) {
            c = true;
          }
          if (!(jt.rhs_end < (begin - ignore - offset) || jt.rhs_begin > (begin - ignore))) {  // NOLINT
            l = true;
          }
          if (!(jt.rhs_end < (end + ignore) || jt.rhs_begin > (end + ignore + offset))) {  // NOLINT
            r = true;
          }
        }
        if (c && l && r) {
          const auto& s = sequences[overlaps.front().lhs_id - lhs_offset];
          // find bounds
          std::uint32_t b = s->inflated_len, e = 0;
          for (const auto& jt : overlaps) {
            auto q = s->InflateData(jt.lhs_begin, jt.lhs_end - jt.lhs_begin);
            if (!jt.strand) {
              biosoup::NucleicAcid q_{"", q};
              q_.ReverseAndComplement();
              q = q_.InflateData();
            }
            auto r = targets[jt.rhs_id]->InflateData(jt.rhs_begin, jt.rhs_end - jt.rhs_begin);  // NOLINT

            EdlibAlignResult result = edlibAlign(
                q.c_str(), q.size(),
                r.c_str(), r.size(),
                edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));  // NOLINT
            if (result.status == EDLIB_STATUS_OK) {
              std::uint32_t q_pos = jt.strand ? jt.lhs_begin : s->inflated_len - jt.lhs_end;  // NOLINT
              std::uint32_t r_pos = jt.rhs_begin;

              for (int a = 0; a <= result.alignmentLength; ++a) {
                if (r_pos == begin || r_pos == end) {
                  b = std::min(b, jt.strand ? q_pos : s->inflated_len - q_pos - 1);  // NOLINT
                  e = std::max(e, jt.strand ? q_pos : s->inflated_len - q_pos - 1);  // NOLINT
                }
                if (a == result.alignmentLength) {
                  break;
                }
                switch (result.alignment[a]) {
                  case 0:
                  case 3: {
                    ++q_pos;
                    ++r_pos;
                    break;
                  }
                  case 1: {
                    ++q_pos;
                    break;
                  }
                  case 2: {
                    ++r_pos;
                    break;
                  }
                  default: break;
                }
              }
            }
            edlibFreeAlignResult(result);
          }

          // print trimmed sequence
          os << ">" << s->name << "_trimmed" << std::endl
             << s->InflateData(b, e - b) << std::endl;

          // print whole sequence
          std::cout << ">" << s->name << std::endl
                    << s->InflateData() << std::endl
                    << ">dummy" << std::endl
                    << std::string(dummy, 'N') << std::endl;
        }
        //   std::cout << sequences[jt.lhs_id - lhs_offset]->name << "\t"
        //             << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
        //             << jt.lhs_begin << "\t"
        //             << jt.lhs_end << "\t"
        //             << (jt.strand ? "+" : "-") << "\t"
        //             << targets[jt.rhs_id - rhs_offset]->name << "\t"
        //             << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
        //             << jt.rhs_begin << "\t"
        //             << jt.rhs_end << "\t"
        //             << jt.score << "\t"
        //             << std::max(
        //                   jt.lhs_end - jt.lhs_begin,
        //                   jt.rhs_end - jt.rhs_begin) << "\t"
        //             << 255
        //             << std::endl;
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

      if (is_ava && biosoup::NucleicAcid::num_objects >= num_targets) {
        break;
      }
    }

    sparser->Reset();
    biosoup::NucleicAcid::num_objects = num_targets;
  }

  os.close();

  std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;

  return 0;
}

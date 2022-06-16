// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"

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
  {"threads", required_argument, nullptr, 't'},
  {"weighted", optional_argument, nullptr, 'W'},
  {"beginend", optional_argument, nullptr, 'S'},
  {"beginendk", optional_argument, nullptr, 'K'},
  {"beginendw", optional_argument, nullptr, 'I'},
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
      "usage: ram [options ...]\n"
      "\n"
      "  # default output is stdout\n"
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
      "    --weighted\n"
      "      use weighted minimizer sampling\n"
      "    --beginend\n"
      "      if number n (greater then 0) provided, first n and last n bases of the sequences will be minimized with kmer length beginendk and window size beginendw\n"  
      "    --beginendk\n"
      "      kmer length used with beginend option for beginning and end of sequences\n"
      "    --beginendw\n"
      "      window size used with beginend option for beginning and end of sequences\n"
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
  double weightedMinimizerSampling = 0;
  std::uint32_t beginAndEndSequenceLength = 0;
  std::uint32_t beginAndEndSequenceK = 0;
  std::uint32_t beginAndEndSequenceW = 0;

  std::vector<std::string> input_paths;

  const char* optstr = "k:w:f:t:h";
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
      case 't': num_threads = std::atoi(optarg); break;
      case 'W': weightedMinimizerSampling = std::atof(optarg); break;
      case 'S': beginAndEndSequenceLength = std::atoi(optarg); break;
      case 'K': beginAndEndSequenceK = std::atoi(optarg); break;
      case 'I': beginAndEndSequenceW = std::atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  input_paths.emplace_back(TEST_DATA);

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

  std::uint32_t exactOverlapLhsBegin  = 100000;
  std::uint32_t exactOverlapLhsEnd    = 150000;
  std::uint32_t exactOverlapRhsBegin  = 0;
  std::uint32_t exactOverlapRhsEnd    = 50000;

  std::ofstream overlapsBetweenWrongReads;
  overlapsBetweenWrongReads.open ("overlapsBetweenWrongReads.txt");
  std::ofstream exactOverlaps;
  exactOverlaps.open ("exactOverlaps.txt");
  std::ofstream correctOverlaps;
  correctOverlaps.open ("correctOverlaps.txt");
  std::ofstream reaminingOverlaps;
  reaminingOverlaps.open ("reaminingOverlaps.txt");
  std::ofstream overlapsWithRightOverhang;
  overlapsWithRightOverhang.open ("overlapsWithRightOverhang.txt");
  std::ofstream overlapsWithLeftOverhang;
  overlapsWithLeftOverhang.open ("overlapsWithLeftOverhang.txt");
  std::ofstream overlapsWithDoubleOverhang;
  overlapsWithDoubleOverhang.open ("overlapsWithDoubleOverhang.txt");

  std::uint32_t totalNumberOfOverlaps               = 0;
  std::uint32_t numberOfOverlapsBetweenWrongReads   = 0; 
  std::uint32_t numberOfExactOverlaps               = 0;
  std::uint32_t numberOfCorrectOverlaps             = 0;
  std::uint32_t numberOfOverlapsWithRightOverhang   = 0;
  std::uint32_t numberOfOverlapsWithLeftOverhang    = 0;
  std::uint32_t numberOfOverlapsWithDoubleOverhang  = 0;

  std::set<std::string> foundOverlapsBetweenReads;

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

    minimizer_engine.Minimize(targets.begin(), targets.end(), minhash, weightedMinimizerSampling, beginAndEndSequenceLength, beginAndEndSequenceK, beginAndEndSequenceW);
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

      std::uint64_t rhs_offset = targets.front()->id;
      std::uint64_t lhs_offset = sequences.front()->id;
      for (auto& it : futures) {
        for (const auto& jt : it.get()) {
          totalNumberOfOverlaps++;

          foundOverlapsBetweenReads.insert(sequences[jt.lhs_id - lhs_offset]->name + targets[jt.rhs_id - rhs_offset]->name);
          
          /*
          std::cout << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                    << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                    << jt.lhs_begin << "\t"
                    << jt.lhs_end << "\t"
                    << (jt.strand ? "+" : "-") << "\t"
                    << targets[jt.rhs_id - rhs_offset]->name << "\t"
                    << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                    << jt.rhs_begin << "\t"
                    << jt.rhs_end << "\t"
                    << jt.score << "\t"
                    << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                    << 255
                    << std::endl;
          */

          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 != stoi(targets[jt.rhs_id - rhs_offset]->name)) {
            
            numberOfOverlapsBetweenWrongReads++;
          
            overlapsBetweenWrongReads << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          } else

          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 == stoi(targets[jt.rhs_id - rhs_offset]->name) &&
              jt.lhs_begin == exactOverlapLhsBegin && jt.lhs_end == exactOverlapLhsEnd &&
              jt.rhs_begin == exactOverlapRhsBegin && jt.rhs_end == exactOverlapRhsEnd) {
            
            numberOfExactOverlaps++;
          
            exactOverlaps << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          } else

          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 == stoi(targets[jt.rhs_id - rhs_offset]->name) &&
              jt.lhs_begin >= exactOverlapLhsBegin && jt.lhs_end <= exactOverlapLhsEnd &&
              jt.rhs_begin >= exactOverlapRhsBegin && jt.rhs_end <= exactOverlapRhsEnd &&
              jt.lhs_begin - exactOverlapLhsBegin == jt.rhs_begin - exactOverlapRhsBegin &&
              exactOverlapLhsEnd - jt.lhs_end == exactOverlapRhsEnd - jt.rhs_end) {
            
            numberOfCorrectOverlaps++;
          
            correctOverlaps << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          } else

          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 == stoi(targets[jt.rhs_id - rhs_offset]->name) &&
              jt.lhs_begin < exactOverlapLhsBegin && jt.lhs_end <= exactOverlapLhsEnd &&
              jt.rhs_begin >= exactOverlapRhsBegin && jt.rhs_end <= exactOverlapRhsEnd) {
            
            numberOfOverlapsWithRightOverhang++;
          
            overlapsWithRightOverhang << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          } else
          
          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 == stoi(targets[jt.rhs_id - rhs_offset]->name) &&
              jt.lhs_begin >= exactOverlapLhsBegin && jt.lhs_end <= exactOverlapLhsEnd &&
              jt.rhs_begin >= exactOverlapRhsBegin && jt.rhs_end > exactOverlapRhsEnd) {
            
            numberOfOverlapsWithLeftOverhang++;
          
            overlapsWithLeftOverhang << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          } else

          if (stoi(sequences[jt.lhs_id - lhs_offset]->name) + 1 == stoi(targets[jt.rhs_id - rhs_offset]->name) &&
              jt.lhs_begin < exactOverlapLhsBegin && jt.lhs_end <= exactOverlapLhsEnd &&
              jt.rhs_begin >= exactOverlapRhsBegin && jt.rhs_end > exactOverlapRhsEnd) {
            
            numberOfOverlapsWithDoubleOverhang++;
          
            overlapsWithDoubleOverhang << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          }

          else {

            reaminingOverlaps << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                      << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                      << jt.lhs_begin << "\t"
                      << jt.lhs_end << "\t"
                      << (jt.strand ? "+" : "-") << "\t"
                      << targets[jt.rhs_id - rhs_offset]->name << "\t"
                      << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                      << jt.rhs_begin << "\t"
                      << jt.rhs_end << "\t"
                      << jt.score << "\t"
                      << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                      << 255
                      << std::endl;
          }

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

      if (is_ava && biosoup::NucleicAcid::num_objects >= num_targets) {
        break;
      }
    }

    sparser->Reset();
    biosoup::NucleicAcid::num_objects = num_targets;
  }

  std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;

  overlapsBetweenWrongReads.close();
  exactOverlaps.close();
  correctOverlaps.close();
  reaminingOverlaps.close();
  overlapsWithRightOverhang.close();
  overlapsWithLeftOverhang.close();
  overlapsWithDoubleOverhang.close();

  std::cout << "Total number of overlaps: " << totalNumberOfOverlaps << std::endl;
  std::cout << std::endl;
  std::cout << "Number of exact overlaps: " << numberOfExactOverlaps << std::endl;
  std::cout << "Number of correct overlaps: " << numberOfCorrectOverlaps << std::endl;
  std::cout << "Number of overlaps with right overhang: " << numberOfOverlapsWithRightOverhang << std::endl;
  std::cout << "Number of overlaps with left overhang: " << numberOfOverlapsWithLeftOverhang << std::endl;
  std::cout << "Number of overlaps with overhang on both sides: " << numberOfOverlapsWithDoubleOverhang << std::endl;

  std::cout << "Missing overlaps between reads: " << std::endl;
  for (int i = 1, j = 2; i < 100; i++, j++) {
    if (foundOverlapsBetweenReads.find(std::to_string(i) + std::to_string(j)) == foundOverlapsBetweenReads.end()) {
        std::cout << "Read " << i << " -> Read " << j << std::endl;
    }
  }

  std::cout << "Number of overlaps between wrong reads: " << numberOfOverlapsBetweenWrongReads << std::endl;

  return 0;
}

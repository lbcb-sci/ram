#include <iostream>
#include <string>
#include <vector>

#include <bioparser/bioparser.hpp>

#include "ram_minimizer.hpp"
#include "ram_aligner.hpp"

class Fast {
public:
    std::string name;
    std::string sequence;
    std::string quality;

    Fast(
            const char *name, uint32_t name_length,
            const char *sequence, uint32_t sequence_length) :
            name{std::string(name, name_length)},
            sequence{std::string(sequence, sequence_length)} {}

    Fast(
            const char *name, uint32_t name_length,
            const char *sequence, uint32_t sequence_length,
            const char *quality, uint32_t quality_length) :
            name{std::string(name, name_length)},
            sequence{std::string(sequence, sequence_length)},
            quality{std::string(quality, quality_length)} {}
};

std::vector<std::unique_ptr<Fast>> parse_fasta(std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;

    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq(std::string fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;

    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024;  // 500 MB

    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }
    return fastq_objects;
}

uint64_t get_complement_letter_value(char letter) {
    switch (letter) {
        case 'C':
            return 'G';
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        default:
            return  0;
    }
}

int get_ceil_index(std::vector<Match> matches, std::vector<uint64_t> &tail_indices, uint64_t left, uint64_t right, uint64_t key) {
    while (right - left > 1) {
      int middle = left + (right - left)/2;
      if (matches[tail_indices[middle]].read_start >= key) {
        right = middle;
      } else {
        left = middle;
      }
    }
    return right;
}

std::vector<Match> find_lis_matches(std::vector<Match> &matches) {

  std::vector<uint64_t> tail_indices(matches.size(), 0);
  std::vector<uint64_t> prev_indices(matches.size(), -1);

  uint64_t len = 1;
  for (uint64_t i = 1; i < matches.size(); i++) {
      if (matches[i].read_start < matches[tail_indices[0]].read_start) {
          tail_indices[0] = i;
      } else if (matches[i].read_start > matches[tail_indices[len-1]].read_start) {
          prev_indices[i] = tail_indices[len-1];
          tail_indices[len++] = i;
      } else {
          uint64_t pos = get_ceil_index(matches, tail_indices, -1, len-1, matches[i].read_start);
          prev_indices[i] = tail_indices[pos-1];
          tail_indices[pos] = i;
      }
  }

  std::vector<Match> best_matches;

  for (int i = tail_indices[len-1]; i >= 0; i = prev_indices[i]) {
    best_matches.push_back(matches[i]);
  }
  return best_matches;
}

std::vector<Match> find_best_matches(std::vector<Match> &matches) {
  sort(matches.begin(), matches.end());

  std::vector<Match> forward_matches;
  std::vector<Match> reversed_matches;

  for (auto const &match: matches) {
    if (match.strand) {
      reversed_matches.push_back(match);
    } else {
      forward_matches.push_back(match);
    }
  }

  reverse(reversed_matches.begin(), reversed_matches.end());
  std::vector<Match> best_matches_reversed;
  if (reversed_matches.size() > 0) {
    best_matches_reversed = find_lis_matches(reversed_matches);
  }

  std::vector<Match> best_matches_forward;
  if (forward_matches.size() > 0) {
    best_matches_forward = find_lis_matches(forward_matches);
  }

  if (best_matches_forward.size() > best_matches_reversed.size()) {
    return best_matches_forward;
  } else {
    return best_matches_reversed;
  }
}

std::vector<Match> get_matches(std::vector<Minimizer> &first_sequence, std::vector<Minimizer> &second_sequence) {
  std::vector<Match> matches;
  uint64_t first_minimizers_index = 0;
  uint64_t second_minimizers_index = 0;

  while (first_minimizers_index < first_sequence.size() && second_minimizers_index < second_sequence.size()) {
    Minimizer first_minimizer = first_sequence[first_minimizers_index];
    Minimizer second_minimizer = second_sequence[second_minimizers_index];

    if (first_minimizer.value == second_minimizer.value) {
      uint64_t first_minimizers_inner_index = first_minimizers_index;
      uint64_t second_minimizers_inner_index = second_minimizers_index;

      while (first_minimizers_inner_index < first_sequence.size()) {
        Minimizer first_minimizer_inner = first_sequence[first_minimizers_inner_index];

        if (first_minimizer.value != first_minimizer_inner.value) {
          break;
        }

        second_minimizers_inner_index = second_minimizers_index;
        while (second_minimizers_inner_index < second_sequence.size()) {
          Minimizer second_minimizer_inner = second_sequence[second_minimizers_inner_index];
          if (second_minimizer_inner.value == first_minimizer_inner.value) {
            matches.emplace_back(first_minimizer_inner, second_minimizer_inner);
            second_minimizers_inner_index += 1;
          } else {
            break;
          }
        }
        first_minimizers_inner_index += 1;
      }

      first_minimizers_index = first_minimizers_inner_index;
      second_minimizers_index = second_minimizers_inner_index;
    } else {
      if (first_minimizer.value < second_minimizer.value) {
        first_minimizers_index += 1;
      } else {
        second_minimizers_index += 1;
      }
    }
  }
  return matches;
}

std::vector<Cluster> create_clusters(std::vector<Match> matches, uint64_t kmer_length) {

  std::vector<Cluster> clusters;

  if (matches.size() <= 0) {
    return clusters;
  }

  uint64_t start_ref = matches[0].ref_start;
  uint64_t len_ref = 0;
  uint64_t start_read = matches[0].read_start;
  uint64_t len_read = 0;

  for (uint64_t i = 1; i < matches.size(); i++) {
    Match match = matches[i];
    Match previous_match = matches[i-1];

    uint64_t ref_len = match.ref_start - previous_match.ref_start;
    uint64_t read_len = match.strand ? previous_match.read_start - match.read_start : match.read_start - previous_match.read_start;

    if (read_len < 15 && ref_len > 100) {
      uint64_t adjusted_start_read = match.strand ? start_read-len_read : start_read;
      clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, match.strand);

      start_ref = match.ref_start;
      start_read = match.read_start;
      len_ref = 0;
      len_read = 0;
    } else {
      len_ref += ref_len;
      len_read += read_len;
    }
  }
  if (len_ref > 0) {
    uint64_t adjusted_start_read = matches[0].strand ? start_read-len_read : start_read;
    clusters.emplace_back(start_ref, len_ref+kmer_length, adjusted_start_read, len_read+kmer_length, matches[0].strand);
  }

  return clusters;
}

int main(int argc, char *argv[]) {

    uint64_t kmer_length = 15;
    uint64_t window_length = 5;

    std::vector<std::unique_ptr<Fast>> reads_containees;
    reads_containees = parse_fasta(argv[1]);
    std::string read_containee = reads_containees.front()->sequence;

    std::vector<std::unique_ptr<Fast>> reads_containeds;
    reads_containeds = parse_fasta(argv[2]);
    std::string read_contained = reads_containeds.front()->sequence;

    std::cout << "Containee size: " << read_containee.length() << std::endl;
    std::cout << "Contained size: " << read_contained.length() << std::endl;

    std::vector<Minimizer> minimizers_contained;
    minimizers_contained = create_minimizers(read_contained.c_str(), read_contained.length(), kmer_length, window_length);

    std::cout << "Minimizers contained: " << minimizers_contained.size() << std::endl;

    std::vector<Minimizer> start_minimizers;
    std::vector<Minimizer> end_minimizers;

    for (auto const &minimizer: minimizers_contained) {
      if(minimizer.location < 1000) {
        start_minimizers.push_back(minimizer);
      }
      if(minimizer.location > (read_contained.length()-1000)) {
        end_minimizers.push_back(minimizer);
      }
    }

    std::cout << "Start minimizers size: " << start_minimizers.size() << std::endl;
    std::cout << "End minimizers size: " << end_minimizers.size() << std::endl;

    std::vector<Minimizer> minimizers_containee;
    minimizers_containee = create_minimizers(read_containee.c_str(), read_containee.length(), kmer_length, window_length);

    std::vector<Match> matches_start = get_matches(start_minimizers, minimizers_containee);
    std::vector<Match> result_matches_start = find_best_matches(matches_start);
    std::vector<Cluster> clusters_start = create_clusters(result_matches_start, kmer_length);

    for (auto const &match: result_matches_start) {
      std::cout << match.to_string() << std::endl;
    }

    for (auto const &cluster: clusters_start) {
      std::cout << cluster.to_string() << std::endl;
    }

    std::vector<Match> matches_end = get_matches(end_minimizers, minimizers_containee);
    std::vector<Match> result_matches_end = find_best_matches(matches_end);
    std::vector<Cluster> clusters_end = create_clusters(result_matches_end, kmer_length);

    for (auto const &match: result_matches_end) {
      std::cout << match.to_string() << std::endl;
    }

    for (auto const &cluster: clusters_end) {
      std::cout << cluster.to_string() << std::endl;
    }


    // todo alignment
    std::string cigar;
    uint64_t target_begin = 0;

    Cluster cluster = clusters_start[0];

    std::string t_sub = read_containee.substr(cluster.ref_start, cluster.ref_length);
    std::string q_sub = read_contained.substr(cluster.read_start, cluster.read_length);

    char q_sub_reverse[q_sub.size()];
    for (uint64_t i = 0; i < q_sub.size(); i++) {
      q_sub_reverse[(q_sub.size()-1)-i] = get_complement_letter_value(q_sub[i]);
    }

    const char* t_sub_char = t_sub.c_str();
    
    pairwise_alignment(&q_sub_reverse[0], q_sub.size(), t_sub_char, t_sub.size(), 2, -1, -2, cigar, target_begin);

    cigar = std::string(cigar.rbegin(), cigar.rend());

    std::cout << cigar << std::endl;

    return 0;
}

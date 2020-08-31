// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <deque>
#include <stdexcept>
#include <fstream>
#include <set>
#include <stack>
#include <string.h>
#include <chrono>
#include "ksw2.h"
#include <iostream>
#include <math.h>

using namespace std;
using namespace std::chrono;

namespace {

static std::uint64_t First(const std::pair<std::uint64_t, std::uint64_t>& pr) {
  return pr.first;
}

static std::uint64_t Second(const std::pair<std::uint64_t, std::uint64_t>& pr) {
  return pr.second;
}

}  // namespace

namespace ram {

#define MAXN 3001
#define STRANDS 2

const std::vector<std::uint64_t> kCoder = {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255};

struct OverlapCluster {
  biosoup::Overlap overlap;
  bool valid;
};

struct AnchorCluster {
  int64_t read_start;
  int64_t read_end;
  int64_t ref_start;
  int64_t ref_end;
  int covered_bases;
  int index;

  AnchorCluster(int64_t _read_start, int64_t _read_end, int64_t _ref_start, int64_t _ref_end, int _covered_bases, int _indexJ) :
    read_start(_read_start), read_end(_read_end), ref_start(_ref_start), ref_end(_ref_end), covered_bases(_covered_bases), index(_indexJ) {}

  friend bool operator <(const AnchorCluster& a, const AnchorCluster& b) {
    if (a.read_start != b.read_start) {
      return a.read_start < b.read_start;
    }
    if (a.ref_start != b.ref_start) {
      return a.ref_start < b.ref_start;
    }
    if (a.read_end != b.read_end) {
      return a.read_end < b.read_end;
    }
    if (a.ref_end != b.ref_end) {
      return a.ref_end < b.ref_end;
    }
    return a.covered_bases < b.covered_bases;
  }
};

// ANCDHOR ALIGNMENT PART

uint8_t seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

template<typename T>
std::vector<T> GenerateSimpleMatchMatrix(T match, T mismatch, size_t alphabet_size) {
  std::vector<T> matrix(alphabet_size * alphabet_size, mismatch);
  for (size_t i=0; i<(alphabet_size - 1); i++) {
    matrix[i*alphabet_size + i] = match;
    matrix[i*alphabet_size + alphabet_size - 1] = 0;
    matrix[(alphabet_size - 1) * alphabet_size + i] = 0;
  }
  return matrix;
}

std::vector<int8_t> ConvertSeqAlphabet(const int8_t* seq, size_t seqlen, const uint8_t* conv_table) {
  std::vector<int8_t> ret(seqlen + 33);
  for (size_t i=0; i<seqlen; i++) {
    ret[i] = (int8_t) conv_table[(uint8_t) seq[i]];
  }
  return ret;
}

void AlignAnchors(const char *_tseq, const char *_qseq, int tlen, int qlen, int sc_mch, int sc_mis, int gapo, int gape, std::vector<std::pair<int, char>> &cigar, bool type = 1) {
  void *km = 0;
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));

  auto mat = GenerateSimpleMatchMatrix<int8_t>((int8_t) 5, (int8_t) -4, 5);

  ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1;
  ez.max = 0;
  ez.mqe = ez.mte = KSW_NEG_INF;
  ez.n_cigar = 0;

  auto qseq = ConvertSeqAlphabet((const int8_t*) _qseq, qlen, &seq_nt4_table[0]);
  auto tseq = ConvertSeqAlphabet((const int8_t*) _tseq, tlen, &seq_nt4_table[0]);
  
  if (type) {
    ksw_exts2_sse(km, qlen, (const uint8_t*) &qseq[0], tlen, (const uint8_t*) &tseq[0], 5, &mat[0], 4, 2, 32, 9, 200, 1600, &ez);
  } else {
    ksw_extd2_sse(km, qlen, (const uint8_t*) &qseq[0], tlen, (const uint8_t*) &tseq[0], 5, &mat[0], 4, 2, 13, 1, -1, -1, 0, &ez);
  }
  
  for (int i = 0; i < ez.n_cigar; ++i) {
    int count = ez.cigar[i]>>4;
    char op = "MIDN"[ez.cigar[i]&0xf];
    if (op == 'N') {
      op = 'D';
    }
    cigar.emplace_back(count, op);
  }
  free(ez.cigar);
}

// ADJUST CIGAR PART

void FindRefOffsets(const char *ref, std::uint32_t reference_size, char first_base, char second_base, char third_base, char fourth_base, std::vector<int> *left_offset, std::vector<int> *right_offset, int64_t start_position, int number_of_bases) {
  int limit = 10;
  for(int i = -limit; i < limit; i++) {
    if(i < 0 || i >= reference_size) {
      continue;
    }
    if(std::toupper(ref[start_position+i]) == first_base && std::toupper(ref[start_position+i+1]) == second_base) {
      left_offset->emplace_back(i);
    }
  }

  for(int i = -limit; i < limit; i++) {
    if(i < 0 || i >= reference_size) {
      continue;
    }
    if(std::toupper(ref[start_position+number_of_bases+i-2]) == third_base && std::toupper(ref[start_position+number_of_bases+i-1]) == fourth_base) {
      right_offset->emplace_back(i);
    }
  }
}

int FindReadLeftOffset(const char *query, int left_offset_ref, int64_t start_position_read, std::stack<std::pair<int, char> > *cigar_stack) {
  int counter = -left_offset_ref;
  int read_offset = 0;

  while (!cigar_stack->empty()) {
    std::pair<int, char> c = cigar_stack->top();
    cigar_stack->pop();
    int count = c.first;

    if(c.second == 'N' || c.second == 'S') {
      cigar_stack->push(c);
      break;
    }

    if(c.second == '=' || c.second == 'X' || c.second == 'M') {
      if(count > counter) {
        read_offset += counter;
        std::pair<int, char> new_c(count - counter, c.second);
        cigar_stack->push(new_c);
        break;
      } else if(count == counter) {
        read_offset += counter;
        break;
      } else {
        read_offset += count;
        counter -= count;
      }
    }

    if(c.second == 'D') {
      if(count > counter) {
        std::pair<int, char> new_c(count - counter, c.second);
        cigar_stack->push(new_c);
        break;
      } else if(count == counter) {
        break;
      } else {
        counter -= count;
      }
    }

    if(c.second == 'I') {
      read_offset += count;
    }
  }

  return read_offset;
}

int FindReadRightOffset(const char *query, int right_offset_ref, int64_t start_position_read, std::deque<std::pair<int, char>> *cigar_queue) {
  int counter = right_offset_ref;
  int read_offset = 0;

  while (!cigar_queue->empty()) {
    std::pair<int, char> c = cigar_queue->front();
    cigar_queue->pop_front();
    int count = c.first;

    if(c.second == 'N' || c.second == 'S') {
      cigar_queue->push_front(c);
      break;
    }

    if(c.second == '=' || c.second == 'X' || c.second == 'M') {
      if(count > counter) {
        read_offset += counter;
        std::pair<int, char> new_c = std::pair<int, char>(count - counter, c.second);
        cigar_queue->push_front(new_c);
        break;
      } else if(count == counter) {
        read_offset += counter;
        break;
      } else {
        read_offset += count;
        counter -= count;
      }
    }

    if(c.second == 'D') {
      if(count > counter) {
        std::pair<int, char> new_c = std::pair<int, char>(count - counter, c.second);
        cigar_queue->push_front(new_c);
        break;
      } else if(count == counter) {
        break;
      } else {
        counter -= count;
      }
    }

    if(c.second == 'I') {
      read_offset += count;
    }
  }

  return read_offset;
}

void AdjustEnds(int left_offset_ref, int right_offset_ref, const char *query, const char *ref, int64_t *start_position_ref, int64_t *start_position_read, int number_of_bases, std::stack<std::pair<int, char>> *cigar_stack, std::deque<std::pair<int, char>> *cigar_queue, bool type) {
  int left_offset_read = left_offset_ref > 0 ? 0 : FindReadLeftOffset(query, left_offset_ref, *start_position_read, cigar_stack);
  int right_offset_read = right_offset_ref < 0 ? 0 : FindReadRightOffset(query, right_offset_ref, *start_position_read + number_of_bases, cigar_queue);

  if (left_offset_ref >= 0 && right_offset_ref <= 0) {
    if(left_offset_ref > 0) {
      std::pair<int, char> c_left = std::pair<int, char>(left_offset_ref, 'D');
      cigar_stack->push(c_left);
    }
    std::pair<int, char> c_gap = std::pair<int, char>(number_of_bases + right_offset_ref + (-left_offset_ref), 'N');
    cigar_stack->push(c_gap);
    if(-right_offset_ref > 0) {
      std::pair<int, char> c_right = std::pair<int, char>(-right_offset_ref, 'D');
      cigar_stack->push(c_right);
    }

    *start_position_ref += number_of_bases;
  } else if (left_offset_ref <= 0 && right_offset_ref >= 0) {
    if(left_offset_read > 0) {
      std::pair<int, char> c_left = std::pair<int, char>(left_offset_read, 'I');
      cigar_stack->push(c_left);
    }
    std::pair<int, char> c_gap = std::pair<int, char>(number_of_bases + right_offset_ref + (-left_offset_ref), 'N');
    cigar_stack->push(c_gap);
    if(right_offset_read > 0) {
      std::pair<int, char> c_right = std::pair<int, char>(right_offset_read, 'I');
      cigar_stack->push(c_right);
    }

    *start_position_ref += (number_of_bases + right_offset_ref);
    *start_position_read += (right_offset_read);
  } else if (left_offset_ref >= 0 && right_offset_ref >= 0) {
    std::string ref_string;
    for(int64_t i = 0; i < left_offset_ref; i++) {
      ref_string.push_back(ref[*start_position_ref + i]);
    }
    std::string read_string;
    for(int64_t i = 0; i < right_offset_read; i++) {
      read_string.push_back(query[*start_position_read + i]);
    }
    const char* ref_sub = ref_string.c_str();
    const char* read_sub = read_string.c_str();

    std::vector<std::pair<int, char>> first_cigar;
    AlignAnchors(ref_sub, read_sub, left_offset_ref, right_offset_read, 5, -2, 4, 2, first_cigar, 0);

    for (int k = 0; k < first_cigar.size(); k++) {
      cigar_stack->push(first_cigar[k]);
    }

    std::pair<int, char> c_gap = std::pair<int, char>(number_of_bases + right_offset_ref + (-left_offset_ref), 'N');
    cigar_stack->push(c_gap);

    *start_position_ref += (number_of_bases + right_offset_ref);
    *start_position_read += (right_offset_read);
  } else if (left_offset_ref <= 0 && right_offset_ref <= 0) {
    std::string ref_string;
    for(int64_t i = right_offset_ref; i < 0; i++) {
      ref_string.push_back(ref[*start_position_ref + number_of_bases - i]);
    }
    std::string read_string;
    for(int64_t i = -left_offset_read; i < 0; i++) {
      read_string.push_back(query[*start_position_read + left_offset_read]);
    }
    const char* ref_sub = ref_string.c_str();
    const char* read_sub = read_string.c_str();

    std::vector<std::pair<int, char>> first_cigar;
    AlignAnchors(ref_sub, read_sub, -right_offset_ref, left_offset_read, 5, -2, 4, 2, first_cigar, 0);

    std::pair<int, char> c_gap = std::pair<int, char>(number_of_bases + right_offset_ref + (-left_offset_ref), 'N');
    cigar_stack->push(c_gap);

    for (int k = 0; k < first_cigar.size(); k++) {
      cigar_stack->push(first_cigar[k]);
    }

    *start_position_ref += (number_of_bases);
  }
}

double AlignEdges(const char *query, const char *ref, int leftRef, int rightRef, int64_t start_position_read, int64_t start_position_ref, int number_of_bases,  std::stack<std::pair<int, char>> cigar_stack, std::deque<std::pair<int, char>> cigar_queue) {
  std::string read_container;
  std::string ref_container;

  int counter = leftRef;

  while(counter > 0) {
    ref_container.insert(0, 1, ref[start_position_ref+counter-1]);
    counter -= 1;
  }

  int offset = std::min(0, leftRef);
  counter = offset;
  int window = 10;

  while(counter > offset - window) {
    ref_container.insert(0, 1, ref[start_position_ref+counter-1]);
    counter -= 1;
  }

  counter = -offset + window;
  int read_offset = 0;

  while (!cigar_stack.empty()) {
    std::pair<int, char> c = cigar_stack.top();
    cigar_stack.pop();
    int count = c.first;

    if(c.second == 'N' || c.second == 'S') {
      break;
    }

    if(c.second == '=' || c.second == 'X' || c.second == 'M') {
      if(count > counter) {
        read_offset += counter;
        std::pair<int, char> new_c(count - counter, c.second);
        break;
      } else if(count == counter) {
        read_offset += counter;
        break;
      } else {
        read_offset += count;
        counter -= count;
      }
    }

    if(c.second == 'D') {
      if(count > counter) {
        std::pair<int, char> new_c(count - counter, c.second);
        break;
      } else if(count == counter) {
        break;
      } else {
        counter -= count;
      }
    }

    if(c.second == 'I') {
      read_offset += count;
    }
  }

  counter = read_offset;

  while(counter > 0) {
    read_container.push_back(query[start_position_read-counter]);
    counter -= 1;
  }

  std::string read_container_right;
  std::string ref_container_right;

  counter = rightRef;

  while(counter < 0) {
    ref_container_right.push_back(ref[start_position_ref+number_of_bases+counter]);
    counter += 1;
  }

  int offset_right = std::max(0, rightRef);
  counter = offset_right;

  while(counter < offset_right + window) {
    ref_container_right.push_back(ref[start_position_ref+number_of_bases+counter]);
    counter += 1;
  }

  counter = offset_right + window;
  int read_offset_right = 0;

  while (!cigar_queue.empty()) {
    std::pair<int, char> c = cigar_queue.front();
    cigar_queue.pop_front();
            
    int count = c.first;

    if(c.second == 'N' || c.second == 'S') {
      cigar_queue.push_front(c);
      break;
    }

    if(c.second == '=' || c.second == 'X' || c.second == 'M') {
      if(count > counter) {
        read_offset_right += counter;
        break;
      } else if(count == counter) {
        read_offset_right += counter;
        break;
      } else {
        read_offset_right += count;
        counter -= count;
      }
    }

    if(c.second == 'D') {
      if(count > counter) {
        std::pair<int, char> new_c(count - counter, c.second);
        cigar_queue.push_front(new_c);
        break;
      } else if(count == counter) {
        break;
      } else {
        counter -= count;
      }
    }

    if(c.second == 'I') {
      read_offset_right += count;
    }
  }

  counter = 0;

  while(counter < read_offset_right) {
    read_container_right.push_back(query[start_position_read+counter]);
    counter += 1;
  }

  std::string ref_cut = ref_container + ref_container_right;
  std::string read_cut = read_container + read_container_right;

  if(std::abs((long int)(read_cut.size() - ref_cut.size())) > (long int) 60) {
    return -1;
  }

  std::vector<std::pair<int, char>> result_cigar;
  AlignAnchors(ref_cut.c_str(), read_cut.c_str(), ref_cut.size(), read_cut.size(), 5, -2, 4, 2, result_cigar, 0);

  int matches = 0;
  int length = 0;
  
  int ref_ln = 0;
  int read_ln = 0;

  for (int k = 0; k < result_cigar.size(); k++) {
    if (result_cigar[k].second == 'M' || result_cigar[k].second == '=' || result_cigar[k].second == 'X') {
      for (int l = 0; l < result_cigar[k].first; l++) {
        if(ref_cut[ref_ln + l] == read_cut[read_ln + l]) {
          matches += 1;
        }
      }
      ref_ln += result_cigar[k].first;
      read_ln += result_cigar[k].first;
    } else if(result_cigar[k].second == 'I') {
      read_ln += result_cigar[k].first;
    } else if(result_cigar[k].second == 'D') {
      ref_ln += result_cigar[k].first;
    }
    length += result_cigar[k].first;
  }

  return (double) matches / (double) length;
}

void ExtractCigarString(std::vector<std::pair<int, char>> &cigar_string, std::vector<std::pair<int, char>> &result_cigar, bool should_cut_gap, int expected_left_ref_len, int cutted_gap_len, int current_desired_read_len, int &total_read_len, int &total_ref_len) {
  for (int k = 0; k < cigar_string.size(); k++) {
    int count = cigar_string[k].first;
    char elem = cigar_string[k].second;
            
    switch (elem) {
      case 'D':
        total_ref_len += count;
        if(should_cut_gap && total_ref_len >= expected_left_ref_len) {
          result_cigar.emplace_back(count + cutted_gap_len, elem);
          total_ref_len += cutted_gap_len;
          should_cut_gap = false;
        } else {
          result_cigar.emplace_back(count, elem);
        }
        break;
      case 'I':
        if(count + total_read_len < current_desired_read_len) {
          total_read_len += count;
          result_cigar.emplace_back(count, elem);
        } else {
          int value = current_desired_read_len - total_read_len;
          result_cigar.emplace_back(value, elem);
          total_read_len += value;
          return;
        }
        break;
      case 'M':
        if (count + total_read_len < current_desired_read_len) {
          total_read_len += count;
          total_ref_len += count;
          result_cigar.emplace_back(count, elem);
        } else {
          int value = current_desired_read_len - total_read_len;
          total_read_len += value;
          total_ref_len += value;
          result_cigar.emplace_back(value, elem);
          return;
        }
        break;
      default:
        break;
    }
  }
}

std::vector<std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>> ReverseClusters(std::vector<std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>> &clusters, std::string &read) {
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>> reversed;
    for (int i = clusters.size()-1; i >= 0; i--) {
        std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> current = clusters[i];
        
        std::pair<uint64_t, uint64_t> first_pair = current.first;
        std::pair<uint64_t, uint64_t> second_pair = current.second;
                
        std::pair<uint64_t, uint64_t> first_pair_rev(read.size() - first_pair.second, read.size() - first_pair.first);
        std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> rev(first_pair_rev, second_pair);
        reversed.push_back(rev);        
    }    
    return reversed;
}

std::string GenerateSamFromAnchorClusters(int ref_index, int strand, std::vector<std::unique_ptr<biosoup::Sequence>> &references, const std::unique_ptr<biosoup::Sequence>& read, std::vector<std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>> cluster) {
  if(ref_index < 0) {
    return "";
  }

  if(cluster.size() == 0) {
    return "";
  }
     
  if (strand == 0) {
    read->ReverseAndComplement();
    cluster = ReverseClusters(cluster, read->data);
  }
                            
  std::string strand_string = "0";
  if (strand == 0) {
    strand_string = "16";
  }

  uint64_t start_position_of_the_reference = (cluster.front().second.first + 1);
  std::string output_result = read->name + "\t" + strand_string + "\t" + references[ref_index]->name + "\t" + std::to_string(start_position_of_the_reference) + "\t0\t";
        
  uint64_t front_clipping = cluster[0].first.first;
  uint64_t back_clipping = 0;
  
  if(cluster.back().first.second < read->data.size()) {
    back_clipping = read->data.size() - cluster.back().first.second;
  }
      
  uint64_t reference_current_position = cluster.front().second.first;
      
  std::vector<std::pair<int, char>> result_cigar;
      
  if (front_clipping > 0) {
    result_cigar.emplace_back(front_clipping, 'S');
  }
              
  for(uint32_t j = 0; j < cluster.size()-1; j++) {
    std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> anchor_left = cluster[j];
    std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> anchor_right = cluster[j+1];
          
    uint64_t read_gap_len = 0;
    if(anchor_right.first.first > anchor_left.first.second) {
      // TODO must check ends of clusters;
      read_gap_len = anchor_right.first.first - anchor_left.first.second;
    }
    
    uint64_t ref_gap_len = 0;
    if(anchor_right.second.first > anchor_left.second.second) {
      ref_gap_len = anchor_right.second.first - anchor_left.second.second;
    }
          
    bool should_cut_gap = ref_gap_len > read_gap_len * 2 && ref_gap_len > 200;
          
    int cutted_gap_len = 0;
    int expected_left_ref_len = 0;

    uint64_t buffer = read_gap_len;
    if(should_cut_gap) {
      if (read_gap_len < 15) {
        buffer += 15;
      }
      cutted_gap_len = (anchor_right.second.first - buffer) - (anchor_left.second.second + buffer);
      expected_left_ref_len = (anchor_left.second.second + buffer) - reference_current_position;
    }
          
    std::string read_seq = read->data.substr(anchor_left.first.first, (anchor_right.first.second - anchor_left.first.first));
    std::string ref_seq = should_cut_gap ?
              ((references[ref_index]->data.substr(reference_current_position, ((anchor_left.second.second + buffer) - reference_current_position)))  +  (references[ref_index]->data.substr(anchor_right.second.first - buffer, (anchor_right.second.second - (anchor_right.second.first - buffer)))))
              : references[ref_index]->data.substr(reference_current_position, (anchor_right.second.second - reference_current_position));
    
    std::vector<std::pair<int, char>> cigar_string;

    if(read_seq.size() > 20000 || ref_seq.size() > 20000) {
      // TODO: fix this where gap is not included in result
      return "";
    }

    AlignAnchors(ref_seq.c_str(), read_seq.c_str(), ref_seq.size(), read_seq.size(), 5, -2, 4, 2, cigar_string);
    int current_desired_read_len = anchor_right.first.first - anchor_left.first.first;
    std::vector<std::pair<int, char>> tmp_cigar;

    int ref_len = 0;
    int read_len = 0;
    ExtractCigarString(cigar_string, result_cigar, should_cut_gap, expected_left_ref_len, cutted_gap_len, current_desired_read_len, read_len, ref_len);
    reference_current_position += ref_len;
  }
      
  std::string read_seq = read->data.substr(cluster.back().first.first, (cluster.back().first.second - cluster.back().first.first));
  std::string ref_seq = references[ref_index]->data.substr(reference_current_position, (cluster.back().second.second - cluster.back().second.first));
  std::vector<std::pair<int, char>> cigar_string;
  AlignAnchors(ref_seq.c_str(), read_seq.c_str(), ref_seq.size(), read_seq.size(), 5, -2, 4, 2, cigar_string);


  if(read_seq.size() > 20000 || ref_seq.size() > 20000) {
    return "";
  }

  int current_desired_read_len = cluster.back().first.second - cluster.back().first.first;
      
  int ref_len = 0;
  int read_len = 0;
  ExtractCigarString(cigar_string, result_cigar, false, 0, 0, current_desired_read_len, read_len, ref_len);

  if (cluster.back().first.second < read->data.size()) {
    result_cigar.emplace_back(read->data.size() - cluster.back().first.second, 'S');
  }
      
  int previous_count = 0;
  char previous_elem = 'x';
  std::string result = "";
  std::vector<std::pair<int, char>> result_cigar2;

  for (int k = 0; k < result_cigar.size(); k++) {
    int count = result_cigar[k].first;
    char elem = result_cigar[k].second;
                
    if (count == 0) {
      continue;
    }
    
    if (elem == 'D' && count > 10) {
      elem = 'N';
    }
    
    if (previous_elem == elem) {
      previous_count += count;
    } else {
      if(previous_elem != 'x') {
        result_cigar2.emplace_back(previous_count, previous_elem);
      }
      previous_count = count;
      previous_elem = elem;
    }
  }
      
  result_cigar2.emplace_back(previous_count, previous_elem);
      
  int64_t start_position_ref = cluster.front().second.first;
  int64_t start_position_read = cluster.front().first.first;

  std::stack<std::pair<int, char>> cigar_stack;
  std::deque<std::pair<int, char>> cigar_queue;

  for(int k = 0; k < result_cigar2.size(); k++) {
    cigar_queue.push_back(result_cigar2[k]);
  }

  while (!cigar_queue.empty()) {
    std::pair<int, char> cigar_op = cigar_queue.front();
    cigar_queue.pop_front();
    int64_t number_of_bases = cigar_op.first;

    if (cigar_op.second == 'S') {
      cigar_stack.push(cigar_op);
      continue;
    }

    if (cigar_op.second == 'N') {
      std::vector<int> leftRef;
      std::vector<int> rightRef;
      FindRefOffsets(references[ref_index]->data.c_str(), references[ref_index]->data.size(), 'G', 'T', 'A', 'G', &leftRef, &rightRef, start_position_ref, number_of_bases);

      double rezult = -1;
      int rezL = 0;
      int rezR = 0;

      for(int left: leftRef) {
        for(int right: rightRef) {
          double rez = AlignEdges(read->data.c_str(), references[ref_index]->data.c_str(), left, right, start_position_read, start_position_ref, number_of_bases, cigar_stack, cigar_queue);
          if(rez > rezult) {
            rezult = rez;
            rezL = left;
            rezR = right;
          }
        }
      }
      
      std::vector<int> leftRef2;
      std::vector<int> rightRef2;

      FindRefOffsets(references[ref_index]->data.c_str(), references[ref_index]->data.size(), 'C', 'T', 'A', 'C', &leftRef2, &rightRef2, start_position_ref, number_of_bases);

      for(int left: leftRef2) {
        for(int right: rightRef2) {
          double rez = AlignEdges(read->data.c_str(), references[ref_index]->data.c_str(), left, right, start_position_read, start_position_ref, number_of_bases, cigar_stack, cigar_queue);
          if(rez > rezult) {
            rezult = rez;
            rezL = left;
            rezR = right;
          }
        }
      }

      if(rezL != 0 || rezR != 0) {
        AdjustEnds(rezL, rezR, read->data.c_str(), references[ref_index]->data.c_str(), &start_position_ref, &start_position_read, number_of_bases, &cigar_stack, &cigar_queue, 1);
      } else {
        cigar_stack.push(cigar_op);

        if(cigar_op.second != 'I') {
          start_position_ref += cigar_op.first;
        }
        if(cigar_op.second != 'D' && cigar_op.second != 'N') {
          start_position_read += cigar_op.first;
        }
      }
    } else {
      cigar_stack.push(cigar_op);
      if(cigar_op.second != 'I') {
        start_position_ref += cigar_op.first;
      }
      if(cigar_op.second != 'D' && cigar_op.second != 'N') {
        start_position_read += cigar_op.first;
      }
    }
  }
  
  std::stack<std::pair<int, char>> tmp_stack;
  std::pair<int, char> previous_op = cigar_stack.top();
  cigar_stack.pop();

  while(!cigar_stack.empty()) {
    std::pair<int, char> tmp_op = cigar_stack.top();
    cigar_stack.pop();
    if(tmp_op.second == previous_op.second) {
      previous_op = std::pair<int, char>(previous_op.first + tmp_op.first, previous_op.second);
    } else {
      tmp_stack.push(previous_op);
      previous_op = tmp_op;
    }
  }
  tmp_stack.push(previous_op);
  result_cigar2.clear();

  while(!tmp_stack.empty()) {
    std::pair<int, char> c = tmp_stack.top();
    tmp_stack.pop();
    result_cigar2.emplace_back(c);
  }
      
  for (int k = 0; k < result_cigar2.size(); k++) {
    result.append(std::to_string(result_cigar2[k].first));
    result.push_back(result_cigar2[k].second);
  }

  output_result += (result + "\t*\t0\t0\t" + read->data + "\t*\n");
  return output_result;
}


// KNAPSACK PART

void CalculateDP(int **dp, int **backtrack, std::vector<AnchorCluster> *clusters) {
  for (int strand = 0; strand < STRANDS; strand++) {
    if (clusters[strand].size() == 0) {
      continue;
    }
    for (int n = clusters[strand].size(), x = n; x >= 0; x--) {
      int index = x - 1;
      int64_t xLeft = (index == -1) ? 0 : clusters[strand][index].read_end;
      int64_t yLeft = (index == -1) ? 0 : clusters[strand][index].ref_end;
      int &ref = dp[strand][x];
      for (int i = x; i < n; i++) {
        if (clusters[strand][i].read_start < xLeft || clusters[strand][i].ref_start < yLeft) {
          continue;
        }
        int tmp = clusters[strand][i].covered_bases + dp[strand][i + 1];
        if (tmp > ref) {
          backtrack[strand][x] = i + 1;
          ref = tmp;
        }
      }
    }
  }
}

int RNAFilterClusters(std::vector<AnchorCluster> *clusters, std::vector<OverlapCluster> &cluster_data, int resulting_strand) {

  int **dp = (int **) calloc(STRANDS, sizeof(int));
  int **backtrack = (int **) calloc(STRANDS, sizeof(int));

  for (int strand = 0; strand < STRANDS; strand++) {
    dp[strand] = (int *) calloc(MAXN, sizeof(int));
    backtrack[strand] = (int *) calloc(MAXN, sizeof(int));
    for (int i = 0; i < MAXN; i++) {
      backtrack[strand][i] = -1;
    }
  }

  CalculateDP(dp, backtrack, clusters);
  int strand = !resulting_strand;

  std::set<int> done;

  int start = 0;
  int curr = backtrack[strand][start];

  while (true) {
    if (curr == -1) {
      break;
    }
    if (done.count(curr) > 0) {
      return 0;
    }
    done.insert(curr);
    int currClusterIndex = curr - 1;
    if (curr == backtrack[strand][curr]) {
      return 0;
    }

    int i = clusters[strand][currClusterIndex].index;

    cluster_data[i].valid = true;
    if (backtrack[strand][curr] == -1) {
      break;
    }
    curr = backtrack[strand][curr];
  }
  return 0;
}


// MINIMIZER ENGINE PART

MinimizerEngine::MinimizerEngine(
    std::uint32_t kmer_len,
    std::uint32_t window_len,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::uint64_t wildcard_mask)
    : k_(std::min(std::max(kmer_len, 1U), 32U)),
      w_(window_len),
      thread_pool_(thread_pool ?
          thread_pool :
          std::make_shared<thread_pool::ThreadPool>(1)),
      wildcard_mask_(wildcard_mask),
      occurrence_(-1),
      minimizers_(1U << std::min(14U, 2 * k_)),
      index_(minimizers_.size()) {}

void MinimizerEngine::Minimize(
    std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator begin,
    std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator end) {

  for (auto& it : minimizers_) {
    it.clear();
  }
  for (auto& it : index_) {
    it.clear();
  }

  if (begin >= end) {
    return;
  }

  {
    std::uint64_t bin_mask = minimizers_.size() - 1;

    std::vector<std::future<std::vector<uint128_t>>> futures;
    for (auto it = begin; it != end; ++it) {
      futures.emplace_back(thread_pool_->Submit(
          [&] (std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator it)  // NOLINT
              -> std::vector<uint128_t> {
            return Minimize(*it);
          },
          it));
    }
    for (auto& it : futures) {
      for (const auto& jt : it.get()) {
        minimizers_[jt.first & bin_mask].emplace_back(jt);
      }
    }
  }

  {
    std::vector<std::future<void>> futures;
    for (std::uint32_t i = 0; i < minimizers_.size(); ++i) {
      if (minimizers_[i].empty()) {
        continue;
      }

      futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t bin) -> void{
            RadixSort(
                minimizers_[bin].begin(),
                minimizers_[bin].end(),
                k_ * 2,
                ::First);

            for (std::uint64_t i = 0, c = 0; i < minimizers_[bin].size(); ++i) {
              if (i > 0 && minimizers_[bin][i - 1].first != minimizers_[bin][i].first) {  // NOLINT
                index_[bin].emplace(
                    minimizers_[bin][i - 1].first,
                    std::make_pair(i - c, c));
                c = 0;
              }
              if (i == minimizers_[bin].size() - 1) {
                index_[bin].emplace(
                    minimizers_[bin][i].first,
                    std::make_pair(i - c, c + 1));
              }
              ++c;
            }
          },
          i));
    }
    for (const auto& it : futures) {
      it.wait();
    }
  }
}

void MinimizerEngine::Filter(double frequency) {
  if (!(0 <= frequency && frequency <= 1)) {
    throw std::invalid_argument(
        "[ram::MinimizerEngine::Filter] error: invalid frequency");
  }

  if (frequency == 0) {
    occurrence_ = -1;
    return;
  }

  std::vector<std::uint32_t> occurrences;
  for (const auto& it : index_) {
    for (const auto& jt : it) {
      occurrences.emplace_back(jt.second.second);
    }
  }

  if (occurrences.empty()) {
    occurrence_ = -1;
    return;
  }

  std::nth_element(
      occurrences.begin(),
      occurrences.begin() + (1 - frequency) * occurrences.size(),
      occurrences.end());
  occurrence_ = occurrences[(1 - frequency) * occurrences.size()] + 1;
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::Sequence>& sequence,
    bool avoid_equal,
    bool avoid_symmetric,
    bool micromize) const {

  auto sketch = Minimize(sequence, micromize);
  if (sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  std::uint64_t bin_mask = minimizers_.size() - 1;
  std::vector<uint128_t> matches;
  for (const auto& it : sketch) {
    std::uint32_t bin = it.first & bin_mask;
    auto match = index_[bin].find(it.first);
    if (match == index_[bin].end() || match->second.second > occurrence_) {
      continue;
    }

    auto jt = minimizers_[bin].begin() + match->second.first;
    auto end = jt + match->second.second;
    for (; jt != end; ++jt) {
      std::uint64_t rhs_id = jt->second >> 32;
      if (avoid_equal && sequence->id == rhs_id) {
        continue;
      }
      if (avoid_symmetric && sequence->id > rhs_id) {
        continue;
      }

      std::uint64_t strand = (it.second & 1) == (jt->second & 1);
      std::uint64_t lhs_begin = it.second << 32 >> 33;
      std::uint64_t rhs_begin = jt->second << 32 >> 33;

      std::uint64_t diff = !strand ?
          rhs_begin + lhs_begin :
          rhs_begin - lhs_begin + (3ULL << 30);  // TODO(rvaser): check this

      matches.emplace_back(
          (((rhs_id << 1) | strand) << 32) | diff,
          (lhs_begin << 32) | rhs_begin);
    }
  }
  return Chain(sequence->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::Sequence>& lhs,
    const std::unique_ptr<biosoup::Sequence>& rhs,
    bool micromize) const {

  auto lhs_sketch = Minimize(lhs, micromize);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs);
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(lhs_sketch.begin(), lhs_sketch.end(), k_ * 2, ::First);
  RadixSort(rhs_sketch.begin(), rhs_sketch.end(), k_ * 2, ::First);

  std::vector<uint128_t> matches;
  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].first < rhs_sketch[j].first) {
        break;
      } else if (lhs_sketch[i].first == rhs_sketch[j].first) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].first != rhs_sketch[k].first) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].second & 1) == (rhs_sketch[k].second & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].second << 32 >> 33;
          std::uint64_t rhs_pos = rhs_sketch[k].second << 32 >> 33;
          std::uint64_t diagonal = !strand ?
              rhs_pos + lhs_pos :
              rhs_pos - lhs_pos + (3ULL << 30);

          matches.emplace_back(
              (((rhs->id << 1) | strand) << 32) | diagonal,
              (lhs_pos << 32) | rhs_pos);
        }
        break;
      } else {
        ++j;
      }
    }
  }

  return Chain(lhs->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id,
    std::vector<uint128_t>&& matches) const {

  RadixSort(matches.begin(), matches.end(), 64, ::First);
  matches.emplace_back(-1, -1);  // stop dummy

  std::vector<uint128_t> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].first - matches[j].first > 500) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].first - matches[j].first > 500) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < 4) {
      continue;
    }

    RadixSort(matches.begin() + j, matches.begin() + i, 64, ::Second);

    std::uint64_t strand = matches[j].first >> 32 & 1;

    std::vector<std::uint64_t> indices;
    if (strand) {  // same strand
      indices = LongestSubsequence(  // increasing
          matches.begin() + j,
          matches.begin() + i,
          std::less<std::uint64_t>());
    } else {  // different strand
      indices = LongestSubsequence(  // decreasing
          matches.begin() + j,
          matches.begin() + i,
          std::greater<std::uint64_t>());
    }

    if (indices.size() < 4) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if ((matches[j + indices[k]].second >> 32) -
          (matches[j + indices[k - 1]].second >> 32) > 10000ULL) {
        if (k - l < 4) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches[j + indices[m]].second >> 32;
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[j + indices[m]].second << 32 >> 32;
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < 100UL) {
          l = k;
          continue;
        }

        dst.emplace_back(
            lhs_id << 1 | strand,
            matches[j + indices[l]].second >> 32,  // lhs_begin
            k_ + (matches[j + indices[k - 1]].second >> 32),  // lhs_end
            matches[j].first >> 33,  // rhs_id
            strand ?  // rhs_begin
                matches[j + indices[l]].second << 32 >> 32 :
                matches[j + indices[k - 1]].second << 32 >> 32,
            k_ + (strand ?  // rhs_end
                matches[j + indices[k - 1]].second << 32 >> 32 :
                matches[j + indices[l]].second << 32 >> 32),
            std::min(lhs_matches, rhs_matches));  // score

        l = k;
      }
    }
  }
  return dst;
}

std::vector<MinimizerEngine::uint128_t> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::Sequence>& sequence,
    bool micromize) const {

  if (sequence->data.size() < k_) {
    return std::vector<uint128_t>{};
  }

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;

  auto hash = [&] (std::uint64_t key) -> std::uint64_t {
    key = ((~key) + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
  };

  std::deque<uint128_t> window;
  auto window_add = [&] (std::uint64_t minimizer, std::uint64_t location) -> void {  // NOLINT
    while (!window.empty() && window.back().first > minimizer) {
      window.pop_back();
    }
    window.emplace_back(minimizer, location);
  };
  auto window_update = [&] (std::uint32_t position) -> void {
    while (!window.empty() && (window.front().second << 32 >> 33) < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = sequence->id << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<uint128_t> dst;

  for (std::uint32_t i = 0; i < sequence->data.size(); ++i) {
    std::uint64_t c = kCoder[sequence->data[i]];
    if (c == 255ULL) {
      throw std::invalid_argument(
          "[ram::MinimizerEngine::Minimize] error: invalid character");
    }
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= k_ - 1U) {
      if ((minimizer & wildcard_mask_) < (reverse_minimizer & wildcard_mask_)) {
        window_add(hash(minimizer & wildcard_mask_), (i - (k_ - 1U)) << 1 | 0);
      } else if ((minimizer & wildcard_mask_) > (reverse_minimizer & wildcard_mask_)) {
        window_add(hash(reverse_minimizer & wildcard_mask_), (i - (k_ - 1U)) << 1 | 1);
      }
    }
    if (i >= (k_ - 1U) + (w_ - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->first != window.front().first) {
          break;
        }
        if (it->second & is_stored) {
          continue;
        }
        dst.emplace_back(it->first, id | it->second);
        it->second |= is_stored;
      }
      window_update(i - (k_ - 1U) - (w_ - 1U) + 1);
    }
  }

  if (micromize) {
    RadixSort(dst.begin(), dst.end(), k_ * 2, ::First);
    dst.resize(sequence->data.size() / k_);
  }

  return dst;
}

template<typename T>
void MinimizerEngine::RadixSort(
    std::vector<uint128_t>::iterator begin,
    std::vector<uint128_t>::iterator end,
    std::uint8_t max_bits,
    T compare) {  //  unary comparison function

  if (begin >= end) {
    return;
  }

  std::vector<MinimizerEngine::uint128_t> dst(end - begin);
  auto dst_begin = dst.begin();
  auto dst_end = dst.end();

  std::uint64_t buckets[0x100]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};
    for (auto it = begin; it != end; ++it) {
      ++counts[compare(*it) >> shift & 0xFF];
    }
    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it = begin; it != end; ++it) {
      *(dst_begin + buckets[compare(*it) >> shift & 0xFF]++) = *it;
    }
    std::swap(dst_begin, begin);
    std::swap(dst_end, end);
  }

  if (shift / 8 & 1) {  // copy the sorted array for odd cases
    for (; begin != end; ++begin, ++dst_begin) {
      *dst_begin = *begin;
    }
  }
}

template<typename T>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<uint128_t>::const_iterator begin,
    std::vector<uint128_t>::const_iterator end,
    T compare) {  // binary comparison function

  if (begin >= end) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(end - begin + 1, 0);
  std::vector<std::uint64_t> predecessor(end - begin, 0);

  std::uint64_t longest = 0;
  for (auto it = begin; it != end; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if (((begin + minimal[mid])->second >> 32) < (it->second >> 32) &&
          compare(
              (begin + minimal[mid])->second << 32 >> 32,
              it->second << 32 >> 32)) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - begin] = minimal[lo - 1];
    minimal[lo] = it - begin;
    longest = std::max(longest, lo);
  }

  std::vector<std::uint64_t> dst;
  for (std::uint64_t i = 0, j = minimal[longest]; i < longest; ++i) {
    dst.emplace_back(j);
    j = predecessor[j];
  }
  std::reverse(dst.begin(), dst.end());

  return dst;
}

// RNA MAPPING PART

bool CompareOverlaps(biosoup::Overlap o1, biosoup::Overlap o2) { 
    if(o1.rhs_id == o2.rhs_id) {
      return (o1.rhs_begin < o2.rhs_begin); 
    } else {
      return o1.rhs_id < o2.rhs_id;
    }
} 

bool CompareOverlaps2(biosoup::Overlap o1, biosoup::Overlap o2) { 
  return o1.lhs_begin < o2.lhs_begin;
}

std::vector<biosoup::Overlap> MinimizerEngine::ChainRNA(std::uint64_t lhs_id, std::vector<uint128_t> matches) const {
  RadixSort(matches.begin(), matches.end(), 64, ::First);
  matches.emplace_back(-1, -1);

  for(int i = 0; i < matches.size(); i++) {
    std::uint32_t lhs_pos = matches[i].second >> 32;
    std::uint32_t rhs_pos = matches[i].second << 32 >> 32;
    std::uint32_t diagonal = matches[i].first << 32 >> 32;
    std::uint64_t strand = matches[i].first >> 32 & 1;
  }

  std::vector<uint128_t> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].first - matches[j].first > 500) {
      if (i - j >= 2) {
        if (!intervals.empty() && intervals.back().second > j) {
          intervals.back().second = i;
        } else {
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].first - matches[j].first > 500) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;

  for(int i = 0; i < intervals.size(); i++) {

    uint64_t first = intervals[i].first;
    uint64_t second = intervals[i].second;

    RadixSort(matches.begin() + first, matches.begin() + second, 64, ::Second);

    std::uint64_t strand = matches[first].first >> 32 & 1;

    std::vector<std::uint64_t> indices;
    if (strand) {
      indices = LongestSubsequence(
          matches.begin() + first,
          matches.begin() + second,
          std::less<std::uint64_t>());
    } else {
      indices = LongestSubsequence(
          matches.begin() + first,
          matches.begin() + second,
          std::greater<std::uint64_t>());
    }

    indices.emplace_back(matches.size() - 1 - first);

    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      uint64_t condition1 = (matches[first + indices[k]].second << 32 >> 32) - (matches[first + indices[k - 1]].second << 32 >> 32);
      uint64_t condition2 = (matches[first + indices[k - 1]].second << 32 >> 32) - (matches[first + indices[k]].second << 32 >> 32);
      if ((condition1 > 30ULL && strand == 1) || (condition2 > 30ULL && strand == 0)) {
        if (k - l < 2) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches[first + indices[m]].second >> 32;
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[first + indices[m]].second << 32 >> 32;
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < 15UL) {
          l = k;
          continue;
        }

        dst.emplace_back(
            lhs_id << 1 | strand,
            matches[first + indices[l]].second >> 32,  // lhs_begin
            k_ + (matches[first + indices[k - 1]].second >> 32),  // lhs_end
            matches[first].first >> 33,  // rhs_id
            strand ?  // rhs_begin
                matches[first + indices[l]].second << 32 >> 32 :
                matches[first + indices[k - 1]].second << 32 >> 32,
            k_ + (strand ?  // rhs_end
                matches[first + indices[k - 1]].second << 32 >> 32 :
                matches[first + indices[l]].second << 32 >> 32),
            std::min(lhs_matches, rhs_matches));  // score

        l = k;
      }
    }
  }

  if (dst.size() == 0) {
    return dst;
  }

  sort(dst.begin(), dst.end(), CompareOverlaps); 

  std::vector<std::vector<biosoup::Overlap>> vector_of_overlaps;
  std::vector<biosoup::Overlap> vector_o;
  vector_o.push_back(dst[0]);
  biosoup::Overlap previous = dst[0];

  for(int i = 1; i < dst.size(); i++) {
    biosoup::Overlap current = dst[i];

    if(previous.rhs_begin + 60000 < current.rhs_begin || previous.rhs_id != current.rhs_id) {
      vector_of_overlaps.push_back(vector_o);
      vector_o.clear();
      vector_o.push_back(current);
    } else {
      vector_o.push_back(current);
    }
    previous = current;
  }

  vector_of_overlaps.push_back(vector_o);

  int max_size = 0;
  int max_index = 0;

  for(int i = 0; i < vector_of_overlaps.size(); i++) {
    int size = 0;
    std::vector<biosoup::Overlap> tmp = vector_of_overlaps[i];

    for(int j = 0; j < tmp.size(); j++) {
      size += (tmp[j].rhs_end - tmp[j].rhs_begin);
    }

    if(size > max_size) {
      max_size = size;
      max_index = i;
    }
  }

  std::vector<biosoup::Overlap> best = vector_of_overlaps[max_index];

  int this_size = 0;
  for(int j = 0; j < best.size(); j++) {
    this_size += best[j].lhs_end - best[j].lhs_begin;
  }

  return best;
}

std::string MinimizerEngine::MapRNA(
    const std::unique_ptr<biosoup::Sequence>& sequence,
    bool avoid_equal,
    bool avoid_symmetric,
    bool should_print,
    std::vector<std::unique_ptr<biosoup::Sequence>> &targets,
    bool micromize) const {

  int calculated_rhs_id = -1;
  auto sketch = Minimize(sequence, micromize);
  
  if (sketch.empty()) {
    return "";
  }

  std::uint64_t bin_mask = minimizers_.size() - 1;
  std::vector<uint128_t> matches;
  for (const auto& it : sketch) {
    std::uint32_t bin = it.first & bin_mask;
    auto match = index_[bin].find(it.first);
    if (match == index_[bin].end() || match->second.second > occurrence_) {
      continue;
    }

    auto jt = minimizers_[bin].begin() + match->second.first;
    auto end = jt + match->second.second;
    for (; jt != end; ++jt) {
      std::uint64_t rhs_id = jt->second >> 32;
      if (avoid_equal && sequence->id == rhs_id) {
        continue;
      }
      if (avoid_symmetric && sequence->id > rhs_id) {
        continue;
      }

      std::uint64_t strand = (it.second & 1) == (jt->second & 1);
      std::uint64_t lhs_begin = it.second << 32 >> 33;
      std::uint64_t rhs_begin = jt->second << 32 >> 33;

      std::uint64_t diff = !strand ?
          rhs_begin + lhs_begin :
          rhs_begin - lhs_begin + (3ULL << 30);  // TODO(rvaser): check this

      matches.emplace_back(
          (((rhs_id << 1) | strand) << 32) | diff,
          (lhs_begin << 32) | rhs_begin);
    }
  }

  std::vector<biosoup::Overlap> best = ChainRNA(sequence->id, std::move(matches));

  if(best.size() == 0) {
    return "";
  }

  calculated_rhs_id = best.front().rhs_id;
  for (int l = 1; l < best.size(); ++l) {
    if(calculated_rhs_id != best[l].rhs_id) {
      return "";
    }
  }

  std::uint32_t left_end = best.front().rhs_begin;
  std::uint32_t right_end = best.back().rhs_end;

  int window = 10000;
  uint64_t left_edge = std::max((int32_t) 0, (int32_t) left_end - window);
  uint64_t right_edge = std::min((uint64_t) targets[calculated_rhs_id]->data.size() , (uint64_t) right_end + window);
  std::string ref_substr = targets[calculated_rhs_id]->data.substr(left_edge, right_edge - left_edge);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(1);
  std::uint64_t wildcard_mask = -1;
  ram::MinimizerEngine minimizer_engine{9, 4, thread_pool, wildcard_mask};

  std::unique_ptr<biosoup::Sequence> ref_seq = std::unique_ptr<biosoup::Sequence>(new biosoup::Sequence("ref", ref_substr));

  auto lhs_sketch = minimizer_engine.Minimize(sequence, false);
  if (lhs_sketch.empty()) {
    return "";
  }

  auto rhs_sketch = minimizer_engine.Minimize(ref_seq, false);
  if (rhs_sketch.empty()) {
    return "";
  }

  RadixSort(lhs_sketch.begin(), lhs_sketch.end(), 9 * 2, ::First);
  RadixSort(rhs_sketch.begin(), rhs_sketch.end(), 9 * 2, ::First);

  std::vector<uint128_t> matches4;
  std::vector<uint128_t> matches3;

  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].first < rhs_sketch[j].first) {
        break;
      } else if (lhs_sketch[i].first == rhs_sketch[j].first) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].first != rhs_sketch[k].first) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].second & 1) == (rhs_sketch[k].second & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].second << 32 >> 33;
          std::uint64_t rhs_pos = rhs_sketch[k].second << 32 >> 33;
          std::uint64_t diagonal = !strand ?
              rhs_pos + lhs_pos :
              rhs_pos - lhs_pos + (3ULL << 30);


          if(strand == 1) {
            matches4.emplace_back(
              (((sequence->id << 1) | strand) << 32) | diagonal,
              (lhs_pos << 32) | rhs_pos);
          } else {
            matches3.emplace_back(
              (((sequence->id << 1) | strand) << 32) | diagonal,
              (lhs_pos << 32) | rhs_pos);
          }
        }
        break;
      } else {
        ++j;
      }
    }
  }

// CHECK VALID
    // return Chain(sequence->id, std::move(matches));

  std::vector<uint128_t> matches2;
  int resulting_strand = 0;

  if (matches4.size() > matches3.size()) {
    matches2 = matches4;
    resulting_strand = 1;
  } else {
    matches2 = matches3;
  }

  RadixSort(matches2.begin(), matches2.end(), 64, ::First);
  matches2.emplace_back(-1, -1);

  std::vector<uint128_t> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches2.size(); ++i) {
    if (matches2[i].first - matches2[j].first > 200) {
      if (i - j >= 2) {
        if (!intervals.empty() && intervals.back().second > j) {
          intervals.back().second = i;
        } else {
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches2[i].first - matches2[j].first > 200) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;

  for(int i = 0; i < intervals.size(); i++) {

    uint64_t first = intervals[i].first;
    uint64_t second = intervals[i].second;

    RadixSort(matches2.begin() + first, matches2.begin() + second, 64, ::Second);

    std::uint64_t strand = matches2[first].first >> 32 & 1;

    for (auto it = matches2.begin() + first; it != matches2.begin() + second; ++it) {
      uint128_t match = *it;
      std::uint32_t lhs_pos = match.second >> 32;
      std::uint32_t rhs_pos = match.second << 32 >> 32;
      std::uint32_t diagonal = match.first << 32 >> 32;
    }

    std::vector<std::uint64_t> indices;
    if (strand) {
      indices = LongestSubsequence(
          matches2.begin() + first,
          matches2.begin() + second,
          std::less<std::uint64_t>());
    } else {
      indices = LongestSubsequence(
          matches2.begin() + first,
          matches2.begin() + second,
          std::greater<std::uint64_t>());
    }

    for(int l = 0; l < indices.size(); l++) {
      uint128_t match = matches2[first + indices[l]];
      std::uint32_t lhs_pos = match.second >> 32;
      std::uint32_t rhs_pos = match.second << 32 >> 32;
    }

    indices.emplace_back(matches2.size() - 1 - first);
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      uint64_t condition1 = (matches2[first + indices[k]].second << 32 >> 32) - (matches2[first + indices[k - 1]].second << 32 >> 32);
      uint64_t condition2 = (matches2[first + indices[k - 1]].second << 32 >> 32) - (matches2[first + indices[k]].second << 32 >> 32);

      if ((condition1 > 30ULL && strand == 1) || (condition2 > 30ULL && strand == 0)) {
        if (k - l < 2) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches2[first + indices[m]].second >> 32;
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + 9;

          std::uint32_t rhs_pos = matches2[first + indices[m]].second << 32 >> 32;
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + 9 - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + 9;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < 15UL) {
          l = k;
          continue;
        }

        dst.emplace_back(
            ref_seq->id << 1 | strand,
            matches2[first + indices[l]].second >> 32,  // lhs_begin
            9 + (matches2[first + indices[k - 1]].second >> 32),  // lhs_end
            matches2[first].first >> 33,  // rhs_id
            strand ?  // rhs_begin
                matches2[first + indices[l]].second << 32 >> 32 :
                matches2[first + indices[k - 1]].second << 32 >> 32,
            9 + (strand ?  // rhs_end
                matches2[first + indices[k - 1]].second << 32 >> 32 :
                matches2[first + indices[l]].second << 32 >> 32),
            std::min(lhs_matches, rhs_matches));  // score

        l = k;
      }
    }
  }

  if(dst.size() == 0) {
    return "";
  }

  sort(dst.begin(), dst.end(), CompareOverlaps2);

  std::vector<AnchorCluster> clusters[STRANDS];
  std::vector<OverlapCluster> cluster_data;

  uint64_t last_end_ref = dst.front().rhs_end;

  for (int32_t l = 0; l < dst.size(); l++) {
    int reverseStrand = 0;
    biosoup::Overlap ovrl = dst[l];
    OverlapCluster cluster;
    cluster.overlap = ovrl;
    cluster.valid = false;
    int coverage = fabs(ovrl.rhs_end - ovrl.rhs_begin) + fabs(ovrl.lhs_end - ovrl.lhs_begin) / 2;

    if (resulting_strand) {
      clusters[0].emplace_back(AnchorCluster(dst[l].lhs_begin, std::max(dst[l].lhs_begin, dst[l].lhs_end - 9), dst[l].rhs_begin, std::max(dst[l].rhs_begin, dst[l].rhs_end - 9), coverage, l));
    } else {
      clusters[1].emplace_back(AnchorCluster(dst[l].lhs_begin, std::max(dst[l].lhs_begin, dst[l].lhs_end - 9), fabs(last_end_ref - dst[l].rhs_end), std::max(fabs(last_end_ref - dst[l].rhs_begin) - 9, fabs(last_end_ref - dst[l].rhs_end)), coverage, l));
    }
        
    cluster_data.push_back(cluster);
  }

  RNAFilterClusters(clusters, cluster_data, resulting_strand);

  std::vector<std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>> cluster_final;

  for (int k = 0; k < cluster_data.size(); ++k) {
    if(cluster_data[k].valid) {
      biosoup::Overlap overlap = cluster_data[k].overlap;
      
      std::pair<uint64_t, uint64_t> first_part(overlap.lhs_begin, overlap.lhs_end);
      std::pair<uint64_t, uint64_t> second_part((overlap.rhs_begin + left_edge), (overlap.rhs_end + left_edge));

      cluster_final.emplace_back(first_part, second_part);
    }
  }

  return GenerateSamFromAnchorClusters(calculated_rhs_id, resulting_strand, targets, sequence, cluster_final);
}

}  // namespace ram

// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"
#include "ram/bloom_filter.hpp"

#include <deque>
#include <iostream>
#include <stdexcept>
#include <map>
#include <set>

namespace ram {

MinimizerEngine::MinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::uint32_t k,
    std::uint32_t w,
    std::uint32_t bandwidth,
    std::uint32_t chain,
    std::uint32_t matches,
    std::uint32_t gap)
    : k_(std::min(std::max(k, 1U), 31U)),
      w_(w),
      bandwidth_(bandwidth),
      chain_(chain),
      matches_(matches),
      gap_(gap),
      occurrence_(-1),
      index_(1U << std::min(14U, 2 * k_)),
      thread_pool_(thread_pool ?
          thread_pool :
          std::make_shared<thread_pool::ThreadPool>(1)) {}

std::uint32_t MinimizerEngine::Index::Find(
    std::uint64_t key,
    const std::uint64_t** dst) const {
  auto it = locator.find(key << 1);
  if (it == locator.end()) {
    return 0;
  }
  if (it->first & 1) {
    *dst = &(it->second);
    return 1;
  }
  *dst = &(origins[it->second >> 32]);
  return static_cast<std::uint32_t>(it->second);
}

void MinimizerEngine::Minimize(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    bool minhash, double weightedMinimizerSampling) {

  for (auto& it : index_) {
    it.origins.clear();
    it.locator.clear();
  }

  if (first >= last) {
    return;
  }

  std::vector<std::vector<Kmer>> minimizers(index_.size());
  {
    std::uint64_t mask = index_.size() - 1;

    while (first != last) {
      std::size_t batch_size = 0;
      std::vector<std::future<std::vector<Kmer>>> futures;
      for (; first != last && batch_size < 50000000; ++first) {
        batch_size += (*first)->inflated_len;
        futures.emplace_back(thread_pool_->Submit(
            [&] (decltype(first) it) -> std::vector<Kmer> {
              return Minimize(*it, minhash, weightedMinimizerSampling);
            },
            first));
      }
      for (auto& it : futures) {
        for (const auto& jt : it.get()) {
          auto& m = minimizers[jt.value & mask];
          if (m.capacity() == m.size()) {
            m.reserve(m.capacity() * 1.5);
          }
          m.emplace_back(jt);
        }
      }
    }
  }

  {
    std::vector<std::future<std::pair<std::size_t, std::size_t>>> futures;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> std::pair<std::size_t, std::size_t> {
            if (minimizers[i].empty()) {
              return std::make_pair(0, 0);
            }

            RadixSort(
                minimizers[i].begin(),
                minimizers[i].end(),
                k_ * 2,
                Kmer::SortByValue);

            minimizers[i].emplace_back(-1, -1);  // stop dummy

            std::size_t num_origins = 0;
            std::size_t num_keys = 0;

            for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {  // NOLINT
              if (minimizers[i][j - 1].value != minimizers[i][j].value) {
                if (c > 1) {
                  num_origins += c;
                }
                ++num_keys;
                c = 0;
              }
            }

            return std::make_pair(num_origins, num_keys);
          },
          i));
    }
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      auto num_entries = futures[i].get();
      if (minimizers[i].empty()) {
        continue;
      }

      index_[i].origins.reserve(num_entries.first);
      index_[i].locator.reserve(num_entries.second);

      for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
        if (minimizers[i][j - 1].value != minimizers[i][j].value) {
          if (c == 1) {
            index_[i].locator.emplace(
                minimizers[i][j - 1].value << 1 | 1,
                minimizers[i][j - 1].origin);
          } else {
            index_[i].locator.emplace(
                minimizers[i][j - 1].value << 1,
                index_[i].origins.size() << 32 | c);
            for (std::uint64_t k = j - c; k < j; ++k) {
              index_[i].origins.emplace_back(minimizers[i][k].origin);
            }
          }
          c = 0;
        }
      }

      std::vector<Kmer>().swap(minimizers[i]);
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
    for (const auto& jt : it.locator) {
      if (jt.first & 1) {
        occurrences.emplace_back(1);
      } else {
        occurrences.emplace_back(static_cast<std::uint32_t>(jt.second));
      }
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
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    bool avoid_equal,
    bool avoid_symmetric,
    bool minhash,
    std::vector<std::uint32_t>* filtered,
    double weightedMinimizerSampling) const {

  auto sketch = Minimize(sequence, minhash, weightedMinimizerSampling);
  if (sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  std::vector<Match> matches;
  auto add_match = [&] (const Kmer& kmer, uint64_t origin) -> void {
    auto id = [] (std::uint64_t origin) -> std::uint32_t {
      return static_cast<std::uint32_t>(origin >> 32);
    };
    auto position = [] (std::uint64_t origin) -> std::uint32_t {
      return static_cast<std::uint32_t>(origin) >> 1;
    };
    auto strand = [] (std::uint64_t origin) -> bool {
      return origin & 1;
    };

    if (avoid_equal && sequence->id == id(origin)) {
      return;
    }
    if (avoid_symmetric && sequence->id > id(origin)) {
      return;
    }

    std::uint64_t rhs_id = id(origin);
    std::uint64_t strand_ = kmer.strand() == strand(origin);
    std::uint64_t lhs_pos = kmer.position();
    std::uint64_t rhs_pos = position(origin);
    std::uint64_t diagonal = !strand_ ?
        rhs_pos + lhs_pos :
        rhs_pos - lhs_pos + (3ULL << 30);

    matches.emplace_back(
        (((rhs_id << 1) | strand_) << 32) | diagonal,
        (lhs_pos << 32) | rhs_pos);
  };

  struct Hit {
    const Kmer* kmer;
    std::uint32_t n;
    const uint64_t* origins;

    Hit(const Kmer* kmer, std::uint32_t n, const uint64_t* origins)
        : kmer(kmer),
          n(n),
          origins(origins) {}

    bool operator<(const Hit& other) const {
      return n < other.n;
    }
  };
  std::vector<Hit> filtered_hits;

  std::uint64_t mask = index_.size() - 1;
  std::uint32_t prev = 0;

  sketch.emplace_back(-1, sequence->inflated_len << 1);  // stop dummy

  for (const auto& kmer : sketch) {
    std::uint32_t i = kmer.value & mask;
    const uint64_t* origins = nullptr;
    auto n = index_[i].Find(kmer.value, &origins);
    if (n > occurrence_) {
      filtered_hits.emplace_back(&kmer, n, origins);
      if (filtered) {
        filtered->emplace_back(kmer.position());
      }
      continue;
    }

    std::size_t rescuees = std::min(
        static_cast<std::size_t>(kmer.position() - prev) / bandwidth_,
        filtered_hits.size());
    if (rescuees) {
      std::partial_sort(
          filtered_hits.begin(),
          filtered_hits.begin() + rescuees,
          filtered_hits.end());
      for (auto it = filtered_hits.begin(); rescuees; rescuees--, ++it) {
        for (; it->n; it->n--, ++it->origins) {
          add_match(*it->kmer, *it->origins);
        }
      }
    }
    filtered_hits.clear();
    prev = kmer.position();

    for (; n; n--, ++origins) {
      add_match(kmer, *origins);
    }
  }

  return Chain(sequence->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs,
    bool minhash) const {

  auto lhs_sketch = Minimize(lhs, minhash);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs);
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(lhs_sketch.begin(), lhs_sketch.end(), k_ * 2, Kmer::SortByValue);
  RadixSort(rhs_sketch.begin(), rhs_sketch.end(), k_ * 2, Kmer::SortByValue);

  std::uint64_t rhs_id = rhs->id;

  std::vector<Match> matches;
  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].value < rhs_sketch[j].value) {
        break;
      } else if (lhs_sketch[i].value == rhs_sketch[j].value) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].value != rhs_sketch[k].value) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].strand() & 1) == (rhs_sketch[k].strand() & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].position();
          std::uint64_t rhs_pos = rhs_sketch[k].position();
          std::uint64_t diagonal = !strand ?
              rhs_pos + lhs_pos :
              rhs_pos - lhs_pos + (3ULL << 30);

          matches.emplace_back(
              (((rhs_id << 1) | strand) << 32) | diagonal,
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
    std::vector<Match>&& matches) const {
  RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);
  matches.emplace_back(-1, -1);  // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].group - matches[j].group > bandwidth_) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].group - matches[j].group > bandwidth_) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < chain_) {
      continue;
    }

    RadixSort(
        matches.begin() + j,
        matches.begin() + i,
        64,
        Match::SortByPositions);

    std::uint64_t strand = matches[j].strand();

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

    if (indices.size() < chain_) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if (matches[j + indices[k]].lhs_position() -
          matches[j + indices[k - 1]].lhs_position() > gap_) {
        if (k - l < chain_) {
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
          std::uint32_t lhs_pos = matches[j + indices[m]].lhs_position();
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[j + indices[m]].rhs_position();
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < matches_) {
          l = k;
          continue;
        }

        dst.emplace_back(
            lhs_id,
            matches[j + indices[l]].lhs_position(),
            k_ + matches[j + indices[k - 1]].lhs_position(),
            matches[j].rhs_id(),
            strand ?
                matches[j + indices[l]].rhs_position() :
                matches[j + indices[k - 1]].rhs_position(),
            k_ + (strand ?
                matches[j + indices[k - 1]].rhs_position() :
                matches[j + indices[l]].rhs_position()),
            std::min(lhs_matches, rhs_matches),
            strand);

        l = k;
      }
    }
  }
  return dst;
}

std::set<std::uint64_t> MinimizerEngine::CountKmer(const std::unique_ptr<biosoup::NucleicAcid>& sequence, double weightedMinimizerSampling) const {

  std::unordered_map<std::uint64_t, std::uint64_t> kmerMap; // map containing kmer - count pairs
  std::uint64_t kmerCount = 0;

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;
  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;

  for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= k_ - 1U) {
      if (minimizer < reverse_minimizer) {
        kmerMap[minimizer]++; // we increase count for the minimizer in the map
      } else if (minimizer > reverse_minimizer) {
        kmerMap[reverse_minimizer]++; // we increase count for the minimizer in the map
      }
      kmerCount++;
    }
  }

  std::uint64_t cutoffCount = (weightedMinimizerSampling/100) * kmerCount;

  // desc sorting based on individual kmer count
  std::vector<std::pair<std::uint64_t, std::uint64_t>> vector;
  auto comparator = [&] (std::pair<std::uint64_t, std::uint64_t> a, std::pair<std::uint64_t, std::uint64_t> b) -> bool {
    return a.second > b.second;
  };
  for (auto& it : kmerMap) {
    vector.push_back(it);
  }
  sort(vector.begin(), vector.end(), comparator);

  // we take all kmers that have count >= cutoffCount (those are the most frequent)
  std::set<std::uint64_t> set;
  for (auto &it : vector) {
    if (it.second >= cutoffCount) {
      set.insert(it.first);
    } else {
      break;
    }   
  }

  return set;
}

void MinimizerEngine::CountKmer(bloom_filter& mostFrequentKmersFilter, const std::unique_ptr<biosoup::NucleicAcid>& sequence, double weightedMinimizerSampling) const {

  std::unordered_map<std::uint64_t, std::uint64_t> kmerMap; // map containing kmer - count pairs
  std::uint64_t kmerCount = 0;

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;
  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;

  for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= k_ - 1U) {
      if (minimizer < reverse_minimizer) {
        kmerMap[minimizer]++; // we increase count for the minimizer in the map
      } else if (minimizer > reverse_minimizer) {
        kmerMap[reverse_minimizer]++; // we increase count for the minimizer in the map
      }
      kmerCount++;
    }
  }

  std::uint64_t cutoffCount = (weightedMinimizerSampling/100) * kmerCount;

  // desc sorting based on individual kmer count
  std::vector<std::pair<std::uint64_t, std::uint64_t>> vector;
  auto comparator = [&] (std::pair<std::uint64_t, std::uint64_t> a, std::pair<std::uint64_t, std::uint64_t> b) -> bool {
    return a.second > b.second;
  };
  for (auto& it : kmerMap) {
    vector.push_back(it);
  }
  sort(vector.begin(), vector.end(), comparator);

  // we take all kmers that have count >= cutoffCount (those are the most frequent)
  for (auto &it : vector) {
    if (it.second >= cutoffCount) {
      mostFrequentKmersFilter.insert(it.first);
    } else {
      break;
    }   
  }

}

std::vector<MinimizerEngine::Kmer> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    bool minhash, double weightedMinimizerSampling) const {
  if (sequence->inflated_len < k_) {
    return std::vector<Kmer>{};
  }

//  std::set<std::uint64_t> mostFrequentKmers;
//  if (weightedMinimizerSampling != 0) {
//    mostFrequentKmers = CountKmer(sequence, weightedMinimizerSampling);
//  }

  bloom_parameters parameters;
  parameters.projected_element_count = 1000;
  parameters.false_positive_probability = 0.0001;
  parameters.random_seed = 0xA5A5A5A5;
  parameters.compute_optimal_parameters();
  bloom_filter mostFrequentKmersFilter(parameters);
  if (weightedMinimizerSampling != 0) {
    CountKmer(mostFrequentKmersFilter, sequence, weightedMinimizerSampling);
  }

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;

  auto murmerhash64 = [&] (std::uint64_t key, std::uint64_t mask) -> std::uint64_t {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccd;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53;
    key ^= key >> 33;
    return key & mask;
  };

  auto applyWeight = [&] (std::uint64_t kmer) -> double {
    std::uint64_t hash = murmerhash64(kmer, UINT64_MAX);
    double x = hash * 1.0 / UINT64_MAX;  

    if (mostFrequentKmersFilter.contains(kmer)) {
      double p2 = x*x;
      double p4 = p2 * p2;
      return -1.0 * (p4 * p4);
    }

    return -1.0 * x;
  };

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

  std::deque<Kmer> window;
  
  auto window_add = [&] (std::uint64_t value, std::uint64_t location) -> void {
    while (!window.empty() && window.back().value > value) {
      window.pop_back();
    }
    window.emplace_back(value, location);
  };

  auto window_add_with_weight = [&] (std::uint64_t value, std::uint64_t location, double weight) -> void {
    while (!window.empty() && window.back().weight < weight) {
      window.pop_back();
    }
    window.emplace_back(value, location);
  };

  auto window_update = [&] (std::uint32_t position) -> void {
    while (!window.empty() && (window.front().position()) < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<Kmer> dst;

  for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
    std::uint64_t c = sequence->Code(i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= k_ - 1U) {
      if (minimizer < reverse_minimizer) {
        if (weightedMinimizerSampling == 0) {
          window_add(hash(minimizer), (i - (k_ - 1U)) << 1 | 0);
        } else {
          window_add_with_weight(hash(minimizer), (i - (k_ - 1U)) << 1 | 0, applyWeight(minimizer));
        }
      } else if (minimizer > reverse_minimizer) {
        if (weightedMinimizerSampling == 0) {
          window_add(hash(reverse_minimizer), (i - (k_ - 1U)) << 1 | 1);
        } else {
          window_add_with_weight(hash(reverse_minimizer), (i - (k_ - 1U)) << 1 | 1, applyWeight(minimizer));
        }
      }
    }
    if (i >= (k_ - 1U) + (w_ - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->value != window.front().value) {
          break;
        }
        if (it->origin & is_stored) {
          continue;
        }
        dst.emplace_back(it->value, id | it->origin);
        it->origin |= is_stored;
      }
      window_update(i - (k_ - 1U) - (w_ - 1U) + 1);
    }
  }

  if (minhash) {
    RadixSort(dst.begin(), dst.end(), k_ * 2, Kmer::SortByValue);
    dst.resize(sequence->inflated_len / k_);
    RadixSort(dst.begin(), dst.end(), 64, Kmer::SortByOrigin);
  }

  return dst;
}

template<typename RandomAccessIterator, typename Compare>
void MinimizerEngine::RadixSort(
    RandomAccessIterator first,
    RandomAccessIterator last,
    std::uint8_t max_bits,
    Compare comp) {  //  unary comparison function
  if (first >= last) {
    return;
  }

  std::vector<typename std::iterator_traits<RandomAccessIterator>::value_type> tmp(last - first);  // NOLINT
  auto begin = tmp.begin();
  auto end = tmp.end();

  std::uint64_t buckets[0x100]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};
    for (auto it = first; it != last; ++it) {
      ++counts[comp(*it) >> shift & 0xFF];
    }
    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it = first; it != last; ++it) {
      *(begin + buckets[comp(*it) >> shift & 0xFF]++) = *it;
    }
    std::swap(begin, first);
    std::swap(end, last);
  }

  if (shift / 8 & 1) {  // copy the sorted array for odd cases
    for (; first != last; ++first, ++begin) {
      *begin = *first;
    }
  }
}

template<typename Compare>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<Match>::const_iterator first,
    std::vector<Match>::const_iterator last,
    Compare comp) {  // binary comparison function
  if (first >= last) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(last - first + 1, 0);
  std::vector<std::uint64_t> predecessor(last - first, 0);

  std::uint64_t longest = 0;
  for (auto it = first; it != last; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if ((first + minimal[mid])->lhs_position() < it->lhs_position() &&
          comp((first + minimal[mid])->rhs_position(), it->rhs_position())) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - first] = minimal[lo - 1];
    minimal[lo] = it - first;
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

}  // namespace ram

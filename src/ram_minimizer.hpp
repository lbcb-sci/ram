
struct Minimizer {
  uint64_t value;
  uint64_t location;
  bool strand;

  Minimizer(uint64_t value, uint64_t location, bool strand) {
    this->value = value;
    this->location = location;
    this->strand = strand;
  }

  bool operator<(Minimizer other) const {
    return value < other.value;
  }

  std::string to_string() const {
    return "Location: " + std::to_string(location) + ", value: " + std::to_string(value) + ", strand: " + std::to_string(strand);
  }
};

struct Match {
  uint64_t ref_start;
  uint64_t read_start;
  bool strand;

  Match(Minimizer read_minimizer, Minimizer ref_minimizer) {
    this->ref_start = ref_minimizer.location;
    this->read_start = read_minimizer.location;
    this->strand = read_minimizer.strand != ref_minimizer.strand;
  }

  bool operator<(Match other) const {
    return ref_start < other.ref_start;
  }

  std::string to_string() const {
    return "ReadStart: " + std::to_string(read_start) + ", refStart: " + std::to_string(ref_start) + ", strand: " + std::to_string(strand);
  }
};

struct Cluster {
  uint64_t ref_start;
  uint64_t read_start;
  uint64_t ref_length;
  uint64_t read_length;
  bool strand;

  Cluster(uint64_t ref_start, uint64_t ref_length, uint64_t read_start, uint64_t read_length, bool strand) {
    this->ref_start = ref_start;
    this->read_start = read_start;
    this->ref_length = ref_length;
    this->read_length = read_length;
    this->strand = strand;
  }

  std::string to_string() const {
    return "ReadStart: " + std::to_string(read_start) + " readLen: " + std::to_string(read_length) + ", refStart: " + std::to_string(ref_start) + ", refLen: " + std::to_string(ref_length) + ", strand: " + std::to_string(strand);
  }
};

std::vector<Minimizer> create_minimizers(const char* sequence, uint64_t sequence_length,
                                                                     uint64_t k,
                                                                     uint64_t window_length);

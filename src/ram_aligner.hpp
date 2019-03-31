
int pairwise_alignment(const char* query, uint64_t query_length,
                       const char* target, uint64_t target_length,
                       int match,
                       int mismatch,
                       int gap,
                       std::string& cigar,
                       uint64_t& target_begin);

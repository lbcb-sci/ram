// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include "bioparser/fasta_parser.hpp"
#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace ram {
namespace test {

class RamMinimizerEngineTest: public ::testing::Test {
 public:
  void SetUp() override {
    biosoup::Sequence::num_objects = 0;
    auto p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(RAM_DATA_PATH);  // NOLINT
    s = p->Parse(-1);
    EXPECT_EQ(2, s.size());
  }

  std::vector<std::unique_ptr<biosoup::Sequence>> s;
};

TEST_F(RamMinimizerEngineTest, Map) {
  MinimizerEngine me{15, 5};
  me.Minimize(s.begin(), s.end());
  me.Filter(0.001);

  auto o = me.Map(s.front(), true, true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(30, o.front().lhs_begin);
  EXPECT_EQ(1869, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(0, o.front().rhs_begin);
  EXPECT_EQ(1893, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = me.Map(s.back(), true, true);
  EXPECT_TRUE(o.empty());

  o = me.Map(s.back(), true, false);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(1, o.front().lhs_id);
  EXPECT_EQ(0, o.front().lhs_begin);
  EXPECT_EQ(1893, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(30, o.front().rhs_begin);
  EXPECT_EQ(1869, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = me.Map(s.front(), false, true);
  EXPECT_EQ(2, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(2, o.front().lhs_begin);
  EXPECT_EQ(1897, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(2, o.front().rhs_begin);
  EXPECT_EQ(1897, o.front().rhs_end);
  EXPECT_EQ(1895, o.front().score);
  EXPECT_TRUE(o.front().strand);

  s.front()->ReverseAndComplement();
  o = me.Map(s.front(), true, true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(31, o.front().lhs_begin);
  EXPECT_EQ(1870, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(0, o.front().rhs_begin);
  EXPECT_EQ(1893, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_FALSE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Pair) {
  MinimizerEngine me{15, 5};
  auto o = me.Map(s.front(), s.back());
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(30, o.front().lhs_begin);
  EXPECT_EQ(1869, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(0, o.front().rhs_begin);
  EXPECT_EQ(1893, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = me.Map(s.back(), s.front());
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(1, o.front().lhs_id);
  EXPECT_EQ(0, o.front().lhs_begin);
  EXPECT_EQ(1893, o.front().lhs_end);
  EXPECT_EQ(0, o.front().rhs_id);
  EXPECT_EQ(30, o.front().rhs_begin);
  EXPECT_EQ(1869, o.front().rhs_end);
  EXPECT_EQ(585, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Filter) {
  MinimizerEngine me{9, 3};
  me.Minimize(s.begin(), s.end());

  me.Filter(0.001);
  auto o = me.Map(s.front(), true, true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(31, o.front().lhs_begin);
  EXPECT_EQ(1888, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(1, o.front().rhs_begin);
  EXPECT_EQ(1914, o.front().rhs_end);
  EXPECT_EQ(994, o.front().score);
  EXPECT_TRUE(o.front().strand);

  me.Filter(0.1);
  o = me.Map(s.front(), true, true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(31, o.front().lhs_begin);
  EXPECT_EQ(1888, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(1, o.front().rhs_begin);
  EXPECT_EQ(1914, o.front().rhs_end);
  EXPECT_EQ(980, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

TEST_F(RamMinimizerEngineTest, Micromize) {
  MinimizerEngine me{15, 5};
  me.Minimize(s.begin(), s.end());
  auto o = me.Map(s.front(), true, true, true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(80, o.front().lhs_begin);
  EXPECT_EQ(1857, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(55, o.front().rhs_begin);
  EXPECT_EQ(1881, o.front().rhs_end);
  EXPECT_EQ(242, o.front().score);
  EXPECT_TRUE(o.front().strand);

  o = me.Map(s.front(), s.back(), true);
  EXPECT_EQ(1, o.size());
  EXPECT_EQ(0, o.front().lhs_id);
  EXPECT_EQ(80, o.front().lhs_begin);
  EXPECT_EQ(1857, o.front().lhs_end);
  EXPECT_EQ(1, o.front().rhs_id);
  EXPECT_EQ(55, o.front().rhs_begin);
  EXPECT_EQ(1881, o.front().rhs_end);
  EXPECT_EQ(242, o.front().score);
  EXPECT_TRUE(o.front().strand);
}

}  // namespace test
}  // namespace ram

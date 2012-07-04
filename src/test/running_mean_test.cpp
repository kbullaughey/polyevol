#include "gtest/gtest.h"
#include "running_mean.h"
#include "error_handling.h"

#include <sstream>
#include <string>

using std::stringstream;
using std::string;

class RunningMeanTest : public ::testing::Test {
protected:
  RunningMeanTest() : some_means(4) { }
  RunningMean some_means;
};

TEST_F(RunningMeanTest, HasCorrectLength) {
  EXPECT_EQ(some_means.size(), 4);
}

TEST_F(RunningMeanTest, StartsOutAllZero) {
  for (int i=0; i < 4; i++) {
    EXPECT_EQ(some_means[i], 0);
  }
}

TEST_F(RunningMeanTest, CanPostSample) {
  some_means.post(0, 0.1);
  EXPECT_EQ(some_means[0], 0.1);
}

TEST_F(RunningMeanTest, ComputesMean) {
  some_means.post(0, 1);
  some_means.post(0, 2);
  EXPECT_EQ(some_means[0], 1.5);
}

TEST_F(RunningMeanTest, ComputesMeanSeveralMeans) {
  for (int i=0; i < 4; i++) {
    for (int j=0; j < 2; j++) {
      some_means.post(i, (double)(i+j));
    }
  }
  EXPECT_EQ(some_means[0], 0.5);
  EXPECT_EQ(some_means[1], 1.5);
  EXPECT_EQ(some_means[2], 2.5);
  EXPECT_EQ(some_means[3], 3.5);
}

TEST_F(RunningMeanTest, CountsCorrectly) {
  some_means.post(0, 5);
  some_means.post(0, 10);
  EXPECT_EQ(some_means.count(0), 2); 
}

TEST_F(RunningMeanTest, PrintsToStream) {
  stringstream output;
  some_means.post(0, -1);
  some_means.post(0, -7);
  some_means.post(1, 2);
  some_means.post(1, 4);
  some_means.post(1, 6);
  some_means.post(3, 1000);
  output << "answer:" << some_means;
  string expected("answer: -4,2 4,3 0,0 1000,1");
  EXPECT_EQ(output.str(), expected);
}

TEST_F(RunningMeanTest, PrecisionIsCorrect) {
  stringstream output;
  for (int i=0; i < 8; i++) {
    some_means.post(0, 0);
  }
  some_means.post(0,1);
  output << "answer:" << some_means;
  string expected("answer: 0.11111111,9 0,0 0,0 0,0");
  EXPECT_EQ(output.str(), expected);
}

TEST_F(RunningMeanTest, PostIndexThrowsException) {
  EXPECT_THROW(some_means.post(10,1), SimError);
}

TEST_F(RunningMeanTest, CountIndexThrowsException) {
  EXPECT_THROW(some_means.count(10), SimError);
}

TEST_F(RunningMeanTest, SubscriptIndexThrowsException) {
  EXPECT_THROW(some_means[10], SimError);
}

TEST_F(RunningMeanTest, PostNegativeIndexThrowsException) {
  EXPECT_THROW(some_means.post(-10,1), SimError);
}

TEST_F(RunningMeanTest, CountNegativeIndexThrowsException) {
  EXPECT_THROW(some_means.count(-10), SimError);
}

TEST_F(RunningMeanTest, SubscriptNegativeIndexThrowsException) {
  EXPECT_THROW(some_means[-10], SimError);
}


/* END */

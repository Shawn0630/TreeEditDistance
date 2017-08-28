#include "RNA.h"
#include "Tree.h"
#include "FileManage.h"
#include "TreeComparison.h"
#include "gtest/gtest.h"  

#include <string>
#include <iostream>
using namespace std;

class TreeComparisonTest : public testing::Test {
 protected:
  virtual void SetUp() {

  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(TreeComparisonTest, getTreeEditDistance) {

  string s1 = "";
  string s2 = "";

  RNA r1, r2;
  r1.setRNAName("A");
  r2.setRNAName("B"); 

  SimiMatrix matrix;
  string simiFileName = "ss_distance";

  FileManage* file = FileManage::getInstance();
  file->setSimiFileName(simiFileName);
  matrix = file->readSimiFromFile();

  Tree* t1 = r1.buildTree();
  Tree* t2 = r2.buildTree();

  s1 = "(F(A)(E(C(B))(D)))";
  s2 = "(F(A)(C(E(B)(D))))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(6);

  TreeComparison tc1(t1, t2, matrix);
  tc1.strategyComputation();
  ASSERT_EQ(tc1.getTreeDistance(), 2);  

} 
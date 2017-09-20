#include "RNA.h"
#include "Tree.h"
#include "FileManage.h"
#include "TreeComparison.h"
#include "TreeMap.h"
#include "gtest/gtest.h"  

#include <string>
#include <iostream>
#include <fstream>
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

  ofstream ou;
  
 

  string s1 = "";
  string s2 = "";

  RNA r1, r2;
  TreeComparison tc;

  Tree *t1, *t2;

  SimiMatrix matrix;
  TreeMap* map;
  string simiFileName = "ss_distance";
  FileManage* file = FileManage::getInstance();
  file->setSimiFileName(simiFileName);
  matrix = file->readSimiFromFile();
  tc.setCostModel(matrix);

  r1.setRNAName("A");
  r2.setRNAName("B"); 

  cout << "Test case #1" << endl;
  s1 = "(F(A)(E(C(B))(D)))";
  s2 = "(F(A)(C(E(B)(D))))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(6);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case1.txt");
  ou.open("tree_case1.txt");
  
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    2);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    2);
  cout << "Operation count: " << map->getCount() << endl;
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    2);
  cout << "Operation count: " << map->getCount() << endl;     
  cout << "Pass" << endl;


  cout << "Test case #2" << endl;
  s1 = "(A(B)(C(E)(F))(D))";
  s2 = "(A(B)(C(E))(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(5);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case2.txt");
  ou.open("tree_case2.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    1); 
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    1);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    1); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    1); 
  cout << map->toString();
  cout << "Operation count: " << map->getCount() << endl;
  cout << "Pass" << endl;


  cout << "Test case #3" << endl;
  s1 = "(A(B)(C(E)(F))(D))";
  s2 = "(A(B)(C)(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case3.txt");
  ou.open("tree_case3.txt");
  map = tc.getTreeMap();
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    2); 
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    2); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    2); 
  cout << map->toString();
  cout << "Operation count: " << map->getCount() << endl;
  cout << "Pass" << endl;


  cout << "Test case #4" << endl;
  s1 = "(A(S)(C(E)(F))(D))";
  s2 = "(A(B)(C)(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case4.txt");
  ou.open("tree_case4.txt");
  map = tc.getTreeMap();
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    3);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    3);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    3);
  cout << map->toString();
  cout << "Operation count: " << map->getCount() << endl;
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    3);       
  cout << map->toString();
  cout << "Operation count: " << map->getCount() << endl;
  cout << "Pass" << endl;


  cout << "Test case #5" << endl;
  s1 = "(A(B)(C(E)(F))(D))";
  s2 = "(A(B)(C))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(3);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case5.txt");
  ou.open("tree_case5.txt");
  map = tc.getTreeMap();
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    3);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    3);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    3);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    3);
  cout << map->toString();
  cout << "Operation count: " << map->getCount() << endl;
  cout << "Pass" << endl; 


  cout << "Test case #6" << endl;
  s1 = "(A(B(I)(J(U)))(C(D)(E(Q(N)(M))))(F(S)))";
  s2 = "(A)";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(13);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(1);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case6.txt");
  ou.open("tree_case6.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    12);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    12);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    12);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    12);
  cout << "Pass" << endl; 


  cout << "Test case #7" << endl;
  s1 = "(A(B(I)(J(U)))(C(D)(E(Q(N)(M))))(F(S)))";
  s2 = "(K)";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(13);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(1);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case7.txt");
  ou.open("tree_case7.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    13);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    13);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    13);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    13);
  cout << "Pass" << endl;


  cout << "Test case #8" << endl;
  s1 = "(K)";
  s2 = "(A(D)(E(L)(F(G)(H))))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(1);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(7);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case8.txt");
  ou.open("tree_case8.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    7);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    7);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    7);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    7);
  cout << "Pass" << endl;


  cout << "Test case #9" << endl;
  s1 = "(A(B)(C(D(F)(G(H)(I)))(E)))";
  s2 = "(G(H)(I))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(9);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(3);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case9.txt");
  ou.open("tree_case9.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    6);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    6);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    6);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    6);
  cout << "Pass" << endl;


  cout << "Test case #10" << endl;
  s1 = "(A(B(D)(E))(C))";
  s2 = "(F(G(H)(I))(K))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(5);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(5);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case10.txt");
  ou.open("tree_case10.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    5);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    5);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    5);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    5);
  cout << "Pass" << endl;


  cout << "Test case #11" << endl;
  s1 = "(A(B(C)(D(E)(F))))";
  s2 = "(B(C)(E)(F))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case11.txt");
  ou.open("tree_case11.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    2);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    2);
  cout << "Pass" << endl;



  cout << "Test case #12" << endl;
  s1 = "(A(B(C)(D(E)(F))))";
  s2 = "(B(C)(E)(F))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case12.txt");
  ou.open("tree_case12.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    2);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    2);
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    2);
  cout << "Pass" << endl;

  cout << "Test case #13" << endl;
  s1 = "(C(E)(F))";
  s2 = "(A(B)(C)(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(3);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case13.txt");
  ou.open("tree_case13.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    4);
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    4);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    4); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    4);    
  cout << "Pass" << endl;


  cout << "Test case #14" << endl;
  s1 = "(A(B)(C(E)(F))(D))";
  s2 = "(A(B)(C(E))(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(6);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(5);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case14.txt");
  ou.open("tree_case14.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    1); 
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    1);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    1); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    1); 
  cout << "Pass" << endl;


  cout << "Test case #15" << endl;
  s1 = "(A(A(A)(A(A)(A)))(A(A)(A(A))(A)))";
  s2 = "(A(A(A)(A)(A))(A(A(A)(A)(A))))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(11);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(10);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case15.txt");
  ou.open("tree_case15.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    3); 
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    3);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    3); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    3); 
  cout << "Pass" << endl;


  cout << "Test case #16" << endl;
  s1 = "(A(B)(C)(D(E))(F)(G)(H))";
  s2 = "(A(B(C))(D))";
  r1.setPreOrderSequence(s1);
  r1.setTreeSize(8);
  r2.setPreOrderSequence(s2);
  r2.setTreeSize(4);
  t1 = r1.buildTree();
  t2 = r2.buildTree();
  tc.setTreeA(t1);
  tc.setTreeB(t2);
  tc.init("out_case16.txt");
  ou.open("tree_case16.txt");
  if(DEBUG) {
    ou << "Tree A" << endl;
    ou << t1->toString() << endl;
    ou << "Tree B" << endl;
    ou << t2->toString() << endl;
  }
  ou.close();
  tc.strategyComputation();
  ASSERT_EQ(\
    tc.getTreeDistance(),\
    6); 
  map = tc.getTreeMap();
  cout << map->toString();
  ASSERT_EQ(\
    map->getCount(),\
    6);
  ASSERT_EQ(\
    tc.getTreeDistance_LL(),\
    6); 
  ASSERT_EQ(\
    tc.getTreeDistance_RR(),\
    6); 
  cout << "Pass" << endl;


  
} 
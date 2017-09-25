#include "TreeComparison.h"
#include "Tree.h"

#include <iostream>
#include <string.h>
#include <float.h>
#include <fstream>
#include <vector>
#include <stack>
#include <climits>

using namespace std;

TreeComparison::TreeComparison() {
}

TreeComparison::TreeComparison(Tree* A, Tree* B, SimiMatrix costModel) {
	A_ = A;
	B_ = B;
	costModel_ = costModel;

	treeSizeA = A_->getTreeSize();
	treeSizeB = B_->getTreeSize();

	int maxSize = treeSizeA < treeSizeB? treeSizeB + 1 : treeSizeA + 1;
	fn = new int[maxSize + 1];
	ft = new int[maxSize + 1];


	fn_ft_length = maxSize + 1;

  counter = 0;
  treeDist = 0;
  
  //map = new TreeMap(A_, B_, costModel_);

	Free = new int*[treeSizeA];
	LeftA = new int*[treeSizeA];
	LeftB = new int*[treeSizeA];
	RightA = new int*[treeSizeA];
	RightB = new int*[treeSizeA];
	AllA = new int*[treeSizeA];
	AllB = new int*[treeSizeA];
	LeftAStrategies = new Strategy*[treeSizeA];
	LeftBStrategies = new Strategy*[treeSizeA];
	RightAStrategies = new Strategy*[treeSizeA];
	RightBStrategies = new Strategy*[treeSizeA];
	AllAStrategies = new Strategy*[treeSizeA];
	AllBStrategies = new Strategy*[treeSizeA];
	FreeStrategies = new Strategy*[treeSizeA];
	delta = new float*[treeSizeA];
	hasVisited = new bool*[treeSizeA];

	for(int i = 0; i < treeSizeA; i++) {
		Free[i] = new int[treeSizeB];
		LeftA[i] = new int[treeSizeB];
		LeftB[i] = new int[treeSizeB];
		RightA[i] = new int[treeSizeB];
		RightB[i] = new int[treeSizeB];
		AllA[i] = new int[treeSizeB];
		AllB[i] = new int[treeSizeB];
		LeftAStrategies[i] = new Strategy[treeSizeB];
		LeftBStrategies[i] = new Strategy[treeSizeB];
		RightAStrategies[i] = new Strategy[treeSizeB];
		RightBStrategies[i] = new Strategy[treeSizeB];
		AllAStrategies[i] = new Strategy[treeSizeB];
		AllBStrategies[i] = new Strategy[treeSizeB];
		FreeStrategies[i] = new Strategy[treeSizeB];
		delta[i] = new float[treeSizeB];
		hasVisited[i] = new bool[treeSizeB];
	}

	s = new float*[maxSize - 1];
	t = new float*[maxSize - 1];
	q = new float[maxSize - 1];
	for(int i = 0; i < maxSize - 1; i++) {
		s[i] = new float[maxSize - 1];
		t[i] = new float[maxSize - 1];
	}

	for(int i = 0; i < treeSizeA; i++) {
		for(int j = 0; j < treeSizeB; j++) {
			Free[i][j] = -1;
			LeftA[i][j] = -1;
			LeftB[i][j] = -1;
			RightA[i][j] = -1;
			RightB[i][j] = -1;
			AllA[i][j] = -1;
			AllB[i][j] = -1;
			delta[i][j] = 0.0f;
			hasVisited[i][j] = false;
		}
	}

	ou.open("out.txt");

  computeSumInsAndDelCost(A_);
  computeSumInsAndDelCost(B_);

  cA_ = new CompressedTree(A_);
  cB_ = new CompressedTree(B_);

  computeSumInsAndDelCost_compressed(cA_);
  computeSumInsAndDelCost_compressed(cB_);

  if(DEBUG) {
    ou << "CompressedTree A" << endl;
    ou << cA_->toString() << endl;
    ou << "CompressedTree B" << endl;
    ou << cB_->toString() << endl;
  }
  deltaInit();

};

void TreeComparison::setTreeA(Tree* A) {
  A_ = A;
};

void TreeComparison::setTreeB(Tree* B) {
  B_ = B;
};

void TreeComparison::setRNAA(RNA rA) {
  rA_ = rA;
} ;
void TreeComparison::setRNAB(RNA rB) {
  rB_ = rB;
};

void TreeComparison::setCostModel(SimiMatrix costModel) {
  costModel_ = costModel;
}

void TreeComparison::init(string fileName) {
  if(A_ == NULL || B_ == NULL) return;
  treeSizeA = A_->getTreeSize();
  treeSizeB = B_->getTreeSize();

  //map = new TreeMap(A_, B_, costModel_);
  int maxSize = treeSizeA < treeSizeB? treeSizeB + 1 : treeSizeA + 1;
  fn = new int[maxSize + 1];
  ft = new int[maxSize + 1];

  fn_ft_length = maxSize + 1;


  counter = 0;
  treeDist = 0;

  Free = new int*[treeSizeA];
  LeftA = new int*[treeSizeA];
  LeftB = new int*[treeSizeA];
  RightA = new int*[treeSizeA];
  RightB = new int*[treeSizeA];
  AllA = new int*[treeSizeA];
  AllB = new int*[treeSizeA];
  LeftAStrategies = new Strategy*[treeSizeA];
  LeftBStrategies = new Strategy*[treeSizeA];
  RightAStrategies = new Strategy*[treeSizeA];
  RightBStrategies = new Strategy*[treeSizeA];
  AllAStrategies = new Strategy*[treeSizeA];
  AllBStrategies = new Strategy*[treeSizeA];
  FreeStrategies = new Strategy*[treeSizeA];
  delta = new float*[treeSizeA];
  hasVisited = new bool*[treeSizeA];

  for(int i = 0; i < treeSizeA; i++) {
    Free[i] = new int[treeSizeB];
    LeftA[i] = new int[treeSizeB];
    LeftB[i] = new int[treeSizeB];
    RightA[i] = new int[treeSizeB];
    RightB[i] = new int[treeSizeB];
    AllA[i] = new int[treeSizeB];
    AllB[i] = new int[treeSizeB];
    LeftAStrategies[i] = new Strategy[treeSizeB];
    LeftBStrategies[i] = new Strategy[treeSizeB];
    RightAStrategies[i] = new Strategy[treeSizeB];
    RightBStrategies[i] = new Strategy[treeSizeB];
    AllAStrategies[i] = new Strategy[treeSizeB];
    AllBStrategies[i] = new Strategy[treeSizeB];
    FreeStrategies[i] = new Strategy[treeSizeB];
    delta[i] = new float[treeSizeB];
    hasVisited[i] = new bool[treeSizeB];
  }

  s = new float*[maxSize - 1];
  t = new float*[maxSize - 1];
  q = new float[maxSize - 1];
  for(int i = 0; i < maxSize - 1; i++) {
    s[i] = new float[maxSize - 1];
    t[i] = new float[maxSize - 1];
  }

  for(int i = 0; i < treeSizeA; i++) {
    for(int j = 0; j < treeSizeB; j++) {
      Free[i][j] = -1;
      LeftA[i][j] = -1;
      LeftB[i][j] = -1;
      RightA[i][j] = -1;
      RightB[i][j] = -1;
      AllA[i][j] = -1;
      AllB[i][j] = -1;
      delta[i][j] = 0.0f;
      hasVisited[i][j] = false;
    }
  }

  ou.open(fileName);

  computeSumInsAndDelCost(A_);
  computeSumInsAndDelCost(B_);

  cA_ = new CompressedTree(A_);
  cB_ = new CompressedTree(B_);

  compressedTreeSizeA = cA_->getTreeSize();
  compressedTreeSizeB = cB_->getTreeSize();

  computeSumInsAndDelCost_compressed(cA_);
  computeSumInsAndDelCost_compressed(cB_);

  if(DEBUG) {
    ou << "CompressedTree A" << endl;
    ou << cA_->toString() << endl;
    ou << "CompressedTree B" << endl;
    ou << cB_->toString() << endl;
  }
  deltaInit();

};


void TreeComparison::deltaInit() {
	int treeSizeA = A_->getTreeSize();
	int treeSizeB = B_->getTreeSize();
	for(int i = 0; i < treeSizeA; i++) {
		Node* a = (*A_)[i];
		for(int j = 0; j < treeSizeB; j++) {
			Node* b = (*B_)[j];
			if(a->getSubTreeSize() == 1 && b->getSubTreeSize() == 1) {
				delta[i][j] = 0.0f;
			} else if(a->getSubTreeSize() == 1) {
				delta[i][j] = B_->preL_to_sumInsCost[j] - costModel_.ins(b->getLabel());
			} else if(b->getSubTreeSize() == 1) {
				delta[i][j] = A_->preL_to_sumDelCost[i] - costModel_.del(a->getLabel());
			}
		}
	}
};


void TreeComparison::computeSumInsAndDelCost_compressed(CompressedTree* tree) {
  int treeSize = tree->getTreeSize();
  for(int i = 0; i < treeSize; i++) {
    int postR_in_preL_in_compressed = treeSize - i - 1;// postOrder from bottom to up
    Node* node_in_compressed = (*tree)[postR_in_preL_in_compressed];
    Node* parent_in_compressed = node_in_compressed->getParent();

    vector<int> preLs_in_original = tree->compressed_to_original[postR_in_preL_in_compressed];
    vector<char> preLs_in_original_label = tree->compressed_to_original_label[postR_in_preL_in_compressed];
    for(int j = 0; j < preLs_in_original.size(); j++) {
      tree->preL_to_sumInsCost[postR_in_preL_in_compressed] += costModel_.ins(preLs_in_original_label[j]);
      tree->preL_to_sumDelCost[postR_in_preL_in_compressed] += costModel_.del(preLs_in_original_label[j]);
      tree->preL_to_InsCost[postR_in_preL_in_compressed] += costModel_.ins(preLs_in_original_label[j]);
      tree->preL_to_DelCost[postR_in_preL_in_compressed] += costModel_.del(preLs_in_original_label[j]);
    }

/*  cout << "Before" << endl;
    cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
    cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;
    cout << "ins " << node->getLabel() << " " << costModel_.ins(node->getLabel()) << endl;
    cout << "del " << node->getLabel() << " " << costModel_.del(node->getLabel()) << endl;
    cout << "After" << endl;
*/
    //tree->preL_to_sumInsCost[postR_in_preL_in_compressed] += costModel_.ins(node_in_compressed->getLabel());
    //tree->preL_to_sumDelCost[postR_in_preL_in_compressed] += costModel_.del(node_in_compressed->getLabel());

/*  cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
    cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;
*/
    if(parent_in_compressed != NULL) {
/*    cout << "Update Parent Before" << endl;
      cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
      cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;
*/
      tree->preL_to_sumInsCost[parent_in_compressed->getID()] += tree->preL_to_sumInsCost[postR_in_preL_in_compressed];
      tree->preL_to_sumDelCost[parent_in_compressed->getID()] += tree->preL_to_sumDelCost[postR_in_preL_in_compressed];
/*    cout << "Update Parent After" << endl;
      cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
      cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;
*/
    }
  }

};


void TreeComparison::computeSumInsAndDelCost(Tree* tree) {
	int treeSize = tree->getTreeSize();
	for(int i = 0; i < treeSize; i++) {
		int nodeForSum = treeSize - i - 1;// postOrder from bottom to up
		Node* node = (*tree)[nodeForSum];
		Node* parent = node->getParent();
	/*	cout << "Before" << endl;
		cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
		cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;
		cout << "ins " << node->getLabel() << " " << costModel_.ins(node->getLabel()) << endl;
		cout << "del " << node->getLabel() << " " << costModel_.del(node->getLabel()) << endl;
		cout << "After" << endl;*/
		tree->preL_to_sumInsCost[nodeForSum] += costModel_.ins(node->getLabel());
		tree->preL_to_sumDelCost[nodeForSum] += costModel_.del(node->getLabel());
/*		cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
		cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;*/
		if(parent != NULL) {
			/*cout << "Update Parent Before" << endl;
			cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
			cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;*/
			tree->preL_to_sumInsCost[parent->getID()] += tree->preL_to_sumInsCost[node->getID()];
			tree->preL_to_sumDelCost[parent->getID()] += tree->preL_to_sumDelCost[node->getID()];
/*			cout << "Update Parent After" << endl;
			cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
			cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;*/
		}
	}

};


Strategy** TreeComparison::APTED_ComputeOptStrategy_postL() {
	Strategy** S = new Strategy*[treeSizeA];
	float** cost1_L = new float*[treeSizeA];//in postL order
  float** cost1_R = new float*[treeSizeA];//in postL order
  float** cost1_I = new float*[treeSizeA];//in postL order
  int pathIDOffset = treeSizeA;

	for(int i = 0; i < treeSizeA; i++) {
		S[i] = new Strategy[treeSizeB];
		cost1_L[i] = new float[treeSizeB];
		cost1_R[i] = new float[treeSizeB];
		cost1_I[i] = new float[treeSizeB];
	}

	for(int i = 0; i < treeSizeA; i++) {
		for(int j = 0; j < treeSizeB; j++) {
			Strategy s;
			s.setLeaf(-1);
			s.setTreeToDecompose(0);
			S[i][j] = s;
		}
	}

	float* cost2_L = new float[treeSizeB];//in postL order
  float* cost2_R = new float[treeSizeB];//in postL order
  float* cost2_I = new float[treeSizeB];//in postL order
  int* cost2_path = new int[treeSizeB];//in postL order
  float* leafRow = new float[treeSizeB];

	float minCost = 0x7fffffffffffffffL;
	int strategyPath = -1;

	int size_i, size_j;

	stack<float*> rowsToReuse_L;
  stack<float*> rowsToReuse_R;
  stack<float*> rowsToReuse_I;

	for(int i = 0 ; i < treeSizeA; i++) {

		int i_in_preL = (A_)->postL_to_preL[i];
		size_i = (*A_)[i_in_preL]->getSubTreeSize();

		bool is_i_leaf = size_i == 1? true : false;
		Node* parent_i = (*A_)[i_in_preL]->getParent();
		int parent_i_preL;
		int parent_i_postL;

		int strategyLeftIndex = i_in_preL;
		int strategyRightIndex;
		int strategy_parent_i_LeftIndex, strategy_parent_i_RightIndex;

		int cost_L_LeftIndex, cost_L_RightIndex;
		int cost_R_LeftIndex, cost_R_RightIndex;
		int cost_I_LeftIndex, cost_I_RightIndex;

		int cost_L_parent_i_LeftIndex, cost_L_parent_i_RightIndex;
		int cost_R_parent_i_LeftIndex, cost_R_parent_i_RightIndex;
		int cost_I_parent_i_LeftIndex, cost_I_parent_i_RightIndex;

		if(parent_i != NULL) {
			parent_i_preL = parent_i->getID();
			parent_i_postL = (A_)->preL_to_postL[parent_i_preL];
		}


    int leftPath_i = (A_)->preL_to_lid[i_in_preL];
    int rightPath_i = (A_)->preL_to_rid[i_in_preL];
    int i_leftmost_forest = (*A_)[i_in_preL]->getLeftmostForestNum();
    int i_rightmost_forest = (*A_)[i_in_preL]->getRightmostForestNum();
    int i_special_forest = (*A_)[i_in_preL]->getSpecialForestNum();

    if(is_i_leaf) {
      if(DEBUG) {
      	ou << to_string(i_in_preL) << " is a leaf" << endl;
      }
      cost1_L[i] = leafRow;
      cost1_R[i] = leafRow;
      cost1_I[i] = leafRow;
      for(int j = 0; j < treeSizeB; j++) {
        strategyRightIndex = (B_)->postL_to_preL[j];
        if(DEBUG) {
        	ou << "set S[" << to_string(strategyLeftIndex) << ", " << to_string(strategyRightIndex) << "] = " << to_string(i_in_preL) << endl;
        }
        S[strategyLeftIndex][strategyRightIndex].setLeaf(i_in_preL);//decomposit the left tree
        S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(0);
        }
      }
      cost_L_LeftIndex = i;
      cost_R_LeftIndex = i;
      cost_I_LeftIndex = i;

      if(parent_i != NULL && cost1_L[parent_i_postL] == NULL) {
        if (rowsToReuse_L.empty()) {
         	cost1_L[parent_i_postL] = new float[treeSizeB];
          cost1_R[parent_i_postL] = new float[treeSizeB];
          cost1_I[parent_i_postL] = new float[treeSizeB];
        } else {
          cost1_L[parent_i_postL] = rowsToReuse_L.top();
          rowsToReuse_L.pop();
          cost1_R[parent_i_postL] = rowsToReuse_R.top();
          rowsToReuse_R.pop();
          cost1_I[parent_i_postL] = rowsToReuse_I.top();
          rowsToReuse_I.pop();
        }
      }

      if (parent_i != NULL) {
        cost_L_parent_i_LeftIndex = parent_i_postL;
        cost_R_parent_i_LeftIndex = parent_i_postL;
        cost_I_parent_i_LeftIndex = parent_i_postL;
        strategy_parent_i_LeftIndex = parent_i_preL;
      }

      fill_n(cost2_L, treeSizeB, 0L);
      fill_n(cost2_R, treeSizeB, 0L);
      fill_n(cost2_I, treeSizeB, 0L);
      fill_n(cost2_path, treeSizeB, 0);
      fill_n(leafRow, treeSizeB, 0);

      for(int j = 0; j < treeSizeB; j++) {
 			  int j_in_preL = (B_)->postL_to_preL[j];

 			  if(DEBUG) {
 				 ou << "compute " << to_string(i_in_preL) << ", " << to_string(j_in_preL) << endl;
 			  }

 			  strategyRightIndex = j_in_preL;
 			  Node* parent_j = (*B_)[j_in_preL]->getParent();
 			  int parent_j_preL, parent_j_postL;
        if (parent_j != NULL) {
        	parent_j_preL = parent_j->getID();
          parent_j_postL = (B_)->preL_to_postL[parent_j_preL];
        }
        if(DEBUG) {
        	ou << "parent i = "; 
        	if(parent_i == NULL) {
        		ou << "NULL";
        	} else {
        		ou << to_string(parent_i_preL);
        	}
        	ou << " parent j = ";
        	if(parent_j == NULL) {
        		ou << "NULL" << endl;
        	} else {
        		ou << to_string(parent_j_preL) << endl;
        	}
        }
        size_j = (*B_)[j_in_preL]->getSubTreeSize();
        bool is_j_leaf = size_j == 1? true : false;
        if (is_j_leaf) {
          cost2_L[j] = 0L;
          cost2_R[j] = 0L;
          cost2_I[j] = 0L;
          cost2_path[j] = j_in_preL;
        }
        minCost = 0x7fffffffffffffffL;
        int pathLeaf = -1;
        int treeToDecompose = 0;
        float tmpCost = 0x7fffffffffffffffL;

        if (size_i == 1) {
        	tmpCost = (float) (*B_)[j_in_preL]->getLeftmostForestNum();
        	if(DEBUG) {
					 ou << "left leaf(j) in " << to_string(j_in_preL) << " = " << to_string((B_)->preL_to_lid[j_in_preL]) << " tmpCost = " << to_string(tmpCost) << endl;
				  }
        	if(tmpCost < minCost) {
        		minCost = tmpCost;
        		pathLeaf = (B_)->preL_to_lid[j_in_preL];
        		treeToDecompose = 1;
        	}
          tmpCost = (float) (*B_)[j_in_preL]->getRightmostForestNum();
          if(DEBUG) {
					 ou << "right leaf(j) in " << to_string(j_in_preL) << " = " << to_string((B_)->preL_to_rid[j_in_preL]) << " tmpCost = " << to_string(tmpCost) << endl;
				  }
          if(tmpCost < minCost) {
          	minCost = tmpCost;
          	pathLeaf = (B_)->preL_to_rid[j_in_preL];
          	treeToDecompose = 1;
          }
        }
        if(size_j == 1) {
        	tmpCost = (float) (*A_)[i_in_preL]->getLeftmostForestNum();
        	if(DEBUG) {
					 ou << "left leaf(i) in " << to_string(i_in_preL) << " = " << to_string(leftPath_i) << " tmpCost = " << to_string(tmpCost) << endl;
				  }
        	if(tmpCost < minCost) {
        		minCost = tmpCost;
        		pathLeaf = leftPath_i;
        		treeToDecompose = 0;
        	}
        	tmpCost = (float) (*A_)[i_in_preL]->getRightmostForestNum();
        	if(DEBUG) {
					 ou << "right leaf(i) in " << to_string(i_in_preL) << " = " << to_string(rightPath_i) << " tmpCost = " << to_string(tmpCost) << endl;
				  }
        	if(tmpCost < minCost) {
        		minCost = tmpCost;
        		pathLeaf = rightPath_i;
        		treeToDecompose = 0;
        	}
        }
        if(size_i != 1 && size_j != 1) {
        	cost_L_RightIndex = j;//left path decomposition in i
				  tmpCost = (float) size_i * (float) (*B_)[j_in_preL]->getLeftmostForestNum() + cost1_L[cost_L_LeftIndex][cost_L_RightIndex];
				  if(DEBUG) {
					 ou << "left leaf(i) in " << to_string(i_in_preL) << " = " << to_string(leftPath_i) << " tmpCost = " << to_string(tmpCost) << endl;
				  }
          if (tmpCost < minCost) {
            minCost = tmpCost;
            pathLeaf = leftPath_i;
            treeToDecompose = 0;
            //S[strategyLeftIndex][strategyRightIndex].setLeaf(leftPath_i);
            //S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(0);
         	 }

         	 cost_R_RightIndex = j;//right path decomposition in i
         	 tmpCost = (float) size_i * (float) (*B_)[j_in_preL]->getRightmostForestNum() + cost1_R[cost_R_LeftIndex][cost_R_RightIndex];
         	 if(DEBUG) {
				    ou << "right leaf(i) in " << to_string(i_in_preL) << " = " << to_string(rightPath_i) << " tmpCost = " << to_string(tmpCost) << endl;
				   }
           if (tmpCost < minCost) {
            minCost = tmpCost;
            pathLeaf = rightPath_i;
            treeToDecompose = 0;
            //S[strategyLeftIndex][strategyRightIndex].setLeaf(rightPath_i);
            //S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(0);
           }

          	cost_I_RightIndex = j;//special path decomposition in i
          	tmpCost = (float) size_i * (float) (*B_)[j_in_preL]->getSpecialForestNum() + cost1_I[cost_I_LeftIndex][cost_I_RightIndex];
          	if(DEBUG) {
					    ou << "special leaf(i) in " << to_string(i_in_preL) << " = " << to_string((int)S[strategyLeftIndex][strategyRightIndex].getLeaf()) << " tmpCost = " << to_string(tmpCost) << " cost1_I[" << to_string(i_in_preL) << ", " << to_string(j_in_preL) << "] = " << to_string(cost1_I[cost_I_LeftIndex][cost_I_RightIndex]) << endl;
				    }
          	if (tmpCost < minCost) {
            	minCost = tmpCost;
            	strategyRightIndex = j_in_preL;
            	pathLeaf = (int)S[strategyLeftIndex][strategyRightIndex].getLeaf();
            	treeToDecompose = 0;
            	//S[strategyLeftIndex][strategyRightIndex].setLeaf((int)S[strategyLeftIndex][strategyRightIndex].getLeaf() + 1);
            	//S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(0);
          	}
          	//left path decomposition in j
          	tmpCost = (float) size_j * (float) i_leftmost_forest + cost2_L[j];
          	if(DEBUG) {
					     ou << "left leaf(j) in " << to_string(j_in_preL) << " = " << to_string((B_)->preL_to_lid[j_in_preL]) << " tmpCost = " << to_string(tmpCost) << endl;
				    }
          	if (tmpCost < minCost) {
            	minCost = tmpCost;
            	pathLeaf = (B_)->preL_to_lid[j_in_preL];
            	treeToDecompose = 1;
            	//S[strategyLeftIndex][strategyRightIndex].setLeaf((B_)->preL_to_lid[j_in_preL]);
            	//S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(1);
          	}

          	//right path decompostion in j
          	tmpCost = (float) size_j * (float) i_rightmost_forest + cost2_R[j];
          	if(DEBUG) {
					     ou << "right leaf(j) in " << to_string(j_in_preL) << " = " << to_string((B_)->preL_to_rid[j_in_preL]) << " tmpCost = " << to_string(tmpCost) << endl;
				    }
          	if (tmpCost < minCost) {
            	minCost = tmpCost;
            	pathLeaf = (B_)->preL_to_rid[j_in_preL];
            	treeToDecompose = 1;
            	//S[strategyLeftIndex][strategyRightIndex].setLeaf((B_)->preL_to_rid[j_in_preL]);
            	//S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(1);
          	}

          	//special path decompostion in j
          	tmpCost = (float) size_j * (float) i_special_forest + cost2_I[j];
          	if(DEBUG) {
					     ou << "special leaf(j) in " << to_string(j_in_preL) << " = " << to_string(cost2_path[j]) << " tmpCost = " << to_string(tmpCost) << " cost2_I[" << to_string(j_in_preL) << "] = " << cost2_I[j] << endl;
				    }
          	if (tmpCost < minCost) {
            	minCost = tmpCost;
            	pathLeaf = cost2_path[j];
            	treeToDecompose = 1;
            	//S[strategyLeftIndex][strategyRightIndex].setLeaf(cost2_path[j] + pathIDOffset + 1);
            	//S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(1);
          	}
          }

        	if (parent_i != NULL) {
        		cost_R_parent_i_RightIndex = j;
        		cost1_R[cost_R_parent_i_LeftIndex][cost_R_parent_i_RightIndex] += minCost;
        		tmpCost = -minCost + cost1_I[i][j];
        		/*if(DEBUG) {
          			ou << "minCost = " << to_string(minCost) << " cost1_I[" << to_string(i_in_preL) << ", " << to_string(j_in_preL) << "] = " << to_string(cost1_I[i][j]) << endl;
          		}*/
        		
        		if (tmpCost < cost1_I[parent_i_postL][j]) {
        			cost_I_parent_i_RightIndex = j;
        			if(DEBUG) {
          				ou << "cost1_I[" << to_string(parent_i_preL) << ", " << to_string(parent_j_preL) << "] = " << to_string(tmpCost) << endl;
          			}
            		cost1_I[cost_I_parent_i_LeftIndex][cost_I_parent_i_RightIndex] = tmpCost;
            		strategy_parent_i_RightIndex = j_in_preL;
            		strategyRightIndex = j_in_preL;
            		S[strategy_parent_i_LeftIndex][strategy_parent_i_RightIndex] = S[strategyLeftIndex][strategyRightIndex];
            		/*if(DEBUG) {
            			ou << "S[" << to_string(strategy_parent_i_LeftIndex) << ", " << to_string(strategy_parent_i_RightIndex) << "] = S[" << to_string(strategyLeftIndex) << ", " << to_string(strategyRightIndex) << "] = " << to_string(S[strategyLeftIndex][strategyRightIndex].getLeaf()) << endl;
            		}*/
          		}

          		vector<Node*> children = parent_i->getChildren();
          		bool is_i_leftmost_child = children[0]->getID() == i_in_preL? true : false;
          		bool is_i_rightmost_child = children[children.size() - 1]->getID() == i_in_preL? true : false;
          		/*if(DEBUG) {
          			ou << "is_i_leftmost_child = " << to_string(is_i_leftmost_child) << endl;
          			ou << "is_i_rightmost_child = " << to_string(is_i_rightmost_child) << endl;
          		}*/
          		if (is_i_rightmost_child) {
          			cost_I_parent_i_RightIndex = j;
          			cost_R_parent_i_RightIndex = j;
          			if(DEBUG) {
          				ou << "cost1_I[" << to_string(parent_i_preL) << ", " << to_string(j_in_preL) << "](" << to_string(cost1_I[parent_i_postL][j]) << ") += ";
          			}
            		cost1_I[cost_I_parent_i_LeftIndex][cost_I_parent_i_RightIndex] += cost1_R[cost_R_parent_i_LeftIndex][cost_R_parent_i_RightIndex];
            		if(DEBUG) {				
          				ou << "cost1_R[" << to_string(i_in_preL) << ", " << to_string(j_in_preL) << "](" << to_string(cost1_R[parent_i_postL][j]) << ") = ";
          				ou << to_string(cost1_I[parent_i_postL][j]) << endl;
          			}
            		cost_R_parent_i_RightIndex = j;
            		cost_R_RightIndex = j;
            		cost1_R[cost_R_parent_i_LeftIndex][cost_R_parent_i_RightIndex] += cost1_R[cost_R_LeftIndex][cost_R_RightIndex] - minCost;
          		}
          		if (is_i_leftmost_child) {
          			cost_L_parent_i_RightIndex = j;
          			cost_L_RightIndex = j;
          			cost1_L[cost_L_parent_i_LeftIndex][cost_L_parent_i_RightIndex] += cost1_L[cost_L_LeftIndex][cost_L_RightIndex];
          		} else {
          			cost_L_parent_i_RightIndex = j;
          			cost1_L[cost_L_parent_i_LeftIndex][cost_L_parent_i_RightIndex] += minCost;
          		}
        	}


        	if (parent_j != NULL) {
          	cost2_R[parent_j_postL] += minCost;
          	tmpCost = -minCost + cost2_I[j];
          	/*if(DEBUG) {
          		ou << "minCost = " << to_string(minCost) << " cost2_I[" << to_string(j_in_preL) << "] = " << to_string(cost2_I[j]) << endl;
          	}
            */
          	if (tmpCost < cost2_I[parent_j_postL]) {
          		if(DEBUG) {
          			ou << "cost2_I[" << to_string(parent_j_preL) << "] = " << to_string(tmpCost) << endl;
          		}
            	cost2_I[parent_j_postL] = tmpCost;
            	cost2_path[parent_j_postL] = cost2_path[j];
          	}
          	vector<Node*> children = parent_j->getChildren();
          	bool is_j_leftmost_child = children[0]->getID() == j_in_preL? true : false;
          	bool is_j_rightmost_child = children[children.size() - 1]->getID() == j_in_preL? true : false;
          /*if(DEBUG) {
          		ou << "is_j_leftmost_child = " << to_string(is_j_leftmost_child) << endl;
          		ou << "is_j_rightmost_child = " << to_string(is_j_rightmost_child) << endl;
          	}*/
          	if (is_j_rightmost_child) {
          		if(DEBUG) {
          			ou << "cost2_I[" << to_string(parent_j_preL) << "](" << to_string(cost2_I[parent_j_postL]) << ") += "; 
          		}
            	cost2_I[parent_j_postL] += cost2_R[parent_j_postL];
            	if(DEBUG) {
          			ou << "cost2_R[" << to_string(parent_j_preL) << "](" << to_string(cost2_R[parent_j_postL]) << ") = ";
          			ou << to_string(cost2_I[parent_j_postL]) << endl;
          		}
            	cost2_R[parent_j_postL] += cost2_R[j] - minCost;
          	}
          	if (is_j_leftmost_child) {
          	/*if(DEBUG) {
          		ou << "cost2_L[" << to_string(parent_j_postL) << "] += cost2_L[" << to_string(j) << "] = " << cost2_L[j] << endl;
          		ou << "cost2_L[" << to_string(parent_j_postL) << "] = " << cost2_L[parent_j_postL] << endl;
          		}*/
            	cost2_L[parent_j_postL] += cost2_L[j];
          	} else {
          	/*if(DEBUG) {
          		ou << "cost2_L[" << to_string(parent_j_postL) << "] += "  << to_string(minCost) << endl;
          		}*/
            	cost2_L[parent_j_postL] += minCost;
          	}
        	}
        	if(DEBUG) {
        		ou << "S[" << to_string(strategyLeftIndex) << ", " << to_string(strategyRightIndex) << "] = " << to_string(pathLeaf) << endl;
        	}
        	S[strategyLeftIndex][strategyRightIndex].setLeaf(pathLeaf);
        	S[strategyLeftIndex][strategyRightIndex].setTreeToDecompose(treeToDecompose);
		    }

        bool is_i_in_preL_leaf = (*A_)[i_in_preL]->getSubTreeSize() == 1? true : false;

        if (!is_i_in_preL_leaf) {
        	fill_n(cost1_L[i], treeSizeB, 0);
        	fill_n(cost1_R[i], treeSizeB, 0);
        	fill_n(cost1_I[i], treeSizeB, 0);
        	rowsToReuse_L.push(cost1_L[i]);
        	rowsToReuse_R.push(cost1_R[i]);
        	rowsToReuse_I.push(cost1_I[i]);
      	}
	}
	return S;
}

float TreeComparison::gted(Node* a, Node* b) {
	
  /*if(DEBUG) {
    ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
    ou << "hasVisited[" << to_string(a->getID()) << ", " << to_string(b->getID()) << "] = " << to_string(hasVisited[a->getID()][b->getID()]) << endl;
  }*/
	//if(hasVisited[a->getID()][b->getID()] == true) return delta[a->getID()][b->getID()] + costModel_.ren(a->getLabel(), b->getLabel());
  //if(hasVisited[a->getID()][b->getID()] == true) return delta[a->getID()][b->getID()];
	//hasVisited[a->getID()][b->getID()] = true;
	int treeSizeA = a->getSubTreeSize();
	int treeSizeB = b->getSubTreeSize();
	/*if(DEBUG) {
		ou << "treeSizeA = " << to_string(treeSizeA) << endl;
		ou << "treeSizeB = " << to_string(treeSizeB) << endl;
	}*/
	
	if ((treeSizeA == 1 || treeSizeB == 1)) {
    return spf1(a, treeSizeA, b, treeSizeB);
  }


	int pathLeaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
  int direction = FreeStrategies[a->getID()][b->getID()].getDirection();
	int treeToDecompose = FreeStrategies[a->getID()][b->getID()].getTreeToDecompose();
	Node* currentPathNode = treeToDecompose == 0? (*A_)[pathLeaf] : (*B_)[pathLeaf];


	if(treeToDecompose == 0) { // decompose tree A
		Node* parent = currentPathNode->getParent();
		int pathType = getPathType(A_, a, pathLeaf);// 0 left 1 right 2 special
		/*if(DEBUG) {
			ou << "getPathType A(" << to_string(a->getID()) << " ," << to_string(pathLeaf) << ") = " << to_string(pathType) << endl;
		}*/
    /*if(pathType == 0) {
       return spfL(a, b, pathLeaf, false);
    }
    else if(pathType == 1) {
       return spfR(a, b, pathLeaf, false);
    }
    else {*/
		  while(parent != NULL && parent->getID() >= a->getID()) {
        vector<Node*> children = parent->getChildren();
        for(int i = 0; i < children.size(); i++) {
          Node* child = children[i];
          /*if(DEBUG) {
            ou << "A child = " << to_string(child->getID()) << " currentPathNode = " << to_string(currentPathNode->getID()) << " parent = " << to_string(parent->getID()) << endl;
          }*/
          if(child->getID() != currentPathNode->getID()) {
            /*if(DEBUG) {
          	 ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
          	 ou << "create problem in A " << "gted(" << to_string(child->getID()) << ", " << to_string(b->getID()) << ")" << endl;
            }*/
            gted(child, b);
          }
        }
        currentPathNode = parent;
        parent = currentPathNode->getParent();
      }
      /*if(DEBUG) {
        ou << "swap = " << "false " << "pathType = " << to_string(pathType) << endl; 
      }*/
      if (pathType == 0) {
        return spfL(a, b, pathLeaf, false);
      }
      else if (pathType == 1) {
        return spfR(a, b, pathLeaf, false);
      }
      else {
        if(direction == 1) return spfA_RL(a, b, pathLeaf, pathType, NULL, false, false);//direction = right first add left then right
        else return spfA_LR(a, b, pathLeaf, pathType, NULL, false, false);//direction = left first add right then left
        //return spfA(a, b, pathLeaf, pathType, false);
    }	
  } 

	else if(treeToDecompose == 1) {
		Node* parent = currentPathNode->getParent();
		int pathType = getPathType(B_, b, pathLeaf);
		/*if(DEBUG) {
			ou << "getPathType B (" << to_string(b->getID()) << " ," << to_string(pathLeaf) << ") = " << to_string(pathType) << endl;
		}*/
   /* if(pathType == 0) {
       return spfL(a, b, pathLeaf, false);
    }
    else if(pathType == 1) {
       return spfR(a, b, pathLeaf, false);
    }
    else {*/
		  while(parent != NULL && parent->getID() >= b->getID()) {
			  vector<Node*> children = parent->getChildren();
			  for(int i = 0; i < children.size(); i++) {
				  Node* child = children[i];
				  /*if(DEBUG) {
            ou << "A child = " << to_string(child->getID()) << " currentPathNode = " << to_string(currentPathNode->getID()) << " parent = " << to_string(parent->getID()) << endl;
          }*/
				  if(child->getID() != currentPathNode->getID()) {
					 /*if(DEBUG) {
						  ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
          	 ou << "create problem in B " << "gted(" << to_string(a->getID()) << ", " << to_string(child->getID()) << ")" << endl;
            }*/
					 gted(a, child);
				  }
			  }
        currentPathNode = parent;
        parent = currentPathNode->getParent();
		  }
		/*if(DEBUG) {
      ou << "swap = " << "true " << "pathType = " << to_string(pathType) << endl; 
    }*/
    if(pathType == 0) {
      return spfL(b, a, pathLeaf, true);
    }
    else if(pathType == 1) {
      return spfR(b, a, pathLeaf, true);
    } else {
      if(direction == 1) return spfA_RL(b, a, pathLeaf, pathType, NULL, true, false);//direction = right first add left then right
      else return spfA_LR(b, a, pathLeaf, pathType, NULL, true, false);//direction = left first add right then left
		  //return spfA(b, a, pathLeaf, pathType, true);
    }
	}

};


float TreeComparison::gted_ND(Node* a, Node* b) {
  
  /*if(DEBUG) {
    ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
    ou << "hasVisited[" << to_string(a->getID()) << ", " << to_string(b->getID()) << "] = " << to_string(hasVisited[a->getID()][b->getID()]) << endl;
  }*/
  //if(hasVisited[a->getID()][b->getID()] == true) return delta[a->getID()][b->getID()] + costModel_.ren(a->getLabel(), b->getLabel());
  //if(hasVisited[a->getID()][b->getID()] == true) return delta[a->getID()][b->getID()];
  //hasVisited[a->getID()][b->getID()] = true;
  int treeSizeA = a->getSubTreeSize();
  int treeSizeB = b->getSubTreeSize();
  /*if(DEBUG) {
    ou << "treeSizeA = " << to_string(treeSizeA) << endl;
    ou << "treeSizeB = " << to_string(treeSizeB) << endl;
  }*/
  
  if ((treeSizeA == 1 || treeSizeB == 1)) {
    return spf1(a, treeSizeA, b, treeSizeB);
  }


  int pathLeaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
  int direction = FreeStrategies[a->getID()][b->getID()].getDirection();
  int treeToDecompose = FreeStrategies[a->getID()][b->getID()].getTreeToDecompose();
  Node* currentPathNode = treeToDecompose == 0? (*A_)[pathLeaf] : (*B_)[pathLeaf];


  if(treeToDecompose == 0) { // decompose tree A
    Node* parent = currentPathNode->getParent();
    int pathType = getPathType(A_, a, pathLeaf);// 0 left 1 right 2 special
    /*if(DEBUG) {
      ou << "getPathType A(" << to_string(a->getID()) << " ," << to_string(pathLeaf) << ") = " << to_string(pathType) << endl;
    }*/
    /*if(pathType == 0) {
       return spfL(a, b, pathLeaf, false);
    }
    else if(pathType == 1) {
       return spfR(a, b, pathLeaf, false);
    }
    else {*/
      while(parent != NULL && parent->getID() >= a->getID()) {
        vector<Node*> children = parent->getChildren();
        for(int i = 0; i < children.size(); i++) {
          Node* child = children[i];
          /*if(DEBUG) {
            ou << "A child = " << to_string(child->getID()) << " currentPathNode = " << to_string(currentPathNode->getID()) << " parent = " << to_string(parent->getID()) << endl;
          }*/
          if(child->getID() != currentPathNode->getID()) {
            /*if(DEBUG) {
             ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
             ou << "create problem in A " << "gted(" << to_string(child->getID()) << ", " << to_string(b->getID()) << ")" << endl;
            }*/
            gted_ND(child, b);
          }
        }
        currentPathNode = parent;
        parent = currentPathNode->getParent();
      }
      /*if(DEBUG) {
        ou << "swap = " << "false " << "pathType = " << to_string(pathType) << endl; 
      }*/
      if (pathType == 0) {
        return spfL(a, b, pathLeaf, false);
      }
      else if (pathType == 1) {
        return spfR(a, b, pathLeaf, false);
      }
      else {
        //if(direction == 1) return spfA_RL(a, b, pathLeaf, pathType, NULL, false, false);//direction = right first add left then right
        //else return spfA_LR(a, b, pathLeaf, pathType, NULL, false, false);//direction = left first add right then left
        return spfA(a, b, pathLeaf, pathType, false);
    } 
  } 

  else if(treeToDecompose == 1) {
    Node* parent = currentPathNode->getParent();
    int pathType = getPathType(B_, b, pathLeaf);
    /*if(DEBUG) {
      ou << "getPathType B (" << to_string(b->getID()) << " ," << to_string(pathLeaf) << ") = " << to_string(pathType) << endl;
    }*/
   /* if(pathType == 0) {
       return spfL(a, b, pathLeaf, false);
    }
    else if(pathType == 1) {
       return spfR(a, b, pathLeaf, false);
    }
    else {*/
      while(parent != NULL && parent->getID() >= b->getID()) {
        vector<Node*> children = parent->getChildren();
        for(int i = 0; i < children.size(); i++) {
          Node* child = children[i];
          /*if(DEBUG) {
            ou << "A child = " << to_string(child->getID()) << " currentPathNode = " << to_string(currentPathNode->getID()) << " parent = " << to_string(parent->getID()) << endl;
          }*/
          if(child->getID() != currentPathNode->getID()) {
           /*if(DEBUG) {
              ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
             ou << "create problem in B " << "gted(" << to_string(a->getID()) << ", " << to_string(child->getID()) << ")" << endl;
            }*/
           gted_ND(a, child);
          }
        }
        currentPathNode = parent;
        parent = currentPathNode->getParent();
      }
    /*if(DEBUG) {
      ou << "swap = " << "true " << "pathType = " << to_string(pathType) << endl; 
    }*/
    if(pathType == 0) {
      return spfL(b, a, pathLeaf, true);
    }
    else if(pathType == 1) {
      return spfR(b, a, pathLeaf, true);
    } else {
      //if(direction == 1) return spfA_RL(b, a, pathLeaf, pathType, NULL, true, false);//direction = right first add left then right
      //else return spfA_LR(b, a, pathLeaf, pathType, NULL, true, false);//direction = left first add right then left
      return spfA(b, a, pathLeaf, pathType, true);
    }
  }

};

float TreeComparison::spf1(Node* a, int treeSizeA, Node* b, int treeSizeB) {
  Tree* F, *G;
  F = A_;
  G = B_;

  if(DEBUG) {
    cout << "spf1(" << a->getID() << ", " << b->getID() << ") counter = " << counter  << endl;
  }
  if (treeSizeA == 1 && treeSizeB == 1) {
    float maxCost = costModel_.del(a->getLabel()) + costModel_.ins(b->getLabel());
    float renCost = costModel_.ren(a->getLabel(), b->getLabel());
    counter++;
    /*if(DEBUG) {
      if(renCost < maxCost) ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << renCost << endl;
      else ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << maxCost << endl;
    }*/
    return renCost < maxCost ? renCost : maxCost;
  }

  if (treeSizeA == 1) {
    float cost = G->preL_to_sumInsCost[b->getID()];
    float maxCost = cost + costModel_.del(a->getLabel());
    float minRenMinusIns = cost;
    float nodeRenMinusIns = 0;
    /*if(DEBUG) {
      ou << "cost = " << to_string(cost) << endl;
    }*/
    for (int i = b->getID(); i < b->getID() + treeSizeB; i++) {
      counter++;
      Node* node = (*G)[i];
      nodeRenMinusIns = costModel_.ren(a->getLabel(), node->getLabel()) - costModel_.ins(node->getLabel());
      /*if(DEBUG) {
        ou << "change " << a->getLabel() << " -> " << node->getLabel();
        ou << " from insert " << node->getLabel() << " nodeRenMinusIns = " << to_string(nodeRenMinusIns) << endl;
      }*/
      if (nodeRenMinusIns < minRenMinusIns) {
        minRenMinusIns = nodeRenMinusIns;
      }
    }
    cost += minRenMinusIns;
    /*if(DEBUG) {
      if(cost < maxCost) ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << cost << endl;
      else ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << maxCost << endl;
    }*/
    return cost < maxCost ? cost : maxCost;
  }

  if (treeSizeB == 1) {
    float cost = F->preL_to_sumDelCost[a->getID()];
    float maxCost = cost + costModel_.ins(b->getLabel());
    float minRenMinusDel = cost;
    float nodeRenMinusDel = 0;
    if(DEBUG) {
      ou << "cost = " << to_string(cost) << endl;
    }
    for (int i = a->getID(); i < a->getID() + treeSizeA; i++) {
      counter++;
      Node* node = (*F)[i];
      nodeRenMinusDel = costModel_.ren(node->getLabel(), b->getLabel()) - costModel_.del(node->getLabel());
      /*if(DEBUG) {
        ou << "change " << node->getLabel() << " -> " << b->getLabel();
        ou << " from insert " << node->getLabel() << " nodeRenMinusDel = " << to_string(nodeRenMinusDel) << endl;
      }*/
      if (nodeRenMinusDel < minRenMinusDel) {
        minRenMinusDel = nodeRenMinusDel;
      }
    }
    cost += minRenMinusDel;
    /*if(DEBUG) {
      if(cost < maxCost) ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << cost << endl;
      else ou << "spf1(" << a->getID() << ", " << b->getID() << ") = " << maxCost << endl;
    }*/
    return cost < maxCost ? cost : maxCost;
  }

  return -1;
};


float TreeComparison::spfL(Node* a, Node* b, int leaf, bool swap) {
	Tree *F, *G;
	if(swap) {
		F = B_;
		G = A_;
	} else {
		F = A_;
		G = B_;
	}
  if(DEBUG) {
    cout << "spfL(" << a->getID() << ", " << b->getID() << ") counter = " << counter << endl;
  }
	int* keyRoots = new int[(*G)[b->getID()]->getSubTreeSize()];
	int firstKeyRoot = computeKeyRoots(G, b, G->preL_to_lid[b->getID()], keyRoots, 0);	

	float** forestdist = new float*[(*F)[a->getID()]->getSubTreeSize() + 1];//consider the null
	for(int i = 0; i < (*F)[a->getID()]->getSubTreeSize() + 1; i++) {//consider the null
		forestdist[i] = new float[(*G)[b->getID()]->getSubTreeSize() + 1];
	}
  float dist = 0;
  for (int i = firstKeyRoot - 1; i >= 0; i--) {
    dist = treeEditDist(a, (*G)[keyRoots[i]], forestdist, swap, false);
  }
  return dist;
}

float TreeComparison::spfLL(Node* a, Node* b, int leaf, bool swap) {
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }
  if(DEBUG) {
    cout << "spfLL(" << a->getID() << ", " << b->getID() << ") counter = " << counter << endl;
  }
  int* GkeyRoots = new int[(*G)[b->getID()]->getSubTreeSize()];
  int GfirstKeyRoot = computeKeyRoots(G, b, G->preL_to_lid[b->getID()], GkeyRoots, 0);  

  int* FkeyRoots = new int[(*F)[a->getID()]->getSubTreeSize()];
  int FfirstKeyRoot = computeKeyRoots(F, a, F->preL_to_lid[a->getID()], FkeyRoots, 0);

  float** forestdist = new float*[(*F)[a->getID()]->getSubTreeSize() + 1];//consider the null
  for(int i = 0; i < (*F)[a->getID()]->getSubTreeSize() + 1; i++) {//consider the null
    forestdist[i] = new float[(*G)[b->getID()]->getSubTreeSize() + 1];
  }
  float dist = 0;
  for (int i = FfirstKeyRoot - 1; i >= 0; i--) {
    for (int j = GfirstKeyRoot - 1; j >= 0; j--) {
      dist = treeEditDist((*F)[FkeyRoots[i]], (*G)[GkeyRoots[j]], forestdist, swap, false);
    }
  }
  return dist;
}


float TreeComparison::spfLL_compressed(Node* a, Node* b, int leaf, bool swap) {
  CompressedTree *cF, *cG;
  Tree *F, *G;
  if(swap) {
    cF = cB_;
    F = B_;
    cG = cA_;
    G = A_;
  } else {
    cF = cA_;
    F = A_;
    cG = cB_;
    G = B_;
  }
  if(DEBUG) {
    cout << "spfLL(" << a->getID() << ", " << b->getID() << ") counter = " << counter << endl;
  }
  int* cGkeyRoots = new int[(*cG)[b->getID()]->getSubTreeSize()];
  int cGfirstKeyRoot = computeKeyRoots(cG, b, cG->preL_to_lid[b->getID()], cGkeyRoots, 0);  

  int* cFkeyRoots = new int[(*cF)[a->getID()]->getSubTreeSize()];
  int cFfirstKeyRoot = computeKeyRoots(cF, a, cF->preL_to_lid[a->getID()], cFkeyRoots, 0);

  int a_in_original = cF->compressed_to_original[a->getID()][0];
  int b_in_original = cG->compressed_to_original[b->getID()][0];
  float** forestdist = new float*[(*F)[a->getID()]->getSubTreeSize() + 1];//consider the null
  for(int i = 0; i < (*F)[a_in_original]->getSubTreeSize() + 1; i++) {//consider the null
    forestdist[i] = new float[(*G)[b_in_original]->getSubTreeSize() + 1];
  }
  float dist = 0;
  for (int i = cFfirstKeyRoot - 1; i >= 0; i--) {
    for (int j = cGfirstKeyRoot - 1; j >= 0; j--) {
      vector<int> cFkeyRoots_in_original = F->compressed_to_original[cFkeyRoots[i]];
      vector<int> cGkeyRoots_in_original = G->compressed_to_original[cGkeyRoots[j]];
      for(int Foff = 0; Foff < cFkeyRoots_in_original.size(); Foff++) {
        for(int Goff = 0; Goff < cGkeyRoots_in_original.size(); Goff++) {
          dist = treeEditDist_compressed((*F)[cFkeyRoots[i]], (*G)[cGkeyRoots[j]], forestdist, swap, Foff, Goff, false);
        }
      }
    }
  }
  return dist;
}

float TreeComparison::spfR(Node* a, Node* b, int leaf, bool swap) {
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }
  if(DEBUG) {
    cout << "spfR(" << a->getID() << ", " << b->getID() << ") counter = " << counter << endl;
  }
  int* keyRoots = new int[(*G)[b->getID()]->getSubTreeSize()];
  int firstKeyRoot = computeRevKeyRoots(G, b, G->preL_to_rid[b->getID()], keyRoots, 0);  

  float** forestdist = new float*[(*F)[a->getID()]->getSubTreeSize() + 1];//consider the null
  for(int i = 0; i < (*F)[a->getID()]->getSubTreeSize() + 1; i++) {//consider the null
    forestdist[i] = new float[(*G)[b->getID()]->getSubTreeSize() + 1];
  }

  float dist = 0;
  for (int i = firstKeyRoot - 1; i >= 0; i--) {
    dist = revTreeEditDist(a, (*G)[keyRoots[i]], forestdist, swap, false);
  }
  return dist;
}

float TreeComparison::spfRR(Node* a, Node* b, int leaf, bool swap) {
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }
  if(DEBUG) {
    cout << "spfRR(" << a->getID() << ", " << b->getID() << ") counter = " << counter << endl;
  }
  int* GkeyRoots = new int[(*G)[b->getID()]->getSubTreeSize()];
  int GfirstKeyRoot = computeRevKeyRoots(G, b, G->preL_to_rid[b->getID()], GkeyRoots, 0); 

  int* FkeyRoots = new int[(*F)[a->getID()]->getSubTreeSize()];
  int FfirstKeyRoot = computeRevKeyRoots(F, a, F->preL_to_rid[a->getID()], FkeyRoots, 0);   

  float** forestdist = new float*[(*F)[a->getID()]->getSubTreeSize() + 1];//consider the null
  for(int i = 0; i < (*F)[a->getID()]->getSubTreeSize() + 1; i++) {//consider the null
    forestdist[i] = new float[(*G)[b->getID()]->getSubTreeSize() + 1];
  }

  float dist = 0;
  for (int i = FfirstKeyRoot - 1; i >= 0; i--) {
    for (int j = GfirstKeyRoot - 1; j >= 0; j--){
      dist = revTreeEditDist((*F)[FkeyRoots[i]], (*G)[GkeyRoots[j]], forestdist, swap, false);
    }
  }
  return dist;
}



int TreeComparison::computeKeyRoots(Tree* G, Node* b, int leaf, int* keyRoots, int index) {
	/*if(DEBUG) {
    ou << "computeKeyRoots(" << to_string(b->getID()) << ", " << to_string(leaf) << ", " << to_string(index) << ")" << endl;
  }*/
	keyRoots[index++] = b->getID();

	int pathNode = leaf;
	while(pathNode > b->getID()) {
		Node* parent = (*G)[pathNode]->getParent();
		vector<Node*> children;
		if(parent != NULL)  children = parent->getChildren();
		for(int i = 0; i < children.size(); i++) {
			if(children[i]->getID() != pathNode) {
				index = computeKeyRoots(G, children[i], G->preL_to_lid[children[i]->getID()], keyRoots, index);
			}
		}
		if(parent != NULL)pathNode = parent->getID();
    else pathNode = -1;
	}
	return index;
}

int TreeComparison::computeRevKeyRoots(Tree* G, Node* b, int leaf, int* keyRoots, int index) {
  /*if(DEBUG) {
    ou << "computeRevKeyRoots(" << to_string(b->getID()) << ", " << to_string(leaf) << ", " << to_string(index) << ")" << endl;
  }*/
  keyRoots[index++] = b->getID();

  int pathNode = leaf;
  while(pathNode > b->getID()) {
    Node* parent = (*G)[pathNode]->getParent();
    vector<Node*> children;
    if(parent != NULL)  children = parent->getChildren();
    for(int i = 0; i < children.size(); i++) {
      if(children[i]->getID() != pathNode) {
        index = computeRevKeyRoots(G, children[i], G->preL_to_rid[children[i]->getID()], keyRoots, index);
      }
    }
    if(parent != NULL)pathNode = parent->getID();
    else pathNode = -1;
  }
  return index;
}
// postL all the trees on the left have already computed
// postR all the trees on the right have already computed
// add right and delete right
float TreeComparison::treeEditDist(Node* a, Node* b, float** forestdist, bool swap, bool mark) {
  if(DEBUG) {
    ou << "TreeDistance(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") swap = " << swap << endl;
  }
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }

  int a_in_preL = a->getID();
  int b_in_preL = b->getID();

  int a_in_postL = F->preL_to_postL[a_in_preL];
  int b_in_postL = G->preL_to_postL[b_in_preL];

  int a_leftmost_leaf_in_preL = F->preL_to_lid[a_in_preL];
  int b_leftmost_leaf_in_preL = G->preL_to_lid[b_in_preL];

  if(DEBUG) {
    ou << "a_leftmost_leaf_in_preL = " << a_leftmost_leaf_in_preL << endl;
    ou << "b_leftmost_leaf_in_preL = " << b_leftmost_leaf_in_preL << endl;
  } 

  int a_leftmost_leaf_in_postL = F->preL_to_postL[a_leftmost_leaf_in_preL];
  int b_leftmost_leaf_in_postL = G->preL_to_postL[b_leftmost_leaf_in_preL];

  int aoff = a_leftmost_leaf_in_postL - 1;//not tree-tree distance but tree-tree remove the root
  int boff = b_leftmost_leaf_in_postL - 1;//not tree-tree distance but tree-tree remove the root

  /*if(DEBUG) {
    ou << "aoff = " << to_string(aoff) << endl;
    ou << "boff = " << to_string(boff) << endl;
  }*/

  float da = 0;
  float db = 0;
  float dc = 0;
  float dist = 0;

  forestdist[0][0] = 0;
  for (int a1 = 1; a1 <= a_in_postL - aoff; a1++) {
    int a1_plus_aoff_in_preL = F->postL_to_preL[a1 + aoff];
    forestdist[a1][0] = forestdist[a1 - 1][0] + (swap ? costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel())); // USE COST MODEL - delete a1.
  }
  for (int b1 = 1; b1 <= b_in_postL - boff; b1++) {
    int b1_plus_boff_in_preL = G->postL_to_preL[b1 + boff];
    forestdist[0][b1] = forestdist[0][b1 - 1] + (swap ? costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel()) : costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - insert b1.
  }

  // Fill in the remaining costs.
  for (int a1 = 1; a1 <= a_in_postL - aoff; a1++) {
    for (int b1 = 1; b1 <= b_in_postL - boff; b1++) {
      counter++;
      int a1_plus_aoff_in_preL = F->postL_to_preL[a1 + aoff];
      int b1_plus_boff_in_preL = G->postL_to_preL[b1 + boff];
      /*if(DEBUG) {
        ou << "Compute forestdist(" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << ")" << endl;
      }*/
      float u = (swap ? costModel_.ren((*G)[b1_plus_boff_in_preL]->getLabel(), (*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.ren((*F)[a1_plus_aoff_in_preL]->getLabel(), (*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - rename a1 to g1.
      
      da = forestdist[a1 - 1][b1] + (swap ? costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel())); // USE COST MODEL - delete a1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1 - 1) << ", " << to_string(b1) << "] = " << to_string(forestdist[a1 - 1][b1]) << " ";
        if(swap) ou << "insert F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel());
        else ou << "delete F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel());
        ou << " da = " << to_string(da) << endl;
      }*/
      db = forestdist[a1][b1 - 1] + (swap ? costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel()) : costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - insert b1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1) << ", " << to_string(b1 - 1) << "] = " << to_string(forestdist[a1][b1 - 1]) << " ";
        if(swap) ou << "delete G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel());
        else ou << "insert G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel());
        ou << " db = " << to_string(db) << endl;
      }*/
      // If current subforests are subtrees. 
      int a1_preL = F->postL_to_preL[a1 + aoff];//a1+aoff the actual postL
      int a1_leftmost_leaf = F->preL_to_lid[a1_preL];
      int b1_preL = G->postL_to_preL[b1 + boff];
      int b1_leftmost_leaf = G->preL_to_lid[b1_preL];
      int a1_leftmost_leaf_in_postL = F->preL_to_postL[a1_leftmost_leaf];
      int b1_leftmost_leaf_in_postL = G->preL_to_postL[b1_leftmost_leaf];
      /*if(DEBUG) {
        ou << "a1_preL = " << to_string(a1_preL) << " a1_leftmost_leaf = " << to_string(a1_leftmost_leaf) << " a_leftmost_leaf_in_preL = " << to_string(a_leftmost_leaf_in_preL) << endl;
        ou << "b1_preL = " << to_string(b1_preL) << " b1_leftmost_leaf = " << to_string(b1_leftmost_leaf) << " b_leftmost_leaf_in_preL = " << to_string(b_leftmost_leaf_in_preL) << endl; 
      }*/
      if (a1_leftmost_leaf == a_leftmost_leaf_in_preL && b1_leftmost_leaf == b_leftmost_leaf_in_preL) {//is a tree
         dc = forestdist[a1 - 1][b1 - 1] + u;
         /*if(DEBUG) {
          if(swap) ou << (*G)[b1_plus_boff_in_preL]->getLabel() << " -> " << (*F)[a1_plus_aoff_in_preL]->getLabel();
          else ou << (*F)[a1_plus_aoff_in_preL]->getLabel() << " -> " << (*G)[b1_plus_boff_in_preL]->getLabel();
          ou << " dc = " << to_string(dc) << endl;
         }*/
         if (swap) {
           delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "] = " << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << endl;
           }*/
         } else {
           delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "] = " << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << endl;
           }*/
         }
      } else {
          dc = forestdist[a1_leftmost_leaf_in_postL - 1 - aoff][b1_leftmost_leaf_in_postL - 1 - boff] +
            (swap ? delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] : delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) + u;
          /*if(DEBUG) {
            ou << "dc = forestdist[" << to_string(a1_leftmost_leaf_in_postL - 1 - aoff) << ", " << to_string(b1_leftmost_leaf_in_postL - 1 - boff) << "](" << to_string(forestdist[a1_leftmost_leaf_in_postL - 1 - aoff][b1_leftmost_leaf_in_postL - 1 - boff]) << ")";
            ou << " + ";
            if(swap) ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "](" << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << ")";
            else ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "](" << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << ")";
            ou << " = " << to_string(dc) << endl;
          }*/
      }
      forestdist[a1][b1] = da >= db ? db >= dc ? dc : db : da >= dc ? dc : da;
      if(DEBUG) {
        ou << "(" << a_leftmost_leaf_in_preL << ", " << F->preL_to_preR[a1_preL] << ", " << b_leftmost_leaf_in_preL << ", " << G->preL_to_preR[b1_preL] << ")" << endl;
        ou << "forestdist[" << a1 << ", " << b1 << "] = " << forestdist[a1][b1] << endl;
      }
      dist = forestdist[a1][b1];
    }
  }
  return dist;
}

float TreeComparison::treeEditDist_compressed(Node* a, Node* b, float** forestdist, bool swap, int aoff_original, int boff_original, bool mark) {
  if(DEBUG) {
    ou << "TreeDistance(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") swap = " << swap << endl;
  }
  Tree *F, *G;
  CompressedTree *cF, *cG;
  if(swap) {
    F = B_;
    cF = cB_;
    G = A_;
    cG = cA_;
  } else {
    F = A_;
    cF = cA_;
    G = B_;
    cG = cB_;
  }



  int a_in_compressed_in_preL = a->getID();
  int b_in_compressed_in_preL = b->getID();

  vector<int> a_in_original = cF->compressed_to_original[a_in_compressed_in_preL];
  vecotr<int> b_in_original = cG->compressed_to_original[b_in_compressed_in_preL];

  int a_in_compressed_in_postL = cF->preL_to_postL[a_in_compressed_in_preL];
  int b_in_compressed_in_postL = cG->preL_to_postL[b_in_compressed_in_preL];

  int a_leftmost_leaf_in_compressed_in_preL = cF->preL_to_lid[a_in_compressed_in_preL];
  int b_leftmost_leaf_in_compressed_in_preL = cG->preL_to_lid[b_in_compressed_in_preL];

  int a_leftmost_leaf_in_compressed_in_postL = cF->preL_to_postL[a_leftmost_leaf_in_compressed_in_preL];
  int b_leftmost_leaf_in_compressed_in_postL = cG->preL_to_postL[b_leftmost_leaf_in_compressed_in_preL];

  int aoff = a_leftmost_leaf_in_compressed_in_postL - 1;//not tree-tree distance but tree-tree remove the root
  int boff = b_leftmost_leaf_in_compressed_in_postL - 1;//not tree-tree distance but tree-tree remove the root

  /*if(DEBUG) {
    ou << "aoff = " << to_string(aoff) << endl;
    ou << "boff = " << to_string(boff) << endl;
  }*/

  float da = 0;
  float db = 0;
  float dc = 0;
  float dist = 0;

  forestdist[0][0] = 0;
  int a1;
  for (a1 = 1; a1 < a_in_compressed_in_postL - aoff; a1++) {
    int a1_plus_aoff_in_preL = cF->postL_to_preL[a1 + aoff];
    forestdist[a1][0] = forestdist[a1 - 1][0] + (swap ? cF->preL_to_InsCost[a1_plus_aoff_in_preL]: cF->preL_to_DelCost[a1_plus_aoff_in_preL]); // USE COST MODEL - delete a1.
  }

  for (int i = 0; i <= aoff_original; i++, a1++) {
    int a_in_original_in_preL = a_in_original[a_in_orignal.size() - 1 - i];
    forestdist[a1][0] = forestdist[a1 - 1][0] + (swap ? costModel_.ins((*F)[a_in_original_in_preL]->getLabel()) : costModel_.del((*F)[a_in_original_in_preL]->getLabel()));
  }
  for (int b1 = 1; b1 < b_in_compressed_in_postL - boff; b1++) {
    int b1_plus_boff_in_preL = cG->postL_to_preL[b1 + boff];
    forestdist[0][b1] = forestdist[0][b1 - 1] + (swap ? cG->preL_to_DelCost[b1_plus_boff_in_preL] : cG->preL_to_InsCost[b1_plus_boff_in_preL]); // USE COST MODEL - insert b1.
  }
  for (int i = 0; i <= boff_original; i++, b1++) {
    int b_in_original_in_preL = b_in_original[b_in_orignal.size() - 1 - i];
    forestdist[0][b1] = forestdist[0][b1 - 1] + (swap ? costModel_.del((*G)[b_in_original_in_preL]->getLabel()) : costModel_.del((*G)[b_in_original_in_preL]->getLabel()));
  }

  // Fill in the remaining costs.
  for (int a1 = 1; a1 <= a_in_postL - aoff; a1++) {
    for (int b1 = 1; b1 <= b_in_postL - boff; b1++) {
      counter++;
      int a1_plus_aoff_in_preL = F->postL_to_preL[a1 + aoff];
      int b1_plus_boff_in_preL = G->postL_to_preL[b1 + boff];
      /*if(DEBUG) {
        ou << "Compute forestdist(" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << ")" << endl;
      }*/
      float u = (swap ? costModel_.ren((*G)[b1_plus_boff_in_preL]->getLabel(), (*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.ren((*F)[a1_plus_aoff_in_preL]->getLabel(), (*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - rename a1 to g1.
      
      da = forestdist[a1 - 1][b1] + (swap ? costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel())); // USE COST MODEL - delete a1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1 - 1) << ", " << to_string(b1) << "] = " << to_string(forestdist[a1 - 1][b1]) << " ";
        if(swap) ou << "insert F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel());
        else ou << "delete F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel());
        ou << " da = " << to_string(da) << endl;
      }*/
      db = forestdist[a1][b1 - 1] + (swap ? costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel()) : costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - insert b1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1) << ", " << to_string(b1 - 1) << "] = " << to_string(forestdist[a1][b1 - 1]) << " ";
        if(swap) ou << "delete G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel());
        else ou << "insert G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel());
        ou << " db = " << to_string(db) << endl;
      }*/
      // If current subforests are subtrees. 
      int a1_preL = F->postL_to_preL[a1 + aoff];//a1+aoff the actual postL
      int a1_leftmost_leaf = F->preL_to_lid[a1_preL];
      int b1_preL = G->postL_to_preL[b1 + boff];
      int b1_leftmost_leaf = G->preL_to_lid[b1_preL];
      int a1_leftmost_leaf_in_postL = F->preL_to_postL[a1_leftmost_leaf];
      int b1_leftmost_leaf_in_postL = G->preL_to_postL[b1_leftmost_leaf];
      /*if(DEBUG) {
        ou << "a1_preL = " << to_string(a1_preL) << " a1_leftmost_leaf = " << to_string(a1_leftmost_leaf) << " a_leftmost_leaf_in_preL = " << to_string(a_leftmost_leaf_in_preL) << endl;
        ou << "b1_preL = " << to_string(b1_preL) << " b1_leftmost_leaf = " << to_string(b1_leftmost_leaf) << " b_leftmost_leaf_in_preL = " << to_string(b_leftmost_leaf_in_preL) << endl; 
      }*/
      if (a1_leftmost_leaf == a_leftmost_leaf_in_preL && b1_leftmost_leaf == b_leftmost_leaf_in_preL) {//is a tree
         dc = forestdist[a1 - 1][b1 - 1] + u;
         /*if(DEBUG) {
          if(swap) ou << (*G)[b1_plus_boff_in_preL]->getLabel() << " -> " << (*F)[a1_plus_aoff_in_preL]->getLabel();
          else ou << (*F)[a1_plus_aoff_in_preL]->getLabel() << " -> " << (*G)[b1_plus_boff_in_preL]->getLabel();
          ou << " dc = " << to_string(dc) << endl;
         }*/
         if (swap) {
           delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "] = " << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << endl;
           }*/
         } else {
           delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "] = " << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << endl;
           }*/
         }
      } else {
          dc = forestdist[a1_leftmost_leaf_in_postL - 1 - aoff][b1_leftmost_leaf_in_postL - 1 - boff] +
            (swap ? delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] : delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) + u;
          /*if(DEBUG) {
            ou << "dc = forestdist[" << to_string(a1_leftmost_leaf_in_postL - 1 - aoff) << ", " << to_string(b1_leftmost_leaf_in_postL - 1 - boff) << "](" << to_string(forestdist[a1_leftmost_leaf_in_postL - 1 - aoff][b1_leftmost_leaf_in_postL - 1 - boff]) << ")";
            ou << " + ";
            if(swap) ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "](" << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << ")";
            else ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "](" << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << ")";
            ou << " = " << to_string(dc) << endl;
          }*/
      }
      forestdist[a1][b1] = da >= db ? db >= dc ? dc : db : da >= dc ? dc : da;
      if(DEBUG) {
        ou << "(" << a_leftmost_leaf_in_preL << ", " << F->preL_to_preR[a1_preL] << ", " << b_leftmost_leaf_in_preL << ", " << G->preL_to_preR[b1_preL] << ")" << endl;
        ou << "forestdist[" << a1 << ", " << b1 << "] = " << forestdist[a1][b1] << endl;
      }
      dist = forestdist[a1][b1];
    }
  }
  return dist;
}

// postL all the trees on the left have already computed
// postR all the trees on the right have already computed
// add left and delete left
float TreeComparison::revTreeEditDist(Node* a, Node* b, float** forestdist, bool swap, bool mark) {
  if(DEBUG) {
    ou << "RevTreeDistance(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") swap = " << swap << endl;
  }
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }

  int a_in_preL = a->getID();
  int b_in_preL = b->getID();

  int a_in_postR = F->preL_to_postR[a_in_preL];
  int b_in_postR = G->preL_to_postR[b_in_preL];

  int a_rightmost_leaf_in_preL = F->preL_to_rid[a_in_preL];
  int b_rightmost_leaf_in_preL = G->preL_to_rid[b_in_preL];
  if(DEBUG) {
    ou << "a_rightmost_leaf_in_preL = " << a_rightmost_leaf_in_preL << endl;
    ou << "b_rightmost_leaf_in_preL = " << b_rightmost_leaf_in_preL << endl;
  }

  int a_rightmost_leaf_in_postR = F->preL_to_postR[a_rightmost_leaf_in_preL];
  int b_rightmost_leaf_in_postR = G->preL_to_postR[b_rightmost_leaf_in_preL];

  int aoff = a_rightmost_leaf_in_postR - 1;//not tree-tree distance but tree-tree remove the root
  int boff = b_rightmost_leaf_in_postR - 1;//not tree-tree distance but tree-tree remove the root
  

  /*if(DEBUG) {
    ou << "aoff = " << to_string(aoff) << endl;
    ou << "boff = " << to_string(boff) << endl;
  }*/

  float da = 0;
  float db = 0;
  float dc = 0;
  float dist = 0;

  forestdist[0][0] = 0;
  for (int a1 = 1; a1 <= a_in_postR - aoff; a1++) {
    int a1_plus_aoff_in_preL = F->postR_to_preL[a1 + aoff];
    forestdist[a1][0] = forestdist[a1 - 1][0] + (swap ? costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel())); // USE COST MODEL - delete a1.
  }
  for (int b1 = 1; b1 <= b_in_postR - boff; b1++) {
    int b1_plus_boff_in_preL = G->postR_to_preL[b1 + boff];
    forestdist[0][b1] = forestdist[0][b1 - 1] + (swap ? costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel()) : costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - insert b1.
  }

  /*if(DEBUG) {
    for(int i = 0; i <= a_in_postR - aoff; i++) {
      for(int j = 0; j <= b_in_postR - boff; j++) {
        ou << forestdist[i][j] << " ";
      }
      ou << endl;
    }
  }*/

  // Fill in the remaining costs.
  for (int a1 = 1; a1 <= a_in_postR - aoff; a1++) {
    for (int b1 = 1; b1 <= b_in_postR - boff; b1++) {
      counter++;
      int a1_plus_aoff_in_preL = F->postR_to_preL[a1 + aoff];
      int b1_plus_boff_in_preL = G->postR_to_preL[b1 + boff];
      /*if(DEBUG) {
        ou << "Compute revforestdist(" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << ")" << endl;
      }*/
      float u = (swap ? costModel_.ren((*G)[b1_plus_boff_in_preL]->getLabel(), (*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.ren((*F)[a1_plus_aoff_in_preL]->getLabel(), (*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - rename a1 to g1.
      
      da = forestdist[a1 - 1][b1] + (swap ? costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel()) : costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel())); // USE COST MODEL - delete a1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1 - 1) << ", " << to_string(b1) << "] = " << to_string(forestdist[a1 - 1][b1]) << " ";
        if(swap) ou << "insert F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.ins((*F)[a1_plus_aoff_in_preL]->getLabel());
        else ou << "delete F " << (*F)[a1_plus_aoff_in_preL]->getLabel() << " + " << costModel_.del((*F)[a1_plus_aoff_in_preL]->getLabel());
        ou << " da = " << to_string(da) << endl;
      }*/
      db = forestdist[a1][b1 - 1] + (swap ? costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel()) : costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel())); // USE COST MODEL - insert b1.
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1) << ", " << to_string(b1 - 1) << "] = " << to_string(forestdist[a1][b1 - 1]) << " ";
        if(swap) ou << "delete G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.del((*G)[b1_plus_boff_in_preL]->getLabel());
        else ou << "insert G " << (*G)[b1_plus_boff_in_preL]->getLabel() << " + " << costModel_.ins((*G)[b1_plus_boff_in_preL]->getLabel());
        ou << " db = " << to_string(db) << endl;
      }*/
      // If current subforests are subtrees. 
      int a1_preL = F->postR_to_preL[a1 + aoff];//a1+aoff the actual postL
      int a1_rightmost_leaf = F->preL_to_rid[a1_preL];
      int b1_preL = G->postR_to_preL[b1 + boff];
      int b1_rightmost_leaf = G->preL_to_rid[b1_preL];
      int a1_rightmost_leaf_in_postR = F->preL_to_postR[a1_rightmost_leaf];
      int b1_rightmost_leaf_in_postR = G->preL_to_postR[b1_rightmost_leaf];
      /*if(DEBUG) {
        ou << "a1_preL = " << to_string(a1_preL) << " a1_rightmost_leaf = " << to_string(a1_rightmost_leaf) << " a_rightmost_leaf_in_preL = " << to_string(a_rightmost_leaf_in_preL) << endl;
        ou << "b1_preL = " << to_string(b1_preL) << " b1_rightmost_leaf = " << to_string(b1_rightmost_leaf) << " b_rightmost_leaf_in_preL = " << to_string(b_rightmost_leaf_in_preL) << endl; 
      }*/
      if (a1_rightmost_leaf == a_rightmost_leaf_in_preL && b1_rightmost_leaf == b_rightmost_leaf_in_preL) {//is a tree
         dc = forestdist[a1 - 1][b1 - 1] + u;
         /*if(DEBUG) {
          if(swap) ou << (*G)[b1_plus_boff_in_preL]->getLabel() << " -> " << (*F)[a1_plus_aoff_in_preL]->getLabel();
          else ou << (*F)[a1_plus_aoff_in_preL]->getLabel() << " -> " << (*G)[b1_plus_boff_in_preL]->getLabel();
          ou << " dc = " << to_string(dc) << endl;
         }*/
         if (swap) {
           delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "] = forestdist[" << to_string(a1 - 1) << ", " << to_string(b1 - 1) << "] = " << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << endl;
           }*/
         } else {
           delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL] = forestdist[a1 - 1][b1 - 1];
           /*if(DEBUG) {
            ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "] = forestdist[" << to_string(a1 - 1) << ", " << to_string(b1 - 1) << "] = " << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << endl;
           }*/
         }
      } else {
          dc = forestdist[a1_rightmost_leaf_in_postR - 1 - aoff][b1_rightmost_leaf_in_postR - 1 - boff] +
            (swap ? delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL] : delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) + u;
          /*if(DEBUG) {
            ou << "dc = forestdist[" << to_string(a1_rightmost_leaf_in_postR - 1 - aoff) << ", " << to_string(b1_rightmost_leaf_in_postR - 1 - boff) << "](" << to_string(forestdist[a1_rightmost_leaf_in_postR - 1 - aoff][b1_rightmost_leaf_in_postR - 1 - boff]) << ")";
            ou << " + ";
            if(swap) ou << "delta[" << to_string(b1_plus_boff_in_preL) << ", " << to_string(a1_plus_aoff_in_preL) << "](" << to_string(delta[b1_plus_boff_in_preL][a1_plus_aoff_in_preL]) << ")";
            else ou << "delta[" << to_string(a1_plus_aoff_in_preL) << ", " << to_string(b1_plus_boff_in_preL) << "](" << to_string(delta[a1_plus_aoff_in_preL][b1_plus_boff_in_preL]) << ")";
            ou << " = " << to_string(dc) << endl;
          }*/
      }
      /*if(DEBUG) {
        ou << "da = " << to_string(da) << endl;
        ou << "db = " << to_string(db) << endl;
        ou << "dc = " << to_string(dc) << endl;
      }*/
      forestdist[a1][b1] = da >= db ? db >= dc ? dc : db : da >= dc ? dc : da;
      if(DEBUG) {
        ou << "(" << a1_preL << ", " << F->preL_to_preR[a_rightmost_leaf_in_preL] << ", " << b1_preL << ", " << G->preL_to_preR[b_rightmost_leaf_in_preL] << ")" << endl;
        ou << "forestdist[" << a1 << ", " << b1 << "] = " << forestdist[a1][b1] << endl;
      }
      dist = forestdist[a1][b1];
      /*if(mark) {
        if(dist == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[a1_plus_aoff_in_preL]->getLabel() << endl;//insert
            }
            map->setMap(-1, a1_plus_aoff_in_preL);
          }
          else {
            if(DEBUG) {
              ou << (*F)[a1_plus_aoff_in_preL]->getLabel() << " -> -" << endl;//delete
            }
            map->setMap(a1_plus_aoff_in_preL, -1);
          }
        }
        else if(dist == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[b1_plus_boff_in_preL]->getLabel() << " -> -" << endl;//delete
            }
            map->setMap(b1_plus_boff_in_preL, -1);
          }
          else {
            if(DEBUG) {
              ou << "- -> " << (*G)[b1_plus_boff_in_preL]->getLabel() << endl;//insert
            }
            map->setMap(-1, b1_plus_boff_in_preL);
          }
        }
        else {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[b1_plus_boff_in_preL]->getLabel() << "-> " << (*F)[a1_plus_aoff_in_preL]->getLabel() << endl;
            }
            map->setMap(b1_plus_boff_in_preL, a1_plus_aoff_in_preL);
          }
          else {
            if(DEBUG) {
              ou << (*F)[a1_plus_aoff_in_preL]->getLabel() << " -> " << (*G)[b1_plus_boff_in_preL]->getLabel() << endl;
            }
            map->setMap(a1_plus_aoff_in_preL, b1_plus_boff_in_preL);
          }
        }
      }*/
      /*if(DEBUG) {
        ou << "forestdist[" << to_string(a1) << ", " << to_string(b1) << "] = " << to_string(forestdist[a1][b1]) << endl;
      }*/
    }
  }
  /*if(mark && DEBUG) {
    for(int i = 0; i < F->getTreeSize() + 1; i++) {
      for(int j = 0; j < G->getTreeSize() + 1; j++) {
        ou << forestdist[i][j] << " ";
      }
      ou << endl;
    }

    int i = F->getTreeSize();
    int j = G->getTreeSize();

    while(i != 0 && j != 0) {
      if(DEBUG) {
        ou << "i = " << i << " j = " << j << endl;
      }
      int i_minus1_in_preL = F->postR_to_preL[i - 1];
      int j_minus1_in_preL = G->postR_to_preL[j - 1];
      float da = forestdist[i - 1][j] + (swap? costModel_.ins((*F)[i_minus1_in_preL]->getLabel()) : costModel_.del((*F)[i_minus1_in_preL]->getLabel()));
      float db = forestdist[i][j - 1] + (swap? costModel_.del((*G)[j_minus1_in_preL]->getLabel()) : costModel_.ins((*G)[j_minus1_in_preL]->getLabel()));
      float dc = forestdist[i - 1][j - 1] + (swap? costModel_.ren((*G)[j_minus1_in_preL]->getLabel(), (*F)[i_minus1_in_preL]->getLabel()) : costModel_.ren((*F)[i_minus1_in_preL]->getLabel(), (*G)[j_minus1_in_preL]->getLabel()));
      if(DEBUG) {
        ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
        ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl; 
        ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
        ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
        ou << "forestdist[" << i - 1 << ", " << j - 1 << "] = " << forestdist[i - 1][j - 1] << endl;
      }
      if(da == forestdist[i][j]) {
        if(swap){
          map->setMap(-1, i_minus1_in_preL);
        } else {
          map->setMap(i_minus1_in_preL, -1);
        }
        i = i - 1;
      } else if(db == forestdist[i][j]) {
        if(swap) {
          map->setMap(j_minus1_in_preL, -1);
        } else {
          map->setMap(-1, j_minus1_in_preL);
        }
        j = j - 1;
      } else if(dc == forestdist[i][j]) {
        if(swap) {
          map->setMap(j_minus1_in_preL, i_minus1_in_preL);
        } else {
          map->setMap(i_minus1_in_preL, j_minus1_in_preL);
        }
        i = i - 1;
        j = j - 1;
      }
    }
    while(i != 0) {
      if(DEBUG) {
        ou << "i = " << i << " j = 0" << endl;
      }
      int i_minus1_in_preL = F->postR_to_preL[i - 1];
      if(swap){
        map->setMap(-1, i_minus1_in_preL);
      } else {
        map->setMap(i_minus1_in_preL, -1);
      }
      i = i - 1;
    }
    while(j != 0) {
      if(DEBUG) {
        ou << "i = 0 j = " << j << endl;
      }
      int j_minus1_in_preL = G->postR_to_preL[j - 1];
      if(swap) {
        map->setMap(j_minus1_in_preL, -1);
      } else {
        map->setMap(-1, j_minus1_in_preL);
      }
      j = j - 1;
    }

  }*/
  return dist;
}

/*
  the original one
*/
float TreeComparison::spfA(Node* a, Node* b, int leaf, int pathType, bool swap) {
	Tree *F, *G;
	if(swap) {
		F = B_;
		G = A_;
	} else {
		F = A_;
		G = B_;
	}
	int endF = a->getID(); 
	int endG = b->getID();
	int sizeF = a->getSubTreeSize();
	int sizeG = b->getSubTreeSize();
	int endF_in_preR = F->preL_to_preR[endF];
	int endG_in_preR = G->preL_to_preR[endG];
	int endPathNode = leaf;
  int endPathNode_in_preR = F->preL_to_preR[endPathNode];
	int startPathNode = -1;
	int lFFirst, lFLast, lF;
	int rFFirst, rFLast, rF;
	int lGFirst, lGLast;
	int rGFirst, rGLast;

  int FcurrentForestSize = 0;
  float FcurrentForestCost = 0;
  int FtmpForestSize = 0;
  float FtmpForestCost = 0;

  float dist = 0;
  int lF_prev = endPathNode;

	
	//loop A
	while(endPathNode >= endF) {
		endPathNode_in_preR = F->preL_to_preR[endPathNode];
		int startPathNode_in_preR = startPathNode == -1? 0x7fffffff : F->preL_to_preR[startPathNode];

		int parent_of_endPathNode_preL = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : (*F)[endPathNode]->getParent()->getID();
		int parent_of_endPathNode_preR = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : F->preL_to_preR[parent_of_endPathNode_preL];

		bool hasLeftPart;
		bool hasRightPart;

		

		if(startPathNode - endPathNode > 1) {
			hasLeftPart = true;
		} else {
			hasLeftPart = false;
		}
		if(startPathNode >= 0 && startPathNode_in_preR - endPathNode_in_preR > 1) {
			hasRightPart = true;
		} else {
			hasRightPart = false;
		}

    if(DEBUG) {
      ou << "spfA(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") " << swap << " " << endl;
    }

    if(DEBUG) {
      cout << "spfA(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") " << swap << " " << endl;
    }


    // right path left decomposition
		if(pathType == 1 || pathType == 2 && hasLeftPart) {
			// lFFirst and LFLast is important in this condition.
			// rGFirst and rGLast is to set in this stage.
			if(startPathNode == -1) {
				lFFirst = endPathNode;//the first node is the node on the path
				rFFirst = endPathNode_in_preR;// the first node is the node on the path
			} else {
				lFFirst = startPathNode - 1;//the first node is the node to the left of the path
				rFFirst = startPathNode_in_preR;//rFFirst set to the node on the path
			}
			if(!hasRightPart) {
				rFLast = endPathNode_in_preR;
			}
			lFLast = hasRightPart? endPathNode + 1 : endPathNode;

			rGLast = G->preL_to_preR[endG];
			rGFirst = (rGLast + sizeG) - 1; // get the leftmost child in G

			fn[fn_ft_length - 1] = -1;
			/*if(DEBUG) {
				ou << "initial fn and ft endG = " << to_string(endG) << " endG + sizeG = " << to_string(endG + sizeG) << endl;
			}*/
      for (int i = endG; i < endG + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }
      
      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
			//loop B
			for(int rG = rGFirst; rG >= rGLast; rG--) {
				/*if(DEBUG) {
					ou << "new Round B" << endl;
				}*/
				int rG_in_preL = (G)->preR_to_preL[rG];
				Node* parent = (*G)[rG_in_preL]->getParent();
				int parent_of_rG_in_preL = parent == NULL? 0x7fffffff : parent->getID();
				int parent_of_rG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent_of_rG_in_preL];
				lGFirst = G->preR_to_preL[rG];// lGFirst is set to rGFirst;
				
				int rGminus1_in_preL = rG <= endG_in_preR? 0x7fffffff : G->preR_to_preL[rG - 1];// rG should greater than endG_in_preR cause rG is the inner node of subtree enG
				int rGminus1_in_preR = rG <= endG_in_preR? 0x7fffffff : rG - 1;
				
				if (pathType == 1){
          if (lGFirst == endG || rGminus1_in_preL != parent_of_rG_in_preL) {// parent not exist or not the rightmost child
            lGLast = lGFirst;//lGLast is set to lGFirst
          } else {
            lGLast = parent_of_rG_in_preL + 1;//lGLast is set to the leftmost child of rG's parent
          }
        } else {
          lGLast = lGFirst == endG? lGFirst : endG + 1;//lGLast is set to the leftmost child of the whole tree
        }
			
        updateFnArray(G->preL_to_ln[lGFirst], lGFirst, endG); //stores the counter in D loop fn[ln] stores the start point
        updateFtArray(G->preL_to_ln[lGFirst], lGFirst); 
				
        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;
        //loop C
				for(int lF = lFFirst; lF >= lFLast; lF--) {
					/*if(DEBUG) {
						ou << "new round C" << endl;
					}*/
          rF = startPathNode_in_preR;
          if (lF == lFLast && !hasRightPart) {
              rF = rFLast;
          }
					int lG = lGFirst;
					int lF_in_preR = (F)->preL_to_preR[lF];
          counter++;

					if(DEBUG) {
          	ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
          	ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
					}

          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
					int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          int GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 

					bool fForestIsTree = lF_in_preR == rF;
					int lFSubtreeSize = (*F)[lF]->getSubTreeSize();
					bool lFIsLeftSiblingOfCurrentPathNode = lF + lFSubtreeSize == startPathNode;
					bool lFIsConsecutiveNodeOfCurrentPathNode = startPathNode - lF == 1;

					int case1SLeftIndex, case1SRightIndex;//S[lF + 1, lG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];
          			
          int case2SLeftIndex, case2SRightIndex;//S[lF, lG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;

          case1SLeftIndex =  lF + 1;//fixed
          case2SLeftIndex = lF;//fixed

          case1TRightIndex = rG;//fixed

          case3TRightIndex = rG;//fixed

          float case1 = 0, case2 = 0, case3 = 0;
          int case1_case, case2_case, case3_case;

          float minCost = 0;

          case1_case = 1;
          case3_case = 1;

          if (fForestIsTree) { // F_{lF,rF} is a tree.
            if (lFSubtreeSize == 1) { // F_{lF,rF} is a single node.
              // F_{lF, rF} - lF = null
              case1_case = 3;
            } else if (lFIsConsecutiveNodeOfCurrentPathNode) { // F_{lF,rF}-lF is the path node subtree.
              // F_{lF, rF} - rF = tree
              case1_case = 2;
            }
            
            case3 = 0;//F_{lF, rF} - F(lF) = null
            /*if(DEBUG) {
          		ou << "case3 = 0" << endl; 
          	}*/
            
            case3_case = 2;
          } else {
            if (lFIsConsecutiveNodeOfCurrentPathNode) {
              // F_{lF, rF} - lF = tree
              case1_case = 2;
            }
            case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[lF] : (F)->preL_to_sumDelCost[lF]); // USE COST MODEL - Delete F_{lF,rF}-F_lF.
            /*if(DEBUG) {
              ou << "case3_case FcurrentForest - F(lF) = " << endl;
              ou << to_string(FcurrentForestCost) << " - ";
              if(swap) ou << to_string((F)->preL_to_sumInsCost[lF]) << endl;
              else ou << to_string((F)->preL_to_sumDelCost[lF]) << endl;
            }*/
              			
        		if (lFIsLeftSiblingOfCurrentPathNode) {
              case3_case = 3;
            }
          }

          if(case3_case == 1) {
            case3SLeftIndex = lF + lFSubtreeSize;
          }

          switch(case1_case) {
            case 1: 
             	case1SRightIndex = lG;
              case1 = s[case1SLeftIndex][case1SRightIndex]; 
              /*if(DEBUG) {
             		ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = ";
                ou << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl;
             	}*/
             	break;
            
            case 2: 
              case1TLeftIndex = lG;
              case1 = t[case1TLeftIndex][case1TRightIndex]; 
              /*if(DEBUG) {
              	ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = ";
                ou << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl; 
              }*/
              break;
            
            case 3: 
              case1 = GcurrentForestCost; 
              /*if(DEBUG) {
              	ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << endl; 
              }*/
              break; // USE COST MODEL - Insert G_{lG,rG}.
          }

          case1 += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
          minCost = case1;

          if (GcurrentForestSize == 1) { // G_{lG,rG} is a single node.
            case2 = FcurrentForestCost; // USE COST MODEL - Delete F_{lF,rF}.
            /*if(DEBUG) {
              ou << "case2_case1 FcurrentForestCost = " << to_string(FcurrentForestCost) << endl;
            }*/
          } else { // G_{lG,rG} is a tree.
            case2 = q[lF];
            /*if(DEBUG) {
            	ou << "case2_case2 q[" << to_string(lF) << "] = " << to_string(q[lF]) << endl;
            }*/
          }
          case2 += (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.del((*G)[lG]->getLabel()));
          /*if(DEBUG) {
            ou << "case2 += ";
            if(swap) ou << "delete " << (*G)[lG]->getLabel() << endl;
            else ou << "insert " << (*G)[lG]->getLabel() << endl;
          }*/
          if(case2 < minCost) minCost = case2;

          if (case3 < minCost) {
            case3 += swap? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
            	if(swap) ou << "case3_case3 delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
            	else ou << "case3_case3 delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
            }*/
            if(case3 < minCost) {
            	case3 += swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
              /*if(DEBUG) {
                if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }*/
            } 
            if(case3 < minCost) {
            	minCost = case3;
            }
          }
          if(DEBUG) {
            ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
          }
          /*if(DEBUG) {
            ou << "case1 = " << case1 << endl;
            ou << "case2 = " << case2 << endl;
            ou << "case3 = " << case3 << endl;
          }*/
          dist = minCost;
          s[lF][lG] = minCost;
					lG = ft[lG];
					
          //loop D
					while (lG >= lGLast) {
            counter++;
						if(DEBUG) {
							ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
							ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
						}

						GcurrentForestSize++;
						GcurrentForestCost += (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
						minCost = 0;

						switch(case1_case) {
              case 1:
                case1SRightIndex = lG;
                case1 = s[case1SLeftIndex][case1SRightIndex] + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                	ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 2: 
                case1TLeftIndex = lG;
                case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                	ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 3: 
                case1 = GcurrentForestCost + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
                /*if(DEBUG) {
                	ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Insert G_{lG,rG} and elete lF, leftmost root node in F_{lF,rF}.
            }
            minCost = case1;

            case2SRightIndex = fn[lG];
            case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
            /*if(DEBUG) {
              ou << "case2 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
              if(swap) ou << "delete " << (*G)[lG]->getLabel() << endl;
              else ou << "insert " << (*G)[lG]->getLabel() << endl;
            }*/
            
            if(case2 < minCost) {
              minCost = case2;
            }


            case3 = swap ? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
              if(swap) {
              	ou << "case3 = delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
              } else {
              	ou << "case3 = delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
              }
            }*/
            
            if (case3 < minCost) {
              switch(case3_case) {
                case 1: 
                  case3SRightIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += s[case3SLeftIndex][case3SRightIndex]; 
                  /*if(DEBUG) {
                    ou << "case3 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl;
                  }*/
                  break;
                
                case 2: 
                  case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 
                  /*if(DEBUG) {
                    ou << "case3 += " << "GcurrentForestCost - G(lG) = ";
                    ou << to_string(GcurrentForestCost) << " - ";
                    if(swap) ou << (G)->preL_to_sumDelCost[lG] << endl;
                    else ou << (G)->preL_to_sumInsCost[lG] << endl;
                  }*/
                  break; // USE COST MODEL - Insert G_{lG,rG}-G_lG.
                
                case 3: 
                  case3TLeftIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += t[case3TLeftIndex][case3TRightIndex]; 
                	/*if(DEBUG) {
                		ou << "case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl;
                	}*/

                  break;
              }
              
              if (case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel())); // USE COST MODEL - Rename the leftmost root nodes in F_{lF,rF} and G_{lG,rG}.
                /*if(DEBUG) {
                  if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                  else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
                }*/
                if (case3 < minCost) {
                  minCost = case3;
                }
              }
            }
            if(DEBUG) {
              ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
            }
            s[lF][lG] = minCost;
            dist = minCost;
						lG = ft[lG];
					}
					lF_prev = lF;
				}
				/*if(DEBUG) {
					ou << "rGminus1_in_preR = " << to_string(rGminus1_in_preR) << " rG = " << to_string(rG) << " parent_of_rG_in_preL = " << to_string(parent_of_rG_in_preL) << " parent_of_rG_in_preR = " << to_string(parent_of_rG_in_preR) << endl;
				}*/
				if(rGminus1_in_preR == parent_of_rG_in_preR && rGminus1_in_preR != 0x7fffffff) {
					if (!hasRightPart) {
            if (hasLeftPart) {
              if(swap) {
                delta[parent_of_rG_in_preL][endPathNode] = s[lFLast + 1][rGminus1_in_preL + 1];
              	if(DEBUG) {
              		ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
              	}
              } else {
                delta[endPathNode][parent_of_rG_in_preL] = s[lFLast + 1][rGminus1_in_preL + 1];
              	if(DEBUG) {
              		ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
              	}
              }
            }

            if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {
              if (swap) {
                delta[parent_of_rG_in_preL][parent_of_endPathNode_preL] = s[lFLast][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              } else {
                delta[parent_of_endPathNode_preL][parent_of_rG_in_preL] = s[lFLast][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              }
            }
          }

          for (int lF = lFFirst; lF >= lFLast; lF--) {
            q[lF] = s[lF][parent_of_rG_in_preL + 1];
            if(DEBUG) {
              ou << "q[" << to_string(lF) << "] = " << "s[" << to_string(lF) << ", " << to_string(parent_of_rG_in_preL + 1) << "]" << endl;
            }
          }
			  }

			  for (int lG = lGFirst; lG >= lGLast; lG = ft[lG]) {
          t[lG][rG] = s[lFLast][lG];
				  if(DEBUG) {
					 ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(lG) << "]" << endl;
				  }
        }
		  }
	  }

    // left path right decomposit
	  if (pathType == 0 || pathType == 2 && hasRightPart || pathType == 2 && !hasLeftPart && !hasRightPart) {
		  if (startPathNode == -1) {
        lFFirst = endPathNode;
        rFFirst = F->preL_to_preR[endPathNode];
      } else {
        rFFirst = F->preL_to_preR[startPathNode] - 1;//the node right to the node on the path
        lFFirst = endPathNode + 1;//lFirst is set to the node on the path
      }

      lFLast = endPathNode;
      rFLast = F->preL_to_preR[endPathNode];

      lGLast = endG;
      lGFirst = (lGLast + sizeG) - 1;

      fn[fn_ft_length - 1] = -1;
  	  /*if(DEBUG) {
			   ou << "initial fn and ft endG_in_preR = " << to_string(endG_in_preR) << " endG_in_preR + sizeG = " << to_string(endG_in_preR + sizeG) << endl;
		  }*/
      
      for (int i = endG_in_preR; i < endG_in_preR + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }

      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
      //loop B'
      for (int lG = lGFirst; lG >= lGLast; lG--) {
        /*if(DEBUG) {
				  ou << "new Round B'" << endl;
			  }*/
        Node* parent = (*G)[lG]->getParent();
        int parent_of_lG_in_preL = parent == NULL? 0x7fffffff: parent->getID();
        int parent_of_lG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent->getID()];// not exist -1;
			  rGFirst = G->preL_to_preR[lG];
			  int lG_in_preR = G->preL_to_preR[lG];

			  int lGminus1_in_preL = lG <= endG? 0x7fffffff : lG - 1;
			  int lGminus1_in_preR = lG <= endG? 0x7fffffff : G->preL_to_preR[lG - 1];

			  if (pathType == 0) {
          if (lG == endG || lGminus1_in_preL != parent_of_lG_in_preL) {//parent not exists or not the leftmost child.
            rGLast = rGFirst;
          } else {
            rGLast = parent_of_lG_in_preR + 1;
          }
        } else {// left and right
          rGLast = rGFirst == endG_in_preR ? rGFirst : endG_in_preR;
        }

        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;

			/*if(DEBUG) {
				ou << "updateFnArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ", " << to_string(endG_in_preR) << ")" << endl;
			}*/
			  updateFnArray(G->preR_to_ln[rGFirst], rGFirst, endG_in_preR);
			
			/*if(DEBUG) {
				ou << "updateFtArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ")" << endl;
			}*/
        updateFtArray(G->preR_to_ln[rGFirst], rGFirst);
          	
      /*if(DEBUG) {
       ou << "endG_in_preR = " << to_string(endG_in_preR) << endl;
       ou << "start from rG = " << to_string(rGFirst) << endl;
       ou << "FN" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << fn[i] << " ";
      }
      ou << endl;
      ou << "FT" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << ft[i] << " ";
      }
      ou << endl;
      }*/
        lF = lF_prev;
      // loop C'
        for(int rF = rFFirst; rF >= rFLast; rF--) {
          /*if(DEBUG) {
            ou << "new Round C'" << endl;
          }*/
          //lF = startPathNode;
          if (rF == rFLast) {
              lF = lFLast;
          }
          int rG = rGFirst;
          int rG_in_preL = (G)->preR_to_preL[rG];
          			
          if(rF == rFLast) lF = F->preR_to_preL[rFLast]; 
          counter++;
          if(DEBUG) {
            ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter  = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
				  }

          int rF_in_preL = (F)->preR_to_preL[rF];
          			
          bool FForestIsTree = lF == rF_in_preL;
          int rFSubtreeSize = (*F)[rF_in_preL]->getSubTreeSize();
          			  			
          int case1SLeftIndex, case1SRightIndex;//S[rF + 1, rG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];

          int case2SLeftIndex, case2SRightIndex;//S[rF, rG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;
        
          float case1 = 0, case2 = 0, case3 = 0;
        	int case1_case, case2_case, case3_case;

          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
         	int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          float GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); // USE COST MODEL - reset to subtree insertion cost.

          float minCost = 0;


          case1SLeftIndex = rF + 1;//fixed
          case1SRightIndex = rG;

          case1TLeftIndex = lG;//fixed
          case1TRightIndex = rG;

          case2SLeftIndex = rF;//fixed

          case3TLeftIndex = lG;//fixed


          bool rFIsConsecutiveNodeOfCurrentPathNode;
          bool rFIsRightSiblingOfCurrentPathNode;

          case1_case = 1;
          case3_case = 1;//otherwise

          if (startPathNode > 0) {
            rFIsConsecutiveNodeOfCurrentPathNode = startPathNode_in_preR - rF == 1;
            rFIsRightSiblingOfCurrentPathNode = rF + rFSubtreeSize == startPathNode_in_preR;
          } else {
            rFIsConsecutiveNodeOfCurrentPathNode = false;//consecutiveNode use T;
            rFIsRightSiblingOfCurrentPathNode = false;//
          }

          if(FForestIsTree) {
            if(rFSubtreeSize == 1) {
          	  // F_{lF, rF} - rF = null
          	  //case1 = GcurrentForestCost;//sumG to be computed 
          	  case1_case = 3;
            } else if(rFIsConsecutiveNodeOfCurrentPathNode) {
          	  // F_{lF, rF} - rF = tree
          	  //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]; to be computed
          	  case1_case = 2;
            }
            case3 = 0;
            /*if(DEBUG) {
          	  ou << "case3 = 0" << endl; 
            }*/
            case3_case = 2;// F_{lF, rF} - F(rF) = null
            } else {
              if (rFIsConsecutiveNodeOfCurrentPathNode) {// F_{lF, rF} - rF = the subforest to the left of the path
          	   //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]
          	   case1_case = 2;
              } else {//otherwise
          	    //case1 = s[case1SLeftIndex][case1SRightIndex];//S[rF + 1, rG];// have calculate
          	    case1_case = 1;
              }
              case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[rF_in_preL] : (F)->preL_to_sumDelCost[rF_in_preL]);// the first case in G should be G_{lG, rG} - l(rG) = null // F_{lF, rF} - F(rF), G_{lG, rG} - G(rG)
              /*if(DEBUG) {
                ou << "case3_case FcurrentForest - F(rF)" << endl;
              }*/
              if (rFIsRightSiblingOfCurrentPathNode) {
          	    case3_case = 3; // use T
              }
            }

            if (case3_case == 1) {
              case3SLeftIndex = rF + rFSubtreeSize;//delete the whole rightmost tree//otherwise
            }
          
            switch(case1_case) {
              case 1:
                case1 = s[case1SLeftIndex][case1SRightIndex];
                /*if(DEBUG) {
            	    ou << "case1_case1 = s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl; 
                }*/
              break;
          
              case 2:
                case1 = t[case1TLeftIndex][case1TRightIndex];
                /*if(DEBUG) {
            	    ou << "case1_case2 = t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl;
                } */
                break;
              case 3:
                case1 = GcurrentForestCost;
                /*if(DEBUG) {
            	    ou << "case1_case3 = GcurrentForestCost = " << GcurrentForestCost << endl;
                }*/
                break;  				
            }
            case1 += (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel()));
            /*if(DEBUG) {
              ou << "case1 += ";
              if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
              else ou << "delete " << (*F)[rF]->getLabel() << endl;
            }*/
            minCost = case1;
          

            if (GcurrentForestSize == 1) {// the first case in G should be a node or a tree
              case2 = FcurrentForestCost;
              /*if(DEBUG) {
                ou << "case2_case1 = " << to_string(FcurrentForestCost) << endl;
              }*/
            } else {
              case2 = q[rF];
              /*if(DEBUG) {
                ou << "case2_case2 = q[" << to_string(rF) << "] = " << to_string(q[rF]) << endl;
              }*/
            }
            case2 += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
            /*if(DEBUG) {
              ou << "case2 += ";
              if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
              else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
            }*/
            if(case2 < minCost) {
              minCost = case2;
            }
        
            if(case3 < minCost) { 
              case3 += swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];// F(rF) - rF
              /*if(DEBUG) {
                if(swap) ou << "case3_case3 delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
                else ou << "case3_case3 delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) << endl;
              }*/
          
              if(case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
                //case3 += costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());
                /*if(DEBUG) {
                  if(swap) ou << "case3_case3 += " << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel())) << endl;
                  else ou << "case3_case3 += " << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel())) << endl;
                }*/
              }
              if(case3 < minCost) {
                minCost = case3;
              }
            }
        /*if(DEBUG) {
          ou << "case1 = " << to_string(case1) << endl;
          ou << "case2 = " << to_string(case2) << endl;
          ou << "case3 = " << to_string(case3) << endl;
        }*/
        s[rF][rG] = minCost;
        dist = minCost;
        if(DEBUG) {
           ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(s[rF][rG]) << endl;
        }
        rG = ft[rG];	
        
        // loop D'
        while(rG >= rGLast) {// every G is a subforest not a subtree
          counter++;
          rG_in_preL = (G)->preR_to_preL[rG];
        	if(DEBUG) {
          	ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
          	ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
					}

					GcurrentForestSize++;
					GcurrentForestCost += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));

					switch (case1_case) {
            case 1:
              case1SRightIndex = rG;
              case1 = s[case1SLeftIndex][case1SRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/ 
              break; 
            case 2: 
              case1TRightIndex = rG;
              case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
            case 3: 
              case1 = GcurrentForestCost + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case3 GcurrentForestCost = " << GcurrentForestCost << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
          }
          
          minCost = case1;

          case2SRightIndex = fn[rG];
          case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));//G is not a tree or a node for sure in D loop
          /*if(DEBUG) {
            ou << "case2_case3 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
            if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
            else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
          }*/
          
          
          if(case2 < minCost) {
            minCost = case2;
          }

          case3 = swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];//F_{rF} - rF, G_{rG} - rG
          /*if(DEBUG) {
            if(swap) {
              ou << "case3_case delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
            } else {
              ou << "case3_case delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) <<  endl;
            }
          }*/
          
          if(case3 < minCost) {
            switch(case3_case) {
              case 1: 
              	case3SRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
              	case3 += s[case3SLeftIndex][case3SRightIndex];
              	/*if(DEBUG) {
              	  ou << "case3_case1 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl; 
              	}*/
              	break;
              
              case 2: 
              	case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[rG_in_preL] : (G)->preL_to_sumInsCost[rG_in_preL]);
              	/*if(DEBUG) {
              		ou << "case3_case2 += " << "GcurrentForestCost - G(rG) = " << to_string(GcurrentForestCost) << " - ";
                  if(swap) ou << (G)->preL_to_sumDelCost[rG_in_preL] << endl;
                  else ou << (G)->preL_to_sumInsCost[rG_in_preL] << endl;
              	}*/
              	break;
              
              case 3: 
              	case3TRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
              	case3 += t[case3TLeftIndex][case3TRightIndex];
              	/*if(DEBUG) {
              		ou << "case3_case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl; 
              	}*/
              	break;
            }
              
            if(case3 < minCost) {
              case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
              /*if(DEBUG) {
                ou << "case3 += ";
                if(swap) ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
                else ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }*/
              if(case3 < minCost) {
              	minCost = case3;
              }
            }
          }
          /*if(DEBUG) {
            ou << "case1 = " << to_string(case1) << endl;
            ou << "case2 = " << to_string(case2) << endl;
            ou << "case3 = " << to_string(case3) << endl; 
          }*/
            if(DEBUG) {
              ou << "s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(minCost) << endl;
            }
            s[rF][rG] = minCost;
            dist = minCost;
            rG = ft[rG];
          }
        }

        if(lGminus1_in_preL == parent_of_lG_in_preL && lGminus1_in_preL != 0x7fffffff) { // lG is the leftmost child of its parent
    
          if(hasRightPart) {
            if(swap) {
              delta[parent_of_lG_in_preL][endPathNode] = s[rFLast + 1][lGminus1_in_preR + 1];
              if(DEBUG) {
          	   ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
               ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
              }
            } else {
              delta[endPathNode][parent_of_lG_in_preL] = s[rFLast + 1][lGminus1_in_preR + 1];
              if(DEBUG) {
                ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
              }
            }
          }
      
          if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {//no left and right
            if(swap) {
              delta[parent_of_lG_in_preL][parent_of_endPathNode_preL] = s[rFLast][lGminus1_in_preR + 1];
              if(DEBUG) {
            	 ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
               ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
              }
            } else {
              delta[parent_of_endPathNode_preL][parent_of_lG_in_preL] = s[rFLast][lGminus1_in_preR + 1];
              if(DEBUG) {
          	    ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
              }
            }
          }

          for (int rF = rFFirst; rF >= rFLast; rF--) {
            q[rF] = s[rF][parent_of_lG_in_preR + 1];
            if(DEBUG) {
              ou << "q[" << to_string(rF) << "] = " << "s[" << to_string(rF) << ", " << to_string(parent_of_lG_in_preR + 1) << "] = ";
              ou << to_string(s[rF][parent_of_lG_in_preR + 1]) << endl;	
            }
          }
        }
        for (int rG = rGFirst; rG >= rGLast; rG = ft[rG]) {
          t[lG][rG] = s[rFLast][rG];
          if(DEBUG) {
            ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(rG) << "] = ";
            ou << to_string(s[rFLast][rG]) << endl;
          }
        }
      }
    }
    rF = endPathNode_in_preR;//in D' loop
    startPathNode = endPathNode;
    endPathNode = (*F)[endPathNode] ->getParent() == NULL? -1 : (*F)[endPathNode] ->getParent()->getID();	
    endPathNode_in_preR = F->preL_to_preR[endPathNode];
  }
  return dist;
};



/*
  first add node to the left of the path then the right node -> first right decompose then left decompse
*/
float TreeComparison::spfA_LR(Node* a, Node* b, int leaf, int pathType, float*** forestdist, bool swap, bool mark) {
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }
  int endF = a->getID(); 
  int endG = b->getID();
  int sizeF = a->getSubTreeSize();
  int sizeG = b->getSubTreeSize();
  int endF_in_preR = F->preL_to_preR[endF];
  int endG_in_preR = G->preL_to_preR[endG];
  int endPathNode = leaf;
  int endPathNode_in_preR = F->preL_to_preR[endPathNode];
  int startPathNode = -1;
  int lFFirst, lFLast, lF;
  int rFFirst, rFLast, rF;
  int lGFirst, lGLast;
  int rGFirst, rGLast;

  int FcurrentForestSize = 0;
  float FcurrentForestCost = 0;
  int FtmpForestSize = 0;
  float FtmpForestCost = 0;

  float dist = 0;

  int lF_prev = endPathNode;

  if(DEBUG) {
    ou << "spfA_LR(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") " << swap << " " << endl;
  }

  if(DEBUG) {
    cout << "spfA_LR(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") counter = " << counter << " " << endl;
  }
  
  //loop A
  while(endPathNode >= endF) {
    endPathNode_in_preR = F->preL_to_preR[endPathNode];
    int startPathNode_in_preR = startPathNode == -1? 0x7fffffff : F->preL_to_preR[startPathNode];

    int parent_of_endPathNode_preL = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : (*F)[endPathNode]->getParent()->getID();
    int parent_of_endPathNode_preR = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : F->preL_to_preR[parent_of_endPathNode_preL];

    bool hasLeftPart;
    bool hasRightPart;

    

    if(startPathNode - endPathNode > 1) {
      hasLeftPart = true;
    } else {
      hasLeftPart = false;
    }
    if(startPathNode >= 0 && startPathNode_in_preR - endPathNode_in_preR > 1) {
      hasRightPart = true;
    } else {
      hasRightPart = false;
    }

    // right path left decomposition
    if(pathType == 1 || pathType == 2 && hasLeftPart) {
      // lFFirst and LFLast is important in this condition.
      // rGFirst and rGLast is to set in this stage.
      if(startPathNode == -1) {
        lFFirst = endPathNode;//the first node is the node on the path
        rFFirst = endPathNode_in_preR;// the first node is the node on the path
      } else {
        lFFirst = startPathNode - 1;//the first node is the node to the left of the path
        rFFirst = startPathNode_in_preR;//rFFirst set to the node on the path
      }
      if(!hasRightPart) {
        rFLast = endPathNode_in_preR;
      }
      lFLast = hasRightPart? endPathNode + 1 : endPathNode;

      rGLast = G->preL_to_preR[endG];
      rGFirst = (rGLast + sizeG) - 1; // get the leftmost child in G

      fn[fn_ft_length - 1] = -1;
      /*if(DEBUG) {
        ou << "initial fn and ft endG = " << to_string(endG) << " endG + sizeG = " << to_string(endG + sizeG) << endl;
      }*/
      for (int i = endG; i < endG + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }
      
      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
      //loop B
      for(int rG = rGFirst; rG >= rGLast; rG--) {
        /*if(DEBUG) {
          ou << "new Round B" << endl;
        }*/
        int rG_in_preL = (G)->preR_to_preL[rG];
        Node* parent = (*G)[rG_in_preL]->getParent();
        int parent_of_rG_in_preL = parent == NULL? 0x7fffffff : parent->getID();
        int parent_of_rG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent_of_rG_in_preL];
        lGFirst = G->preR_to_preL[rG];// lGFirst is set to rGFirst;
        
        int rGminus1_in_preL = rG <= endG_in_preR? 0x7fffffff : G->preR_to_preL[rG - 1];// rG should greater than endG_in_preR cause rG is the inner node of subtree enG
        int rGminus1_in_preR = rG <= endG_in_preR? 0x7fffffff : rG - 1;
        
        if (pathType == 1){
          if (lGFirst == endG || rGminus1_in_preL != parent_of_rG_in_preL) {// parent not exist or not the rightmost child
            lGLast = lGFirst;//lGLast is set to lGFirst
          } else {
            lGLast = parent_of_rG_in_preL + 1;//lGLast is set to the leftmost child of rG's parent
          }
        } else {
          lGLast = lGFirst == endG? lGFirst : endG + 1;//lGLast is set to the leftmost child of the whole tree
        }
      
        updateFnArray(G->preL_to_ln[lGFirst], lGFirst, endG); //stores the counter in D loop fn[ln] stores the start point
        updateFtArray(G->preL_to_ln[lGFirst], lGFirst); 
        
        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;
        //loop C
        for(int lF = lFFirst; lF >= lFLast; lF--) {
          /*if(DEBUG) {
            ou << "new round C" << endl;
          }*/
          rF = startPathNode_in_preR;
          if (lF == lFLast && !hasRightPart) {
              rF = rFLast;
          }
          int lG = lGFirst;
          int lF_in_preR = (F)->preL_to_preR[lF];
          counter++;

          if(DEBUG && mark) {
            ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
          }

          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
          if(mark) {
            for(int i = 0; i <= (*G)[endG]->getSubTreeSize(); i++) {
              forestdist[lF_in_preR][i][(*G)[endG]->getSubTreeSize()] = FcurrentForestCost;
              forestdist[lF_in_preR][(*G)[endG]->getSubTreeSize()][i] = FcurrentForestCost;
            }   
            if(DEBUG) {
              ou << "forestdist[" << lF_in_preR << ", " << (*G)[endG]->getSubTreeSize() << ", " << (*G)[endG]->getSubTreeSize() << ") = " << forestdist[lF_in_preR][(*G)[endG]->getSubTreeSize()][(*G)[endG]->getSubTreeSize()] << endl;
            }
          }
          int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          int GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
            if(DEBUG) {
              ou << "forestdist[" << (*F)[endF]->getSubTreeSize() << ", " << lG << ", " << rG << "] = " << forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] << endl;
            }
          }

          bool fForestIsTree = lF_in_preR == rF;
          int lFSubtreeSize = (*F)[lF]->getSubTreeSize();
          bool lFIsLeftSiblingOfCurrentPathNode = lF + lFSubtreeSize == startPathNode;
          bool lFIsConsecutiveNodeOfCurrentPathNode = startPathNode - lF == 1;

          int case1SLeftIndex, case1SRightIndex;//S[lF + 1, lG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];
                
          int case2SLeftIndex, case2SRightIndex;//S[lF, lG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;

          case1SLeftIndex =  lF + 1;//fixed
          case2SLeftIndex = lF;//fixed

          case1TRightIndex = rG;//fixed

          case3TRightIndex = rG;//fixed

          float case1 = 0, case2 = 0, case3 = 0;
          int case1_case, case2_case, case3_case;

          float minCost = 0;

          case1_case = 1;
          case3_case = 1;

          if (fForestIsTree) { // F_{lF,rF} is a tree.
            if (lFSubtreeSize == 1) { // F_{lF,rF} is a single node.
              // F_{lF, rF} - lF = null
              case1_case = 3;
            } else if (lFIsConsecutiveNodeOfCurrentPathNode) { // F_{lF,rF}-lF is the path node subtree.
              // F_{lF, rF} - rF = tree
              case1_case = 2;
            }
            
            case3 = 0;//F_{lF, rF} - F(lF) = null
            /*if(DEBUG) {
              ou << "case3 = 0" << endl; 
            }*/
            
            case3_case = 2;
          } else {
            if (lFIsConsecutiveNodeOfCurrentPathNode) {
              // F_{lF, rF} - lF = tree
              case1_case = 2;
            }
            case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[lF] : (F)->preL_to_sumDelCost[lF]); // USE COST MODEL - Delete F_{lF,rF}-F_lF.
            /*if(DEBUG) {
              ou << "case3_case FcurrentForest - F(lF) = " << endl;
              ou << to_string(FcurrentForestCost) << " - ";
              if(swap) ou << to_string((F)->preL_to_sumInsCost[lF]) << endl;
              else ou << to_string((F)->preL_to_sumDelCost[lF]) << endl;
            }*/
                    
            if (lFIsLeftSiblingOfCurrentPathNode) {
              case3_case = 3;
            }
          }

          if(case3_case == 1) {
            case3SLeftIndex = lF + lFSubtreeSize;
          }

          switch(case1_case) {
            case 1: 
              case1SRightIndex = lG;
              case1 = s[case1SLeftIndex][case1SRightIndex]; 
              /*if(DEBUG) {
                ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = ";
                ou << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl;
              }*/
              break;
            
            case 2: 
              case1TLeftIndex = lG;
              case1 = t[case1TLeftIndex][case1TRightIndex]; 
              /*if(DEBUG) {
                ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = ";
                ou << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl; 
              }*/
              break;
            
            case 3: 
              case1 = GcurrentForestCost; 
              /*if(DEBUG) {
                ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << endl; 
              }*/
              break; // USE COST MODEL - Insert G_{lG,rG}.
          }

          case1 += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
          minCost = case1;
          
          if (GcurrentForestSize == 1) { // G_{lG,rG} is a single node.
            case2 = FcurrentForestCost; // USE COST MODEL - Delete F_{lF,rF}.
            /*if(DEBUG) {
              ou << "case2_case1 FcurrentForestCost = " << to_string(FcurrentForestCost) << endl;
            }*/
          } else { // G_{lG,rG} is a tree.
            case2 = q[lF];
            /*if(DEBUG) {
              ou << "case2_case2 q[" << to_string(lF) << "] = " << to_string(q[lF]) << endl;
            }*/
          }
          case2 += (swap ? costModel_.ins((*G)[lG]->getLabel()) : costModel_.del((*G)[lG]->getLabel()));
          if(case2 < minCost) {
            minCost = case2;
          }

          if (case3 < minCost) {
            case3 += swap? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
              if(swap) ou << "case3_case3 delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
              else ou << "case3_case3 delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
            }*/
            if(case3 < minCost) {
              case3 += swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
              /*if(DEBUG) {
                if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }*/
            } 
            if(case3 < minCost) {
              minCost = case3;
            }
          }
          if(DEBUG) {
            ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
          }
          dist = minCost;
          s[lF][lG] = minCost;
          if(mark) {
            forestdist[lF_in_preR][lG][rG] = minCost;
            if(DEBUG) {
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }
          }
          /*if(mark) {
            if(case1 == minCost) {
              if(swap) {
                if(DEBUG) {
                  ou << "- -> " << (*F)[lF]->getLabel() << endl;
                }
                map->setMap(-1, lF);
              } else {
                if(DEBUG) {
                  ou << (*F)[lF]->getLabel() << " -> " << "-" << endl;
                }
                map->setMap(lF, -1);
              }
            } else if(case2 == minCost) {
              if(swap) {
                if(DEBUG) {
                  ou << (*G)[lG]->getLabel() << " -> " << "-" << endl;
                }
                map->setMap(lG, -1);
              } else {
                if(DEBUG) {
                  ou << "- -> " << (*G)[lG]->getLabel() << endl;
                }
                map->setMap(-1, lG);
              }
            } else if(case3 == minCost) {
              if(swap) {
                if(DEBUG) {
                  ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                } 
                map->setMap(lG, lF);
              } else {
                if(DEBUG) {
                  ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
                } 
                map->setMap(lF, lG);
              }
            }
          }*/
          lG = ft[lG];
          
          //loop D
          while (lG >= lGLast) {
            counter++;
            if(DEBUG && mark) {
              ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
              ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
            }

            GcurrentForestSize++;
            GcurrentForestCost += (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
            if(mark) {
              forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
              if(DEBUG) {
                ou << "forestdist[" << (*F)[endF]->getSubTreeSize() << ", " << lG << ", " << rG << "] = " << forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] << endl;
              }
            }
            minCost = 0;

            switch(case1_case) {
              case 1:
                case1SRightIndex = lG;
                case1 = s[case1SLeftIndex][case1SRightIndex] + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                  ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 2: 
                case1TLeftIndex = lG;
                case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                  ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 3: 
                case1 = GcurrentForestCost + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
                /*if(DEBUG) {
                  ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Insert G_{lG,rG} and elete lF, leftmost root node in F_{lF,rF}.
            }
            minCost = case1;
            /*if(DEBUG) {
              if(swap) ou << "G: " << "- ->" << (*F)[lF]->getLabel() << endl;
              else ou << "F: " << (*F)[lF]->getLabel() << " -> " << "-" << endl;
            }*/

            case2SRightIndex = fn[lG];
            case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
            /*if(DEBUG) {
              ou << "case2 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
              if(swap) ou << "delete " << (*G)[lG]->getLabel() << endl;
              else ou << "insert " << (*G)[lG]->getLabel() << endl;
            }*/
            
            if(case2 < minCost) {
              minCost = case2;
              /*if(DEBUG) {
                if(swap) ou << "F: " << (*G)[lG]->getLabel() << " -> " << "-" << endl;
                else ou << "G: " << "- ->" << (*G)[lG]->getLabel() << endl;
              }*/
            }


            case3 = swap ? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
              if(swap) {
                ou << "case3 = delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
              } else {
                ou << "case3 = delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
              }
            }*/
            
            if (case3 < minCost) {
              switch(case3_case) {
                case 1: 
                  case3SRightIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += s[case3SLeftIndex][case3SRightIndex]; 
                  /*if(DEBUG) {
                    ou << "case3 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl;
                  }*/
                  break;
                
                case 2: 
                  case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 
                  /*if(DEBUG) {
                    ou << "case3 += " << "GcurrentForestCost - G(lG) = ";
                    ou << to_string(GcurrentForestCost) << " - ";
                    if(swap) ou << (G)->preL_to_sumDelCost[lG] << endl;
                    else ou << (G)->preL_to_sumInsCost[lG] << endl;
                  }*/
                  break; // USE COST MODEL - Insert G_{lG,rG}-G_lG.
                
                case 3: 
                  case3TLeftIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += t[case3TLeftIndex][case3TRightIndex]; 
                  /*if(DEBUG) {
                    ou << "case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl;
                  }*/

                  break;
              }
              
              if (case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel())); // USE COST MODEL - Rename the leftmost root nodes in F_{lF,rF} and G_{lG,rG}.
                /*if(DEBUG) {
                  if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                  else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
                }*/
                if (case3 < minCost) {
                  minCost = case3;
                }
              }
            }
            if(DEBUG) {
              ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
            }
            s[lF][lG] = minCost;
            if(mark) {
              forestdist[lF_in_preR][lG][rG] = minCost;
              if(DEBUG) {
                ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
              }
            }
            /*if(mark) {
              if(case1 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << "- -> " << (*F)[lF]->getLabel() << endl;
                  }
                  map->setMap(-1, lF);
                } else {
                  if(DEBUG) {
                    ou << (*F)[lF]->getLabel() << " -> " << "-" << endl;
                  }
                  map->setMap(lF, -1);
                }
              } else if(case2 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << (*G)[lG]->getLabel() << " -> " << "-" << endl;
                  }
                  map->setMap(lG, -1);
                } else {
                  if(DEBUG) {
                    ou << "- -> " << (*G)[lG]->getLabel() << endl;
                  }
                  map->setMap(-1, lG);
                }
              } else if(case3 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                  }
                  map->setMap(lG, lF);
                } else {
                  if(DEBUG) {
                    ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
                  }
                  map->setMap(lF, lG);
                }
              }
            }*/
            dist = minCost;
            lG = ft[lG];
          }
          lF_prev = lF;
        }
        /*if(DEBUG) {
          ou << "rGminus1_in_preR = " << to_string(rGminus1_in_preR) << " rG = " << to_string(rG) << " parent_of_rG_in_preL = " << to_string(parent_of_rG_in_preL) << " parent_of_rG_in_preR = " << to_string(parent_of_rG_in_preR) << endl;
        }*/
        if(rGminus1_in_preR == parent_of_rG_in_preR && rGminus1_in_preR != 0x7fffffff) {
          if (!hasRightPart) {
            if (hasLeftPart) {
              if(swap) {
                delta[parent_of_rG_in_preL][endPathNode] = s[lFLast + 1][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              } else {
                delta[endPathNode][parent_of_rG_in_preL] = s[lFLast + 1][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              }
            }

            if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {
              if (swap) {
                delta[parent_of_rG_in_preL][parent_of_endPathNode_preL] = s[lFLast][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              } else {
                delta[parent_of_endPathNode_preL][parent_of_rG_in_preL] = s[lFLast][rGminus1_in_preL + 1];
                if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }
              }
            }
          }

          for (int lF = lFFirst; lF >= lFLast; lF--) {
            q[lF] = s[lF][parent_of_rG_in_preL + 1];
            if(DEBUG) {
              ou << "q[" << to_string(lF) << "] = " << "s[" << to_string(lF) << ", " << to_string(parent_of_rG_in_preL + 1) << "]" << endl;
            }
          }
        }

        for (int lG = lGFirst; lG >= lGLast; lG = ft[lG]) {
          t[lG][rG] = s[lFLast][lG];
          if(DEBUG) {
           ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(lG) << "]" << endl;
          }
        }
      }
    }

    // left path right decompose
    if (pathType == 0 || pathType == 2 && hasRightPart || pathType == 2 && !hasLeftPart && !hasRightPart) {
      if (startPathNode == -1) {
        lFFirst = endPathNode;
        rFFirst = F->preL_to_preR[endPathNode];
      } else {
        rFFirst = F->preL_to_preR[startPathNode] - 1;//the node right to the node on the path
        lFFirst = endPathNode + 1;//lFirst is set to the node on the path
      }

      lFLast = endPathNode;
      rFLast = F->preL_to_preR[endPathNode];

      lGLast = endG;
      lGFirst = (lGLast + sizeG) - 1;

      fn[fn_ft_length - 1] = -1;
      /*if(DEBUG) {
         ou << "initial fn and ft endG_in_preR = " << to_string(endG_in_preR) << " endG_in_preR + sizeG = " << to_string(endG_in_preR + sizeG) << endl;
      }*/
      
      for (int i = endG_in_preR; i < endG_in_preR + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }

      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
      //loop B'
      for (int lG = lGFirst; lG >= lGLast; lG--) {
        /*if(DEBUG) {
          ou << "new Round B'" << endl;
        }*/
        Node* parent = (*G)[lG]->getParent();
        int parent_of_lG_in_preL = parent == NULL? 0x7fffffff: parent->getID();
        int parent_of_lG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent->getID()];// not exist -1;
        rGFirst = G->preL_to_preR[lG];
        int lG_in_preR = G->preL_to_preR[lG];

        int lGminus1_in_preL = lG <= endG? 0x7fffffff : lG - 1;
        int lGminus1_in_preR = lG <= endG? 0x7fffffff : G->preL_to_preR[lG - 1];

        if (pathType == 0) {
          if (lG == endG || lGminus1_in_preL != parent_of_lG_in_preL) {//parent not exists or not the leftmost child.
            rGLast = rGFirst;
          } else {
            rGLast = parent_of_lG_in_preR + 1;
          }
        } else {// left and right
          if(endPathNode == endF) {
            if(lG == endG || lGminus1_in_preL != parent_of_lG_in_preL) rGLast = rGFirst;
            else rGLast = parent_of_lG_in_preR + 1;
          }
          else rGLast = rGFirst == endG_in_preR ? rGFirst : endG_in_preR;
        }

        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;

      /*if(DEBUG) {
        ou << "updateFnArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ", " << to_string(endG_in_preR) << ")" << endl;
      }*/
        updateFnArray(G->preR_to_ln[rGFirst], rGFirst, endG_in_preR);
      
      /*if(DEBUG) {
        ou << "updateFtArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ")" << endl;
      }*/
        updateFtArray(G->preR_to_ln[rGFirst], rGFirst);
            
      /*if(DEBUG) {
       ou << "endG_in_preR = " << to_string(endG_in_preR) << endl;
       ou << "start from rG = " << to_string(rGFirst) << endl;
       ou << "FN" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << fn[i] << " ";
      }
      ou << endl;
      ou << "FT" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << ft[i] << " ";
      }
      ou << endl;
      }*/
        lF = lF_prev;
      // loop C'
        for(int rF = rFFirst; rF >= rFLast; rF--) {
          /*if(DEBUG) {
            ou << "new Round C'" << endl;
          }*/
          //lF = startPathNode;
          if (rF == rFLast) {
              lF = lFLast;
          }
          int rG = rGFirst;
          int rG_in_preL = (G)->preR_to_preL[rG];
                
          //if(rF == rFLast) lF = F->preR_to_preL[rFLast]; 
          counter++;
          if(DEBUG && mark) {
            ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter  = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
          }

          int rF_in_preL = (F)->preR_to_preL[rF];
                
          bool FForestIsTree = lF == rF_in_preL;
          int rFSubtreeSize = (*F)[rF_in_preL]->getSubTreeSize();
                        
          int case1SLeftIndex, case1SRightIndex;//S[rF + 1, rG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];

          int case2SLeftIndex, case2SRightIndex;//S[rF, rG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;
        
          float case1 = 0, case2 = 0, case3 = 0;
          int case1_case, case2_case, case3_case;

          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
          if(mark) {
            for(int i = 0; i <= (*G)[endG]->getSubTreeSize(); i++) {
              forestdist[rF][i][(*G)[endG]->getSubTreeSize()] = FcurrentForestCost;
              forestdist[rF][(*G)[endG]->getSubTreeSize()][i] = FcurrentForestCost;
            }
            if(DEBUG) {
              ou << "forestdist[" << rF << ", " << (*G)[endG]->getSubTreeSize() << ", " << (*G)[endG]->getSubTreeSize() << "] = " << forestdist[rF][(*G)[endG]->getSubTreeSize()][(*G)[endG]->getSubTreeSize()] << endl;
            }
          }
          int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          float GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); // USE COST MODEL - reset to subtree insertion cost.
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
            if(DEBUG) {
              ou << "forestdist[" << (*F)[endF]->getSubTreeSize() << ", " << lG << ", " << rG << "] = " << forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] << endl;
            }
          }

          float minCost = 0;


          case1SLeftIndex = rF + 1;//fixed
          case1SRightIndex = rG;

          case1TLeftIndex = lG;//fixed
          case1TRightIndex = rG;

          case2SLeftIndex = rF;//fixed

          case3TLeftIndex = lG;//fixed


          bool rFIsConsecutiveNodeOfCurrentPathNode;
          bool rFIsRightSiblingOfCurrentPathNode;

          case1_case = 1;
          case3_case = 1;//otherwise

          if (startPathNode > 0) {
            rFIsConsecutiveNodeOfCurrentPathNode = startPathNode_in_preR - rF == 1;
            rFIsRightSiblingOfCurrentPathNode = rF + rFSubtreeSize == startPathNode_in_preR;
          } else {
            rFIsConsecutiveNodeOfCurrentPathNode = false;//consecutiveNode use T;
            rFIsRightSiblingOfCurrentPathNode = false;//
          }

          if(FForestIsTree) {
            if(rFSubtreeSize == 1) {
              // F_{lF, rF} - rF = null
              //case1 = GcurrentForestCost;//sumG to be computed 
              case1_case = 3;
            } else if(rFIsConsecutiveNodeOfCurrentPathNode) {
              // F_{lF, rF} - rF = tree
              //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]; to be computed
              case1_case = 2;
            }
            case3 = 0;
            /*if(DEBUG) {
              ou << "case3 = 0" << endl; 
            }*/
            case3_case = 2;// F_{lF, rF} - F(rF) = null
            } else {
              if (rFIsConsecutiveNodeOfCurrentPathNode) {// F_{lF, rF} - rF = the subforest to the left of the path
               //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]
               case1_case = 2;
              } else {//otherwise
                //case1 = s[case1SLeftIndex][case1SRightIndex];//S[rF + 1, rG];// have calculate
                case1_case = 1;
              }
              case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[rF_in_preL] : (F)->preL_to_sumDelCost[rF_in_preL]);// the first case in G should be G_{lG, rG} - l(rG) = null // F_{lF, rF} - F(rF), G_{lG, rG} - G(rG)
              /*if(DEBUG) {
                ou << "case3_case FcurrentForest - F(rF)" << endl;
              }*/
              if (rFIsRightSiblingOfCurrentPathNode) {
                case3_case = 3; // use T
              }
            }

            if (case3_case == 1) {
              case3SLeftIndex = rF + rFSubtreeSize;//delete the whole rightmost tree//otherwise
            }
          
            switch(case1_case) {
              case 1:
                case1 = s[case1SLeftIndex][case1SRightIndex];
                /*if(DEBUG) {
                  ou << "case1_case1 = s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl; 
                }*/
              break;
          
              case 2:
                case1 = t[case1TLeftIndex][case1TRightIndex];
                /*if(DEBUG) {
                  ou << "case1_case2 = t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl;
                }*/ 
                break;
              case 3:
                case1 = GcurrentForestCost;
                /*if(DEBUG) {
                  ou << "case1_case3 = GcurrentForestCost = " << GcurrentForestCost << endl;
                }*/
                break;          
            }
            case1 += (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel()));
            /*if(DEBUG) {
              ou << "case1 += ";
              if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
              else ou << "delete " << (*F)[rF]->getLabel() << endl;
            }*/
            minCost = case1;

            if (GcurrentForestSize == 1) {// the first case in G should be a node or a tree
              case2 = FcurrentForestCost;
              /*if(DEBUG) {
                ou << "case2_case1 = " << to_string(FcurrentForestCost) << endl;
              }*/
            } else {
              case2 = q[rF];
              /*if(DEBUG) {
                ou << "case2_case2 = q[" << to_string(rF) << "] = " << to_string(q[rF]) << endl;
              }*/
            }

            case2 += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
            /*if(DEBUG) {
              ou << "case2 += ";
              if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
              else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
            }*/
            if(case2 < minCost) {
              minCost = case2;
            }
        
            if(case3 < minCost) { 
              case3 += swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];// F(rF) - rF
              /*if(DEBUG) {
                if(swap) ou << "case3_case3 delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
                else ou << "case3_case3 delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) << endl;
              }*/
          
              if(case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
                //case3 += costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());
                /*if(DEBUG) {
                  if(swap) ou << "case3_case3 += " << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel())) << endl;
                  else ou << "case3_case3 += " << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel())) << endl;
                }*/
              }
              if(case3 < minCost) {
                minCost = case3;
              }
            }
        /*if(DEBUG) {
          ou << "case1 = " << to_string(case1) << endl;
          ou << "case2 = " << to_string(case2) << endl;
          ou << "case3 = " << to_string(case3) << endl;
        }*/
        s[rF][rG] = minCost;
        if(mark) {
          forestdist[rF][lG][rG] = minCost;
          if(DEBUG) {
            ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
          }
        }
        /*if(mark) {
          if(case1 == minCost) {
            if(swap) {
              if(DEBUG) {
                ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
            } else {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << "-" << endl;
              }
              map->setMap(rF_in_preL, -1);
            }
          } else if(case2 == minCost) {
            if(swap) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " << "-" << endl;
              }
              map->setMap(rG_in_preL, -1);
            } else {
              if(DEBUG) {
                ou << "- -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rG_in_preL);
            }
          } else if(case3 == minCost) {
            if(swap) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(rG_in_preL, rF_in_preL);
            } else {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, rG_in_preL);
            }
          }
        }*/
        dist = minCost;
        if(DEBUG) {
           ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(s[rF][rG]) << endl;
        }
        rG = ft[rG];  
        
        // loop D'
        while(rG >= rGLast) {// every G is a subforest not a subtree
          counter++;
          rG_in_preL = (G)->preR_to_preL[rG];
          if(DEBUG && mark) {
            ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
          }

          GcurrentForestSize++;
          GcurrentForestCost += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
            if(DEBUG) {
              ou << "forestdist[" << (*F)[endF]->getSubTreeSize() << ", " << lG << ", " << rG << "] = " << forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] << endl;
            }
          }

          switch (case1_case) {
            case 1:
              case1SRightIndex = rG;
              case1 = s[case1SLeftIndex][case1SRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/ 
              break; 
            case 2: 
              case1TRightIndex = rG;
              case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
            case 3: 
              case1 = GcurrentForestCost + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case3 GcurrentForestCost = " << GcurrentForestCost << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
          }
          
          minCost = case1;

          case2SRightIndex = fn[rG];
          case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));//G is not a tree or a node for sure in D loop
          /*if(DEBUG) {
            ou << "case2_case3 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
            if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
            else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
          }*/
          
          
          if(case2 < minCost) {
            minCost = case2;
          }

          case3 = swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];//F_{rF} - rF, G_{rG} - rG
          /*if(DEBUG) {
            if(swap) {
              ou << "case3_case delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
            } else {
              ou << "case3_case delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) <<  endl;
            }
          }*/
          
          if(case3 < minCost) {
            switch(case3_case) {
              case 1: 
                case3SRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
                case3 += s[case3SLeftIndex][case3SRightIndex];
                /*if(DEBUG) {
                  ou << "case3_case1 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl; 
                }*/
                break;
              
              case 2: 
                case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[rG_in_preL] : (G)->preL_to_sumInsCost[rG_in_preL]);
                /*if(DEBUG) {
                  ou << "case3_case2 += " << "GcurrentForestCost - G(rG) = " << to_string(GcurrentForestCost) << " - ";
                  if(swap) ou << (G)->preL_to_sumDelCost[rG_in_preL] << endl;
                  else ou << (G)->preL_to_sumInsCost[rG_in_preL] << endl;
                }*/
                break;
              
              case 3: 
                case3TRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
                case3 += t[case3TLeftIndex][case3TRightIndex];
                /*if(DEBUG) {
                  ou << "case3_case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl; 
                }*/
                break;
            }
              
            if(case3 < minCost) {
              case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
              /*if(DEBUG) {
                ou << "case3 += ";
                if(swap) ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
                else ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }*/
              if(case3 < minCost) {
                minCost = case3;
              }
            }
          }
          /*if(DEBUG) {
            ou << "case1 = " << to_string(case1) << endl;
            ou << "case2 = " << to_string(case2) << endl;
            ou << "case3 = " << to_string(case3) << endl; 
          }*/
            if(DEBUG) {
              ou << "s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(minCost) << endl;
            }
            s[rF][rG] = minCost;
            if(mark) {
              forestdist[rF][lG][rG] = minCost;
              if(DEBUG) {
                ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
              }
            }
            /*if(mark) {
              if(case1 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
                  }
                  map->setMap(-1, rF_in_preL);
                } else {
                  if(DEBUG) {
                    ou << (*F)[rF_in_preL]->getLabel() << " -> " << "-" << endl;
                  }
                  map->setMap(rF_in_preL, -1);
                }
              } else if(case2 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << (*G)[rG_in_preL]->getLabel() << " -> " << "-" << endl;
                  }
                  map->setMap(rG_in_preL, -1);
                } else {
                  if(DEBUG) {
                    ou << "- -> " << (*G)[rG_in_preL]->getLabel() << endl;
                  }
                  map->setMap(-1, rG_in_preL);
                }
              } else if(case3 == minCost) {
                if(swap) {
                  if(DEBUG) {
                    ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
                  }
                  map->setMap(rG_in_preL, rF_in_preL);
                } else {
                  if(DEBUG) {
                    ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
                  } 
                  map->setMap(rF_in_preL, rG_in_preL);
                }
              }
            }*/
            dist = minCost;
            rG = ft[rG];
          }
        }

        if(lGminus1_in_preL == parent_of_lG_in_preL && lGminus1_in_preL != 0x7fffffff) { // lG is the leftmost child of its parent
    
          if(hasRightPart) {
            if(swap) {
              delta[parent_of_lG_in_preL][endPathNode] = s[rFLast + 1][lGminus1_in_preR + 1];
              if(DEBUG) {
               ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
               ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
              }
            } else {
              delta[endPathNode][parent_of_lG_in_preL] = s[rFLast + 1][lGminus1_in_preR + 1];
              if(DEBUG) {
                ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
              }
            }
          }
      
          if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {//no left and right
            if(swap) {
              delta[parent_of_lG_in_preL][parent_of_endPathNode_preL] = s[rFLast][lGminus1_in_preR + 1];
              if(DEBUG) {
               ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
               ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
              }
            } else {
              delta[parent_of_endPathNode_preL][parent_of_lG_in_preL] = s[rFLast][lGminus1_in_preR + 1];
              if(DEBUG) {
                ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
              }
            }
          }

          for (int rF = rFFirst; rF >= rFLast; rF--) {
            q[rF] = s[rF][parent_of_lG_in_preR + 1];
            if(DEBUG) {
              ou << "q[" << to_string(rF) << "] = " << "s[" << to_string(rF) << ", " << to_string(parent_of_lG_in_preR + 1) << "] = ";
              ou << to_string(s[rF][parent_of_lG_in_preR + 1]) << endl; 
            }
          }
        }
        for (int rG = rGFirst; rG >= rGLast; rG = ft[rG]) {
          t[lG][rG] = s[rFLast][rG];
          if(DEBUG) {
            ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(rG) << "] = ";
            ou << to_string(s[rFLast][rG]) << endl;
          }
        }
      }
    }
    rF = endPathNode_in_preR;//in D' loop
    startPathNode = endPathNode;
    endPathNode = (*F)[endPathNode] ->getParent() == NULL? -1 : (*F)[endPathNode] ->getParent()->getID(); 
    endPathNode_in_preR = F->preL_to_preR[endPathNode];
  }

  /*if(mark) {
    int lF = 0;
    int rF = 0;
    int lG = 0;
    int rG = 0;
    int direction = 0;//0 for right 1 for left;
    int leaf = FreeStrategies[0][0].getLeaf();
    Node* parent = (*F)[leaf]->getParent();
    int* favouriteChild = new int[F->getTreeSize()];
    int count = 0;
    while(parent != NULL) {
      favouriteChild[count++] = leaf;
      leaf = parent->getID();
      parent = (*F)[leaf]->getParent();
    }

    if(DEBUG) {
      ou << "favouriteChild: " << endl;
      for(int i = count - 1; i >= 0; i--) {
        ou << favouriteChild[i] << " ";
      }
      ou << endl;
    }


    int i = count - 1;
    int prev_fav_child = 0;
    while(lF <= F->getTreeSize() && rF <= F->getTreeSize() && lG <= G->getTreeSize() && rG <= G->getTreeSize() ) {
      if(DEBUG) {
        ou << "(" << lF << ", " << rF << ", " << lG << ", " << rG << ")" << endl;
      }
      int rF_in_preL = F->preR_to_preL[rF];
      int rG_in_preL = G->preR_to_preL[rG];
      int lF_in_preR = F->preL_to_preR[lF];
      int favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];

      float da = 0;
      float db = 0;
      float dc = 0;

      bool isFTree = rF_in_preL == lF;
      bool isGTree = rG_in_preL == lG;

      if(rF == favouriteChild_in_preR) {
        direction = 1;
      }
      else if(lF == favouriteChild[i]) {
        i--;
        direction = 0;
      }
      
      if(isGTree && isFTree) {
        da = forestdist[rF + 1][lG][rG] + (swap? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel()));
        db = forestdist[rF][lG + 1][rG + 1] + (swap? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
        dc = forestdist[rF + 1][lG + 1][rG + 1] + (swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel()));
        if(DEBUG) {
          ou << "isGTree && isFTree" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << rF << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[rF][lG + 1][rG + 1] << endl;
          ou << "forestdist[" << rF + 1 << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[rF + 1][lG + 1][rG + 1] << endl;
          ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
        }

        if(forestdist[rF][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(-1, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rF_in_preL, -1);
          }
          rF++;
          lF++;
        }

        else if(forestdist[rF][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> -" << endl;
            }
            map->setMap(lG, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(-1, lG);
          }
          lG++;
          rG++;
        }

        else if(forestdist[rF][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(lG, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(rF_in_preL, lG);
          }
          lG++;
          rG++;
          lF++;
          rF++;
        }

      } else if(isFTree && direction == 0) {
        da = forestdist[rF + 1][lG][rG] + (swap? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel()));
        db = forestdist[rF][lG][rG + 1] + (swap? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
        int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
        dc = forestdist[rF + 1][lG][rG + size_of_rG] + delta[rF_in_preL][rG_in_preL] +(swap? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));

        if(DEBUG) {
          ou << "isFTree && direction == 0" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
        }
        if(forestdist[rF][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(-1, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rF_in_preL, -1);
          }
          lF++;
          rF++;
        }

        else if(forestdist[rF][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[rG_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rG_in_preL, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[rG_in_preL] << endl; 
            }
            map->setMap(-1, rG_in_preL);
          }
          rG++;
        }

        else if(forestdist[rF][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(rG_in_preL, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
            }
            map->setMap(rF_in_preL, rG_in_preL);
          }
          lF++;
          rF++;
          rG = rG + size_of_rG;
        }
      } else if(isFTree && direction == 1) {
        da = forestdist[rF + 1][lG][rG] + (swap? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel()));
        db = forestdist[rF][lG + 1][rG] + (swap? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
        int size_of_lG = (*G)[lG]->getSubTreeSize();
        dc = forestdist[rF + 1][lG + size_of_lG][rG] + delta[rF_in_preL][lG] + (swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel()));

        if(DEBUG) {
          ou << "isFTree && direction == 1" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
        }

        if(forestdist[lF_in_preR][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(-1, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rF_in_preL, -1);
          }
          lF++;
          rF++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> -" << endl;
            }
            map->setMap(lG, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(-1, lG);
          }
          lG++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(lG, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(rF_in_preL, lG);
          }
          lF++;
          rF++;
          lG = lG + size_of_lG;
        }
      } else if(isGTree && direction == 0) {
        da = forestdist[rF + 1][lG][rG] + (swap? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel()));
        db = forestdist[rF][lG + 1][rG + 1] + (swap? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
        int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
        dc = forestdist[rF + size_of_rF][lG + 1][rG + 1] + delta[rF_in_preL][lG] + (swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel()));

        if(DEBUG) {
          ou << "isGTree && direction == 0" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
        }

        if(forestdist[rF][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(-1, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rF_in_preL, -1);
          }
          rF++;
        }
        else if(forestdist[rF][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> -" << endl;
            }
            map->setMap(lG, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(-1, lG);
          }
          lG++;
          rG++;
        }
        else if(forestdist[rF][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(lG, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(rF_in_preL, lG);
          }
          rF = rF + size_of_rF;
          lG++;
          rG++;
        }
      } else if(isGTree && direction == 1) {
        da = forestdist[lF_in_preR - 1][lG][rG] + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
        db = forestdist[rF][lG + 1][rG + 1] + (swap? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
        int size_of_lF = (*F)[lF]->getSubTreeSize();
        dc = forestdist[lF_in_preR - size_of_lF][lG + 1][rG + 1] + delta[lF][lG] + (swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel()));

        if(DEBUG) {
          ou << "isGTree && direction == 1" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
        }

        if(forestdist[lF_in_preR][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[lF]->getLabel() << endl; 
            }
            map->setMap(-1, lF);
          } else {
            if(DEBUG) {
              ou << (*F)[lF]->getLabel() << " -> -" << endl; 
            }
            map->setMap(lF, -1);
          }
          lF++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> -" << endl;
            }
            map->setMap(lG, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(-1, lG);
          }
          lG++;
          rG++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
            }
            map->setMap(lG, lF);
          } else {
            if(DEBUG) {
              ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(lF, lG);
          }
          lF = lF + size_of_lF;
          lG++;
          rG++;
        }
      }
      //F and G are forest and direction = right
      else if(direction == 0) {
        da = forestdist[rF + 1][lG][rG] + (swap? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel()));
        db = forestdist[rF][lG][rG + 1] + (swap? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
        int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
        int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
        dc = forestdist[rF + size_of_rF][lG][rG + size_of_rG] + delta[rF_in_preL][rG_in_preL] + (swap? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));

        if(DEBUG) {
          ou << "Both forest && direction == 0" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
        }

        if(forestdist[rF][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(-1, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rF_in_preL, -1);
          }
          rF++;
        }

        else if(forestdist[rF][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[rG_in_preL]->getLabel() << " -> -" << endl;
            }
            map->setMap(rG_in_preL, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[rG_in_preL] << endl;
            }
            map->setMap(-1, rG_in_preL);
          }
          rG++;
        }

        else if(forestdist[rF][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
            }
            map->setMap(rG_in_preL, rF_in_preL);
          } else {
            if(DEBUG) {
              ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
            }
            map->setMap(rF_in_preL, rG_in_preL);
          }
          rF = rF + size_of_rF;
          rG = rG + size_of_rG;
        }
      }
      //F and G are forest and direction = left
      else if(direction == 1) {
        da = forestdist[lF_in_preR - 1][lG][rG] + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
        db = forestdist[lF_in_preR][lG + 1][rG] + (swap? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
        int size_of_lF = (*F)[lF]->getSubTreeSize();
        int size_of_lG = (*G)[lG]->getSubTreeSize();
        dc = forestdist[lF_in_preR - size_of_lF][lG + size_of_lG][rG] + delta[lF][lG] + (swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel()));

        if(DEBUG) {
          ou << "Both forest && direction == 1" << endl;
          ou << "da = " << da << endl;
          ou << "db = " << db << endl;
          ou << "dc = " << dc << endl; 
          ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
        }

        if(forestdist[lF_in_preR][lG][rG] == da) {
          if(swap) {
            if(DEBUG) {
              ou << "- -> " << (*F)[lF]->getLabel() << endl;
            }
            map->setMap(-1, lF);
          } else {
            if(DEBUG) {
              ou << (*F)[lF]->getLabel() << " -> -" << endl;
            }
            map->setMap(lF, -1);
          }
          lF++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == db) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> -" << endl;
            }
            map->setMap(lG, -1);
          } else {
            if(DEBUG) {
              ou << "- -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(-1, lG);
          }
          lG++;
        }

        else if(forestdist[lF_in_preR][lG][rG] == dc) {
          if(swap) {
            if(DEBUG) {
              ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
            }
            map->setMap(lG, lF);
          } else {
            if(DEBUG) {
              ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
            }
            map->setMap(lF, lG);
          }
          lF = lF + size_of_lF;
          lG = lG + size_of_lG; 
        }
      }
    }
  }*/
  return dist;
};


/*
  first add right then add left, left decompose
*/
float TreeComparison::spfA_RL(Node* a, Node* b, int leaf, int pathType, float*** forestdist, bool swap, bool mark) {
  Tree *F, *G;
  if(swap) {
    F = B_;
    G = A_;
  } else {
    F = A_;
    G = B_;
  }
  int endF = a->getID(); 
  int endG = b->getID();
  int sizeF = a->getSubTreeSize();
  int sizeG = b->getSubTreeSize();
  int endF_in_preR = F->preL_to_preR[endF];
  int endG_in_preR = G->preL_to_preR[endG];
  int endPathNode = leaf;
  int endPathNode_in_preR = F->preL_to_preR[endPathNode];
  int startPathNode = -1;
  int lFFirst, lFLast, lF;
  int rFFirst, rFLast, rF;
  int lGFirst, lGLast;
  int rGFirst, rGLast;

  int FcurrentForestSize = 0;
  float FcurrentForestCost = 0;
  int FtmpForestSize = 0;
  float FtmpForestCost = 0;

  float dist = 0;

  int rF_prev = endPathNode_in_preR;

  if(DEBUG) {
    ou << "spfA_RL(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") " << swap << " " << endl;;
  }

  if(DEBUG) {
    cout << "spfA_RL(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") swap = " << swap << " counter = " << counter << " " << endl;
  }
  
  //loop A
  while(endPathNode >= endF) {
    endPathNode_in_preR = F->preL_to_preR[endPathNode];
    int startPathNode_in_preR = startPathNode == -1? 0x7fffffff : F->preL_to_preR[startPathNode];

    int parent_of_endPathNode_preL = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : (*F)[endPathNode]->getParent()->getID();
    int parent_of_endPathNode_preR = (*F)[endPathNode]->getParent() == NULL? 0x7fffffff : F->preL_to_preR[parent_of_endPathNode_preL];

    bool hasLeftPart;
    bool hasRightPart;

    
    if(startPathNode - endPathNode > 1) {
      hasLeftPart = true;
    } else {
      hasLeftPart = false;
    }
    if(startPathNode >= 0 && startPathNode_in_preR - endPathNode_in_preR > 1) {
      hasRightPart = true;
    } else {
      hasRightPart = false;
    }

    /*if(DEBUG) {
      ou << "startPathNode = " << startPathNode << " endPathNode = " << endPathNode << endl;
      ou << "startPathNode_in_preR = " << startPathNode_in_preR << " endPathNode_in_preR = " << endPathNode_in_preR << endl;
      ou << "hasLeftPart = " << hasLeftPart << " hasRightPart = " << hasRightPart << endl;
    }*/

    // add right -> right decompose left path 
    if(pathType == 0 || pathType == 2 && hasRightPart) {
      if (startPathNode == -1) {
        lFFirst = endPathNode;
        rFFirst = F->preL_to_preR[endPathNode];
      } else {
        rFFirst = F->preL_to_preR[startPathNode] - 1;//the node right to the node on the path
        lFFirst = endPathNode + 1;//lFirst is set to the node on the path
      }

      if(!hasLeftPart) {
        lFLast = endPathNode;
      }
      rFLast = hasLeftPart? endPathNode_in_preR + 1 : endPathNode_in_preR;

      lGLast = endG;
      lGFirst = (lGLast + sizeG) - 1;

      fn[fn_ft_length - 1] = -1;
      /*if(DEBUG) {
         ou << "initial fn and ft endG_in_preR = " << to_string(endG_in_preR) << " endG_in_preR + sizeG = " << to_string(endG_in_preR + sizeG) << endl;
      }*/
      
      for (int i = endG_in_preR; i < endG_in_preR + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }

      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
      //loop B'
      for (int lG = lGFirst; lG >= lGLast; lG--) {
        /*if(DEBUG) {
          ou << "new Round B'" << endl;
        }*/
        Node* parent = (*G)[lG]->getParent();
        int parent_of_lG_in_preL = parent == NULL? 0x7fffffff: parent->getID();
        int parent_of_lG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent->getID()];// not exist -1;
        rGFirst = G->preL_to_preR[lG];
        int lG_in_preR = G->preL_to_preR[lG];

        int lGminus1_in_preL = lG <= endG? 0x7fffffff : lG - 1;
        int lGminus1_in_preR = lG <= endG? 0x7fffffff : G->preL_to_preR[lG - 1];

        if (pathType == 0) {
          if (lG == endG || lGminus1_in_preL != parent_of_lG_in_preL) {//parent not exists or not the leftmost child.
            rGLast = rGFirst;
          } else {
            rGLast = parent_of_lG_in_preR + 1;
          }
        } else {// left and right
          rGLast = rGFirst == endG_in_preR ? rGFirst : endG_in_preR;
        }

        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;

      /*if(DEBUG) {
        ou << "updateFnArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ", " << to_string(endG_in_preR) << ")" << endl;
      }*/
        updateFnArray(G->preR_to_ln[rGFirst], rGFirst, endG_in_preR);
      
      /*if(DEBUG) {
        ou << "updateFtArray(" << to_string(G->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ")" << endl;
      }*/
        updateFtArray(G->preR_to_ln[rGFirst], rGFirst);
            
      /*if(DEBUG) {
       ou << "endG_in_preR = " << to_string(endG_in_preR) << endl;
       ou << "start from rG = " << to_string(rGFirst) << endl;
       ou << "FN" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << fn[i] << " ";
      }
      ou << endl;
      ou << "FT" << endl;
      for(int i = endG; i < endG + sizeG + 1; i++) {
        ou << ft[i] << " ";
      }
      ou << endl;
      }*/
      // loop C'
        for(int rF = rFFirst; rF >= rFLast; rF--) {
          /*if(DEBUG) {
            ou << "new Round C'" << endl;
          }*/

          lF = startPathNode;
          if (rF == rFLast && !hasLeftPart) {
              lF = lFLast;
          }
          int rG = rGFirst;
          int rG_in_preL = (G)->preR_to_preL[rG];
                
          //if(rF == rFLast) lF = F->preR_to_preL[rFLast]; 
          counter++;
          if(DEBUG && mark) {
            ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter  = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
          }

          int rF_in_preL = (F)->preR_to_preL[rF];
                
          bool FForestIsTree = lF == rF_in_preL;
          int rFSubtreeSize = (*F)[rF_in_preL]->getSubTreeSize();
                        
          int case1SLeftIndex, case1SRightIndex;//S[rF + 1, rG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];

          int case2SLeftIndex, case2SRightIndex;//S[rF, rG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;
        
          float case1 = 0, case2 = 0, case3 = 0;
          int case1_case, case2_case, case3_case;

          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[rF_in_preL]->getLabel()) : costModel_.del((*F)[rF_in_preL]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
          if(mark) {
            for(int i = 0; i <= (*G)[endG]->getSubTreeSize(); i++) {
              forestdist[rF_in_preL][i][(*G)[endG]->getSubTreeSize()] = FcurrentForestCost;
              forestdist[rF_in_preL][(*G)[endG]->getSubTreeSize()][i] = FcurrentForestCost;
            }
          }
          int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          float GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); // USE COST MODEL - reset to subtree insertion cost.
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
          }

          float minCost = 0;


          case1SLeftIndex = rF + 1;//fixed
          case1SRightIndex = rG;

          case1TLeftIndex = lG;//fixed
          case1TRightIndex = rG;

          case2SLeftIndex = rF;//fixed

          case3TLeftIndex = lG;//fixed


          bool rFIsConsecutiveNodeOfCurrentPathNode;
          bool rFIsRightSiblingOfCurrentPathNode;

          case1_case = 1;
          case3_case = 1;//otherwise

          if (startPathNode > 0) {
            rFIsConsecutiveNodeOfCurrentPathNode = startPathNode_in_preR - rF == 1;
            rFIsRightSiblingOfCurrentPathNode = rF + rFSubtreeSize == startPathNode_in_preR;
          } else {
            rFIsConsecutiveNodeOfCurrentPathNode = false;//consecutiveNode use T;
            rFIsRightSiblingOfCurrentPathNode = false;//
          }

          if(FForestIsTree) {
            if(rFSubtreeSize == 1) {
              // F_{lF, rF} - rF = null
              //case1 = GcurrentForestCost;//sumG to be computed 
              case1_case = 3;
            } else if(rFIsConsecutiveNodeOfCurrentPathNode) {
              // F_{lF, rF} - rF = tree
              //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]; to be computed
              case1_case = 2;
            }
            case3 = 0;
            /*if(DEBUG) {
              ou << "case3 = 0" << endl; 
            }*/
            case3_case = 2;// F_{lF, rF} - F(rF) = null
            } else {
              if (rFIsConsecutiveNodeOfCurrentPathNode) {// F_{lF, rF} - rF = the subforest to the left of the path
               //case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]
               case1_case = 2;
              } else {//otherwise
                //case1 = s[case1SLeftIndex][case1SRightIndex];//S[rF + 1, rG];// have calculate
                case1_case = 1;
              }
              case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[rF_in_preL] : (F)->preL_to_sumDelCost[rF_in_preL]);// the first case in G should be G_{lG, rG} - l(rG) = null // F_{lF, rF} - F(rF), G_{lG, rG} - G(rG)
              /*if(DEBUG) {
                ou << "case3_case FcurrentForest - F(rF)" << endl;
              }*/
              if (rFIsRightSiblingOfCurrentPathNode) {
                case3_case = 3; // use T
              }
            }

            if (case3_case == 1) {
              case3SLeftIndex = rF + rFSubtreeSize;//delete the whole rightmost tree//otherwise
              /*if(DEBUG) {
                ou << "case3SLeftIndex = rF(" << rF << ") + rFSubtreeSize(" << rFSubtreeSize << ") = " << case3SLeftIndex << endl;
                ou << "rF_in_preL = " << rF_in_preL << endl;
              }*/
            }
          
            switch(case1_case) {
              case 1:
                case1 = s[case1SLeftIndex][case1SRightIndex];
                /*if(DEBUG) {
                  ou << "case1_case1 = s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl; 
                }*/
              break;
          
              case 2:
                case1 = t[case1TLeftIndex][case1TRightIndex];
                /*if(DEBUG) {
                  ou << "case1_case2 = t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl;
                }*/ 
                break;
              case 3:
                case1 = GcurrentForestCost;
                /*if(DEBUG) {
                  ou << "case1_case3 = GcurrentForestCost = " << GcurrentForestCost << endl;
                }*/
                break;          
            }
            case1 += (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel()));
            /*if(DEBUG) {
              ou << "case1 += ";
              if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
              else ou << "delete " << (*F)[rF]->getLabel() << endl;
            }*/
            minCost = case1;
          
            if (GcurrentForestSize == 1) {// the first case in G should be a node or a tree
              case2 = FcurrentForestCost;
              /*if(DEBUG) {
                ou << "case2_case1 = " << to_string(FcurrentForestCost) << endl;
              }*/
            } else {
              case2 = q[rF];
              /*if(DEBUG) {
                ou << "case2_case2 = q[" << to_string(rF) << "] = " << to_string(q[rF]) << endl;
              }*/
            }

            case2 += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
            /*if(DEBUG) {
              ou << "case2 += ";
              if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
              else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
            }*/
            if(case2 < minCost) {
              minCost = case2;
            }
        
            if(case3 < minCost) { 
              case3 += swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];// F(rF) - rF
              /*if(DEBUG) {
                if(swap) ou << "case3_case3 delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
                else ou << "case3_case3 delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) << endl;
              }*/
          
              if(case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
                //case3 += costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());
                /*if(DEBUG) {
                  if(swap) ou << "case3_case3 += " << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel())) << endl;
                  else ou << "case3_case3 += " << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << " = " << to_string(costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel())) << endl;
                }*/
              }
              if(case3 < minCost) {
                minCost = case3;
              }
            }
          /*if(DEBUG) {
            ou << "case1 = " << case1 << endl;
            ou << "case2 = " << case2 << endl;
            ou << "case3 = " << case3 << endl;
          }*/
        s[rF][rG] = minCost;
        if(mark) {
          forestdist[rF_in_preL][lG][rG] = minCost;
          if(DEBUG) {
            ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
          }
        }
        dist = minCost;
        if(DEBUG) {
           ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(s[rF][rG]) << endl;
        }
        rG = ft[rG];  
        
        // loop D'
        while(rG >= rGLast) {// every G is a subforest not a subtree
          counter++;
          rG_in_preL = (G)->preR_to_preL[rG];
          if(DEBUG && mark) {
            ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
          }

          GcurrentForestSize++;
          GcurrentForestCost += (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
          }

          switch (case1_case) {
            case 1:
              case1SRightIndex = rG;
              case1 = s[case1SLeftIndex][case1SRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/ 
              break; 
            case 2: 
              case1TRightIndex = rG;
              case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
            case 3: 
              case1 = GcurrentForestCost + (swap ? costModel_.ins((*F)[rF]->getLabel()) : costModel_.del((*F)[rF]->getLabel())); 
              /*if(DEBUG) {
                ou << "case1_case3 GcurrentForestCost = " << GcurrentForestCost << " + ";
                if(swap) ou << "insert " << (*F)[rF]->getLabel() << endl;
                else ou << "delete " << (*F)[rF]->getLabel() << endl;
              }*/
              break; 
          }
          
          minCost = case1;

          case2SRightIndex = fn[rG];
          case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[rG_in_preL]->getLabel()) : costModel_.ins((*G)[rG_in_preL]->getLabel()));//G is not a tree or a node for sure in D loop
          /*if(DEBUG) {
            ou << "case2_case3 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
            if(swap) ou << "delete " << (*G)[rG_in_preL]->getLabel() << endl;
            else ou << "insert " << (*G)[rG_in_preL]->getLabel() << endl;
          }*/
          
          
          if(case2 < minCost) {
            minCost = case2;
          }

          case3 = swap ? delta[rG_in_preL][rF_in_preL] : delta[rF_in_preL][rG_in_preL];//F_{rF} - rF, G_{rG} - rG
          /*if(DEBUG) {
            if(swap) {
              ou << "case3_case delta[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "] = " << to_string(delta[rG_in_preL][rF_in_preL]) << endl;
            } else {
              ou << "case3_case delta[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "] = " << to_string(delta[rF_in_preL][rG_in_preL]) <<  endl;
            }
          }*/
          
          if(case3 < minCost) {
            switch(case3_case) {
              case 1: 
                case3SRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
                case3 += s[case3SLeftIndex][case3SRightIndex];
                /*if(DEBUG) {
                  ou << "case3_case1 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl; 
                }*/
                break;
              
              case 2: 
                case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[rG_in_preL] : (G)->preL_to_sumInsCost[rG_in_preL]);
                /*if(DEBUG) {
                  ou << "case3_case2 += " << "GcurrentForestCost - G(rG) = " << to_string(GcurrentForestCost) << " - ";
                  if(swap) ou << (G)->preL_to_sumDelCost[rG_in_preL] << endl;
                  else ou << (G)->preL_to_sumInsCost[rG_in_preL] << endl;
                }*/
                break;
              
              case 3: 
                case3TRightIndex = fn[(rG + (*G)[rG_in_preL]->getSubTreeSize()) - 1];
                case3 += t[case3TLeftIndex][case3TRightIndex];
                /*if(DEBUG) {
                  ou << "case3_case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl; 
                }*/
                break;
            }
              
            if(case3 < minCost) {
              case3 += (swap ? costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel()) : costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel()));
              /*if(DEBUG) {
                ou << "case3 += ";
                if(swap) ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
                else ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }*/
              if(case3 < minCost) {
                minCost = case3;
              }
            }
          }
          if(DEBUG) {
            ou << "s[" << to_string(rF) << ", " << to_string(rG) << "] = " << to_string(minCost) << endl;
          }
          /*if(DEBUG) {
            ou << "case1 = " << case1 << endl;
            ou << "case2 = " << case2 << endl;
            ou << "case3 = " << case3 << endl;
          }*/
          s[rF][rG] = minCost;
          if(mark) {
            forestdist[rF_in_preL][lG][rG] = minCost;
            if(DEBUG) {
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
            }
          }
          dist = minCost;
          rG = ft[rG];
          }
          rF_prev = rF;
          /*if(DEBUG) {
            ou << "rF_prev = " << rF_prev << endl;
          }*/
        }

        if(lGminus1_in_preL == parent_of_lG_in_preL && lGminus1_in_preL != 0x7fffffff) { // lG is the leftmost child of its parent
          if (!hasLeftPart) {
            if(hasRightPart) {
              if(swap) {
                delta[parent_of_lG_in_preL][endPathNode] = s[rFLast + 1][lGminus1_in_preR + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(rFLast + 1) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                  ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
                }*/
              } else {
                delta[endPathNode][parent_of_lG_in_preL] = s[rFLast + 1][lGminus1_in_preR + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast + 1) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                  ou << to_string(s[rFLast + 1][lGminus1_in_preR + 1]) << endl;
                }*/
              }
            }
      
            if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {//no left and right
              if(swap) {
                delta[parent_of_lG_in_preL][parent_of_endPathNode_preL] = s[rFLast][lGminus1_in_preR + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_lG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                  ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
                }*/
              } else {
                delta[parent_of_endPathNode_preL][parent_of_lG_in_preL] = s[rFLast][lGminus1_in_preR + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "] = "; //rightmosts child of p(lG)
                  ou << to_string(s[rFLast][lGminus1_in_preR + 1]) << endl;
                }*/
              }
            }
          }

          for (int rF = rFFirst; rF >= rFLast; rF--) {
            q[rF] = s[rF][parent_of_lG_in_preR + 1];
            /*if(DEBUG) {
              ou << "q[" << to_string(rF) << "] = " << "s[" << to_string(rF) << ", " << to_string(parent_of_lG_in_preR + 1) << "] = ";
              ou << to_string(s[rF][parent_of_lG_in_preR + 1]) << endl; 
            }*/
          }
        }
        for (int rG = rGFirst; rG >= rGLast; rG = ft[rG]) {
          t[lG][rG] = s[rFLast][rG];
          /*if(DEBUG) {
            ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(rFLast) << ", " << to_string(rG) << "] = ";
            ou << to_string(s[rFLast][rG]) << endl;
          }*/
        }
      }
    }

    // add left right path decompse
    if (pathType == 1 || pathType == 2 && hasLeftPart || pathType == 2 && !hasLeftPart && !hasRightPart) {
      // lFFirst and LFLast is important in this condition.
      // rGFirst and rGLast is to set in this stage.
      if(startPathNode == -1) {
        lFFirst = endPathNode;//the first node is the node on the path
        rFFirst = endPathNode_in_preR;// the first node is the node on the path
      } else {
        lFFirst = startPathNode - 1;//the first node is the node to the left of the path
        rFFirst = startPathNode_in_preR;//rFFirst set to the node on the path
      }

      lFLast = endPathNode;
      rFLast = F->preL_to_preR[endPathNode];

      rGLast = G->preL_to_preR[endG];
      rGFirst = (rGLast + sizeG) - 1; // get the leftmost child in G

      fn[fn_ft_length - 1] = -1;
      /*if(DEBUG) {
        ou << "initial fn and ft endG = " << to_string(endG) << " endG + sizeG = " << to_string(endG + sizeG) << endl;
      }*/
      for (int i = endG; i < endG + sizeG; i++) {
        fn[i] = -1;
        ft[i] = -1;
      }
      
      FtmpForestSize = FcurrentForestSize;
      FtmpForestCost = FcurrentForestCost;
      //loop B
      for(int rG = rGFirst; rG >= rGLast; rG--) {
        /*if(DEBUG) {
          ou << "new Round B" << endl;
        }*/
        int rG_in_preL = (G)->preR_to_preL[rG];
        Node* parent = (*G)[rG_in_preL]->getParent();
        int parent_of_rG_in_preL = parent == NULL? 0x7fffffff : parent->getID();
        int parent_of_rG_in_preR = parent == NULL? 0x7fffffff : G->preL_to_preR[parent_of_rG_in_preL];
        lGFirst = G->preR_to_preL[rG];// lGFirst is set to rGFirst;
        
        int rGminus1_in_preL = rG <= endG_in_preR? 0x7fffffff : G->preR_to_preL[rG - 1];// rG should greater than endG_in_preR cause rG is the inner node of subtree enG
        int rGminus1_in_preR = rG <= endG_in_preR? 0x7fffffff : rG - 1;
        
        if (pathType == 1){
          if (lGFirst == endG || rGminus1_in_preL != parent_of_rG_in_preL) {// parent not exist or not the rightmost child
            lGLast = lGFirst;//lGLast is set to lGFirst
          } else {
            lGLast = parent_of_rG_in_preL + 1;//lGLast is set to the leftmost child of rG's parent
          }
        } else {
          if(endPathNode == endF){
            if(lGFirst == endG || rGminus1_in_preL != parent_of_rG_in_preL) lGLast = lGFirst;
            else lGLast = parent_of_rG_in_preL + 1;
          }
          else lGLast = lGFirst == endG? lGFirst : endG + 1;//lGLast is set to the leftmost child of the whole tree
        }
      
        updateFnArray(G->preL_to_ln[lGFirst], lGFirst, endG); //stores the counter in D loop fn[ln] stores the start point
        updateFtArray(G->preL_to_ln[lGFirst], lGFirst); 
        
        FcurrentForestSize = FtmpForestSize;
        FcurrentForestCost = FtmpForestCost;
        rF = rF_prev;
        //loop C
        for(int lF = lFFirst; lF >= lFLast; lF--) {
          /*if(DEBUG) {
            ou << "new round C" << endl;
          }*/
          //rF = startPathNode_in_preR;
          if (lF == lFLast) {
              rF = rFLast;
          }
          /*if(DEBUG) {
            ou << "hasRightPart = " << hasRightPart << endl;
            ou << "rFLast = " << rFLast << endl;
            ou << "lFLast = " << lFLast << endl;
            ou << "rF_prev = " << rF_prev << endl;
          }*/
          int lG = lGFirst;
          int lF_in_preR = (F)->preL_to_preR[lF];
          counter++;

          if(DEBUG && mark) {
            ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
            ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
          }


          FcurrentForestSize++;
          FcurrentForestCost += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); // USE COST MODEL - sum up deletion cost of a forest.
          if(mark) {
            for(int i = 0; i <= (*G)[endG]->getSubTreeSize(); i++) {
              forestdist[lF][i][(*G)[endG]->getSubTreeSize()] = FcurrentForestCost;
              forestdist[lF][(*G)[endG]->getSubTreeSize()][i] = FcurrentForestCost;
            }
          }
          int GcurrentForestSize = (*G)[lG]->getSubTreeSize();
          int GcurrentForestCost = (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 
          if(mark) {
            forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
          }

          bool fForestIsTree = lF_in_preR == rF;
          int lFSubtreeSize = (*F)[lF]->getSubTreeSize();
          bool lFIsLeftSiblingOfCurrentPathNode = lF + lFSubtreeSize == startPathNode;
          bool lFIsConsecutiveNodeOfCurrentPathNode = startPathNode - lF == 1;

          int case1SLeftIndex, case1SRightIndex;//S[lF + 1, lG];
          int case1TLeftIndex, case1TRightIndex;//T[lG, rG];
                
          int case2SLeftIndex, case2SRightIndex;//S[lF, lG];

          int case3SLeftIndex, case3SRightIndex;
          int case3TLeftIndex, case3TRightIndex;

          case1SLeftIndex =  lF + 1;//fixed
          case2SLeftIndex = lF;//fixed

          case1TRightIndex = rG;//fixed

          case3TRightIndex = rG;//fixed

          float case1 = 0, case2 = 0, case3 = 0;
          int case1_case, case2_case, case3_case;

          float minCost = 0;

          case1_case = 1;
          case3_case = 1;

          if (fForestIsTree) { // F_{lF,rF} is a tree.
            if (lFSubtreeSize == 1) { // F_{lF,rF} is a single node.
              // F_{lF, rF} - lF = null
              case1_case = 3;
            } else if (lFIsConsecutiveNodeOfCurrentPathNode) { // F_{lF,rF}-lF is the path node subtree.
              // F_{lF, rF} - rF = tree
              case1_case = 2;
            }
            
            case3 = 0;//F_{lF, rF} - F(lF) = null
            /*if(DEBUG) {
              ou << "case3 = 0" << endl; 
            }*/
            
            case3_case = 2;
          } else {
            if (lFIsConsecutiveNodeOfCurrentPathNode) {
              // F_{lF, rF} - lF = tree
              case1_case = 2;
            }
            case3 = FcurrentForestCost - (swap ? (F)->preL_to_sumInsCost[lF] : (F)->preL_to_sumDelCost[lF]); // USE COST MODEL - Delete F_{lF,rF}-F_lF.
            /*if(DEBUG) {
              ou << "case3_case FcurrentForest - F(lF) = " << endl;
              ou << to_string(FcurrentForestCost) << " - ";
              if(swap) ou << to_string((F)->preL_to_sumInsCost[lF]) << endl;
              else ou << to_string((F)->preL_to_sumDelCost[lF]) << endl;
            }*/
                    
            if (lFIsLeftSiblingOfCurrentPathNode) {
              case3_case = 3;
            }
          }

          if(case3_case == 1) {
            case3SLeftIndex = lF + lFSubtreeSize;
            /*if(DEBUG) {
              ou << "case3SLeftIndex = lF(" << lF << ") + lFSubtreeSize(" << lFSubtreeSize << ") = " << case3SLeftIndex << endl;
            }*/
          }

          switch(case1_case) {
            case 1: 
              case1SRightIndex = lG;
              case1 = s[case1SLeftIndex][case1SRightIndex]; 
              /*if(DEBUG) {
                ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = ";
                ou << to_string(s[case1SLeftIndex][case1SRightIndex]) << endl;
              }*/
              break;
            
            case 2: 
              case1TLeftIndex = lG;
              case1 = t[case1TLeftIndex][case1TRightIndex]; 
              /*if(DEBUG) {
                ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = ";
                ou << to_string(t[case1TLeftIndex][case1TRightIndex]) << endl; 
              }*/
              break;
            
            case 3: 
              case1 = GcurrentForestCost; 
              /*if(DEBUG) {
                ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << endl; 
              }*/
              break; // USE COST MODEL - Insert G_{lG,rG}.
          }

          case1 += (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
          minCost = case1;

          if (GcurrentForestSize == 1) { // G_{lG,rG} is a single node.
            case2 = FcurrentForestCost; // USE COST MODEL - Delete F_{lF,rF}.
            /*if(DEBUG) {
              ou << "case2_case1 FcurrentForestCost = " << to_string(FcurrentForestCost) << endl;
            }*/
          } else { // G_{lG,rG} is a tree.
            case2 = q[lF];
            /*if(DEBUG) {
              ou << "case2_case2 q[" << to_string(lF) << "] = " << to_string(q[lF]) << endl;
            }*/
          }
          case2 += (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
          if(case2 < minCost) {
            minCost = case2;
          }

          if (case3 < minCost) {
            case3 += swap? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
              if(swap) ou << "case3_case3 delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
              else ou << "case3_case3 delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
            }*/
            if(case3 < minCost) {
              case3 += swap? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
              /*if(DEBUG) {
                if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }*/
            } 
            if(case3 < minCost) {
              minCost = case3;
            }
          }
          if(DEBUG) {
            ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
          }
          /*if(DEBUG) {
            ou << "case1 = " << case1 << endl;
            ou << "case2 = " << case2 << endl;
            ou << "case3 = " << case3 << endl;
          }*/
          dist = minCost;
          s[lF][lG] = minCost;
          if(mark) {
            forestdist[lF][lG][rG] = minCost;
            if(DEBUG) {
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }
          }
          lG = ft[lG];
          
          //loop D
          while (lG >= lGLast) {
            counter++;
            if(DEBUG && mark) {
              ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ") counter = " << to_string(counter) << endl;
              ou << "Save to s[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
            }

            GcurrentForestSize++;
            GcurrentForestCost += (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));
            if(mark) {
              forestdist[(*F)[endF]->getSubTreeSize()][lG][rG] = GcurrentForestCost;
            }
            minCost = 0;

            switch(case1_case) {
              case 1:
                case1SRightIndex = lG;
                case1 = s[case1SLeftIndex][case1SRightIndex] + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                  ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "] = " << to_string(s[case1SLeftIndex][case1SRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 2: 
                case1TLeftIndex = lG;
                case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel())); 
                /*if(DEBUG) {
                  ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "] = " << to_string(t[case1TLeftIndex][case1TRightIndex]) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
              
              case 3: 
                case1 = GcurrentForestCost + (swap? costModel_.ins((*F)[lF]->getLabel()) : costModel_.del((*F)[lF]->getLabel()));
                /*if(DEBUG) {
                  ou << "case1_case3 GcurrentForestCost = " << to_string(GcurrentForestCost) << " + ";
                  if(swap) ou << "insert " << (*F)[lF]->getLabel() << endl;
                  else ou << "delete " << (*F)[lF]->getLabel() << endl;
                }*/
                break; // USE COST MODEL - Insert G_{lG,rG} and elete lF, leftmost root node in F_{lF,rF}.
            }
            minCost = case1;

            case2SRightIndex = fn[lG];
            case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*G)[lG]->getLabel()) : costModel_.ins((*G)[lG]->getLabel()));

            /*if(DEBUG) {
              ou << "case2 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "] = " << to_string(s[case2SLeftIndex][case2SRightIndex]) << " + ";
              if(swap) ou << "delete " << (*G)[lG]->getLabel() << endl;
              else ou << "insert " << (*G)[lG]->getLabel() << endl;
            }*/
            
            if(case2 < minCost) {
              minCost = case2;
            }


            case3 = swap ? delta[lG][lF] : delta[lF][lG];
            /*if(DEBUG) {
              if(swap) {
                ou << "case3 = delta[" << to_string(lG) << ", " << to_string(lF) << "] = " << to_string(delta[lG][lF]) << endl;
              } else {
                ou << "case3 = delta[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(delta[lF][lG]) << endl;
              }
            }*/
            
            if (case3 < minCost) {
              switch(case3_case) {
                case 1: 
                  case3SRightIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += s[case3SLeftIndex][case3SRightIndex]; 
                  /*if(DEBUG) {
                    ou << "case3 += s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "] = " << to_string(s[case3SLeftIndex][case3SRightIndex]) << endl;
                  }*/
                  break;
                
                case 2: 
                  case3 += GcurrentForestCost - (swap ? (G)->preL_to_sumDelCost[lG] : (G)->preL_to_sumInsCost[lG]); 
                  /*if(DEBUG) {
                    ou << "case3 += " << "GcurrentForestCost - G(lG) = ";
                    ou << to_string(GcurrentForestCost) << " - ";
                    if(swap) ou << (G)->preL_to_sumDelCost[lG] << endl;
                    else ou << (G)->preL_to_sumInsCost[lG] << endl;
                  }*/
                  break; // USE COST MODEL - Insert G_{lG,rG}-G_lG.
                
                case 3: 
                  case3TLeftIndex = fn[lG + (*G)[lG]->getSubTreeSize() - 1];
                  case3 += t[case3TLeftIndex][case3TRightIndex]; 
                  /*if(DEBUG) {
                    ou << "case3 += t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "] = " << to_string(t[case3TLeftIndex][case3TRightIndex]) << endl;
                  }*/

                  break;
              }
              
              if (case3 < minCost) {
                case3 += (swap ? costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel()) : costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel())); // USE COST MODEL - Rename the leftmost root nodes in F_{lF,rF} and G_{lG,rG}.
                /*if(DEBUG) {
                  if(swap) ou << "case3 += " << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
                  else ou << "case3 += " << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
                }*/
                if (case3 < minCost) {
                  minCost = case3;
                }
              }
            }
            if(DEBUG) {
              ou << "s[" << to_string(lF) << ", " << to_string(lG) << "] = " << to_string(minCost) << endl;
            }

            /*if(DEBUG) {
              ou << "case1 = " << case1 << endl;
              ou << "case2 = " << case2 << endl;
              ou << "case3 = " << case3 << endl;
            }*/
            s[lF][lG] = minCost;
            if(mark) {
              forestdist[lF][lG][rG] = minCost;
              if(DEBUG) {
                ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
              }
            }
            dist = minCost;
            lG = ft[lG];
          }
        }
        /*if(DEBUG) {
          ou << "rGminus1_in_preR = " << to_string(rGminus1_in_preR) << " rG = " << to_string(rG) << " parent_of_rG_in_preL = " << to_string(parent_of_rG_in_preL) << " parent_of_rG_in_preR = " << to_string(parent_of_rG_in_preR) << endl;
        }*/
        if(rGminus1_in_preR == parent_of_rG_in_preR && rGminus1_in_preR != 0x7fffffff) {
            if (hasLeftPart) {
              if(swap) {
                delta[parent_of_rG_in_preL][endPathNode] = s[lFLast + 1][rGminus1_in_preL + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(endPathNode) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }*/
              } else {
                delta[endPathNode][parent_of_rG_in_preL] = s[lFLast + 1][rGminus1_in_preL + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(endPathNode) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }*/
              }
            }

            if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {
              if (swap) {
                delta[parent_of_rG_in_preL][parent_of_endPathNode_preL] = s[lFLast][rGminus1_in_preL + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_rG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }*/
              } else {
                delta[parent_of_endPathNode_preL][parent_of_rG_in_preL] = s[lFLast][rGminus1_in_preL + 1];
                /*if(DEBUG) {
                  ou << "save to delta[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                }*/
              }
            }

          for (int lF = lFFirst; lF >= lFLast; lF--) {
            q[lF] = s[lF][parent_of_rG_in_preL + 1];
            /*if(DEBUG) {
              ou << "q[" << to_string(lF) << "] = " << "s[" << to_string(lF) << ", " << to_string(parent_of_rG_in_preL + 1) << "]" << endl;
            }*/
          }
        }

        for (int lG = lGFirst; lG >= lGLast; lG = ft[lG]) {
          t[lG][rG] = s[lFLast][lG];
          /*if(DEBUG) {
           ou << "t[" << to_string(lG) << ", " << to_string(rG) << "] = " << "s[" << to_string(lFLast) << ", " << to_string(lG) << "]" << endl;
          }*/
        }
      }
    }
    rF = endPathNode_in_preR;//in D' loop
    startPathNode = endPathNode;
    endPathNode = (*F)[endPathNode] ->getParent() == NULL? -1 : (*F)[endPathNode] ->getParent()->getID(); 
    endPathNode_in_preR = F->preL_to_preR[endPathNode];
  }
  return dist;
};


void TreeComparison::updateFnArray(int lnForNode, int node, int currentSubtreePreL) {
    if (lnForNode >= currentSubtreePreL) {
      fn[node] = fn[lnForNode];//the last leaf node whose next leaf is lnForNode
      fn[lnForNode] = node;// fn[lnfornode] points to the start point
     /* if(DEBUG) {
      	ou << "fn[" << to_string(node) << "] = fn[" << to_string(lnForNode) << "] = " << to_string(fn[lnForNode]) << endl; 
      	ou << "fn[" << to_string(lnForNode) << "] = " << to_string(node) << endl;
      }*/
    } else {
      fn[node] = fn[fn_ft_length - 1];
      fn[fn_ft_length - 1] = node;
      /*if(DEBUG) {
      	ou << "O fn[" << to_string(node) << "] = fn[" << to_string(fn_ft_length - 1) << "] = " << fn[fn_ft_length - 1] << endl; 
      	ou << "O fn[" << to_string(fn_ft_length - 1) << "] = " << to_string(node) << endl;
      }*/
    }
}

void TreeComparison::updateFtArray(int lnForNode, int node) {
    ft[node] = lnForNode;
    /*if(DEBUG) {
    	ou << "ft[" << to_string(node) << "] = " << to_string(lnForNode) << endl;
    }*/
    if(fn[node] > -1) {
      ft[fn[node]] = node;
     /* if(DEBUG) {
      	ou << "ft[fn[" << to_string(node) << "]] = " << to_string(node) << endl;  
      }*/
    }
}


int TreeComparison::getPathType(Tree* tree, Node* node, int leaf) {
	if(tree->preL_to_lid[node->getID()] == leaf) return 0;
	else if(tree->preL_to_rid[node->getID()] == leaf) return 1;
	else return 2;
}

void TreeComparison::strategyComputation() {
	vector<Node*> preA = A_->getPreL();
	vector<Node*> preB = B_->getPreL();

	free(preA[0], preB[0]);
	//Strategy** S = APTED_ComputeOptStrategy_postL();
	if(DEBUG) {
		ou << "RESULT" << endl;
		for(int i = 0; i < treeSizeA; i++) {
			for(int j = 0; j < treeSizeB; j++) {
				if(&FreeStrategies[i][j] != NULL) {
					ou << FreeStrategies[i][j].getLeaf() << " in ";
					if(FreeStrategies[i][j].getTreeToDecompose() == 0) ou << "A ";
					else ou << "B ";
				}
			}
			ou << endl;
		}



	/*	ou << "RESULT 2" << endl;

		for(int i = 0; i < treeSizeA; i++) {
			for(int j = 0; j < treeSizeB; j++) {
				if(&S[i][j] != NULL) {
					ou << S[i][j].getLeaf() << " in ";
					if(S[i][j].getTreeToDecompose() == 0) ou << "A ";
					else ou << "B ";
				}
			}
			ou << endl;
		}*/

    ou << "Free" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << Free[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "LeftA" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << LeftA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "RightA" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << RightA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "AllA" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << AllA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "LeftB" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << LeftB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "RightB" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << RightB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;


    ou << "AllB" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << AllB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

	}
};

void TreeComparison::strategyComputation_compressed() {
  vector<Node*> preA = cA_->getPreL();
  vector<Node*> preB = cB_->getPreL();

  free(preA[0], preB[0]);
  //Strategy** S = APTED_ComputeOptStrategy_postL();
  if(DEBUG) {
    ou << "RESULT" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        if(&FreeStrategies[i][j] != NULL) {
          ou << FreeStrategies[i][j].getLeaf() << " in ";
          if(FreeStrategies[i][j].getTreeToDecompose() == 0) ou << "A ";
          else ou << "B ";
        }
      }
      ou << endl;
    }



  /*  ou << "RESULT 2" << endl;

    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        if(&S[i][j] != NULL) {
          ou << S[i][j].getLeaf() << " in ";
          if(S[i][j].getTreeToDecompose() == 0) ou << "A ";
          else ou << "B ";
        }
      }
      ou << endl;
    }*/

    ou << "Free" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << Free[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "LeftA" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << LeftA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "RightA" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << RightA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "AllA" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << AllA[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "LeftB" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << LeftB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

    ou << "RightB" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << RightB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;


    ou << "AllB" << endl;
    for(int i = 0; i < compressedTreeSizeA; i++) {
      for(int j = 0; j < compressedTreeSizeB; j++) {
        ou << AllB[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;

  }
};

float TreeComparison::getTreeDistance_compressed(void) {
  vector<Node*> preA = cA_->getPreL();
  vector<Node*> preB = cB_->getPreL();
  counter = 0;
  treeDist = gted_compressed(preA[0], preB[0]);
  if(DEBUG) {
    ou << "delta Result" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << delta[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;
  }
  if(DEBUG) {
    cout << "counter = " << counter << " Free[0][0] = " << Free[0][0] << endl;
  }
  //ou.close();
  return treeDist;
};



float TreeComparison::getTreeDistance(void) {
  vector<Node*> preA = A_->getPreL();
  vector<Node*> preB = B_->getPreL();
  counter = 0;
  treeDist = gted(preA[0], preB[0]);
  if(DEBUG) {
    ou << "delta Result" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << delta[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;
  }
  if(DEBUG) {
    cout << "counter = " << counter << " Free[0][0] = " << Free[0][0] << endl;
  }
  //ou.close();
  return treeDist;
};

float TreeComparison::getTreeDistance_ND(void) {
  vector<Node*> preA = A_->getPreL();
  vector<Node*> preB = B_->getPreL();
  counter = 0;
  treeDist = gted_ND(preA[0], preB[0]);
  if(DEBUG) {
    ou << "delta Result" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << delta[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;
  }
  if(DEBUG) {
    cout << "counter = " << counter << " Free[0][0] = " << Free[0][0] << endl;
  }
  //ou.close();
  return treeDist;
};

float TreeComparison::getTreeDistance_RR(void) {
  vector<Node*> preA = A_->getPreL();
  vector<Node*> preB = B_->getPreL();
  counter = 0;
  map->init();
  for(int i = 0; i < treeSizeA; i++) {
    for(int j = 0; j < treeSizeB; j++) {
      delta[i][j] = 0.0f;
    }
  }
  treeDist = spfRR(preA[0], preB[0], A_->preL_to_rid[preA[0]->getID()], false);
  if(DEBUG) {
    ou << "delta Result RR" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << delta[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;
  }
  if(DEBUG) {
    cout << "counter = " << counter << " Free[0][0] = " << Free[0][0] << endl;
  }
  ou.close();
  return treeDist;
};

float TreeComparison::getTreeDistance_LL(void) {
  vector<Node*> preA = A_->getPreL();
  vector<Node*> preB = B_->getPreL();
  counter = 0;
  map->init();
  for(int i = 0; i < treeSizeA; i++) {
    for(int j = 0; j < treeSizeB; j++) {
      delta[i][j] = 0.0f;
    }
  }
  treeDist = spfLL(preA[0], preB[0], A_->preL_to_lid[preA[0]->getID()], false);
  if(DEBUG) {
    ou << "delta Result LL" << endl;
    for(int i = 0; i < treeSizeA; i++) {
      for(int j = 0; j < treeSizeB; j++) {
        ou << delta[i][j] << " ";
      }
      ou << endl;
    }
    ou << endl;
  }
  if(DEBUG) {
    cout << "counter = " << counter << " Free[0][0] = " << Free[0][0] << endl;
  }
  //ou.close();
  return treeDist;
};


int TreeComparison::free(Node* a, Node* b) {
	if(Free[a->getID()][b->getID()] != -1) return Free[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	vector<Node*> childrenB = b->getChildren();
	Strategy freeS;
	int min = INT_MAX;
	int freeSumA = 0;
	vector<int> childrenSizeSumA;
	int freeSumB = 0;
	vector<int> childrenSizeSumB;


	if(childrenA.empty()) {
		//int left = b->getRightmostForestNum();
		//int right = b->getLeftmostForestNum();
    int num = b->getSubTreeSize();
		/*if(min > right) {
			min = right;
			freeS.setDirection(0);
			freeS.setTreeToDecompose(0);
			freeS.setKeyNode(a->getID());
			freeS.setLeaf(a->getID());
		} else if(min > left){
			min = left;
			freeS.setDirection(1);
			freeS.setTreeToDecompose(0);
			freeS.setKeyNode(a->getID());
			freeS.setLeaf(a->getID());
		}*/

    min = num;
    freeS.setDirection(0);
    freeS.setTreeToDecompose(0);
    freeS.setKeyNode(a->getID());
    freeS.setLeaf(a->getID());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(a->getID()) << " in Tree A #Subproblem: " << to_string(min) << " Direction: Left";
		}
	}
	if(childrenB.empty()) {
		//int left = a->getRightmostForestNum();
		//int right = a->getLeftmostForestNum();
    int num = a->getSubTreeSize();
		/*if(min > right) {
			min = right;
			freeS.setDirection(0);
			freeS.setTreeToDecompose(1);
			freeS.setKeyNode(b->getID());
			freeS.setLeaf(b->getID());
		} else if(min > left){
			min = left;
			freeS.setDirection(1);
			freeS.setTreeToDecompose(1);
			freeS.setKeyNode(b->getID());
			freeS.setLeaf(b->getID());
		}*/

    min = num;
    freeS.setDirection(0);
    freeS.setTreeToDecompose(1);
    freeS.setKeyNode(b->getID());
    freeS.setLeaf(b->getID());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(b->getID()) << " in Tree B #Subproblem: " << to_string(min) << " Direction: Left";
		}
	} 

	int prev = 0;
	if(!childrenA.empty()) {
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
	}
	prev = 0;
	if(!childrenB.empty()) {
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
	}

  if(DEBUG) {
    ou << "free[" << a->getID() << "][" << b->getID() << "] freeSumA = " << freeSumA << " freeSumB = " << freeSumB << endl;
  }

	if(!childrenA.empty()) {
		int aleftmost = freeSumA - free(childrenA[0], b) + leftA(childrenA[0], b) + b->getLeftmostForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenA[0]->getID()) << "(leftmost) in Tree A #Subproblem: " << to_string(aleftmost) << " Direction: Right" << endl;
		}
		if(min > aleftmost) {
			min = aleftmost;
			freeS.setKeyNode(childrenA[0]->getID());
			freeS.setTreeToDecompose(0);
			freeS.setDirection(0);
			freeS.setLeaf(LeftAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
		}
    /*if(DEBUG) {
      ou << "free(" << childrenA[0]->getID() << ", " << b->getID() << ")(" << Free[childrenA[0]->getID()][b->getID()] << ") leftA(" << childrenA[0]->getID() <<  ", " << b->getID() << ")(" << LeftA[childrenA[0]->getID()][b->getID()] << ") " << " b->getLeftmostForestNum() = " << b->getLeftmostForestNum() << " a->getSubTreeSize() - childrenA[0]->getSubTreeSize()" << to_string(a->getSubTreeSize() - childrenA[0]->getSubTreeSize()) << endl;
      ou << "aleftmost = " << aleftmost << endl;
    }*/
		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getRightmostForestNum() * (1 + childrenSizeSumA[i - 1])  + b->getSpecialForestNum() * (childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]);//direction left = first left decompose then right decompose
			int right = b->getLeftmostForestNum() * (1 + childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]) + b->getSpecialForestNum() * (childrenSizeSumA[i - 1]);//direction right = first right decompse then left decompose
			/*if(DEBUG) {
        ou << "1 + childrenSizeSumA[" << i << "] = " << 1 + childrenSizeSumA[i - 1] << " childrenSizeSumA[" << childrenA.size() - 1 << "] - childrenSizeSumA[" << i << "] = " << childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i] << endl;
        ou << "left = " << left << endl;
        ou << "1 + childrenSizeSumA[" << childrenA.size() - 1 << "] - childrenSizeSumA[" << i << "] = " << 1 + childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i] << " childrenSizeSumA[" << i - 1 << "] = " << childrenSizeSumA[i - 1] << endl;
        ou << "right = " << right << endl;
      }*/
      //int sum = prefix + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
      int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				freeS.setKeyNode(childrenA[i]->getID());
				freeS.setTreeToDecompose(0);
				freeS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)freeS.setDirection(0);
				else freeS.setDirection(1);
			}
			if(DEBUG) {
				ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
				ou << "If select " << to_string(childrenA[i]->getID()) << " in Tree A #Subproblem: " << to_string(sum) << " Direction: ";
				if(left > right) ou << "Right" << endl;
				else ou << "Left" << endl;
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + rightA(childrenA[childrenA.size() - 1], b) + b->getRightmostForestNum() * (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		/*if(DEBUG) {
      ou << "arightmost = " << arightmost << endl;
    }*/
    if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenA[childrenA.size() - 1]->getID()) << "(rightmost) in Tree A #Subproblem: " << to_string(arightmost) << " Direction: Left" << endl;
		}
		if(min > arightmost) {
			min = arightmost;
			freeS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			freeS.setLeaf(RightAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			freeS.setTreeToDecompose(0);
			freeS.setDirection(1);
		}
	}

	if(!childrenB.empty()) {
		int bleftmost = freeSumB - free(a, childrenB[0]) + leftB(a, childrenB[0]) + a->getLeftmostForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenB[0]->getID()) << "(leftmost) in Tree B #Subproblem: " << to_string(bleftmost) << " Direction: Right" << endl;
		}
		if(min > bleftmost) {
			min = bleftmost;
			freeS.setKeyNode(childrenB[0]->getID());
			freeS.setLeaf(LeftBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			freeS.setTreeToDecompose(1);
			freeS.setDirection(0);
		}
    /*if(DEBUG) {
      ou << "bleftmost = " << bleftmost << endl;
    }*/
		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getRightmostForestNum() * (1 + childrenSizeSumB[i - 1]) + a->getSpecialForestNum() * (childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]);
			int right = a->getLeftmostForestNum() * (1 + childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]) + a->getSpecialForestNum() * (childrenSizeSumB[i - 1]);
      /*if(DEBUG) {
        ou << "1 + childrenSizeSumB[" << i << "] = " << 1 + childrenSizeSumB[i - 1] << " childrenSizeSumB[" << childrenB.size() - 1 << "] - childrenSizeSumB[" << i << "] = " << childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i] << endl;
        ou << "left = " << left << endl;
        ou << "1 + childrenSizeSumB[" << childrenB.size() - 1 << "] - childrenSizeSumB[" << i << "] = " << 1 + childrenSizeSumB[childrenA.size() - 1] - childrenSizeSumB[i] << " childrenSizeSumB[" << i - 1 << "] = " << childrenSizeSumB[i - 1] << endl;
        ou << "right = " << right << endl;
      }*/
			//int sum = prefix + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
      int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				freeS.setKeyNode(childrenB[i]->getID());
				freeS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				freeS.setTreeToDecompose(1);
				if(left > right)freeS.setDirection(0);
				else freeS.setDirection(1);
			}
			if(DEBUG) {
				ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
				ou << "If select " << to_string(childrenB[i]->getID()) << " in Tree B #Subproblem: " << to_string(sum) << " Direction: ";
				if(left > right) ou << "Right" << endl;
				else ou << "Left" << endl;
			}
		}

		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + rightB(a, childrenB[childrenB.size() - 1]) + a->getRightmostForestNum() * (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		/*if(DEBUG) {
      ou << "brightmost = " << brightmost << endl;
    }*/
    if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenB[childrenB.size() - 1]->getID()) << "(rightmost) in Tree B #Subproblem: " << to_string(brightmost) << " Direction: Left" << endl;
		}
		if(min > brightmost) {
			min = brightmost;
			freeS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			freeS.setLeaf(RightBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			freeS.setTreeToDecompose(1);
			freeS.setDirection(1);
		}
	}
	Free[a->getID()][b->getID()] = min;
  if(DEBUG) {
    ou << "Free[" << a->getID() << "][" << b->getID() << "] = " << Free[a->getID()][b->getID()] << endl;
  }
	FreeStrategies[a->getID()][b->getID()] = freeS;
	if(DEBUG) {
		ou << "FreeS(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
		ou << FreeStrategies[a->getID()][b->getID()].toString() << endl;
	}
	return min;
};

int TreeComparison::leftA(Node* a, Node* b) {
	if(LeftA[a->getID()][b->getID()] != -1) return LeftA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy leftAS;
	if(childrenA.empty()) {
		min = b->getLeftmostForestNum();
		leftAS.setKeyNode(a->getID());
		leftAS.setLeaf(a->getID());
		leftAS.setTreeToDecompose(0);
		leftAS.setDirection(0);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + leftA(childrenA[0], b) + b->getLeftmostForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		/*if(DEBUG) {
			ou << "aleftmost = " << freeSumA - free(childrenA[0], b) + leftA(childrenA[0], b) << " + " << b->getLeftmostForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize()) << " = " << aleftmost << endl;
		}*/
		if(min >aleftmost) {
			min = aleftmost;
			leftAS.setKeyNode(childrenA[0]->getID());
			leftAS.setLeaf(LeftAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			leftAS.setTreeToDecompose(0);
			leftAS.setDirection(0);
		}

		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			int right = b->getLeftmostForestNum() * (1 + childrenSizeSumA[childrenSizeSumA.size() - 1] - childrenSizeSumA[i]) + b->getSpecialForestNum() * (childrenSizeSumA[i - 1]);
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				leftAS.setKeyNode(childrenA[i]->getID());
				leftAS.setTreeToDecompose(0);
				leftAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)leftAS.setDirection(0);
				else leftAS.setDirection(1);
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + allA(childrenA[childrenA.size() - 1], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		if(min > arightmost) {
			min = arightmost;
			leftAS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			leftAS.setLeaf(AllAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			leftAS.setTreeToDecompose(0);
			leftAS.setDirection(1);
		}
	}
	/*if(DEBUG) {
		ou << "LeftA[" << to_string(a->getID()) << ", " << to_string(b->getID()) << "] = " << to_string(min) << endl;
	}*/
	LeftA[a->getID()][b->getID()] = min;
	LeftAStrategies[a->getID()][b->getID()] = leftAS;
	return min;

};

int TreeComparison::rightA(Node* a, Node* b) {
	if(RightA[a->getID()][b->getID()] != -1) return RightA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy rightAS;
	if(childrenA.empty()) {
		min = b->getRightmostForestNum();
		rightAS.setKeyNode(a->getID());
		rightAS.setLeaf(a->getID());
		rightAS.setTreeToDecompose(0);
		rightAS.setDirection(1);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + allA(childrenA[0], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(min >aleftmost) {
			min = aleftmost;
			rightAS.setKeyNode(childrenA[0]->getID());
			rightAS.setLeaf(AllAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			rightAS.setTreeToDecompose(0);
			rightAS.setDirection(0);
		}

		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getRightmostForestNum() * (1 + childrenSizeSumA[i - 1]) + b->getSpecialForestNum() * (childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]);
			int right = b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				rightAS.setKeyNode(childrenA[i]->getID());
				rightAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				rightAS.setTreeToDecompose(0);
				if(left > right)rightAS.setDirection(0);
				else rightAS.setDirection(1);
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + rightA(childrenA[childrenA.size() - 1], b) + b->getRightmostForestNum()* (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		if(min >= arightmost) {
			min = arightmost;
			rightAS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			rightAS.setLeaf(RightAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			rightAS.setTreeToDecompose(0);
			rightAS.setDirection(1);
		}
	}
	RightA[a->getID()][b->getID()] = min;
	RightAStrategies[a->getID()][b->getID()] = rightAS;
	/*if(DEBUG) {
		ou << "RightA[" << to_string(a->getID()) << ", " << to_string(b->getID()) << "] set to " << to_string(rightAS.getLeaf()) << endl;
	}*/
	return min;
};

int TreeComparison::allA(Node* a, Node* b) {
	if(AllA[a->getID()][b->getID()] != -1) return AllA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy allAS;
	if(childrenA.empty()) {
		min = b->getSpecialForestNum();
		allAS.setKeyNode(a->getID());
		allAS.setLeaf(a->getID());
		allAS.setTreeToDecompose(0);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + allA(childrenA[0], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(min >aleftmost) {
			min = aleftmost;
			allAS.setKeyNode(childrenA[0]->getID());
			allAS.setLeaf(AllAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			allAS.setTreeToDecompose(0);
		}

		for(int i = 1; i < childrenA.size(); i++) {
			int sum = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			if(min > sum) {
				min = sum;
				allAS.setKeyNode(childrenA[i]->getID());
				allAS.setTreeToDecompose(0);
				allAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)allAS.setDirection(0);
				else allAS.setDirection(1);
			}
		}
	}
	AllA[a->getID()][b->getID()] = min;
	AllAStrategies[a->getID()][b->getID()] = allAS;
	return min;
};

int TreeComparison::leftB(Node* a, Node* b) {
	if(LeftB[a->getID()][b->getID()] != -1) return LeftB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy leftBS;
	if(childrenB.empty()) {
		min = a->getLeftmostForestNum();
		leftBS.setKeyNode(b->getID());
		leftBS.setLeaf(b->getID());
		leftBS.setTreeToDecompose(1);
		leftBS.setDirection(0);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + leftB(a, childrenB[0]) + a->getLeftmostForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(min >bleftmost) {
			min = bleftmost;
			leftBS.setKeyNode(childrenB[0]->getID());
			leftBS.setLeaf(LeftBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			leftBS.setTreeToDecompose(1);
			leftBS.setDirection(0);
		}

		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			int right = a->getLeftmostForestNum() * (1 + childrenSizeSumB[childrenSizeSumB.size() - 1] - childrenSizeSumB[i]) + a->getSpecialForestNum() * (childrenSizeSumB[i - 1]);
			int sum;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				leftBS.setKeyNode(childrenB[i]->getID());
				leftBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				leftBS.setTreeToDecompose(1);
				if(left > right)leftBS.setDirection(0);
				else leftBS.setDirection(1);
			}
		}
		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + allB(a, childrenB[childrenB.size() - 1]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		if(min > brightmost) {
			min = brightmost;
			leftBS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			leftBS.setLeaf(AllBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			leftBS.setTreeToDecompose(1);
			leftBS.setDirection(1);
		}
	}
	/*if(DEBUG) {
		ou << "LeftB[" << to_string(a->getID()) << ", " << to_string(b->getID()) << "] = " << to_string(min) << endl;
	}*/
	LeftB[a->getID()][b->getID()] = min;
	LeftBStrategies[a->getID()][b->getID()] = leftBS;
	return min;
};

int TreeComparison::rightB(Node* a, Node* b) {
	if(RightB[a->getID()][b->getID()] != -1) return RightB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy rightBS;
	if(childrenB.empty()) {
		min = a->getRightmostForestNum();
		rightBS.setKeyNode(b->getID());
		rightBS.setLeaf(b->getID());
		rightBS.setTreeToDecompose(1);
		rightBS.setDirection(1);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + allB(a, childrenB[0]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		/*if(DEBUG) {
			ou << "bleftmost = " << bleftmost << endl;
		}*/
		if(min > bleftmost) {
			min = bleftmost;
			rightBS.setKeyNode(childrenB[0]->getID());
			rightBS.setLeaf(AllBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			rightBS.setTreeToDecompose(1);
			rightBS.setDirection(0);
		}

		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getRightmostForestNum() * (1 + childrenSizeSumB[i - 1]) + a->getSpecialForestNum() * (childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]);
			int right = a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				rightBS.setKeyNode(childrenB[i]->getID());
				rightBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				rightBS.setTreeToDecompose(1);
				if(left > right)rightBS.setDirection(0);
				else rightBS.setDirection(1);
			}
			/*if(DEBUG) {
				ou << "childrenB " << childrenB[i]->getID() << " = " << sum << endl; 
			}*/
		}
		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + rightB(a, childrenB[childrenB.size() - 1]) + a->getRightmostForestNum()* (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		/*if(DEBUG) {
			ou << "brightmost = " << freeSumB - free(a, childrenB[childrenB.size() - 1]) + rightB(a, childrenB[childrenB.size() - 1]) << " + " << a->getRightmostForestNum()* (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize()) << " = " << brightmost << endl;
		}*/
		if(min >= brightmost) {
			min = brightmost;
			rightBS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			rightBS.setLeaf(RightBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			rightBS.setTreeToDecompose(1);
			rightBS.setDirection(1);
		}
	}
	RightB[a->getID()][b->getID()] = min;
	RightBStrategies[a->getID()][b->getID()] = rightBS;
	/*if(DEBUG) {
		ou << "RightB(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") = " << to_string(min) << endl;
	}*/
	return min;
};

int TreeComparison::allB(Node* a, Node* b) {
	if(AllB[a->getID()][b->getID()] != -1) return AllB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy allBS;
	if(childrenB.empty()) {
		min = a->getSpecialForestNum();
		allBS.setLeaf(b->getID());
		allBS.setKeyNode(b->getID());
		allBS.setTreeToDecompose(1);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + allB(a, childrenB[0]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(min >bleftmost) {
			min = bleftmost;
			allBS.setKeyNode(childrenB[0]->getID());
			allBS.setLeaf(AllBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			allBS.setTreeToDecompose(1);
		}

		for(int i = 1; i < childrenB.size(); i++) {
			int sum = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			if(min > sum) {
				min = sum;
				allBS.setKeyNode(childrenB[i]->getID());
				allBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				allBS.setTreeToDecompose(1);
				if(left > right)allBS.setDirection(0);
				else allBS.setDirection(1);
			}
		}
	}
	AllB[a->getID()][b->getID()] = min;
	AllBStrategies[a->getID()][b->getID()] = allBS;
	return min;
};

int TreeComparison::getCounter(void) {
  return counter;
};
TreeMap* TreeComparison::getTreeMap(void) {
  map = new TreeMap(A_, B_, costModel_);
  //gteo((*A_)[0], (*B_)[0]);
  gteo_LL((*A_)[0], (*B_)[0]);
  return map;
};

void TreeComparison::gteo_LL(Node* a, Node* b) {
  Tree *F = A_;
  Tree *G = B_;
  if(DEBUG) {
    ou << "gteo_LL(" << a->getID() << ", " << b->getID() << ")" << endl;
  }

  if(a->getSubTreeSize() == 1 || b->getSubTreeSize() == 1) {
    F = A_;
    G = B_;
    ou << "spf1" << endl;
    float min = FLT_MAX;
    int from = 0;
    int to = 0;
    if(a->getSubTreeSize() == 1) {
      int a_in_preL = a->getID();
      for(int i = 0; i < b->getSubTreeSize(); i++) {
        float cost = costModel_.ren((*F)[a_in_preL]->getLabel(), (*G)[i]->getLabel()) - costModel_.del((*F)[a_in_preL]->getLabel());
        if(cost < min) {
          from = a_in_preL;
          to = i;
          min = cost;
        }
      }
    } else if(b->getSubTreeSize() == 1) {
      int b_in_preL = b->getID();
      for(int i = 0; i < a->getSubTreeSize(); i++) {
        float cost = costModel_.ren((*F)[i]->getLabel(), (*G)[b->getID()]->getLabel()) - costModel_.ins((*G)[b_in_preL]->getLabel());
        if(cost < min) {
          from = i;
          to = b_in_preL;
          min = cost;
        }
      }
    }
    map->setMap(from, to);
    return;
  }

  float** forestdist;
  forestdist = new float*[a->getSubTreeSize() + 1];
  for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
    forestdist[i] = new float[b->getSubTreeSize() + 1];
  }
  treeEditDist(a, b, forestdist, false, true);

  int a_leftmost_leaf_in_preL = F->preL_to_lid[a->getID()];
  int b_leftmost_leaf_in_preL = G->preL_to_lid[b->getID()];

  int a_leftmost_leaf_in_postL = F->preL_to_postL[a_leftmost_leaf_in_preL];
  int b_leftmost_leaf_in_postL = G->preL_to_postL[b_leftmost_leaf_in_preL];

  int aoff = a_leftmost_leaf_in_postL - 1;//consider gap
  int boff = b_leftmost_leaf_in_postL - 1;//consider gap

  if(DEBUG) {
    ou << "forestdist" << endl;
    for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
      for(int j = 0; j < b->getSubTreeSize() + 1; j++) {
        ou << forestdist[i][j] << " ";
      }
      ou << endl;
    } 
  }

  int i = F->preL_to_postL[a->getID()] - aoff;
  int j = G->preL_to_postL[b->getID()] - boff;

  int FcurrentForestSize = a->getSubTreeSize();
  int GcurrentForestSize = b->getSubTreeSize();

  if(a->getID() != 0 && b->getID() != 0) {
    i--;
    j--;
    FcurrentForestSize--;
    GcurrentForestSize--;
  }

  while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {

    int i_minus1_in_preL = F->postL_to_preL[i + aoff];
    int j_minus1_in_preL = G->postL_to_preL[j + boff];
    bool isFTree = F->preL_to_lid[i_minus1_in_preL] == F->preL_to_lid[a->getID()];
    bool isGTree = G->preL_to_lid[j_minus1_in_preL] == G->preL_to_lid[b->getID()];
    float da = forestdist[i - 1][j] + costModel_.del((*F)[i_minus1_in_preL]->getLabel());
    float db = forestdist[i][j - 1] + costModel_.ins((*G)[j_minus1_in_preL]->getLabel());
    int size_of_i_minus1_in_preL = (*F)[i_minus1_in_preL]->getSubTreeSize();
    int size_of_j_minus1_in_preL = (*G)[j_minus1_in_preL]->getSubTreeSize();
    float dc = forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] + delta[i_minus1_in_preL][j_minus1_in_preL] + costModel_.ren((*F)[i_minus1_in_preL]->getLabel(), (*G)[j_minus1_in_preL]->getLabel());
       
    if(DEBUG) {
      ou << "i = " << i_minus1_in_preL << " j = " << j_minus1_in_preL << endl;
      ou << "FcurrentForestSize = " << FcurrentForestSize << ", " << "GcurrentForestSize = " << GcurrentForestSize << endl;
    }

    if(DEBUG) {
      ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
      ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
      ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
      ou << "forestdist[" << i - size_of_i_minus1_in_preL << ", " << j - size_of_j_minus1_in_preL << "] = " << forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] << endl;
      ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl;    
    }
    if(da == forestdist[i][j]) {
      if(DEBUG) {
        ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ")" << " -> -" << endl;
      }
      map->setMap(i_minus1_in_preL, -1);
      i = i - 1;
      FcurrentForestSize--;
    } else if(db == forestdist[i][j]) {
      if(DEBUG) {
        ou << "- -> " << (*G)[j_minus1_in_preL] << "(" << j_minus1_in_preL << ")" << endl;
      }
      map->setMap(-1, j_minus1_in_preL);
      j = j - 1;
      GcurrentForestSize--;
    } else if(dc == forestdist[i][j]) {
      if(DEBUG) {
        ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ") -> " << (*G)[j_minus1_in_preL]->getLabel() << "(" << j_minus1_in_preL << ")" << endl;
      }
      map->setMap(i_minus1_in_preL, j_minus1_in_preL);
      if(isFTree && isGTree) {
        i = i - 1;
        j = j - 1;
        FcurrentForestSize--;
        GcurrentForestSize--;
      } else {
        if((*F)[i_minus1_in_preL]->getSubTreeSize() > 1 && (*G)[j_minus1_in_preL]->getSubTreeSize() > 1) {
          gteo_LL((*F)[i_minus1_in_preL], (*G)[j_minus1_in_preL]);
        }
        i = i - size_of_i_minus1_in_preL;
        j = j - size_of_j_minus1_in_preL;
        FcurrentForestSize -= size_of_i_minus1_in_preL;
        GcurrentForestSize -= size_of_j_minus1_in_preL;
      }
    }
  }
  while(FcurrentForestSize > 0) {
    if(DEBUG) {
      ou << "i = " << i << " j = 0" << endl;
    } 
    int i_minus1_in_preL = F->postL_to_preL[i - 1];
    map->setMap(i_minus1_in_preL, -1);
    i = i - 1;
    FcurrentForestSize--;
  }
  while(GcurrentForestSize > 0) {
    if(DEBUG) {
      ou << "i = 0 j = " << j << endl;
    }
    int j_minus1_in_preL = G->postL_to_preL[j - 1];
    map->setMap(-1, j_minus1_in_preL);
    j = j - 1;
    GcurrentForestSize--;
  }
}





/*void TreeComparison::gteo(Node* a, Node* b) {
  int pathLeaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
  int direction = FreeStrategies[a->getID()][b->getID()].getDirection();
  int treeToDecompose = FreeStrategies[a->getID()][b->getID()].getTreeToDecompose();
  Tree *F, *G;

  int pathType;// 0 left 1 right 2 special
  if(DEBUG) {
    ou << "gteo(" << a->getID() << ", " << b->getID() << ")" << endl;
    ou << "treeToDecompose = " << treeToDecompose << " pathLeaf = " << pathLeaf << endl;
  }

  if(a->getSubTreeSize() == 1 || b->getSubTreeSize() == 1) {
    F = A_;
    G = B_;
    ou << "spf1" << endl;
    float min = FLT_MAX;
    int from = 0;
    int to = 0;
    if(a->getSubTreeSize() == 1) {
      int a_in_preL = a->getID();
      for(int i = 0; i < b->getSubTreeSize(); i++) {
        float cost = costModel_.ren((*F)[a_in_preL]->getLabel(), (*G)[i]->getLabel()) - costModel_.del((*F)[a_in_preL]->getLabel());
        if(cost < min) {
          from = a_in_preL;
          to = i;
          min = cost;
        }
      }
    } else if(b->getSubTreeSize() == 1) {
      int b_in_preL = b->getID();
      for(int i = 0; i < a->getSubTreeSize(); i++) {
        float cost = costModel_.ren((*F)[i]->getLabel(), (*G)[b->getID()]->getLabel()) - costModel_.ins((*G)[b_in_preL]->getLabel());
        if(cost < min) {
          from = i;
          to = b_in_preL;
          min = cost;
        }
      }
    }
    map->setMap(from, to);
    return;
  }

  if(treeToDecompose == 0) {
    F = A_;
    G = B_;
    pathType = getPathType(A_, a, pathLeaf);
    if(DEBUG) {
      ou << "pathType = " << pathType << endl;
    }

    if(pathType == 0) {
      float** forestdist;
      forestdist = new float*[a->getSubTreeSize() + 1];
      for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
        forestdist[i] = new float[b->getSubTreeSize() + 1];
      }
      treeEditDist(a, b, forestdist, false, true);



      int a_leftmost_leaf_in_preL = F->preL_to_lid[a->getID()];
      int b_leftmost_leaf_in_preL = G->preL_to_lid[b->getID()];

      int a_leftmost_leaf_in_postL = F->preL_to_postL[a_leftmost_leaf_in_preL];
      int b_leftmost_leaf_in_postL = G->preL_to_postL[b_leftmost_leaf_in_preL];

      int aoff = a_leftmost_leaf_in_postL - 1;//consider gap
      int boff = b_leftmost_leaf_in_postL - 1;//consider gap

      if(DEBUG) {
        ou << "forestdist" << endl;
        for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
          for(int j = 0; j < b->getSubTreeSize() + 1; j++) {
            ou << forestdist[i][j] << " ";
          }
          ou << endl;
        } 
      }

      int i = F->preL_to_postL[a->getID()] - aoff;
      int j = G->preL_to_postL[b->getID()] - boff;

      int FcurrentForestSize = a->getSubTreeSize();
      int GcurrentForestSize = b->getSubTreeSize();

      if(a->getID() != 0 && b->getID() != 0) {
        i--;
        j--;
        FcurrentForestSize--;
        GcurrentForestSize--;
      }

      while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {

        int i_minus1_in_preL = F->postL_to_preL[i + aoff];
        int j_minus1_in_preL = G->postL_to_preL[j + boff];
        bool isFTree = F->preL_to_lid[i_minus1_in_preL] == F->preL_to_lid[a->getID()];
        bool isGTree = G->preL_to_lid[j_minus1_in_preL] == G->preL_to_lid[b->getID()];
        float da = forestdist[i - 1][j] + costModel_.del((*F)[i_minus1_in_preL]->getLabel());
        float db = forestdist[i][j - 1] + costModel_.ins((*G)[j_minus1_in_preL]->getLabel());
        int size_of_i_minus1_in_preL = (*F)[i_minus1_in_preL]->getSubTreeSize();
        int size_of_j_minus1_in_preL = (*G)[j_minus1_in_preL]->getSubTreeSize();
        float dc = forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] + delta[i_minus1_in_preL][j_minus1_in_preL] + costModel_.ren((*F)[i_minus1_in_preL]->getLabel(), (*G)[j_minus1_in_preL]->getLabel());
       
        if(DEBUG) {
          ou << "i = " << i_minus1_in_preL << " j = " << j_minus1_in_preL << endl;
        }

        if(DEBUG) {
          ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
          ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
          ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
          ou << "forestdist[" << i - size_of_i_minus1_in_preL << ", " << j - size_of_j_minus1_in_preL << "] = " << forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] << endl;
          ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl;    
        }
        if(da == forestdist[i][j]) {
          if(DEBUG) {
            ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ")" << " -> -" << endl;
          }
          map->setMap(i_minus1_in_preL, -1);
          i = i - 1;
          FcurrentForestSize--;
        } else if(db == forestdist[i][j]) {
          if(DEBUG) {
            ou << "- -> " << (*G)[j_minus1_in_preL] << "(" << j_minus1_in_preL << ")" << endl;
          }
          map->setMap(-1, j_minus1_in_preL);
          j = j - 1;
          GcurrentForestSize--;
        } else if(dc == forestdist[i][j]) {
          if(DEBUG) {
            ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ") -> " << (*G)[j_minus1_in_preL]->getLabel() << "(" << j_minus1_in_preL << ")" << endl;
          }
          map->setMap(i_minus1_in_preL, j_minus1_in_preL);
          if(isFTree && isGTree) {
            i = i - 1;
            j = j - 1;
            FcurrentForestSize--;
            GcurrentForestSize--;
          } else {
            if((*F)[i_minus1_in_preL]->getSubTreeSize() > 1 && (*G)[j_minus1_in_preL]->getSubTreeSize() > 1) {
              gteo((*F)[i_minus1_in_preL], (*G)[j_minus1_in_preL]);
            }
            i = i - size_of_i_minus1_in_preL;
            j = j - size_of_j_minus1_in_preL;
            FcurrentForestSize -= size_of_i_minus1_in_preL;
            GcurrentForestSize -= size_of_j_minus1_in_preL;
          }
        }
      }
      while(FcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = " << i << " j = 0" << endl;
        } 
        int i_minus1_in_preL = F->postL_to_preL[i - 1];
        map->setMap(i_minus1_in_preL, -1);
        i = i - 1;
        FcurrentForestSize--;
      }
      while(GcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = 0 j = " << j << endl;
        }
        int j_minus1_in_preL = G->postL_to_preL[j - 1];
        map->setMap(-1, j_minus1_in_preL);
        j = j - 1;
        GcurrentForestSize--;
      }
    }
    else if(pathType == 1) {
      float** forestdist;

      forestdist = new float*[a->getSubTreeSize() + 1];
      for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
        forestdist[i] = new float[b->getSubTreeSize() + 1];
      }
      revTreeEditDist(a, b, forestdist, false, true);

      int a_rightmost_leaf_in_preL = F->preL_to_rid[a->getID()];
      int b_rightmost_leaf_in_preL = G->preL_to_rid[b->getID()];

      int a_rightmost_leaf_in_postR = F->preL_to_postR[a_rightmost_leaf_in_preL];
      int b_rightmost_leaf_in_postR = G->preL_to_postR[b_rightmost_leaf_in_preL];

      int aoff = a_rightmost_leaf_in_postR - 1;//consider gap
      int boff = b_rightmost_leaf_in_postR - 1;//consider gap

      int FcurrentForestSize = a->getSubTreeSize();
      int GcurrentForestSize = b->getSubTreeSize();

      if(DEBUG) {
        ou << "forestdist" << endl;
        for(int i = 0; i < a->getSubTreeSize() + 1; i++) {
          for(int j = 0; j < b->getSubTreeSize() + 1; j++) {
            ou << forestdist[i][j] << " ";
          }
          ou << endl;
        } 
      }


      int i = F->preL_to_postR[a->getID()] - aoff;
      int j = G->preL_to_postR[b->getID()] - boff;

      if(a->getID() != 0 && b->getID() != 0) {
        i--;
        j--;
        FcurrentForestSize--;
        GcurrentForestSize--;
      }

      while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {

        int i_minus1_in_preL = F->postR_to_preL[i + aoff];
        int j_minus1_in_preL = G->postR_to_preL[j + boff];
        bool isFTree = F->preL_to_rid[i_minus1_in_preL] == F->preL_to_rid[a->getID()];
        bool isGTree = G->preL_to_rid[j_minus1_in_preL] == G->preL_to_rid[b->getID()];
        float da = forestdist[i - 1][j] + costModel_.del((*F)[i_minus1_in_preL]->getLabel());
        float db = forestdist[i][j - 1] + costModel_.ins((*G)[j_minus1_in_preL]->getLabel());
        int size_of_i_minus1_in_preL = (*F)[i_minus1_in_preL]->getSubTreeSize();
        int size_of_j_minus1_in_preL = (*G)[j_minus1_in_preL]->getSubTreeSize();
        float dc = forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] + delta[i_minus1_in_preL][j_minus1_in_preL] + costModel_.ren((*F)[i_minus1_in_preL]->getLabel(), (*G)[j_minus1_in_preL]->getLabel());
   

        if(DEBUG) {
          ou << "i = " << i_minus1_in_preL << " j = " << j_minus1_in_preL << endl;
        }

        if(DEBUG) {
          ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
          ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
          ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
          ou << "forestdist[" << i - size_of_i_minus1_in_preL << ", " << j - size_of_j_minus1_in_preL << "] = " << forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] << endl;
          ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl;          
        }
        if(da == forestdist[i][j]) {
          if(DEBUG) {
            ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ")" << " -> -" << endl;
          }
          map->setMap(i_minus1_in_preL, -1);
          i = i - 1;
          FcurrentForestSize--;
        } else if(db == forestdist[i][j]) {
          if(DEBUG) {
            ou << "- -> " << (*G)[j_minus1_in_preL] << "(" << j_minus1_in_preL << ")" << endl;
          }
          map->setMap(-1, j_minus1_in_preL);
          j = j - 1;
          GcurrentForestSize--;
        } else if(dc == forestdist[i][j]) {
          if(DEBUG) {
            ou << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ") -> " << (*G)[j_minus1_in_preL]->getLabel() << "(" << j_minus1_in_preL << ")" << endl;
          }
          map->setMap(i_minus1_in_preL, j_minus1_in_preL);
          if(isFTree && isGTree) {
            i = i - 1;
            j = j - 1;
            FcurrentForestSize--;
            GcurrentForestSize--;
          }
          else {
            i = i - size_of_i_minus1_in_preL;
            j = j - size_of_j_minus1_in_preL;
            FcurrentForestSize -= size_of_i_minus1_in_preL;
            GcurrentForestSize -= size_of_j_minus1_in_preL;
            if((*F)[i_minus1_in_preL]->getSubTreeSize() > 1 && (*G)[j_minus1_in_preL]->getSubTreeSize() > 1) {
              gteo((*F)[i_minus1_in_preL], (*G)[j_minus1_in_preL]);
            }
          }
        }
      }
      while(FcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = " << i << " j = 0" << endl;
        } 
        int i_minus1_in_preL = F->postR_to_preL[i - 1];
        map->setMap(i_minus1_in_preL, -1);
        i = i - 1;
      }
      while(GcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = 0 j = " << j << endl;
        }
        int j_minus1_in_preL = G->postR_to_preL[j - 1];
        map->setMap(-1, j_minus1_in_preL);
        j = j - 1;
      }
    }

    else if(pathType == 2) {
      if(direction == 1) {
        float*** forestdist;
        forestdist = new float**[treeSizeA+ 1];
        for(int i = 0; i < treeSizeA + 1; i++) {
          forestdist[i] = new float*[treeSizeB + 1];
          for(int j = 0; j < treeSizeB + 1; j++) {
            forestdist[i][j] = new float[treeSizeB + 1];
          }
        }
        spfA_RL(a, b, pathLeaf, pathType, forestdist, false, true);

        int lF = a->getID();
        int rF = F->preL_to_preR[lF];
        int FcurrentForestSize = a->getSubTreeSize();
        
        int lG = b->getID();
        int rG = G->preL_to_preR[lG];
        int GcurrentForestSize = b->getSubTreeSize();


        int direction = 1;//0 for right 1 for left;
        int leaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
        Node* parent = (*F)[leaf]->getParent();
        int* favouriteChild = new int[F->getTreeSize()];
        int count = 0;

        while(leaf != a->getID()) {
          favouriteChild[count++] = leaf;
          leaf = parent->getID();
          parent = (*F)[leaf]->getParent();
        }
        favouriteChild[count++] = leaf;
        int i = count - 1;
        if(!(a->getID() == 0 && b->getID() == 0)) {
          lF++;
          rF++;
          lG++;
          rG++;
          FcurrentForestSize--;
          GcurrentForestSize--;
          i--;
        }
        int prev_fav_child = favouriteChild[count - 1];
        int prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];

        if(DEBUG) {
          ou << "favouriteChild: " << endl;
          for(int i = count - 1; i >= 0; i--) {
            ou << favouriteChild[i] << " ";
          }
          ou << endl;
        }

        while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
          if(DEBUG) {
            ou << "(" << lF << ", " << rF << ", " << lG << ", " << rG << ")" << endl;
            ou << "FcurrentForestSize = " << FcurrentForestSize << endl;
            ou << "GcurrentForestSize = " << GcurrentForestSize << endl;
          }

          int rF_in_preL = F->preR_to_preL[rF];
          int rG_in_preL = G->preR_to_preL[rG];
          int lF_in_preR = F->preL_to_preR[lF];

          float da = 0;
          float db = 0;
          float dc = 0;

          bool isFTree = rF_in_preL == lF;
          bool isGTree = rG_in_preL == lG;
          int favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];

          if(lF == favouriteChild[i]) {
            direction = 0;
          }
          if(rF == favouriteChild_in_preR) {
            prev_fav_child = favouriteChild[i];
            prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];
            i--;
            direction = 1;
            favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];
          }


          bool hasRightPart = (abs(favouriteChild_in_preR - prev_fav_child_in_preR) > 1);
          bool hasLeftPart = (abs(favouriteChild[i] - prev_fav_child) > 1); 

          if(DEBUG) {
            ou << "favouriteChild = " << favouriteChild[i] << endl;
            ou << "hasRightPart = " << hasRightPart << endl;
            ou << "hasLeftPart = " << hasLeftPart << endl;
          }

          if(isGTree && isFTree) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            //dc = forestdist[lF + 1][lG + 1][rG + 1] + delta[lF][lG] +costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
            dc = delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
            if(DEBUG) {
              ou << "isGTree && isFTree" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              if(hasLeftPart) {
                ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG << "] = " << forestdist[lF + 1][lG][rG] << endl;
              } else {
                ou << "forestdist[" << rF_plus_one_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_plus_one_in_preL][lG][rG] << endl;
              }
              ou << "forestdist[" << lF << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[lF][lG + 1][rG + 1] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl;
              }
              map->setMap(lF, -1);
              rF++;
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              lG++;
              rG++;
              lF++;
              rF++;
              FcurrentForestSize--;
              GcurrentForestSize--;
            }
          }

          else if(isFTree && direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.del((*F)[lF]->getLabel());
            //da = forestdist[lF + 1][lG][rG] + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF][lG][rG + 1] + costModel_.ins((*G)[rG_in_preL]->getLabel());
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = (hasLeftPart? forestdist[lF + 1][lG][rG + size_of_rG] : forestdist[rF_plus_one_in_preL][lG][rG + size_of_rG]) + costModel_.ren((*F)[lF]->getLabel(), (*G)[rG_in_preL]->getLabel());
            //dc = forestdist[lF + 1][lG][rG + size_of_rG] + delta[lF][rG_in_preL] + costModel_.ren((*F)[lF]->getLabel(), (*G)[rG_in_preL]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }
            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl;
              }
              map->setMap(lF, -1);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[rG_in_preL] << endl; 
              }
              map->setMap(-1, rG_in_preL);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(lF, rG_in_preL);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[rG_in_preL]);
              }
              lF++;
              rF++;
              rG = rG + size_of_rG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_rG;
            }
          }

          else if(isFTree && direction == 1) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.del((*F)[lF]->getLabel());
            //da = forestdist[lF + 1][lG][rG] + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = (hasLeftPart? forestdist[lF + 1][lG + size_of_lG][rG] : forestdist[rF_plus_one_in_preL][lG + size_of_lG][rG]) + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());
            //dc = forestdist[lF + 1][lG + size_of_lG][rG] + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl;
              }
              map->setMap(lF, -1);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              } 
              map->setMap(-1, lG);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[lG]);
              }
              lF++;
              rF++;
              lG = lG + size_of_lG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_lG;
            }
          }

          else if(isGTree && direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = forestdist[rF_plus_one_in_preL][lG][rG] + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF_in_preL][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int rF_plus_size_of_rF_in_preL = F->preR_to_preL[rF + size_of_rF];
            dc = forestdist[rF_in_preL - size_of_rF][lG + 1][rG + 1] + delta[rF_in_preL][lG] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());

     
            if(DEBUG) {
              ou << "isGTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
            }

            if(forestdist[rF_in_preL][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              rF++;
              FcurrentForestSize--;
            }
            else if(forestdist[rF_in_preL][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }
            else if(forestdist[rF_in_preL][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, lG);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[lG]);
              }
              rF = rF + size_of_rF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize--;
            }
          }

          else if(isGTree && direction == 1) {
            da = (lF + 1 == favouriteChild[i] ? forestdist[rF_in_preL][lG][rG] : forestdist[lF + 1][lG][rG]) + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            dc = (lF + size_of_lF == favouriteChild[i] ? forestdist[rF_in_preL][lG + 1][rG + 1] : forestdist[lF + size_of_lF][lG + 1][rG + 1]) + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());

            if(DEBUG) {
              ou << "isGTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl; 
              }
              map->setMap(lF, -1);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[lG]);
              }
              lF = lF + size_of_lF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize--;
            }
          }

          else if(direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = forestdist[rF_plus_one_in_preL][lG][rG] + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF_in_preL][lG][rG + 1] + costModel_.ins((*G)[rG_in_preL]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            int rF_plus_size_of_rF_in_preL = F->preR_to_preL[rF + size_of_rF];
            dc = forestdist[rF_plus_size_of_rF_in_preL][lG][rG + size_of_rG] + delta[rF_in_preL][rG_in_preL] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel());

            if(DEBUG) {
              ou << "Both forest && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl;
              ou << "forestdist[" << rF_plus_one_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_plus_one_in_preL][lG][rG] << endl;
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG + 1 << "] = " << forestdist[rF_in_preL][lG][rG + 1] << endl;
              ou << "forestdist[" << rF_plus_size_of_rF_in_preL << ", " << lG << ", " << rG + size_of_rG << "] = " << forestdist[rF_plus_size_of_rF_in_preL][lG][rG + size_of_rG] << endl;
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
            }

            if(forestdist[rF_in_preL][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF_in_preL][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[rG_in_preL] << endl;
              }
              map->setMap(-1, rG_in_preL);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF_in_preL][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, rG_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[rG_in_preL]);
              }
              rF = rF + size_of_rF;
              rG = rG + size_of_rG;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize -= size_of_rG;
            }
          }


          //F and G are forest and direction = left
          else if(direction == 1) {
            da = (lF + 1 == favouriteChild[i] ? forestdist[rF_in_preL][lG][rG] : forestdist[lF + 1][lG][rG]) + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = (lF + size_of_lF == favouriteChild[i] ? forestdist[rF_in_preL][lG + size_of_lG][rG] : forestdist[lF + size_of_lF][lG + size_of_lG][rG]) + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());


            if(DEBUG) {
              ou << "Both forest && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              if(lF + 1 == favouriteChild[i]) ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
              else ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG << "] = " << forestdist[lF + 1][lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG + 1 << ", " << rG << "] = " << forestdist[lF][lG + 1][rG] << endl;
              if(lF + size_of_lF == favouriteChild[i]) ou << "forestdist[" << rF_in_preL << ", " << lG + size_of_lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG + size_of_lG][rG] << endl;
              else ou << "forestdist[" << lF + size_of_lF << ", " << lG + size_of_lG << ", " << rG << "] = " << forestdist[lF + size_of_lF][lG + size_of_lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl;
              }
              map->setMap(lF, -1);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[lG]);
              }
              lF = lF + size_of_lF;
              lG = lG + size_of_lG;
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize -= size_of_lG; 
            }
          }     
        }
      } 




      else if(direction == 0) {
        float*** forestdist; 
        forestdist = new float**[treeSizeA + 1];
        for(int i = 0; i < treeSizeA + 1; i++) {
          forestdist[i] = new float*[treeSizeB + 1];
          for(int j = 0; j < treeSizeB + 1; j++) {
            forestdist[i][j] = new float[treeSizeB + 1];
          }
        }
        spfA_LR(a, b, pathLeaf, pathType, forestdist, false, true);

        int lF = a->getID();
        int rF = F->preL_to_preR[lF];
        int FcurrentForestSize = a->getSubTreeSize();

        int lG = b->getID();
        int rG = G->preL_to_preR[lG];
        int GcurrentForestSize = b->getSubTreeSize();


        int direction = 0;//0 for right 1 for left;
        int leaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
        Node* parent = (*F)[leaf]->getParent();
        int* favouriteChild = new int[F->getTreeSize()];
        int count = 0;
        while(leaf != a->getID()) {
          favouriteChild[count++] = leaf;
          leaf = parent->getID();
          parent = (*F)[leaf]->getParent();
        }
        favouriteChild[count++] = leaf;
        int i = count - 1;
        if(!(a->getID() == 0 && b->getID() == 0)) {
          lF++;
          rF++;
          lG++;
          rG++;
          FcurrentForestSize--;
          GcurrentForestSize--;
          i--;
        } 

        int prev_fav_child = favouriteChild[count - 1];
        int prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];
        if(DEBUG) {
          ou << "favouriteChild: " << endl;
          for(int i = count - 1; i >= 0; i--) {
            ou << favouriteChild[i] << " ";
          }
          ou << endl;
        }


        while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
          if(DEBUG) {
            ou << "(" << lF << ", " << rF << ", " << lG << ", " << rG << ")" << endl;
            ou << "FcurrentForestSize = " << FcurrentForestSize << endl;
            ou << "GcurrentForestSize = " << GcurrentForestSize << endl;
          }
          int rF_in_preL = F->preR_to_preL[rF];
          int rG_in_preL = G->preR_to_preL[rG];
          int lF_in_preR = F->preL_to_preR[lF];


          float da = 0;
          float db = 0;
          float dc = 0;

          bool isFTree = rF_in_preL == lF;
          bool isGTree = rG_in_preL == lG;

          int favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];

          if(rF == favouriteChild_in_preR) {
            direction = 1;
          }
          if(lF == favouriteChild[i]) {
            prev_fav_child = favouriteChild[i];
            prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];
            i--;
            direction = 0;
            favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];
          }


          bool hasRightPart = (abs(favouriteChild_in_preR - prev_fav_child_in_preR) > 1);
          bool hasLeftPart = (abs(favouriteChild[i] - prev_fav_child) > 1); 

      
          if(isGTree && isFTree) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            //dc = forestdist[rF + 1][lG + 1][rG + 1] + delta[rF_in_preL][lG] +costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());
            dc = delta[rF_in_preL][lG] +costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());
            if(DEBUG) {
              ou << "isGTree && isFTree" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_plus_one_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_plus_one_in_preR][lG][rG] << endl;
              ou << "forestdist[" << rF << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[rF][lG + 1][rG + 1] << endl;
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              rF++;
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, lG);
              lG++;
              rG++;
              lF++;
              rF++;
              FcurrentForestSize--;
              GcurrentForestSize--;
            }
          } 

          else if(isFTree && direction == 0) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.del((*F)[rF_in_preL]->getLabel());
            //da = forestdist[rF + 1][lG][rG] + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG][rG + 1] + costModel_.ins((*G)[rG_in_preL]->getLabel());
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = (hasRightPart? forestdist[rF + 1][lG][rG + size_of_rG] : forestdist[lF_plus_one_in_preR][lG][rG + size_of_rG]) + delta[rF_in_preL][rG_in_preL] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel());
            //dc = forestdist[rF + 1][lG][rG + size_of_rG] + delta[rF_in_preL][rG_in_preL] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }
            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[rG_in_preL] << endl; 
              }
              map->setMap(-1, rG_in_preL);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, rG_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[rG_in_preL]);
              }
              lF++;
              rF++;
              rG = rG + size_of_rG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_rG;
            }
          } 

          else if(isFTree && direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.del((*F)[rF_in_preL]->getLabel());
            //da = forestdist[rF + 1][lG][rG] + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = (hasRightPart? forestdist[rF + 1][lG + size_of_lG][rG] :forestdist[lF_plus_one_in_preR][lG + size_of_lG][rG]) + delta[rF_in_preL][lG] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());
            //dc = forestdist[rF + 1][lG + size_of_lG][rG] + delta[rF_in_preL][lG] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              } 
              map->setMap(-1, lG);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, lG);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[lG]);
              }
              lF++;
              rF++;
              lG = lG + size_of_lG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_lG;
            }
          } 

          else if(isGTree && direction == 0) {
            da = (rF + 1 == favouriteChild_in_preR ? forestdist[lF_in_preR][lG][rG] : forestdist[rF + 1][lG][rG]) + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            dc = (rF + size_of_rF == favouriteChild_in_preR ? forestdist[lF_in_preR][lG][rG] : forestdist[rF + size_of_rF][lG + 1][rG + 1]) + delta[rF_in_preL][lG] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[lG]->getLabel());

     
            if(DEBUG) {
              ou << "isGTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              rF++;
              FcurrentForestSize--;
            }
            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }
            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, lG);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[lG]);
              }
              rF = rF + size_of_rF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize--;
            }
          } 

          else if(isGTree && direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = forestdist[lF_plus_one_in_preR][lG][rG] + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int lF_plus_size_of_lF_in_preR = F->preL_to_preR[lF + size_of_lF];
            dc = forestdist[lF_plus_size_of_lF_in_preR][lG + 1][rG + 1] + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());

            if(DEBUG) {
              ou << "isGTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl; 
              }
              map->setMap(lF, -1);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[lG]);
              }
              lF = lF + size_of_lF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize--;
            }
          }
          
          //F and G are forest and direction = right
          else if(direction == 0) {
            da = (rF + 1 == favouriteChild_in_preR ? forestdist[lF_in_preR][lG][rG] : forestdist[rF + 1][lG][rG]) + costModel_.del((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG][rG + 1] + costModel_.ins((*G)[rG_in_preL]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = (rF + size_of_rF == favouriteChild_in_preR ? forestdist[lF_in_preR][lG][rG] : forestdist[rF + size_of_rF][lG][rG + size_of_rG]) + delta[rF_in_preL][rG_in_preL] + costModel_.ren((*F)[rF_in_preL]->getLabel(), (*G)[rG_in_preL]->getLabel());

            if(DEBUG) {
              ou << "Both forest && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> -" << endl;
              }
              map->setMap(rF_in_preL, -1);
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[rG_in_preL] << endl;
              }
              map->setMap(-1, rG_in_preL);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[rF_in_preL]->getLabel() << " -> " << (*G)[rG_in_preL]->getLabel() << endl;
              }
              map->setMap(rF_in_preL, rG_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*F)[rF_in_preL], (*G)[rG_in_preL]);
              }
              rF = rF + size_of_rF;
              rG = rG + size_of_rG;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize -= size_of_rG;
            }
          }
          
          //F and G are forest and direction = left
          else if(direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = forestdist[lF_plus_one_in_preR][lG][rG] + costModel_.del((*F)[lF]->getLabel());
            db = forestdist[lF_in_preR][lG + 1][rG] + costModel_.ins((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            int lF_plus_size_of_lF_in_preR = F->preL_to_preR[lF + size_of_lF];
            dc = forestdist[lF_plus_size_of_lF_in_preR][lG + size_of_lG][rG] + delta[lF][lG] + costModel_.ren((*F)[lF]->getLabel(), (*G)[lG]->getLabel());


            if(DEBUG) {
              ou << "Both forest && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> -" << endl;
              }
              map->setMap(lF, -1);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << "- -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(-1, lG);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*F)[lF]->getLabel() << " -> " << (*G)[lG]->getLabel() << endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*F)[lF], (*G)[lG]);
              }
              lF = lF + size_of_lF;
              lG = lG + size_of_lG; 
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize -= size_of_lG;
            }
          }


        }
      }
    }

  } else if(treeToDecompose == 1) {
    F = B_;
    G = A_;
    pathType = getPathType(B_, b, pathLeaf);
    if(DEBUG) {
      ou << "pathType = " << pathType << endl;
    }

    if(pathType == 0) {
      float** forestdist;
      forestdist = new float*[b->getSubTreeSize() + 1];
      for(int i = 0; i < b->getSubTreeSize() + 1; i++) {
        forestdist[i] = new float[a->getSubTreeSize() + 1];
      }
      treeEditDist(b, a, forestdist, true, true);


      int a_leftmost_leaf_in_preL = F->preL_to_lid[a->getID()];
      int b_leftmost_leaf_in_preL = G->preL_to_lid[b->getID()];

      int a_leftmost_leaf_in_postL = F->preL_to_postL[a_leftmost_leaf_in_preL];
      int b_leftmost_leaf_in_postL = G->preL_to_postL[b_leftmost_leaf_in_preL];

      int aoff = a_leftmost_leaf_in_postL - 1;//consider gap
      int boff = b_leftmost_leaf_in_postL - 1;//consider gap

      int FcurrentForestSize = b->getSubTreeSize();
      int GcurrentForestSize = a->getSubTreeSize();

      if(DEBUG) {
        ou << "forestdist" << endl;
        for(int i = 0; i < b->getSubTreeSize() + 1; i++) {
          for(int j = 0; j < a->getSubTreeSize() + 1; j++) {
            ou << forestdist[i][j] << " ";
          }
          ou << endl;
        } 
      }

      int i = F->preL_to_postL[b->getID()] - boff;
      int j = G->preL_to_postL[a->getID()] - aoff;

      if(b->getID() != 0 && a->getID() != 0) {
        i--;
        j--;
        FcurrentForestSize--;
        GcurrentForestSize--;
      }

      while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
        
        int i_minus1_in_preL = F->postL_to_preL[i + boff];
        int j_minus1_in_preL = G->postL_to_preL[j + aoff];
        bool isFTree = F->preL_to_lid[i_minus1_in_preL] == F->preL_to_lid[b->getID()];
        bool isGTree = G->preL_to_lid[j_minus1_in_preL] == G->preL_to_lid[a->getID()];
        float da = forestdist[i - 1][j] + costModel_.ins((*F)[i_minus1_in_preL]->getLabel());
        float db = forestdist[i][j - 1] + costModel_.del((*G)[j_minus1_in_preL]->getLabel());
        int size_of_i_minus1_in_preL = (*F)[i_minus1_in_preL]->getSubTreeSize();
        int size_of_j_minus1_in_preL = (*G)[j_minus1_in_preL]->getSubTreeSize();

        float dc = forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] + delta[j_minus1_in_preL][i_minus1_in_preL] + costModel_.ren((*G)[j_minus1_in_preL]->getLabel(), (*F)[i_minus1_in_preL]->getLabel());
          
        if(DEBUG) {
          ou << "i = " << i_minus1_in_preL << " j = " << j_minus1_in_preL << endl;
        }

        if(DEBUG) {
          ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
          ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
          ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
          ou << "forestdist[" << i - size_of_i_minus1_in_preL << ", " << j - size_of_j_minus1_in_preL << "] = " << forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] << endl;
          ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl; 
        }
        if(da == forestdist[i][j]) {
          map->setMap(-1, i_minus1_in_preL);
          if(DEBUG) {
            ou << "- -> " << (*F)[i_minus1_in_preL]->getLabel() << endl;
          }
          i = i - 1;
          FcurrentForestSize--;
        } else if(db == forestdist[i][j]) {
          map->setMap(j_minus1_in_preL, -1);
          if(DEBUG) {
            ou << (*G)[j_minus1_in_preL]->getLabel() << " - ->" << endl;
          }
          j = j - 1;
          GcurrentForestSize--;
        } else if(dc == forestdist[i][j]) {
          map->setMap(j_minus1_in_preL, i_minus1_in_preL);
          if(DEBUG) {
            ou << (*G)[j_minus1_in_preL]->getLabel() << "(" << j_minus1_in_preL << ")" << " -> " << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ")" << endl;
          }
          if(isFTree && isGTree) {
            i = i - 1;
            j = j - 1;
            FcurrentForestSize--;
            GcurrentForestSize--;
          }
          else {
            i = i - size_of_i_minus1_in_preL;
            j = j - size_of_j_minus1_in_preL;
            FcurrentForestSize -= size_of_i_minus1_in_preL;
            GcurrentForestSize -= size_of_j_minus1_in_preL;
            if((*F)[i_minus1_in_preL]->getSubTreeSize() > 1 && (*G)[j_minus1_in_preL]->getSubTreeSize() > 1) {
              gteo((*G)[j_minus1_in_preL], (*F)[i_minus1_in_preL]);
            }
          }
        }
      }
      while(FcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = " << i << " j = 0" << endl;
        } 
        int i_minus1_in_preL = F->postL_to_preL[i - 1];
        map->setMap(-1, i_minus1_in_preL);
        i = i - 1;
        FcurrentForestSize--;
      }
      while(GcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = 0 j = " << j << endl;
        }
        int j_minus1_in_preL = G->postL_to_preL[j - 1];
        map->setMap(j_minus1_in_preL, -1);
        j = j - 1;
        GcurrentForestSize--;
      }
    }
    else if(pathType == 1) {
      float** forestdist;

      forestdist = new float*[b->getSubTreeSize() + 1];
      for(int i = 0; i < b->getSubTreeSize() + 1; i++) {
        forestdist[i] = new float[a->getSubTreeSize() + 1];
      }
      revTreeEditDist(b, a, forestdist, true, true);

      int a_rightmost_leaf_in_preL = F->preL_to_rid[a->getID()];
      int b_rightmost_leaf_in_preL = G->preL_to_rid[b->getID()];

      int a_rightmost_leaf_in_postR = F->preL_to_postR[a_rightmost_leaf_in_preL];
      int b_rightmost_leaf_in_postR = G->preL_to_postR[b_rightmost_leaf_in_preL];

      int aoff = a_rightmost_leaf_in_postR - 1;//consider gap
      int boff = b_rightmost_leaf_in_postR - 1;//consider gap

      int FcurrentForestSize = b->getSubTreeSize();
      int GcurrentForestSize = a->getSubTreeSize();

      if(DEBUG) {
        ou << "forestdist" << endl;
        for(int i = 0; i < b->getSubTreeSize() + 1; i++) {
          for(int j = 0; j < a->getSubTreeSize() + 1; j++) {
            ou << forestdist[i][j] << " ";
          }
          ou << endl;
        } 
      }


      int i = F->preL_to_postR[b->getID()] - boff;
      int j = G->preL_to_postR[a->getID()] - aoff;

      if(b->getID() != 0 && a->getID() != 0) {
        i--;
        j--;
        FcurrentForestSize--;
        GcurrentForestSize--;
      }

      while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
 
        int i_minus1_in_preL = F->postR_to_preL[i + boff];
        int j_minus1_in_preL = G->postR_to_preL[j + aoff];
        bool isFTree = F->preL_to_rid[i_minus1_in_preL] == F->preL_to_rid[b->getID()];
        bool isGTree = G->preL_to_rid[j_minus1_in_preL] == G->preL_to_rid[a->getID()];
        float da = forestdist[i - 1][j] + costModel_.ins((*F)[i_minus1_in_preL]->getLabel());
        float db = forestdist[i][j - 1] + costModel_.del((*G)[j_minus1_in_preL]->getLabel());
        int size_of_i_minus1_in_preL = (*F)[i_minus1_in_preL]->getSubTreeSize();
        int size_of_j_minus1_in_preL = (*G)[j_minus1_in_preL]->getSubTreeSize();
        float dc = forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] + delta[j_minus1_in_preL][i_minus1_in_preL] + costModel_.ren((*G)[j_minus1_in_preL]->getLabel(), (*F)[i_minus1_in_preL]->getLabel());


        if(DEBUG) {
          ou << "i = " << i_minus1_in_preL << " j = " << j_minus1_in_preL << endl;
        }

        if(DEBUG) {
          ou << "da = " << da << " db = " << db << " dc = " << dc << endl;
          ou << "forestdist[" << i - 1 << ", " << j << "] = " << forestdist[i - 1][j] << endl;
          ou << "forestdist[" << i << ", " << j - 1 << "] = " << forestdist[i][j - 1] << endl;
          ou << "forestdist[" << i - size_of_i_minus1_in_preL << ", " << j - size_of_j_minus1_in_preL << "] = " << forestdist[i - size_of_i_minus1_in_preL][j - size_of_j_minus1_in_preL] << endl;
          ou << "forestdist[" << i << ", " << j << "] = " << forestdist[i][j] << endl; 
        }
        if(da == forestdist[i][j]) {
          map->setMap(-1, i_minus1_in_preL);
          if(DEBUG) {
            ou << "- -> " << (*F)[i_minus1_in_preL]->getLabel() << endl;
          }
          i = i - 1;
          FcurrentForestSize--;
        } else if(db == forestdist[i][j]) {
          map->setMap(j_minus1_in_preL, -1);
          if(DEBUG) {
            ou << (*G)[j_minus1_in_preL]->getLabel() << " - ->" << endl;
          }
          j = j - 1;
          GcurrentForestSize--;
        } else if(dc == forestdist[i][j]) {
          map->setMap(j_minus1_in_preL, i_minus1_in_preL);
          if(DEBUG) {
            ou << (*G)[j_minus1_in_preL]->getLabel() << "(" << j_minus1_in_preL << ")" << " -> " << (*F)[i_minus1_in_preL]->getLabel() << "(" << i_minus1_in_preL << ")" << endl;
          }
          if(isFTree && isGTree) {
            i = i - 1;
            j = j - 1;
            FcurrentForestSize--;
            GcurrentForestSize--;
          }
          else {
            i = i - size_of_i_minus1_in_preL;
            j = j - size_of_j_minus1_in_preL;
            FcurrentForestSize -= size_of_i_minus1_in_preL;
            GcurrentForestSize -= size_of_j_minus1_in_preL;
            if((*F)[i_minus1_in_preL]->getSubTreeSize() > 1 && (*G)[j_minus1_in_preL]->getSubTreeSize() > 1) {
              gteo((*G)[j_minus1_in_preL], (*F)[i_minus1_in_preL]);
            } 
          }
        }
      }
      while(FcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = " << i << " j = 0" << endl;
        } 
        int i_minus1_in_preL = F->postR_to_preL[i - 1];
        map->setMap(-1, i_minus1_in_preL);
        i = i - 1;
        FcurrentForestSize--;
      }
      while(GcurrentForestSize > 0) {
        if(DEBUG) {
          ou << "i = 0 j = " << j << endl;
        }
        int j_minus1_in_preL = G->postR_to_preL[j - 1];
        map->setMap(j_minus1_in_preL, -1);
        j = j - 1;
        GcurrentForestSize--;
      }
    }

    else if(pathType == 2) {
      
      if(direction == 0) {

        float*** forestdist; 
        forestdist = new float**[treeSizeB + 1];
        for(int i = 0; i < treeSizeB + 1; i++) {
          forestdist[i] = new float*[treeSizeA + 1];
          for(int j = 0; j < treeSizeA + 1; j++) {
            forestdist[i][j] = new float[treeSizeA + 1];
          }
        }
        spfA_LR(b, a, pathLeaf, pathType, forestdist, true, true);

        int lF = b->getID();
        int rF = F->preL_to_preR[lF];
        int FcurrentForestSize = b->getSubTreeSize();

        int lG = a->getID();
        int rG = G->preL_to_preR[lG];
        int GcurrentForestSize = a->getSubTreeSize();


        int direction = 0;//0 for right 1 for left;
        int leaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
        Node* parent = (*F)[leaf]->getParent();
        int* favouriteChild = new int[F->getTreeSize()];
        int count = 0;
        while(leaf != b->getID()) {
          favouriteChild[count++] = leaf;
          leaf = parent->getID();
          parent = (*F)[leaf]->getParent();
        }
        favouriteChild[count++] = leaf;
        int i = count - 1;
        if(!(a->getID() == 0 && b->getID() == 0)) {
          lF++;
          rF++;
          lG++;
          rG++;
          FcurrentForestSize--;
          GcurrentForestSize--;
          i--;
        }

        int prev_fav_child = favouriteChild[count - 1];
        int prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];

        if(DEBUG) {
          ou << "favouriteChild: " << endl;
          for(int i = count - 1; i >= 0; i--) {
            ou << favouriteChild[i] << " ";
          }
          ou << endl;
        }


        while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
          if(DEBUG) {
            ou << "(" << lF << ", " << rF << ", " << lG << ", " << rG << ")" << endl;
            ou << "FcurrentForestSize = " << FcurrentForestSize << endl;
            ou << "GcurrentForestSize = " << GcurrentForestSize << endl;
          }
          int rF_in_preL = F->preR_to_preL[rF];
          int rG_in_preL = G->preR_to_preL[rG];
          int lF_in_preR = F->preL_to_preR[lF];

          float da = 0;
          float db = 0;
          float dc = 0;

          bool isFTree = rF_in_preL == lF;
          bool isGTree = rG_in_preL == lG;
          int favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];

          if(rF == favouriteChild_in_preR) {
            direction = 1;
          }
          if(lF == favouriteChild[i]) {
            prev_fav_child = favouriteChild[i];
            prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];
            i--;
            direction = 0;
            favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];
          }


          bool hasRightPart = (abs(favouriteChild_in_preR - prev_fav_child_in_preR) > 1);
          bool hasLeftPart = (abs(favouriteChild[i] - prev_fav_child) > 1); 
      
          if(isGTree && isFTree) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            //dc = forestdist[rF + 1][lG + 1][rG + 1] + delta[lG][rF_in_preL] +costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());
            dc = delta[lG][rF_in_preL] +costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());
            if(DEBUG) {
              ou << "isGTree && isFTree" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_plus_one_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_plus_one_in_preR][lG][rG] << endl;
              ou << "forestdist[" << rF << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[rF][lG + 1][rG + 1] << endl;
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou <<  " - ->" << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
              rF++;
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << "-> -" << endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(lG, rF_in_preL);
              lG++;
              rG++;
              lF++;
              rF++;
              FcurrentForestSize--;
              GcurrentForestSize--;
            }
          } 

          else if(isFTree && direction == 0) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.ins((*F)[rF_in_preL]->getLabel());
            //da = forestdist[rF + 1][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG][rG + 1] + costModel_.del((*G)[rG_in_preL]->getLabel());
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = (hasRightPart? forestdist[rF + 1][lG][rG + size_of_rG] : forestdist[lF_plus_one_in_preR][lG][rG + size_of_rG]) + delta[rG_in_preL][rF_in_preL] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());
            //dc = forestdist[rF + 1][lG][rG + size_of_rG] + delta[rG_in_preL][rF_in_preL] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }
            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << " - ->" << (*F)[rF_in_preL]->getLabel() <<  endl;
              }
              map->setMap(-1, rF_in_preL);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL] << "-> - " << endl; 
              }
              map->setMap(rG_in_preL, -1);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(rG_in_preL, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*G)[rG_in_preL], (*F)[rF_in_preL]);
              }
              lF++;
              rF++;
              rG = rG + size_of_rG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_rG;
            }
          } 

          else if(isFTree && direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = (hasRightPart? forestdist[rF + 1][lG][rG] : forestdist[lF_plus_one_in_preR][lG][rG]) + costModel_.ins((*F)[rF_in_preL]->getLabel());
            //da = forestdist[rF + 1][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = (hasRightPart? forestdist[rF + 1][lG + size_of_lG][rG] : forestdist[lF_plus_one_in_preR][lG + size_of_lG][rG]) + delta[lG][rF_in_preL] + costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());
            //dc = forestdist[rF + 1][lG + size_of_lG][rG] + delta[lG][rF_in_preL] + costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou <<  "- ->" << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << "-> -" << endl;
              } 
              map->setMap(lG, -1);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(lG, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[rF_in_preL]);
              }
              lF++;
              rF++;
              lG = lG + size_of_lG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_lG;
            }
          } 

          else if(isGTree && direction == 0) {
            da = forestdist[rF + 1][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            dc = forestdist[rF + size_of_rF][lG + 1][rG + 1] + delta[lG][rF_in_preL] + costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());

     
            if(DEBUG) {
              ou << "isGTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou <<  "- ->" << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
              rF++;
              FcurrentForestSize--;
            }
            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() <<  "-> -" << endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }
            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() <<  endl;
              }
              map->setMap(lG, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[rF_in_preL]);
              }
              rF = rF + size_of_rF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize--;
            }
          } 

          else if(isGTree && direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = forestdist[lF_plus_one_in_preR][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[rF][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int lF_plus_size_of_lF_in_preR = F->preL_to_preR[lF + size_of_lF];
            dc = forestdist[lF_plus_size_of_lF_in_preR][lG + 1][rG + 1] + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());

            if(DEBUG) {
              ou << "isGTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou << "- ->" << (*F)[lF]->getLabel() << endl; 
              }
              map->setMap(-1, lF);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << "-> -" << endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(lG, lF);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[lF]);
              }
              lF = lF + size_of_lF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize--;
            }
          }
          
          //F and G are forest and direction = right
          else if(direction == 0) {
            da = forestdist[rF + 1][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF][lG][rG + 1] + costModel_.del((*G)[rG_in_preL]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = forestdist[rF + size_of_rF][lG][rG + size_of_rG] + delta[rG_in_preL][rF_in_preL] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());

            if(DEBUG) {
              ou << "Both forest && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF << ", " << lG << ", " << rG << "] = " << forestdist[rF][lG][rG] << endl;
            }

            if(forestdist[rF][lG][rG] == da) {
              if(DEBUG) {
                ou << " - -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL] << "-> -" << endl;
              }
              map->setMap(rG_in_preL, -1);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " <<  (*F)[rF_in_preL]->getLabel() <<  endl;
              }
              map->setMap(rG_in_preL, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*G)[rG_in_preL], (*F)[rF_in_preL]);
              }
              rF = rF + size_of_rF;
              rG = rG + size_of_rG;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize -= size_of_rG;
            }
          }
          
          //F and G are forest and direction = left
          else if(direction == 1) {
            int lF_plus_one_in_preR = F->preL_to_preR[lF + 1];
            da = forestdist[lF_plus_one_in_preR][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF_in_preR][lG + 1][rG] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            int lF_plus_size_of_lF_in_preR = F->preL_to_preR[lF + size_of_lF];
            dc = forestdist[lF_plus_size_of_lF_in_preR][lG + size_of_lG][rG] + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());


            if(DEBUG) {
              ou << "Both forest && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF_in_preR << ", " << lG << ", " << rG << "] = " << forestdist[lF_in_preR][lG][rG] << endl;
            }

            if(forestdist[lF_in_preR][lG][rG] == da) {
              if(DEBUG) {
                ou << " - ->" << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(-1, lF);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() <<  "-> -" << endl;
              }
              map->setMap(lG, -1);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF_in_preR][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() <<  endl;
              }
              map->setMap(lG, lF);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[lF]);
              }
              lF = lF + size_of_lF;
              lG = lG + size_of_lG; 
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize -= size_of_lG;
            }
          }


        }
      } 

      else if(direction == 1) {

        float*** forestdist;
        forestdist = new float**[treeSizeB + 1];
        for(int i = 0; i < treeSizeB + 1; i++) {
          forestdist[i] = new float*[treeSizeA + 1];
          for(int j = 0; j < treeSizeA + 1; j++) {
            forestdist[i][j] = new float[treeSizeA + 1];
          }
        }
        spfA_RL(b, a, pathLeaf, pathType, forestdist, true, true);

        int lF = b->getID();
        int rF = F->preL_to_preR[lF];
        int FcurrentForestSize = b->getSubTreeSize();

        int lG = a->getID();
        int rG = G->preL_to_preR[lG];
        int GcurrentForestSize = a->getSubTreeSize();


        int direction = 1;//0 for right 1 for left;
        int leaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
        Node* parent = (*F)[leaf]->getParent();
        int* favouriteChild = new int[F->getTreeSize()];
        int count = 0;
        while(leaf != b->getID()) {
          favouriteChild[count++] = leaf;
          leaf = parent->getID();
          parent = (*F)[leaf]->getParent();
        }
        favouriteChild[count++] = leaf;
        int i = count - 1;
        if(!(a->getID() == 0 && b->getID() == 0)) {
          lF++;
          rF++;
          lG++;
          rG++;
          FcurrentForestSize--;
          GcurrentForestSize--;
          i--;
        }

        int prev_fav_child = favouriteChild[count - 1];
        int prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];

        if(DEBUG) {
          ou << "favouriteChild: " << endl;
          for(int i = count - 1; i >= 0; i--) {
            ou << favouriteChild[i] << " ";
          }
          ou << endl;
        }


        while(FcurrentForestSize > 0 && GcurrentForestSize > 0) {
          if(DEBUG) {
            ou << "(" << lF << ", " << rF << ", " << lG << ", " << rG << ")" << endl;
            ou << "FcurrentForestSize = " << FcurrentForestSize << endl;
            ou << "GcurrentForestSize = " << GcurrentForestSize << endl;
          }
          int rF_in_preL = F->preR_to_preL[rF];
          int rG_in_preL = G->preR_to_preL[rG];
          int lF_in_preR = F->preL_to_preR[lF];

          float da = 0;
          float db = 0;
          float dc = 0;

          int favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];

          if(lF == favouriteChild[i]) {
            direction = 0;
          }
          if(rF == favouriteChild_in_preR) {
            prev_fav_child = favouriteChild[i];
            prev_fav_child_in_preR = F->preL_to_preR[prev_fav_child];
            i--;
            direction = 1;
            favouriteChild_in_preR = F->preL_to_preR[favouriteChild[i]];
          }

          bool isFTree = rF_in_preL == lF;
          bool isGTree = rG_in_preL == lG;


          bool hasRightPart = (abs(favouriteChild_in_preR - prev_fav_child_in_preR) > 1);
          bool hasLeftPart = (abs(favouriteChild[i] - prev_fav_child) > 1); 
          if(DEBUG) {
            ou << "hasRightPart = " << hasRightPart << endl;
            ou << "hasLeftPart = " << hasLeftPart << endl;
            ou << "favouriteChild = " << favouriteChild[i] << endl;
            ou << "favouriteChild_in_preR = " << favouriteChild_in_preR << endl;
          }

          if(isGTree && isFTree) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            //dc = forestdist[lF + 1][lG + 1][rG + 1] + delta[lG][lF] +costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());
            dc = delta[lG][lF] +costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());
            if(DEBUG) {
              ou << "isGTree && isFTree" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF_plus_one_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_plus_one_in_preL][lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG + 1 << ", " << rG + 1 << "] = " << forestdist[lF][lG + 1][rG + 1] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << " -> -" << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(-1, lF);
              rF++;
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() <<  "-> -" << endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(lG, lF);
              lG++;
              rG++;
              lF++;
              rF++;
              FcurrentForestSize--;
              GcurrentForestSize--;
            }
          }

          else if(isFTree && direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.ins((*F)[lF]->getLabel());
            //da = forestdist[lF + 1][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF][lG][rG + 1] + costModel_.del((*G)[rG_in_preL]->getLabel());
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            dc = (hasLeftPart ? forestdist[lF + 1][lG][rG + size_of_rG] : forestdist[rF_plus_one_in_preL][lG][rG + size_of_rG]) + delta[rG_in_preL][lF] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[lF]->getLabel());
            //dc = forestdist[lF + 1][lG][rG + size_of_rG] + delta[rG_in_preL][lF] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[lF]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG << "] = " << forestdist[lF + 1][lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG + 1 << "] = " << forestdist[lF][lG][rG + 1] << endl;
              ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG + size_of_rG << "] = " << forestdist[lF + 1][lG][rG + size_of_rG] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }
            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << "- -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(-1, lF);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL] << "-> -" << endl; 
              }
              map->setMap(rG_in_preL, -1);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(rG_in_preL, lF);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*G)[rG_in_preL], (*F)[lF]);
              }
              lF++;
              rF++;
              rG = rG + size_of_rG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_rG;
            }
          }

          else if(isFTree && direction == 1) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = (hasLeftPart ? forestdist[lF + 1][lG][rG] : forestdist[rF_plus_one_in_preL][lG][rG]) + costModel_.ins((*F)[lF]->getLabel());
            //da = forestdist[lF + 1][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = (hasLeftPart ? forestdist[lF + 1][lG + size_of_lG][rG] : forestdist[rF_plus_one_in_preL][lG + size_of_lG][rG]) + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());
            //dc = forestdist[lF + 1][lG + size_of_lG][rG] + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());

            if(DEBUG) {
              ou << "isFTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG << "] = " << forestdist[lF + 1][lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG + 1 << ", " << rG << "] = " << forestdist[lF][lG + 1][rG] << endl;
              ou << "forestdist[" << lF + 1 << ", " << lG + size_of_lG << ", " << rG << "] = " << forestdist[lF + 1][lG + size_of_lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << "- -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(-1, lF);
              lF++;
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> - " << endl;
              } 
              map->setMap(lG, -1);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() <<  endl;
              }
              map->setMap(lF, lG);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[lF]);
              }
              lF++;
              rF++;
              lG = lG + size_of_lG;
              FcurrentForestSize--;
              GcurrentForestSize -= size_of_lG;
            }
          }

          else if(isGTree && direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = forestdist[rF_plus_one_in_preL][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF_in_preL][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int rF_plus_size_of_rF_in_preL = F->preR_to_preL[rF + size_of_rF];
            dc = forestdist[rF_plus_size_of_rF_in_preL][lG + 1][rG + 1] + delta[lG][rF_in_preL] + costModel_.ren((*G)[lG]->getLabel(), (*F)[rF_in_preL]->getLabel());

     
            if(DEBUG) {
              ou << "isGTree && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
            }

            if(forestdist[rF_in_preL][lG][rG] == da) {
              if(DEBUG) {
                ou << " - -> " << (*F)[rF_in_preL]->getLabel() <<  endl;
              }
              map->setMap(-1, rF_in_preL);
              rF++;
              FcurrentForestSize--;
            }
            else if(forestdist[rF_in_preL][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> -" <<  endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }
            else if(forestdist[rF_in_preL][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(lG, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[rF_in_preL]);
              }
              rF = rF + size_of_rF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize--;
            }
          }

          else if(isGTree && direction == 1) {
            da = forestdist[lF + 1][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG + 1] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            dc = forestdist[lF + size_of_lF][lG + 1][rG + 1] + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());

            if(DEBUG) {
              ou << "isGTree && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << " - -> " << (*F)[lF]->getLabel() <<  endl; 
              }
              map->setMap(-1, lF);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> -" <<  endl;
              }
              map->setMap(lG, -1);
              lG++;
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() <<  endl;
              }
              map->setMap(lG, lF);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[lF]);
              }
              lF = lF + size_of_lF;
              lG++;
              rG++;
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize--;
            }
          }

          else if(direction == 0) {
            int rF_plus_one_in_preL = F->preR_to_preL[rF + 1];
            da = forestdist[rF_plus_one_in_preL][lG][rG] + costModel_.ins((*F)[rF_in_preL]->getLabel());
            db = forestdist[rF_in_preL][lG][rG + 1] + costModel_.del((*G)[rG_in_preL]->getLabel());
            int size_of_rF = (*F)[rF_in_preL]->getSubTreeSize();
            int size_of_rG = (*G)[rG_in_preL]->getSubTreeSize();
            int rF_plus_size_of_rF_in_preL = F->preR_to_preL[rF + size_of_rF];
            dc = forestdist[rF_plus_size_of_rF_in_preL][lG][rG + size_of_rG] + delta[rG_in_preL][rF_in_preL] + costModel_.ren((*G)[rG_in_preL]->getLabel(), (*F)[rF_in_preL]->getLabel());

            if(DEBUG) {
              ou << "Both forest && direction == 0" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << rF_in_preL << ", " << lG << ", " << rG << "] = " << forestdist[rF_in_preL][lG][rG] << endl;
            }

            if(forestdist[rF_in_preL][lG][rG] == da) {
              if(DEBUG) {
                ou << " - -> " <<  (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(-1, rF_in_preL);
              rF++;
              FcurrentForestSize--;
            }

            else if(forestdist[rF_in_preL][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL] << " -> -" <<  endl;
              }
              map->setMap(rG_in_preL, -1);
              rG++;
              GcurrentForestSize--;
            }

            else if(forestdist[rF_in_preL][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[rG_in_preL]->getLabel() << " -> " << (*F)[rF_in_preL]->getLabel() << endl;
              }
              map->setMap(rG_in_preL, rF_in_preL);
              if((*F)[rF_in_preL]->getSubTreeSize() > 1 && (*G)[rG_in_preL]->getSubTreeSize() > 1) {
                gteo((*G)[rG_in_preL], (*F)[rF_in_preL]);
              }
              rF = rF + size_of_rF;
              rG = rG + size_of_rG;
              FcurrentForestSize -= size_of_rF;
              GcurrentForestSize -= size_of_rG;
            }
          }


          //F and G are forest and direction = left
          else if(direction == 1) {
            da = forestdist[lF + 1][lG][rG] + costModel_.ins((*F)[lF]->getLabel());
            db = forestdist[lF][lG + 1][rG] + costModel_.del((*G)[lG]->getLabel());
            int size_of_lF = (*F)[lF]->getSubTreeSize();
            int size_of_lG = (*G)[lG]->getSubTreeSize();
            dc = forestdist[lF + size_of_lF][lG + size_of_lG][rG] + delta[lG][lF] + costModel_.ren((*G)[lG]->getLabel(), (*F)[lF]->getLabel());


            if(DEBUG) {
              ou << "Both forest && direction == 1" << endl;
              ou << "da = " << da << endl;
              ou << "db = " << db << endl;
              ou << "dc = " << dc << endl; 
              ou << "forestdist[" << lF + 1 << ", " << lG << ", " << rG << "] = " << forestdist[lF + 1][lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG + 1 << ", " << rG << "] = " << forestdist[lF][lG + 1][rG] << endl;
              ou << "forestdist" << lF + size_of_lF << ", " << lG + size_of_lG << ", " << rG << "] = " << forestdist[lF + size_of_lF][lG + size_of_lG][rG] << endl;
              ou << "forestdist[" << lF << ", " << lG << ", " << rG << "] = " << forestdist[lF][lG][rG] << endl;
            }

            if(forestdist[lF][lG][rG] == da) {
              if(DEBUG) {
                ou << " - -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(-1, lF);
              lF++;
              FcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == db) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> - " << endl;
              }
              map->setMap(lG, -1);
              lG++;
              GcurrentForestSize--;
            }

            else if(forestdist[lF][lG][rG] == dc) {
              if(DEBUG) {
                ou << (*G)[lG]->getLabel() << " -> " << (*F)[lF]->getLabel() << endl;
              }
              map->setMap(lG, lF);
              if((*F)[lF]->getSubTreeSize() > 1 && (*G)[lG]->getSubTreeSize() > 1) {
                gteo((*G)[lG], (*F)[lF]);
              }
              lF = lF + size_of_lF;
              lG = lG + size_of_lG; 
              FcurrentForestSize -= size_of_lF;
              GcurrentForestSize -= size_of_lG;
            }
          }     
        }




      }
    }
  }
};
*/


char** TreeComparison::getResult(void) {
  int maxSize = rB_.getRNASize() + rA_.getRNASize();


  char** result;

  result = new char*[4];
  for(int i = 0; i < 4; i++) {
    result[i] = new char[maxSize];
  }


/*  int i;
  for(int rnaAIndex = 2; rnaAIndex < rA_.RNASize_; rnaAIndex++) {
    if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
      if((*map)[treeAIndex] == -1)
      result[0][i - 2] = '-';
      result[2][i - 2] = rA_.originalSequence[i];
    } else if(rA_.secondaryStructure[i] > i) {
      result[0][i - 2] = '(';
      result[2][i - 2] = rA_.originalSequence[i];
      result[0][rA_.secondaryStructure[i] - 2] = ')';
      result[2][rA_.secondaryStructure[i] - 2] = rA_.originalSequence[rA_.secondaryStructure[i]];
    }
  }
  result[0][i - 2] = '\0';
  result[2][i - 2] = '\0';

  for(i = 2; i < rB_.RNASize_; i++) {
    if(rB_.secondaryStructure[i] == i) {
      result[1][i - 2] = '-';
      result[3][i - 2] = rB_.originalSequence[i];
    } else if(rB_.secondaryStructure[i] > i) {
      result[1][i - 2] = '(';
      result[3][i - 2] = rB_.originalSequence[i];
      result[1][rB_.secondaryStructure[i] - 2] = ')';
      result[3][rB_.secondaryStructure[i] - 2] = rB_.originalSequence[rB_.secondaryStructure[i]];
    }
  }
  result[1][i - 2] = '\0';
  result[3][i - 2] = '\0';
*/

  int rnaAIndex = 2;
  int rnaBIndex = 2;
  int resultIndex = 0;
  
  while(rnaAIndex < rA_.RNASize_ && rnaBIndex < rB_.RNASize_) {
    int treeAIndex = rA_.original_to_tree[rnaAIndex];
    int treeBIndex = rB_.original_to_tree[rnaBIndex];

    if((*map)[treeAIndex] == treeBIndex) {
      result[2][resultIndex] = rA_.originalSequence[rnaAIndex];
      result[3][resultIndex] = rB_.originalSequence[rnaBIndex];
      
      if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
        result[0][resultIndex] = '-';
      } else if(rA_.secondaryStructure[rnaAIndex] > rnaAIndex) {
        result[0][resultIndex] = '(';
      } else if(rA_.secondaryStructure[rnaAIndex] < rnaAIndex) {
        result[0][resultIndex] = ')';
      }

      if(rB_.secondaryStructure[rnaBIndex] == rnaBIndex) {
        result[1][resultIndex] = '-';
      } else if(rB_.secondaryStructure[rnaBIndex] > rnaBIndex) {
        result[1][resultIndex] = '(';
      } else if(rB_.secondaryStructure[rnaBIndex] < rnaBIndex) {
        result[1][resultIndex] = ')';
      }

      rnaAIndex++;
      rnaBIndex++;
    }

    else if((*map)[treeAIndex] == -1) {
      result[2][resultIndex] = rA_.originalSequence[rnaAIndex];
      result[3][resultIndex] = '-';

      if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
        result[0][resultIndex] = '-';
      } else if(rA_.secondaryStructure[rnaAIndex] > rnaAIndex) {
        result[0][resultIndex] = '(';
      } else if(rA_.secondaryStructure[rnaAIndex] < rnaAIndex) {
        result[0][resultIndex] = ')';
      }

      result[1][resultIndex] = '-';

      /*if(rB_.secondaryStructure[rnaBIndex] == rnaBIndex) {
        result[1][resultIndex] = '-';
      } else if(rB_.secondaryStructure[rnaBIndex] > rnaBIndex) {
        result[1][resultIndex] = '(';
      } else if(rB_.secondaryStructure[rnaBIndex] < rnaBIndex) {
        result[1][resultIndex] = ')';
      }*/

      rnaAIndex++;
    }

    else if((*map)(treeBIndex) == -1) {
      result[2][resultIndex] = '-';
      result[3][resultIndex] = rB_.originalSequence[rnaBIndex];

      result[0][resultIndex] = '-';

      /*if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
        result[0][resultIndex] = '-';
      } else if(rA_.secondaryStructure[rnaAIndex] > rnaAIndex) {
        result[0][resultIndex] = '(';
      } else if(rA_.secondaryStructure[rnaAIndex] < rnaAIndex) {
        result[0][resultIndex] = ')';
      }*/

      if(rB_.secondaryStructure[rnaBIndex] == rnaBIndex) {
        result[1][resultIndex] = '-';
      } else if(rB_.secondaryStructure[rnaBIndex] > rnaBIndex) {
        result[1][resultIndex] = '(';
      } else if(rB_.secondaryStructure[rnaBIndex] < rnaBIndex) {
        result[1][resultIndex] = ')';
      }

      rnaBIndex++;
    }
    resultIndex++;
  }
  while(rnaAIndex < rA_.RNASize_) {
     result[2][resultIndex] = rA_.originalSequence[rnaAIndex];
      result[3][resultIndex] = '-';

      if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
        result[0][resultIndex] = '-';
      } else if(rA_.secondaryStructure[rnaAIndex] > rnaAIndex) {
        result[0][resultIndex] = '(';
      } else if(rA_.secondaryStructure[rnaAIndex] < rnaAIndex) {
        result[0][resultIndex] = ')';
      }

      result[1][resultIndex] = '-';

      /*if(rB_.secondaryStructure[rnaBIndex] == rnaBIndex) {
        result[1][resultIndex] = '-';
      } else if(rB_.secondaryStructure[rnaBIndex] > rnaBIndex) {
        result[1][resultIndex] = '(';
      } else if(rB_.secondaryStructure[rnaBIndex] < rnaBIndex) {
        result[1][resultIndex] = ')';
      }*/

      rnaAIndex++;
      resultIndex++;
  }
  while(rnaBIndex < rB_.RNASize_) {
    result[2][resultIndex] = '-';
      result[3][resultIndex] = rB_.originalSequence[rnaBIndex];

      result[0][resultIndex] = '-';

      /*if(rA_.secondaryStructure[rnaAIndex] == rnaAIndex) {
        result[0][resultIndex] = '-';
      } else if(rA_.secondaryStructure[rnaAIndex] > rnaAIndex) {
        result[0][resultIndex] = '(';
      } else if(rA_.secondaryStructure[rnaAIndex] < rnaAIndex) {
        result[0][resultIndex] = ')';
      }*/

      if(rB_.secondaryStructure[rnaBIndex] == rnaBIndex) {
        result[1][resultIndex] = '-';
      } else if(rB_.secondaryStructure[rnaBIndex] > rnaBIndex) {
        result[1][resultIndex] = '(';
      } else if(rB_.secondaryStructure[rnaBIndex] < rnaBIndex) {
        result[1][resultIndex] = ')';
      }

      rnaBIndex++;
      resultIndex++;
  }

  result[0][resultIndex] = '\0';
  result[1][resultIndex] = '\0';
  result[2][resultIndex] = '\0';
  result[3][resultIndex] = '\0';

/*  int rAIndex = 2;
  int rBIndex = 2;
  int tAIndex = 1;
  int tBIndex = 1;

  string s1 = "";
  string s2 = "";

  while(tAIndex < treeSizeA && tBIndex < treeSizeB) {
    //cout << tAIndex << ", " << tBIndex << endl;
    if((*map)[tAIndex] == tBIndex) {
      s1 += (*A_)[tAIndex]->getLabel();
      s2 += (*B_)[tBIndex]->getLabel();
      tAIndex++;
      tBIndex++;
    }

     else if((*map)[tAIndex] == -1) {
      s1 += (*A_)[tAIndex]->getLabel();
      s2 += '-';
      tAIndex++;
    }

    else if((*map)(tBIndex) == -1) {
      s1 += '-';
      s2 += (*B_)[tBIndex]->getLabel();
      tBIndex++;
    }

  }

  while(tAIndex < treeSizeA) {
    s1 += (*A_)[tAIndex]->getLabel();
    s2 += '-';
    tAIndex++;
  }

  while(tBIndex < treeSizeB) {
    s1 += '-';
    s2 += (*B_)[tBIndex]->getLabel();
    tBIndex++;
  }

  cout << s1 << endl;
  cout << s2 << endl;*/

  return result;

};
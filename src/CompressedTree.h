#ifndef COMPRESSEDTREE_H
#define COMPRESSEDTREE_H

#include "Tree.h"
#include "Node.h"

class CompressedTree {
public:
	CompressedTree(Tree* );
	int getTreeSize(void) const;
	Node* operator[](int);


	string toString(void) const;

	friend class TreeComparison;

private:
	Tree* originalTree_;
	int* original_to_compressed;
	int originalTreeSize_;


	int compressedTreeSize_;
	vector<vector<int> > compressed_to_original; 
	vector<vector<char> > compressed_to_original_label;

	vector<Node*> preL;

	int* preL_to_preR;
	int* preR_to_preL;
	int* preL_to_postL;
	int* postL_to_preL;
	int* preL_to_postR;
	int* postR_to_preL;

	int* preL_to_lid;//left-to-right preorder to the leftmost tree leaf's id in left-to-right preorder
	int* preL_to_rid;//left-to-right preorder to the rightmost tree leaf's id in left-to-right preorder
	int* preL_to_ln;//first leaf node to the left of n in left-to-right preorder.
	int* preR_to_ln;//first leaf node to the right of n in right-to-left preorder.


	int* preL_to_sumDelCost;
	int* preL_to_sumInsCost;

	int* preL_to_DelCost;
	int* preL_to_InsCost;

};

#endif
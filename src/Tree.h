#ifndef Tree_H
#define Tree_H

#include "Node.h"
#include <string>
#include <vector>

using namespace std;
//using namespace FPM;

class Tree {

private:
	string label_;
	int treeSize_;
	vector<Node*> preL;
	//vector<Node*> post;
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

	int* keyRoot_L;
	int* keyRoot_R;

	int keyRoot_L_size;
	int keyRoot_R_size;

public:
	Tree(string, int);
	string getLabel(void) const;
	int getTreeSize(void) const;
	vector<Node*> getPreL(void) const;
	Node* operator[](int);
	//vector<Node*> getPost(void) const;
	void pushNodeToPreL(Node*);
	//void pushNodeToPost(Node*);

	friend class RNA;
	friend class TreeComparison;
	friend class CompressedTree;

	string toString()const;

};


#endif
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
public:
	Tree(string, int);
	string getLabel(void) const;
	int getTreeSize(void) const;
	vector<Node*> getPreL(void) const;
	//vector<Node*> getPost(void) const;
	void pushNodeToPreL(Node*);
	//void pushNodeToPost(Node*);

	friend class RNA;

	string toString()const;

};


#endif
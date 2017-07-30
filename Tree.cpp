#include "Tree.h"
#include "Node.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;
//using namespace FPM;


Tree::Tree(string label, int treeSize) {
	label_ = label;
	treeSize_ = treeSize;
	preL_to_preR = new int[treeSize_];
	preR_to_preL = new int[treeSize_];
	preL_to_postL = new int[treeSize_];
	postL_to_preL = new int[treeSize_];
	preL_to_postR = new int[treeSize_];
	postR_to_preL = new int[treeSize_];
	preL_to_lid = new int[treeSize_];
	preL_to_rid = new int[treeSize_];
};

string Tree::getLabel(void) const{
	return label_;
};

int Tree::getTreeSize(void) const {
	return treeSize_;
};


vector<Node*> Tree::getPreL(void) const {
	return preL;
};

Node* Tree::operator[](int i) {
	if(i < 0 || i >= treeSize_) {
		cout << "Overflow" << endl;
		return preL[0];
	}
	return preL[i];
};

/*vector<Node*> Tree::getPost(void) const {
	return post;
};*/
	
void Tree::pushNodeToPreL(Node* node) {
	preL.push_back(node);
};
	
/*void Tree::pushNodeToPost(Node* node) {
	post.push_back(node);
};*/

string Tree::toString() const {
	string res = "";
	res += "preL size = " + to_string(preL.size()) + "\n";
	for(int i = 0; i < preL.size(); i++) {
		res = res + preL[i]->toString() + "\n";
	}
	res += "\n";
	res += "preL_to_preR\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_preR[i]) + " "; 
	}
	res += "\n";
	res += "preL_to_postL\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_postL[i]) + " ";
	}
	res += "\n";
	res += "preL_to_postR\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_postR[i]) + " ";
	}
	res += "\n";
	res += "preL_to_lid\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_lid[i]) + " ";
	}
	res += "\n";
	res += "preL_to_rid\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_rid[i]) + " ";
	}
	return res;
};
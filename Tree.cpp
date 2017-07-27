#include "Tree.h"
#include "Node.h"
#include <vector>
#include <string>

using namespace std;
//using namespace FPM;


Tree::Tree(string label, int treeSize) {
	Tree();
	label_ = label;
	treeSize_ = treeSize;
	preL_to_preR = new int[treeSize_];
	preR_to_preL = new int[treeSize_];
	preL_to_postL = new int[treeSize_];
	postL_to_preL = new int[treeSize_];
	preL_to_postR = new int[treeSize_];
	postR_to_preL = new int[treeSize_];
};

string Tree::getLabel(void) const{
	return label_;
};

int Tree::getTreeSize(void) const {
	return treeSize_;
};


vector<Node*> Tree::getPreL(void) const {
	return pre;
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
	res += "preL\n";
	for(int i = 0; i < pre.size(); i++) {
		res += pre[i] + "\n"
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
	return res;
};
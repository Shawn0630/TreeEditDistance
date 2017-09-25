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
	keyRoot_L_size = 0;
	keyRoot_R_size = 0;

	preL_to_preR = new int[treeSize_];
	preR_to_preL = new int[treeSize_];
	preL_to_postL = new int[treeSize_];
	postL_to_preL = new int[treeSize_];
	preL_to_postR = new int[treeSize_];
	postR_to_preL = new int[treeSize_];
	preL_to_lid = new int[treeSize_];
	preL_to_rid = new int[treeSize_];
	preL_to_ln = new int[treeSize_];
	preR_to_ln = new int[treeSize_];
	
	preL_to_sumDelCost = new int[treeSize_];
	preL_to_sumInsCost = new int[treeSize_];

	keyRoot_L = new int[treeSize_];
	keyRoot_R = new int[treeSize_];


	fill_n(preL_to_preR, treeSize_, 0); 
	fill_n(preR_to_preL, treeSize_, 0); 
	fill_n(preL_to_postL, treeSize_, 0); 
	fill_n(postL_to_preL, treeSize_, 0); 
	fill_n(preL_to_postR, treeSize_, 0); 
	fill_n(postR_to_preL, treeSize_, 0); 
	fill_n(preL_to_lid, treeSize_, 0); 
	fill_n(preL_to_rid, treeSize_, 0); 
	fill_n(preL_to_sumDelCost, treeSize_, 0); 
	fill_n(preL_to_sumInsCost, treeSize_, 0);
	fill_n(keyRoot_L, treeSize_, 0);
	fill_n(keyRoot_R, treeSize_, 0); 

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
		cout << "i = " << i << " treeSize_ = " << treeSize_ << endl;
		cout << "Tree Overflow" << endl;
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
	res += "\n";
	res += "preL_to_sumDelCost\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_sumDelCost[i]) + " ";
	}
	res += "\n";
	res += "preL_to_sumInsCost\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_sumInsCost[i]) + " ";
	}
	res += "\n";
	res += "preL_to_ln\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preL_to_ln[i]) + " ";
	}
	res += "\n";
	res += "preR_to_ln\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(preR_to_ln[i]) + " ";
	}
	res += "\n";
	res += "postL_to_preL\n";
	for(int i = 0; i < treeSize_; i++) {
		res += to_string(postL_to_preL[i]) + " ";
	}
	return res;
};
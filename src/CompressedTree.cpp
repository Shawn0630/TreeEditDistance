#include "CompressedTree.h"
#include "Node.h"

#include <string>
#include <iostream>
using namespace std;


CompressedTree::CompressedTree(Tree *t) {

	originalTree_ = t;
	originalTreeSize_ = t->getTreeSize();
	original_to_compressed = new int[originalTreeSize_];
	int prev_single_path_node_in_original = -1;
	int compressed_index = 0;
	for(int i = 0; i < originalTreeSize_; i++) {
		Node* n = (*t)[i];
		vector<Node*> children = n->getChildren();
		if(children.size() == 1) {
			if(prev_single_path_node_in_original == -1) {
				prev_single_path_node_in_original = i;
				vector<int> v;
				vector<char> v_label;
				v.push_back(i);
				v_label.push_back(n->getLabel());
				compressed_to_original.push_back(v);
				compressed_to_original_label.push_back(v_label);
				original_to_compressed[i] = compressed_index;
			} else {
				compressed_to_original[compressed_index].push_back(i);
				compressed_to_original_label[compressed_index].push_back(n->getLabel());
				original_to_compressed[i] = compressed_index;
			}
		} else {
			if(prev_single_path_node_in_original != -1) {
				prev_single_path_node_in_original = -1;
				compressed_to_original[compressed_index].push_back(i);
				compressed_to_original_label[compressed_index].push_back(n->getLabel());
				original_to_compressed[i] = compressed_index;
				compressed_index++;
			}
			else {
				vector<int> v;
				vector<char> v_label;
				v.push_back(i);
				v_label.push_back(n->getLabel());
				compressed_to_original.push_back(v);
				compressed_to_original_label.push_back(v_label);
				original_to_compressed[i] = compressed_index;
				compressed_index++;
			}
		}
	}
	compressedTreeSize_ = compressed_index;

	preL_to_preR = new int[compressedTreeSize_];
 	preR_to_preL = new int[compressedTreeSize_];
	preL_to_postL = new int[compressedTreeSize_];
	postL_to_preL = new int[compressedTreeSize_];
	preL_to_postR = new int[compressedTreeSize_];
	postR_to_preL = new int[compressedTreeSize_];

	preL_to_lid = new int[compressedTreeSize_];
	preL_to_rid = new int[compressedTreeSize_];
	preL_to_ln = new int[compressedTreeSize_];
	preR_to_ln = new int[compressedTreeSize_];

	preL_to_sumDelCost = new int[compressedTreeSize_];
	preL_to_sumInsCost = new int[compressedTreeSize_];

	preL_to_DelCost = new int[compressedTreeSize_];
	preL_to_InsCost = new int[compressedTreeSize_];

	fill_n(preL_to_preR, compressedTreeSize_, 0); 
	fill_n(preR_to_preL, compressedTreeSize_, 0); 
	fill_n(preL_to_postL, compressedTreeSize_, 0); 
	fill_n(postL_to_preL, compressedTreeSize_, 0); 
	fill_n(preL_to_postR, compressedTreeSize_, 0); 
	fill_n(postR_to_preL, compressedTreeSize_, 0); 
	fill_n(preL_to_lid, compressedTreeSize_, 0); 
	fill_n(preL_to_rid, compressedTreeSize_, 0); 
	fill_n(preL_to_sumDelCost, compressedTreeSize_, 0); 
	fill_n(preL_to_sumInsCost, compressedTreeSize_, 0);
	fill_n(preL_to_DelCost, compressedTreeSize_, 0);
	fill_n(preL_to_InsCost, compressedTreeSize_, 0);



	int compressed_postL = 0;

	for(int i = 0; i < compressedTreeSize_; i++) {
		int original_preL = compressed_to_original[i][0];
		Node* n = new Node(i, (*t)[original_preL]->getLabel());
		preL.push_back(n);
	}

	for(int i = 0; i < originalTreeSize_; i++) {

		int original_preL = t->postL_to_preL[i];
		int compressed_preL = original_to_compressed[original_preL];

		if(compressed_to_original[compressed_preL][0] == original_preL) {
			vector<Node*> children_in_compressed = preL[compressed_preL]->getChildren();
			int size = 0;
			int sum = 0;
			int leftmostForestNum = 0;
			int rightmostForestNum = 0;
			int specialForestNum = 0;
			int leftmostSize = 0;
			int rightmostSize = 0;

			for(int i = 0; i < children_in_compressed.size(); i++) {
				sum += children_in_compressed[i]->getSubTreeSizeSum();
				size += children_in_compressed[i]->getSubTreeSize();
				leftmostForestNum += children_in_compressed[i]->getLeftmostForestNum();
				rightmostForestNum += children_in_compressed[i]->getRightmostForestNum();
			}
			if(!children_in_compressed.empty()) {
				leftmostSize = children_in_compressed[0]->getSubTreeSize();
				rightmostSize = children_in_compressed[children_in_compressed.size() - 1]->getSubTreeSize();
			}
			size += 1;
			sum += size;
			leftmostForestNum += size - leftmostSize;
			rightmostForestNum += size - rightmostSize;
			specialForestNum = size * (size + 3) / 2 - sum;
			preL[compressed_preL]->setSubTreeSize(size);
			preL[compressed_preL]->setLeftmostForestNum(leftmostForestNum);
			preL[compressed_preL]->setRightmostForestNum(rightmostForestNum);
			preL[compressed_preL]->setSpecialForestNum(specialForestNum);
			preL[compressed_preL]->setSubTreeSizeSum(sum);
			if(size == 1) {
				preL_to_lid[compressed_preL] = compressed_preL;
				preL_to_rid[compressed_preL] = compressed_preL;
			} else {
				preL_to_lid[compressed_preL] = preL_to_lid[preL[compressed_preL]->getLeftmostChild()->getID()];
				preL_to_rid[compressed_preL] = preL_to_rid[preL[compressed_preL]->getRightmostChild()->getID()];
			}

			Node* parent_in_original = (*t)[original_preL]->getParent();
			if(parent_in_original != NULL) {
				int parent_in_original_preL = parent_in_original->getID();
				int parent_in_compressed_preL = original_to_compressed[parent_in_original_preL];
				preL[compressed_preL]->setParent(preL[parent_in_compressed_preL]);
				preL[parent_in_compressed_preL]->pushChild(preL[compressed_preL]);
			}
			preL_to_postL[compressed_preL] = compressed_postL;
			postL_to_preL[compressed_postL] = compressed_preL;
			int compressed_postR = compressedTreeSize_ - 1- compressed_preL;
			preL_to_postR[compressed_preL] = compressed_postR;
			postR_to_preL[compressed_postR] = compressed_preL;
			int compressed_preR = compressedTreeSize_ - 1 - compressed_postL;
			preL_to_preR[compressed_preL] = compressed_preR;
			preR_to_preL[compressed_preR] = compressed_preL;
			compressed_postL++;
		}
	}

	int currentLeaf = -1;
	for(int i = 0; i < compressedTreeSize_; i++) {
		preL_to_ln[i] = currentLeaf;
		if(preL[i]->getChildren().size() == 0) currentLeaf = i;

	}

	currentLeaf = -1;
	for(int i = 0; i < compressedTreeSize_; i++) {
	  preR_to_ln[i] = currentLeaf;
      if(preL[preR_to_preL[i]]->getChildren().size() == 0) {
        currentLeaf = i;
      }
	}

};

int CompressedTree::getTreeSize(void) const {
	return compressedTreeSize_;
};

Node* CompressedTree::operator[](int i) {
	if(i < 0 || i >= compressedTreeSize_) {
		cout << "i = " << i << " compressedTreeSize_ = " << compressedTreeSize_ << endl;
		cout << "CompressedTree Overflow" << endl;
		return preL[0];
	}
	return preL[i];
};

string CompressedTree::toString(void) const {
	string res = "";
	res += "compressedTreeSize = " + to_string(compressedTreeSize_);
	res += '\n';
	for(int i = 0; i < preL.size(); i++) {
		res = res + preL[i]->toString() + "\n";
	}
	res += "original_to_compressed\n";
	for(int i = 0; i < originalTreeSize_; i++) {
		res += to_string(original_to_compressed[i]) + " ";
	}
	res += '\n';
	res += "compressed_to_original\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		vector<int> v = compressed_to_original[i];
		vector<char> v_label = compressed_to_original_label[i];
		res += to_string(v.size()) + " node(s) compressed to one node, they are "; 
		for(int j = 0; j < v.size(); j++) {
			res += to_string(v[j]) + "(" + v_label[j] + ") ";
		}
		res += '\n';
	}
	res += "preL_to_preR\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_preR[i]) + " ";
	}
	res += '\n';
	res += "preL_to_postL\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_postL[i]) + " ";
	}
	res += '\n';
	res += "preL_to_postR\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_postR[i]) + " ";
	}
	res += '\n';
	res += "preL_to_lid\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_lid[i]) + " ";
	}
	res += "\n";
	res += "preL_to_rid\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_rid[i]) + " ";
	}
	res += '\n';
	res += "preL_to_ln\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_ln[i]) + " ";
	}
	res += "\n";
	res += "preR_to_ln\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preR_to_ln[i]) + " ";
	}
	res += "\n";
	res += "preL_to_sumDelCost\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_sumDelCost[i]) + " ";
	}
	res += "\n";
	res += "preL_to_sumInsCost\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_sumInsCost[i]) + " ";
	}
	res += "preL_to_DelCost\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_DelCost[i]) + " ";
	}
	res += "\n";
	res += "preL_to_InsCost\n";
	for(int i = 0; i < compressedTreeSize_; i++) {
		res += to_string(preL_to_InsCost[i]) + " ";
	}
	return res;

};
#include "RNA.h"
#include "Node.h"
#include "Constant.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stack>
#include <string>
using namespace std;

RNA::RNA() {
	for(int i = 0; i < maxSize; i++) {
		secondaryStructure[i] = i;
	}
	RNASize_ = 0;
}

RNA::RNA(string RNAName) {
	RNA();
	RNAName_ = RNAName;
};

void RNA::setRNAName(string RNAName) {
	RNAName_ = RNAName;
};
string RNA::getRNAName(void) const {
	return RNAName_;
};

void RNA::setRNASize(int RNASize) {
	RNASize_ = RNASize;
};
int RNA::getRNASize(void) const {
	return RNASize_;
};

int RNA::getPreOrderSequenceSize(void)const {
	return preLSequenceSize_;
};

char& RNA::operator[] (int i) {
	if(i >= maxSize || i < 0) {
		cout << "RNA 1 Overflow" << endl;
		return originalSequence[0];
	}
	return originalSequence[i];

};

int& RNA::operator() (int i) {
	if(i >= maxSize || i < 0) {
		cout << "RNA 2 Overflow" << endl;
		return secondaryStructure[0];
	}
	return secondaryStructure[i];
};

char* RNA::getOriginalSequence(void) {
	char* sequence = new char[RNASize_ + 1];
	for(int i = 1; i <= RNASize_; i++) {
		sequence[i] = originalSequence[i];
	}
	return sequence;
};

int* RNA::getSecondaryStructure(void) {
	int* structure = new int[RNASize_ + 1];
	for(int i = 1; i <= RNASize_; i++) {
		structure[i] = secondaryStructure[i];
	}
	return structure;
};

void RNA::setPreOrderSequence(string preLSequence_) {
	for(int i = 0; i < preLSequence_.length(); i++) {
		preLSequence[i] = preLSequence_[i];
	}
};

int RNA::getTreeSize() const {
	return treeSize_;
};
void RNA::setTreeSize(int treeSize) {
	treeSize_ = treeSize;
};

char* RNA::getPreLSequence(void) {
	int j = 0;
	int treeSize = 0;
	for(int i = 1; i <= RNASize_; i++) {
		if(secondaryStructure[i] == i){
			preLSequence[j++] = '(';
			preLSequence[j++] = originalSequence[i];
			preLSequence[j++] = ')';
			tree_to_original[treeSize] = i;
			original_to_tree[i] = treeSize;
			treeSize++;
		}
		// if the base have pair, process the smaller one 
		else if(secondaryStructure[i] > i){
			preLSequence[j++] = '(';
			// there are altogether 16 cases.
			if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'A'){
				preLSequence[j++] = 'B';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'G'){
				preLSequence[j++] = 'D';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'C'){
				preLSequence[j++] = 'E';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'U'){
				preLSequence[j++] = 'F';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'A'){
				preLSequence[j++] = 'H';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'G'){
				preLSequence[j++] = 'I';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'C'){
				preLSequence[j++] = 'J';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'U'){
				preLSequence[j++] = 'K';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'A'){
				preLSequence[j++] = 'L';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'G'){
				preLSequence[j++] = 'M';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'C'){
				preLSequence[j++] = 'N';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'U'){
				preLSequence[j++] = 'O';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'A'){
				preLSequence[j++] = 'P';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'G'){
				preLSequence[j++] = 'Q';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'C'){
				preLSequence[j++] = 'R';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'U'){
				preLSequence[j++] = 'S';
			}
			tree_to_original[treeSize] = i;
			original_to_tree[i] = treeSize;
			original_to_tree[secondaryStructure[i]] = treeSize;
			treeSize++;
		}
		else {
			preLSequence[j++] = ')';
		}
	}
	preLSequenceSize_ = j;
	treeSize_ = treeSize;
	return preLSequence;
};
/*
	Before building the tree, preOrderSequence and the corresponding treeSize_ should be known.
*/
Tree* RNA::buildTree(void) {
	stack<Node*> nodes;
	Tree* tree = new Tree(RNAName_, treeSize_);
	int preN = 0;
	int cursor = 1;
	int postN = 0;
	int left = 1;
	int currentLeaf = -1;
	Node* node = new Node(preN, preLSequence[cursor]);
	nodes.push(node);
	tree->preL_to_postR[node->getID()] = tree->treeSize_ - 1 - node->getID();
	tree->postR_to_preL[tree->treeSize_ - 1 - node->getID()] = node->getID();
	tree->preL_to_ln[node->getID()] = currentLeaf;
	tree->pushNodeToPreL(node);
	while(left > 0) {
		cursor++;
		if(preLSequence[cursor] == '(') {
			left++;
			Node* node = new Node(++preN, preLSequence[++cursor]);
			tree->preL_to_postR[node->getID()] = tree->treeSize_ - 1 - node->getID();
			tree->postR_to_preL[tree->treeSize_ - 1 - node->getID()] = node->getID();
			tree->preL_to_ln[node->getID()] = currentLeaf;
			node->setParent(nodes.top());
			nodes.top()->pushChild(node);
			nodes.push(node);
			tree->pushNodeToPreL(node);
		} else if(preLSequence[cursor] == ')') {
			left--;
			vector<Node*> children = nodes.top()->getChildren();
			int size = 0;
			int sum = 0;
			int leftmostForestNum = 0;
			int rightmostForestNum = 0;
			int specialForestNum = 0;
			int leftmostSize = 0;
			int rightmostSize = 0;
			for(int i = 0; i < children.size(); i++) {
				sum += children[i]->getSubTreeSizeSum();
				size += children[i]->getSubTreeSize();
				leftmostForestNum += children[i]->getLeftmostForestNum();
				rightmostForestNum += children[i]->getRightmostForestNum();
			}
			if(!children.empty()) {
				leftmostSize = children[0]->getSubTreeSize();
				rightmostSize = children[children.size() - 1]->getSubTreeSize();
			}
			size += 1;
			sum += size;
			leftmostForestNum += size - leftmostSize;
			rightmostForestNum += size - rightmostSize;
			specialForestNum = size * (size + 3) / 2 - sum;
			if(size == 1) currentLeaf = nodes.top()->getID();
			nodes.top()->setSubTreeSize(size);
			nodes.top()->setLeftmostForestNum(leftmostForestNum);
			nodes.top()->setRightmostForestNum(rightmostForestNum);
			nodes.top()->setSpecialForestNum(specialForestNum);
			nodes.top()->setSubTreeSizeSum(sum);
			
			tree->preL_to_postL[nodes.top()->getID()] = postN;						//tree->pushNodeToPost(nodes.top());
			tree->postL_to_preL[postN++] = nodes.top()->getID();
			if(nodes.top()->getSubTreeSize() == 1) {
				tree->preL_to_lid[nodes.top()->getID()] = nodes.top()->getID();
				tree->preL_to_rid[nodes.top()->getID()] = nodes.top()->getID();
			} else {
				tree->preL_to_lid[nodes.top()->getID()] = tree->preL_to_lid[nodes.top()->getLeftmostChild()->getID()];
      			tree->preL_to_rid[nodes.top()->getID()] = tree->preL_to_rid[nodes.top()->getRightmostChild()->getID()];
			}
			tree->preL_to_preR[nodes.top()->getID()] = tree->treeSize_ - 1 - tree->preL_to_postL[nodes.top()->getID()];
			tree->preR_to_preL[tree->treeSize_ - 1 - tree->preL_to_postL[nodes.top()->getID()]] = nodes.top()->getID();
			nodes.pop();
		}
	}

	currentLeaf = -1;
    for(int i = 0; i < tree->getTreeSize(); i++) {
      tree->preR_to_ln[i] = currentLeaf;
      if((*tree)[tree->preR_to_preL[i]]->getChildren().size() == 0) {
        currentLeaf = i;
      }
    }
	return tree;
};


string RNA::toString(void) {
	string res;
	res = res + "RNAName: " + RNAName_ + "\nRNASize: " + to_string(RNASize_) + "\n";
	res = res + "RNA Sequence\n";
	for(int i = 2; i < RNASize_; i++) {
		res += originalSequence[i];
	}
/*	res += "\n";
	res = res + "RNA secondary Structure\n";
	for(int i = 2; i < RNASize_; i++) {
		res += to_string(secondaryStructure[i]) + " ";
	}*/
	return res;
};
	
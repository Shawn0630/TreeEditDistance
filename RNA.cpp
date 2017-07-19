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
	return preOrderSequenceSize_;
};

char& RNA::operator[] (int i) {
	if(i >= maxSize) {
		cout << "Overflow" << endl;
		return originalSequence[0];
	}
	return originalSequence[i];

};

int& RNA::operator() (int i) {
	if(i >= maxSize) {
		cout << "Overflow" << endl;
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

void RNA::setPreOrderSequence(string preOrderSequence_) {
	for(int i = 0; i < preOrderSequence_.length(); i++) {
		preOrderSequence[i] = preOrderSequence_[i];
	}
};

char* RNA::getPreOrderSequence(void) {
	int j = 0;
	for(int i = 1; i <= RNASize_; i++) {
		if(secondaryStructure[i] == i){
			preOrderSequence[j++] = '(';
			preOrderSequence[j++] = originalSequence[i];
			preOrderSequence[j++] = ')';
		}
		// if the base have pair, process the smaller one 
		else if(secondaryStructure[i] > i){
			preOrderSequence[j++] = '(';
			// there are altogether 16 cases.
			if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'A'){
				preOrderSequence[j++] = 'B';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'G'){
				preOrderSequence[j++] = 'D';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'C'){
				preOrderSequence[j++] = 'E';
			}
			else if(originalSequence[i] == 'A' && originalSequence[secondaryStructure[i]] == 'U'){
				preOrderSequence[j++] = 'F';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'A'){
				preOrderSequence[j++] = 'H';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'G'){
				preOrderSequence[j++] = 'I';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'C'){
				preOrderSequence[j++] = 'J';
			}
			else if(originalSequence[i] == 'G' && originalSequence[secondaryStructure[i]] == 'U'){
				preOrderSequence[j++] = 'K';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'A'){
				preOrderSequence[j++] = 'L';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'G'){
				preOrderSequence[j++] = 'M';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'C'){
				preOrderSequence[j++] = 'N';
			}
			else if(originalSequence[i] == 'C' && originalSequence[secondaryStructure[i]] == 'U'){
				preOrderSequence[j++] = 'O';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'A'){
				preOrderSequence[j++] = 'P';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'G'){
				preOrderSequence[j++] = 'Q';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'C'){
				preOrderSequence[j++] = 'R';
			}
			else if(originalSequence[i] == 'U' && originalSequence[secondaryStructure[i]] == 'U'){
				preOrderSequence[j++] = 'S';
			}
		}
		else {
			preOrderSequence[j++] = ')';
		}
	}
	preOrderSequenceSize_ = j;
	return preOrderSequence;
};

Tree* RNA::buildTree(void) {
	stack<Node*> nodes;
	Tree* tree = new Tree(RNAName_);
	int preN = 1;
	int cursor = 1;
	int postN;
	int left = 1;
	Node* node = new Node(preN, preOrderSequence[cursor]);
	nodes.push(node);
	tree->pushNodeToPre(node);
	while(left > 0) {
		cursor++;
		if(preOrderSequence[cursor] == '(') {
			left++;
			Node* node = new Node(++preN, preOrderSequence[++cursor]);
			node->setParent(nodes.top());
			nodes.top()->pushChild(node);
			nodes.push(node);
			tree->pushNodeToPre(node);
		} else if(preOrderSequence[cursor] == ')') {
			left--;
			vector<Node*> children = nodes.top()->getChildren();
			int size = 0;
			int sum = 0;
			int left = 0;
			int right = 0;
			int special = 0;
			int leftmostSize = 0;
			int rightmostSize = 0;
			for(int i = 0; i < children.size(); i++) {
				sum += children[i]->getSubTreeSizeSum();
				size += children[i]->getSubTreeSize();
				left += children[i]->getLeftmostForestNum();
				right += children[i]->getRightmostForestNum();
			}
			if(!children.empty()) {
				leftmostSize = children[0]->getSubTreeSize();
				rightmostSize = children[children.size() - 1]->getSubTreeSize();
			}
			size += 1;
			sum += size;
			left += size - leftmostSize;
			right += size - rightmostSize;
			special = size * (size + 3) / 2 - sum;
			nodes.top()->setSubTreeSize(size);
			nodes.top()->setLeftmostForestNum(left);
			nodes.top()->setRightmostForestNum(right);
			nodes.top()->setSpecialForestNum(special);
			nodes.top()->setSubTreeSizeSum(sum);
			tree->pushNodeToPost(nodes.top());
			nodes.pop();
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
	res += "\n";
	res = res + "RNA secondary Structure\n";
	for(int i = 2; i < RNASize_; i++) {
		res += to_string(secondaryStructure[i]) + " ";
	}
	return res;
};
	
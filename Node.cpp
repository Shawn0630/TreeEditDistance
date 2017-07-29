#include "Node.h"
#include <string>
#include <vector>

using namespace std;
//using namespace FPM;

Node::Node(int id, char label) {
	id_ = id;
	label_ = label;
	parent_ = NULL;
	subTreeSize_ = 0;
	right_ = 0;
	left_ = 0;
	special_ = 0;
}

int Node::getID(void) const {
	return id_;
}

char Node::getLabel(void) const {
	return label_;
}

Node* Node::getParent(void) const {
  return parent_;
}

void Node::setParent(Node* parent) {
	parent_ = parent;
}

void Node::pushChild(Node* child) {
	children_.push_back(child);
} 

int Node::getChildrenNum(void) const {
	return children_.size();
}

vector<Node*> Node::getChildren(void) const {
	return children_;
}

Node* Node::getLeftmostChild(void) const {
	if(children_.size() == 0) return NULL;
	return children_[0];
};

Node* Node::getRightmostChild(void) const {
	if(children_.size() == 0) return NULL;
	return children_[children_.size() - 1];
};


void Node::setSubTreeSize(int subTreeSize) {
	subTreeSize_ = subTreeSize;
};
int Node::getSubTreeSize(void) const {
	return subTreeSize_;
};

void Node::setSubTreeSizeSum(int subTreeSizeSum) {
	subTreeSizeSum_ = subTreeSizeSum;
};
int Node::getSubTreeSizeSum(void) const {
	return subTreeSizeSum_;
};

void Node::setLeftmostForestNum(int left) {
	left_ = left;
};
int Node::getLeftmostForestNum(void) const {
	return left_;
};

void Node::setRightmostForestNum(int right) {
	right_ = right;
};
int Node::getRightmostForestNum(void) const {
	return right_;
};

void Node::setSpecialForestNum(int special) {
	special_ = special;
} ;
int Node::getSpecialForestNum(void) const {
	return special_;
};

string Node::toString(void) const {
	string res;
	res = res + "ID: " + to_string(id_) + " label: " + label_ + "\n";
	if(parent_ == NULL) res = res + "Root" + "\n";
	else res = res + "parent is " + to_string(parent_->getID()) + "\n";
	res = res + "childrenNum: " + to_string(children_.size()) + "\n";
	for(int i = 0; i < children_.size(); i++) {
		res = res + to_string(children_[i]->getID()) + " ";
	}
	res = res + "\n";
	res = res + "subTreeSize: " + to_string(subTreeSize_) + "\n"; 
	res = res + "subTreeSizeSum: " + to_string(subTreeSizeSum_) + "\n";

	res = res + "leftmostForestNum: " + to_string(left_) + "\n";
	res = res + "rightmostForestNum: " + to_string(right_) + "\n";
	res = res + "specialmostForestNum: " + to_string(special_) + "\n";

	return res;
}
#ifndef Node_H
#define Node_H 1

#include <string>
#include <vector>

using namespace std;
//using namespace FPM;

class Node {

private:
	int id_;
	char label_;
	int subTreeSize_;
	int subTreeSizeSum_;
	int right_;
	int left_;
	int special_;
	Node* parent_;
	vector<Node*> children_;
public:
	Node(int, char);
	int getID(void) const;
	char getLabel(void) const;
	void setParent(Node*);
	Node* getParent(void) const;
	void pushChild(Node*);
	int getChildrenNum(void) const;
	void setSubTreeSize(int);
	int getSubTreeSize(void) const;
	void setSubTreeSizeSum(int);
	int getSubTreeSizeSum(void) const;
	/*
	leftmost forest = left path = right decompose
	*/
	void setLeftmostForestNum(int);
	int getLeftmostForestNum(void)const;
	/*
	rightmost forest = right path = left decompose
	*/
	void setRightmostForestNum(int);
	int getRightmostForestNum(void)const;
	void setSpecialForestNum(int);
	int getSpecialForestNum(void)const;
	vector<Node*> getChildren(void) const;
	Node* getLeftmostChild(void) const;
	Node* getRightmostChild(void) const;

	string toString(void) const;
};

#endif
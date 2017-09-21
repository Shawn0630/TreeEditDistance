#ifndef RNA_H
#define RNA_H 1

#include "Constant.h"
#include "Tree.h"

#include <string>

using namespace std;


class RNA {
private:
	string RNAName_;
	int RNASize_; // from 1 to RNASize_(inclusive)
	int treeSize_;
	int preLSequenceSize_;
	char originalSequence[maxSize];
	int secondaryStructure[maxSize];
	char preLSequence[2 * maxSize];
	int tree_to_original[maxSize];
	int original_to_tree[maxSize];
	
public:
	RNA();
	RNA(string);
	void setRNAName(string);
	string getRNAName(void) const;
	void setRNASize(int);
	int getRNASize(void) const;
	int getPreOrderSequenceSize(void) const;
	char& operator[] (int);
	int& operator() (int);

	char* getOriginalSequence(void);
	int* getSecondaryStructure(void);
	char* getPreLSequence(void);

	int getTreeSize() const;
	void setTreeSize(int);

	void setPreOrderSequence(string);

	Tree* buildTree(void);

	string toString(void);

	friend class TreeComparison;

};

#endif
#ifndef RNA_H
#define RNA_H

#include "Constant.h"
#include "Tree.h"

#include <string>

using namespace std;


class RNA {
private:
	string RNAName_;
	int RNASize_; // from 1 to RNASize_(inclusive)
	int preOrderSequenceSize_;
	char originalSequence[maxSize];
	int secondaryStructure[maxSize];
	char preOrderSequence[2 * maxSize];
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
	char* getPreOrderSequence(void);

	void setPreOrderSequence(string);

	Tree* buildTree(void);

	string toString(void);

};

#endif
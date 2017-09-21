#ifndef TREEMAP_H
#define TREEMAP_H 1

#include "SimiMatrix.h"
#include "Tree.h"

#include <string>


using namespace std;

class TreeMap{
private:
	int FTreeSize;
	int GTreeSize;

	Tree* F;
	Tree* G;

	int* FtoG;
	int* GtoF;

	SimiMatrix costModel;
	int counter;
public:

	TreeMap();

	TreeMap(Tree*, Tree*, SimiMatrix);

	void setTreeF(Tree*);
	void setTreeG(Tree*);
	void setCostModel(SimiMatrix);

	void init(void);

	void setMap(int, int);

	int operator[](int);//index F
	int operator()(int);//index G

	string toString(void);

	int getCount();

	//friend class TreeComparison;
};


#endif
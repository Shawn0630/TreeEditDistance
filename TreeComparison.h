#include "Constant.h"
#include "Strategy.h"
#include "Node.h"
#include "Tree.h"
#include "SimiMatrix.h"

#include <fstream>
using namespace std;

class TreeComparison {
private:	
	Tree* A_;
	Tree* B_;
	SimiMatrix costModel_;
	int treeSizeA;
	int treeSizeB;
	int** Free;
	int** LeftA;
	int** LeftB;
	int** RightA;
	int** RightB;
	int** AllA;
	int** AllB;
	Strategy** LeftAStrategies;
	Strategy** RightAStrategies;
	Strategy** LeftBStrategies;
	Strategy** RightBStrategies;
	Strategy** AllAStrategies;
	Strategy** AllBStrategies;
	Strategy** FreeStrategies;

	float** delta;

	int free(Node*, Node*);
	int leftA(Node*, Node*);
	int leftB(Node*, Node*);
	int rightA(Node*, Node*);
	int rightB(Node*, Node*);
	int allA(Node*, Node*);
	int allB(Node*, Node*);

	void deltaInit();
	void gted(Node*, Node*);
	void getPathType(Tree*, Node*);
	void computeSumInsAndDelCost(Tree*);
	void computeTreeDistance();

	ofstream ou;

public:
	TreeComparison();
	TreeComparison(Tree*, Tree*, SimiMatrix);
	void strategyComputation();

};
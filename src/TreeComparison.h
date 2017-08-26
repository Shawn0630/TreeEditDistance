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

	bool** hasVisited;

	float** delta;
	float** s;
	float** t;
	float* q;


	int* fn;
	int* ft;

	int fn_ft_length;

	int free(Node*, Node*);
	int leftA(Node*, Node*);
	int leftB(Node*, Node*);
	int rightA(Node*, Node*);
	int rightB(Node*, Node*);
	int allA(Node*, Node*);
	int allB(Node*, Node*);

	void deltaInit();
	float gted(Node*, Node*);
	int getPathType(Tree*, Node*, int);
	void computeSumInsAndDelCost(Tree*);
	void computeTreeDistance();

	void updateFnArray(int, int, int);
	void updateFtArray(int, int);

	float spfA(Node*, Node*, int, int, bool);
	float spfL(Node*, Node*, int, bool);
	float spfR(Node*, Node*, int, bool);

	int computeKeyRoots(Tree*, Node*, int, int*, int);
	void treeEditDist(Node*, Node*, float**, bool);	

	Strategy** APTED_ComputeOptStrategy_postL();

	ofstream ou;

public:
	TreeComparison();
	TreeComparison(Tree*, Tree*, SimiMatrix);
	void strategyComputation();

};
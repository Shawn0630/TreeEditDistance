#include "Constant.h"
#include "Strategy.h"
#include "Node.h"
#include "Tree.h"

#include <fstream>
using namespace std;

class TreeComparison {
private:	
	Tree* A_;
	Tree* B_;
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

	Path** Paths;

	int free(Node*, Node*);
	int leftA(Node*, Node*);
	int leftB(Node*, Node*);
	int rightA(Node*, Node*);
	int rightB(Node*, Node*);
	int allA(Node*, Node*);
	int allB(Node*, Node*);

	void computePath();

	ofstream ou;

public:
	TreeComparison();
	TreeComparison(Tree*, Tree*);
	void strategyComputation();

};
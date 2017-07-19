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
	int Free[maxSize][maxSize];
	int LeftA[maxSize][maxSize];
	int LeftB[maxSize][maxSize];
	int RightA[maxSize][maxSize];
	int RightB[maxSize][maxSize];
	int AllA[maxSize][maxSize];
	int AllB[maxSize][maxSize];
	Strategy* LeftAStrategies[maxSize][maxSize];
	Strategy* RightAStrategies[maxSize][maxSize];
	Strategy* LeftBStrategies[maxSize][maxSize];
	Strategy* RightBStrategies[maxSize][maxSize];
	Strategy* AllAStrategies[maxSize][maxSize];
	Strategy* AllBStrategies[maxSize][maxSize];
	Strategy* FreeStrategies[maxSize][maxSize];

	int free(Node*, Node*);
	int leftA(Node*, Node*);
	int leftB(Node*, Node*);
	int rightA(Node*, Node*);
	int rightB(Node*, Node*);
	int allA(Node*, Node*);
	int allB(Node*, Node*);

	ofstream ou;

public:
	TreeComparison();
	TreeComparison(Tree*, Tree*);
	void strategyComputation();
	void setTreeA(Tree*);
	void setTreeB(Tree*);

};
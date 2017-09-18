#include "Constant.h"
#include "Strategy.h"
#include "Node.h"
#include "Tree.h"
#include "SimiMatrix.h"
#include "TreeMap.h"

#include <fstream>
#include <map>
using namespace std;

struct Tuple{
	int Findex;
	int Gindex;
	Tuple(int Findex_, int Gindex_): Findex(Findex_), Gindex(Gindex_){};
};

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

	float treeDist;

	TreeMap* map;

	
	int free(Node*, Node*);
	int leftA(Node*, Node*);
	int leftB(Node*, Node*);
	int rightA(Node*, Node*);
	int rightB(Node*, Node*);
	int allA(Node*, Node*);
	int allB(Node*, Node*);

	void deltaInit();
	float gted(Node*, Node*);
	void gteo(Node*, Node*);
	int getPathType(Tree*, Node*, int);
	void computeSumInsAndDelCost(Tree*);
	void computeTreeDistance();

	void updateFnArray(int, int, int);
	void updateFtArray(int, int);

	float spfA_LR(Node*, Node*, int, int, bool);
	float spfA_RL(Node*, Node*, int, int, bool);
	float spfA(Node*, Node*, int, int, bool);
	float spfL(Node*, Node*, int, bool);
	float spfLL(Node*, Node*, int, bool);
	float spfR(Node*, Node*, int, bool);
	float spfRR(Node*, Node*, int, bool);
	float spf1(Node*, int, Node*, int);

	int computeKeyRoots(Tree*, Node*, int, int*, int);
	int computeRevKeyRoots(Tree*, Node*, int, int*, int);
	float treeEditDist(Node*, Node*, float**, bool, bool);
	float revTreeEditDist(Node*, Node*, float**, bool, bool);	

	Strategy** APTED_ComputeOptStrategy_postL();

	ofstream ou;

	int counter;


public:
	TreeComparison(void);
	TreeComparison(Tree*, Tree*, SimiMatrix);
	void setTreeA(Tree*);
	void setTreeB(Tree*);
	void setCostModel(SimiMatrix);
	void init(string);
	void strategyComputation(void);
	float getTreeDistance(void);
	float getTreeDistance_LL(void);
	float getTreeDistance_RR(void);
	int getCounter(void);
	TreeMap* getTreeMap(void);


};
#include <string>

using namespace std;

class Strategy {
private:
	int keynode_;
	int leaf_;
	int treeToDecompose_;// 0 for tree A 1 for tree B
	int direction_; // 0 for right 1 for left

public:
	Strategy();

	void setKeyNode(int);
	int getKeyNode(void)const;

	void setLeaf(int);
	int getLeaf(void)const;

	void setTreeToDecompose(int);
	int getTreeToDecompose(void)const;

	void setDirection(int);
	int getDirection(void)const;

	string toString()const;

};
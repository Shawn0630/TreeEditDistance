#ifndef Tree_H
#define Tree_H

#include "Node.h"
#include <string>
#include <vector>

using namespace std;
//using namespace FPM;

class Tree {

private:
	string label_;
	vector<Node*> pre;
	vector<Node*> post;

public:
	Tree();
	Tree(string);
	string getLabel(void) const;
	vector<Node*> getPre(void) const;
	vector<Node*> getPost(void) const;
	void pushNodeToPre(Node*);
	void pushNodeToPost(Node*);

};


#endif
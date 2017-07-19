#include "Node.h"
#include "Tree.h"
#include <iostream>

using namespace std;


int main() {
	Tree tree("testing");
	Node* r = new Node(0, 'A');
	Node* r_l = new Node(1, 'B');
	Node* r_r = new Node(2, 'C');
	Node* r_l_l = new Node(3, 'D');
	Node* r_l_r = new Node(4, 'E');
	Node* r_r_l = new Node(5, 'F');
	Node* r_r_r = new Node(6, 'G');

	r->pushChild(r_l);
	r->pushChild(r_r);

	r_l->setParent(r);
	r_l->pushChild(r_l_l);
	r_l_l->setParent(r_l);
	r_l->pushChild(r_l_r);
	r_l_r->setParent(r_l);

	r_r->setParent(r);
	r_r->pushChild(r_r_l);
	r_r_l->setParent(r_r);
	r_r->pushChild(r_r_r);
	r_r_r->setParent(r_r);

	cout << "Build Tree " << tree.getLabel() << endl;
	cout << endl;

	tree.pushNodeToPre(r);
	tree.pushNodeToPre(r_l);
	tree.pushNodeToPre(r_l_l);
	tree.pushNodeToPre(r_l_r);
	tree.pushNodeToPre(r_r);
	tree.pushNodeToPre(r_r_l);
	tree.pushNodeToPre(r_r_r);

	vector<Node*> pre = tree.getPre();

	for(int i = 0; i < pre.size(); i++) {
		Node* node = pre[i];
		cout << node->toString() << endl;;
	}

}
#include "Tree.h"
#include "Node.h"
#include <vector>
#include <string>

using namespace std;
//using namespace FPM;

Tree::Tree(){};

Tree::Tree(string label) {
	label_ = label;
};

string Tree::getLabel(void) const{
	return label_;
};

vector<Node*> Tree::getPre(void) const {
	return pre;
};

vector<Node*> Tree::getPost(void) const {
	return post;
};
	
void Tree::pushNodeToPre(Node* node) {
	pre.push_back(node);
};
	
void Tree::pushNodeToPost(Node* node) {
	post.push_back(node);

};
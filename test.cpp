#include <string>
#include <fstream>
#include <iostream>

#include "RNA.h"
#include "FileManage.h"
#include "TreeComparison.h"
using namespace std;

int main() {
	/*string fileName = "rna17.data";
	vector<RNA> RNAs;
	//ofstream ou("out.txt");
	FileManage* file = FileManage::getInstance();
	file->setRNAFileName(fileName);
	RNAs = file->readRNAsFromFile();
	//ou << RNAs[0].toString() << endl;
	//ou << RNAs[1].toString() << endl;
	RNAs[0].getPreOrderSequence();
	RNAs[1].getPreOrderSequence();
	Tree* t1 = RNAs[0].buildTree();
	Tree* t2 = RNAs[1].buildTree();
	TreeComparison tc(t1, t2);
	tc.strategyComputation();*/


	/*if(DEBUG) {
		for(int i = 0; i < RNAs.size(); i++) {
			char* sequence = RNAs[i].getOriginalSequence();
			int* structure = RNAs[i].getSecondaryStructure();
			char* preOrderSequence = RNAs[i].getPreOrderSequence();
			Tree* t = RNAs[i].buildTree();
			vector<Node*> pre = t->getPre();
			vector<Node*> post = t->getPost();
			ou << "Sequence" << endl;
			for(int j = 1; j < RNAs[i].getRNASize(); j++) ou << sequence[j];
			ou << endl;
			ou << "Structure" << endl;
			for(int j = 1; j < RNAs[i].getRNASize(); j++) ou << structure[j] << " ";
			ou << endl;
			ou << "preOrderSequence" << endl;
			for(int j = 0; j < RNAs[i].getPreOrderSequenceSize(); j++) ou << preOrderSequence[j];
			ou << endl;

			ou << "pre size = " << pre.size() << endl;
			for(int j = 0; j < pre.size(); j++) {
				ou << pre[j]->toString() << endl << endl;
			}
			ou << endl;

			ou << "post size = " << post.size() << endl;
			for(int j = 0; j < post.size(); j++) {
				ou << post[j]->toString() << endl << endl;			
			}
			ou << endl;

		}
	}*/

	RNA r1, r2;
	//ofstream ou("out.txt");
	r1.setRNAName("A");
	//r2.setRNAName("B");
	string s1 = "(B(C)(D(F)(G))(E))";
	string s2 = "(B(C(E(F)(G))))";
	//string s1 = "(B(C)(D)(E))";
	//string s2 = "(B(C)(D(E)))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(5);
	//r2.setPreOrderSequence(s2);
	Tree* t1 = r1.buildTree();
	//ou << t1->toString() << endl;
	Tree* t2 = r2.buildTree();
	TreeComparison tc(t1, t2);
	tc.strategyComputation();
}
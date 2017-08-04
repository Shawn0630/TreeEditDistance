#include <string>
#include <fstream>
#include <iostream>

#include "RNA.h"
#include "SimiMatrix.h"
#include "FileManage.h"
#include "TreeComparison.h"

using namespace std;

int main() {
	string fileName = "rna16.data";
	string simiFileName = "ss_simi3";
	SimiMatrix matrix;
	vector<RNA> RNAs;
	//ofstream ou("out.txt");
/*	FileManage* file = FileManage::getInstance();
	file->setRNAFileName(fileName);
	file->setSimiFileName(simiFileName);
	RNAs = file->readRNAsFromFile();
	RNAs[0].getPreLSequence();
	RNAs[1].getPreLSequence();
	Tree* t1 = RNAs[0].buildTree();
	Tree* t2 = RNAs[1].buildTree();
	matrix = file->readSimiFromFile();
	TreeComparison tc(t1, t2, matrix);
	tc.strategyComputation();*/
	/*cout << matrix.toString() << endl;
	cout << matrix['C'] << " " << matrix['A'] << " " << matrix['-'] << endl;*/
	//ou << RNAs[0].toString() << endl;
	//ou << RNAs[1].toString() << endl;
/*	RNAs[0].getPreLSequence();
	RNAs[1].getPreLSequence();
	Tree* t1 = RNAs[0].buildTree();
	Tree* t2 = RNAs[1].buildTree();
	ou << "TreeA" << endl;
	ou << t1->toString() << endl;
	ou << "TreeB" << endl;
	ou << t2->toString() << endl;
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
	FileManage* file = FileManage::getInstance();
	file->setSimiFileName(simiFileName);
	matrix = file->readSimiFromFile();
	ofstream ou("out.txt");
	r1.setRNAName("A");
	r2.setRNAName("B");



	/*
	 	       B                         B
		      / \                       / \
    	   	 C   D                     C   D
    	  	/ \                           / \
    	    E  F                          E   F 
           / \                               / \
       	  G   H                             G   H
         / \                                   / \
        I   J                                 I   J

	*/
/*	string s1 = "(B(C(E(G(I)(J))(H))(F))(D))";
	string s2 = "(B(C)(D(E)(F(G)(H(I)(J)))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(9);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(9);*/


	/*
		 A                               A
	   / | \                            / \
      B  C  E                          B   D
         |                             |
         D                             C
	
	*/

    string s1 = "(A(B)(C(D))(E))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(5);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);

    /*
		   A                                 A
	   / / | \ \                            / \
      B C  D  F G                          B   D
           |                               |
           E                               C
	
	*/

/*    string s1 = "(A(B)(C)(D(E))(F)(G))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(7);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);*/

    /*
		   A                                 A
	   / / | \ \ \                            / \
      B C  D  F G H                          B   D
           |                               |
           E                               C
	
	*/

  /*  string s1 = "(A(B)(C)(D(E))(F)(G)(H))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(8);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);*/

	Tree* t1 = r1.buildTree();
	Tree* t2 = r2.buildTree();
	//ou << t2->toString() << endl;
	TreeComparison tc(t1, t2, matrix);
	tc.strategyComputation();
/*	ou << "TreeA" << endl;
	ou << t1->toString() << endl;
	ou << "TreeB" << endl;
	ou << t2->toString() << endl;*/


}
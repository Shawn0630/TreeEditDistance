#include <string>
#include <fstream>
#include <iostream>

#include "src/RNA.h"
#include "src/Errors.h"
#include "src/SimiMatrix.h"
#include "src/FileManage.h"
#include "src/TreeComparison.h"

using namespace std;

int main() {
	string fileName = "rna24.data";
	string simiFileName = "ss_distance";
	SimiMatrix matrix;
	vector<RNA> RNAs;
	FileManage* file = FileManage::getInstance();
	file->setRNAFileName(fileName);
	file->setSimiFileName(simiFileName);
	RNAs = file->readRNAsFromFile();
	RNAs[0].getPreLSequence();
	RNAs[1].getPreLSequence();
	Tree* t1 = RNAs[0].buildTree();
	Tree* t2 = RNAs[1].buildTree();
	matrix = file->readSimiFromFile();
	TreeComparison tc(t1, t2, matrix);
	tc.strategyComputation();
	float dist = tc.getTreeDistance();
	cout << "The distance is " << dist << endl;
/*	RNA r1, r2;
	FileManage* file = FileManage::getInstance();
	file->setSimiFileName(simiFileName);
	matrix = file->readSimiFromFile();
	if(readSimiFileError) {
		cout << "Error in reading simi file!" << endl;
		exit(1);
	}
	ofstream ou("out2.txt");
	r1.setRNAName("A");
	r2.setRNAName("B");
*/

/*	file->setRNAFileName(fileName);
	RNAs = file->readRNAsFromFile();
	if(readRNAFileError) {
		cout << "Error in reading RNA file!" << endl;
		exit(1);
	}
	RNAs[0].getPreLSequence();
	RNAs[1].getPreLSequence();
	Tree* t1 = RNAs[0].buildTree();
	Tree* t2 = RNAs[1].buildTree();
	ou << "TreeA" << endl;
	ou << t1->toString() << endl;
	ou << "TreeB" << endl;
	ou << t2->toString() << endl;
	TreeComparison tc(t1, t2, matrix);
	tc.strategyComputation();
*/



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
	r2.setTreeSize(9);
*/


/*
		 A                               A
	   / | \                            / \
      B  C  E                          B   D
         |                             |
         D                             C
	
*/
/*    string s1 = "(A(B)(C(D))(E))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(5);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);
*/


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
    r2.setTreeSize(4);
*/


/*
		   A                                   A
	   / / | \ \ \                            / \
      B C  D  F G H                          B   D
           |                                 |
           E                                 C	
*/
/*  string s1 = "(A(B)(C)(D(E))(F)(G)(H))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(8);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);
*/



/*
          A                              A 
    	/ | \                          / | \
       B  C  E                        B  C  E
          |                              |
          D                              D
*/
/*  string s1 = "(A(B)(C(D))(E))";
    string s2 = "(A(B)(C(D))(E))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(5);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(5);
*/



/*
	  	  A		        A                                   
	    /  \   		/ / | \ \ \ \ \                         
       B   D 		B C D  F G H I J                        
       |       		    |                                 
       C       		    E                                 
	
*/
/*  string s1 = "(A(B(C))(D))";
    string s2 = "(A(B)(C)(D(E))(F)(G)(H)(I)(J))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(4);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(10);
*/


/*
	  	  A		        A                                   
	    /  \   		/ / | \ \ \ \ \                         
       B   D 		B C D  F G H I J                        
       |       		    |                                 
       C       		    E                                 
	
*/
/*  string s1 = "(A(B)(C)(D(E))(F)(G)(H)(I)(J))";
    string s2 = "(A(B(C))(D))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(10);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(4);
*/



/*
	  	  A		        A                                   
	    /  \   		/ / | \ \ \ \ \                         
       B    C 		B C D  F G H I J                        
              		    |                                 
              		    E                                 
	
*/
/*  
	string s1 = "(A(B)(C)(D(E))(F)(G)(H)(I)(J))";
    string s2 = "(A(B)(C))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(10);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(3);
*/



/*        	    A                          A
     	  /     |   \                     / \
    	 B      C    D                   B   C
   	   / | \   / \    \                    / | \
  	  E  F G  H   I    J                  D  E  F	  
*/
/* 
    string s1 = "(A(B(E)(F)(G))(C(H)(I))(D(J)))";
    string s2 = "(A(B)(C(D)(E)(F)))";
    r1.setPreOrderSequence(s1);
    r1.setTreeSize(10);
    r2.setPreOrderSequence(s2);
    r2.setTreeSize(6);
*/


/*
	 	         B                          B
		      /  | \                       / \
    	   	 C   E  D                     C   D
    	  	/ \  |                           / \
    	   E  F  F...F                      E   F 
          / \                                  / \
       	 G   H                                G   H
         / \                                     / \
        I     J                                 I   J
		  /   |   \
         B    C    D
           /  | \
          E   F   G
            / | \
           B  C  D
*/
/*	string s1 = "(B(C(E(G(I)(J(B)(C(E)(F(B)(C)(D))(G))(D)))(H))(F))(E(F)(F)(F)(F)(F)(F)(F)(F)(F)(F))(D))";
	string s2 = "(B(C)(D(E)(F(G)(H(I)(J)))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(29);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(9);
*/



/*	              B                               B
	       /      |    \                         / \
	      C       D     E                       C   D
	    / | \   / | \                              
       F  G  H F  G  H                            
      / \       / | \                                 
     I  J      I  J  K                              
                                                       
                                                      
*/
/*	string s1 = "(B(C(F(I)(J))(G)(H))(D(F)(G(I)(J)(K))(H))(E))";
	string s2 = "(B(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(15);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(3);
*/


/*	                                          B
	                                         / \
	       B                                C   D
	     / | \                                 / \
        C  D  E                               E   F
            / | \                                / \ 
           F  G  H                              G   H
                                                   / \
                                                  I   J
*/
/*  string s1 = "(B(C)(D)(E(F)(G)(H)))";
	string s2 = "(B(C)(D(E)(F(G)(H(I)(J)))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(7);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(9);
*/


/*
	      B                                        B
	 /  /   \  \                              /  /    \  \ 
    C  D     E  F                            C  D      E  F
     / | \    \                                 |          \
    G  H  I    G                                G           G
   / \                                                     / \
  K   L                                                   H   I
 */
 /* string s1 = "(B(C)(D(G(K)(L))(H)(I))(E(G))(F))";
	//string s2 = "(B(C)(D(G))(E)(F(G(H)(I))))";
	string s2 = "(B(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(11);
	r2.setPreOrderSequence(s2);
	//r2.setTreeSize(9);
	r2.setTreeSize(3);
*/


/*	      
          B                                        B
	 /  /   \  \                              /  /    \  \ 
    C  D     E  F                            C  D      E  F
   /  / \    \                                  |          \
  G  H  I     G                                 G          G
 / \     \                                                / \
K   L     K                                              H   I
   / \     \
  K   L     K
             \
              L
               \
                L
*/
/*  string s1 = "(B(C(G(K)(L(K)(L))))(D(H)(I(K(K(L(L))))))(E(G))(F))";
	//string s2 = "(B(C)(D(G))(E)(F(G(H)(I))))";
	string s2 = "(B(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(17);
	r2.setPreOrderSequence(s2);
	//r2.setTreeSize(9);
	r2.setTreeSize(3);
*/


/*	 
     B                              B
	 |                             / \
	 C                            C   D
	 |
	 D
	/ \
   E   F
       |
       G
       |
       H
*/
/*  	string s1 = "(B(C(D(E)(F(G(H))))))";
	string s2 = "(B(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(7);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(3);*/


/*	 
     B                               B
	 |                             / | \
	 C                            C  D  F
	 |                               |
	 D                               E
	/ \
   E   F
   |   
   G   
   |   
   H   
*/
/*    string s1 = "(B(C(D(E(G(H)))(F))))";
	string s2 = "(B(C)(D(E))(F))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(7);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(5);*/


/*
d = 2
*/

/*	string s1 = "(F(A)(E(C(B))(D)))";
	string s2 = "(F(A)(C(E(B)(D))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(6);*/


/*
d = 1
*/

/*    string s1 = "(A(B)(C(E)(F))(D))";
	string s2 = "(A(B)(C(E))(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(5);*/

/*
d = 2
*/
/*    string s1 = "(A(B)(C(E)(F))(D))";
	string s2 = "(A(B)(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(4);*/

/*
d = 3
*/

/*    string s1 = "(A(S)(C(E)(F))(D))";
	string s2 = "(A(B)(C)(D))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(4);*/

/*
d = 3
*/
/*    string s1 = "(A(B)(C(E)(F))(D))";
	string s2 = "(A(B)(C))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(3);*/


/*
d = 12
*/
/*    string s1 = "(A(B(I)(J(U)))(C(D)(E(Q(N)(M))))(F(S)))";
	string s2 = "(A)";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(13);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(1);*/

/*
d = 13
*/

/*    string s1 = "(A(B(I)(J(U)))(C(D)(E(Q(N)(M))))(F(S)))";
	string s2 = "(K)";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(13);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(1);*/

/*
d = 7
*/

/*    string s1 = "(K)";
	string s2 = "(A(D)(E(L)(T(G)(H))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(1);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(7);*/


/*
d = 6
*/
/*    string s1 = "(A(B)(C(D(F)(G(H)(I)))(E)))";
	string s2 = "(G(H)(I))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(9);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(3);*/

/*
d = 3 something wrong
*/
/*    string s1 = "(A(A(A)(A(A)(A)))(A(A)(A(A))(A)))";
	string s2 = "(A(A(A)(A)(A))(A(A(A)(A)(A))))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(11);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(10);*/

/*
d = 5
*/
/*    string s1 = "(A(B(D)(E))(C))";
	string s2 = "(F(G(H)(I))(K))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(5);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(5);*/

/*
d = 2
*/

/*    string s1 = "(A(B(C)(D(E)(F))))";
	string s2 = "(B(C)(E)(F))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(4);*/

/*
d = 2
*/
/*  	string s1 = "(A(B(C)(D(E)(F))))";
	string s2 = "(B(C)(E)(F))";
	r1.setPreOrderSequence(s1);
	r1.setTreeSize(6);
	r2.setPreOrderSequence(s2);
	r2.setTreeSize(4);*/

/*	Tree* t1 = r1.buildTree();
	Tree* t2 = r2.buildTree();
	ou << "TreeA" << endl;
	ou << t1->toString() << endl;
	ou << "TreeB" << endl;
	ou << t2->toString() << endl;
	TreeComparison tc(t1, t2, matrix);
	tc.strategyComputation();
	ou << "Counter" << endl;
	ou << to_string(tc.getCounter()) << endl;*/


}

#include "TreeComparison.h"
#include "Tree.h"

#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <climits>

using namespace std;

TreeComparison::TreeComparison() {
}

TreeComparison::TreeComparison(Tree* A, Tree* B) {
	A_ = A;
	B_ = B;
	treeSizeA = A_->getTreeSize();
	treeSizeB = B_->getTreeSize();


	Free = new int*[treeSizeA];
	LeftA = new int*[treeSizeA];
	LeftB = new int*[treeSizeA];
	RightA = new int*[treeSizeA];
	RightB = new int*[treeSizeA];
	AllA = new int*[treeSizeA];
	AllB = new int*[treeSizeA];
	LeftAStrategies = new Strategy*[treeSizeA];
	LeftBStrategies = new Strategy*[treeSizeA];
	RightAStrategies = new Strategy*[treeSizeA];
	RightBStrategies = new Strategy*[treeSizeA];
	AllAStrategies = new Strategy*[treeSizeA];
	AllBStrategies = new Strategy*[treeSizeA];
	FreeStrategies = new Strategy*[treeSizeA];
	for(int i = 0; i < treeSizeA; i++) {
		Free[i] = new int[treeSizeB];
		LeftA[i] = new int[treeSizeB];
		LeftB[i] = new int[treeSizeB];
		RightA[i] = new int[treeSizeB];
		RightB[i] = new int[treeSizeB];
		AllA[i] = new int[treeSizeB];
		AllB[i] = new int[treeSizeB];
		LeftAStrategies[i] = new Strategy[treeSizeB];
		LeftBStrategies[i] = new Strategy[treeSizeB];
		RightAStrategies[i] = new Strategy[treeSizeB];
		RightBStrategies[i] = new Strategy[treeSizeB];
		AllAStrategies[i] = new Strategy[treeSizeB];
		AllBStrategies[i] = new Strategy[treeSizeB];
		FreeStrategies[i] = new Strategy[treeSizeB];
	}

	for(int i = 0; i < treeSizeA; i++) {
		for(int j = 0; j < treeSizeB; j++) {
			Free[i][j] = -1;
			LeftA[i][j] = -1;
			LeftB[i][j] = -1;
			RightA[i][j] = -1;
			RightB[i][j] = -1;
			AllA[i][j] = -1;
			AllB[i][j] = -1;
		}
	}

	ou.open("out.txt");

};


/*void TreeComparison::strategyComputation() {
	vector<Node*> postA = A_->getPost();
	vector<Node*> postB = B_->getPost();

	for(int i = 0; i < postA.size(); i++) {
		vector<Node*> children = postA[i]->getChildren();
		
		for(int j = 0; j < postB.size(); j++) {
			if(DEBUG) {
				ou << "Calculate" << endl;
				ou << postA[i]->toString() << endl;
				ou << postB[j]->toString() << endl;
			}
			if(children.empty()) {
				LeftA[postA[i]->getID()][postB[j]->getID()] = postB[j]->getLeftmostForestNum();
				RightA[postA[i]->getID()][postB[j]->getID()] = postB[j]->getRightmostForestNum();
				AllA[postA[i]->getID()][postB[j]->getID()] = postB[j]->getSpecialForestNum();
				if(DEBUG) {
					ou << "LeftA[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << to_string(postB[j]->getLeftmostForestNum()) << endl;
					ou << "RightA[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << to_string(postB[j]->getRightmostForestNum()) << endl;
					ou << "AllA[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << to_string(postB[j]->getSpecialForestNum()) << endl;
				}
			} else {
				int freeSum = 0;
				vector<int> subTreeSizeSum;
				int prevSubTreeSizeSum = 0;
				for(int z = 0; z < children.size(); z++) {
					freeSum += Free[children[z]->getID()][postB[j]->getID()];
					prevSubTreeSizeSum += children[z]->getSubTreeSize();
					subTreeSizeSum.push_back(prevSubTreeSizeSum);
				}
				if(DEBUG) {
					ou << "Compute LeftA" << endl;
				}
				//Compute LeftA inherited from the last right decomposition.
				Strategy* leftAS = new Strategy();
				leftAS->setKeyNode(children[0]->getID());
				leftAS->setDirection(0);
				leftAS->setTreeToDecompose(0);
				int leftALeftmost = freeSum - Free[children[0]->getID()][postB[j]->getID()] + LeftA[children[0]->getID()][postB[j]->getID()] + postB[j]->getLeftmostForestNum() * (postA[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				int leftAMin = leftALeftmost;
				if(DEBUG) {
					ou << "If select Tree A " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(leftALeftmost) << " Direction: Right" << endl;
				}
				for(int z = 1; z < children.size() - 1; z++) {
					if(DEBUG) {
						ou << "If select Tree A " << to_string(children[z]->getID()) << " #Subproblem: ";
					}
					int prefix = freeSum - Free[children[z]->getID()][postB[j]->getID()] + AllA[children[z]->getID()][postB[j]->getID()];
					int left = postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					int right = postB[j]->getLeftmostForestNum() * (1 + subTreeSizeSum[children.size() - 1] - subTreeSizeSum[z - 1]) + postB[j]->getSpecialForestNum() * (subTreeSizeSum[z - 1]);
					int sum = 0;
					if(left > right) sum = prefix + right;
					else sum = prefix + left;
					if(leftAMin > sum) {
						leftAMin = sum;
						leftAS->setKeyNode(children[z]->getID());
						if(left > right) leftAS->setDirection(0);
						else leftAS->setDirection(1);
					}
					if(DEBUG) {
						ou << to_string(sum) << " Direction: ";
						if(left > right) ou << "Right" << endl;
						else ou << "Left" << endl;
					}
				}
				int leftARightmost = freeSum - Free[children[children.size() - 1]->getID()][postB[j]->getID()] + AllA[children[children.size() - 1]->getID()][postB[j]->getID()] + postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[children.size() - 1]->getSubTreeSize());
				if(leftAMin > leftARightmost) {
					leftAS->setKeyNode(children[children.size() - 1]->getID());
					leftAS->setDirection(1);
				}
				if(DEBUG) {
					ou << "If select Tree A " << to_string(children[children.size() - 1]->getID()) << " #Subproblem: " << to_string(leftARightmost) << " Direction: Left" << endl;
				}
				LeftA[postA[i]->getID()][postB[j]->getID()] = leftAMin;
				LeftAStrategies[postA[i]->getID()][postB[j]->getID()] = leftAS;
				if(DEBUG) {
					ou << "LeftAStrategies[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << leftAS->toString() << endl;
				}


				if(DEBUG) {
					ou << "Compute RightA" << endl;
				}
				//Compute RightA inherited from the last left decomposition.
				Strategy* rightAS = new Strategy();
				rightAS->setKeyNode(children[0]->getID());
				rightAS->setTreeToDecompose(0);
				rightAS->setDirection(1);
				int rightARightmost = freeSum - Free[children[children.size() - 1]->getID()][postB[j]->getID()] + RightA[children[children.size() - 1]->getID()][postB[j]->getID()] + postB[j]->getRightmostForestNum() * (postA[i]->getSubTreeSize() - children[children.size() - 1]->getSubTreeSize());
				int rightAMin = rightARightmost;
				if(DEBUG) {
					ou << "If select Tree A " << to_string(children[children.size() - 1]->getID()) << " #Subproblem: " << to_string(rightARightmost) << " Direction: Left" << endl;
				}
				for(int z = 1; z < children.size() - 1; z++) {
					if(DEBUG) {
						ou << "If select Tree A " << to_string(children[z]->getID()) << " #Subproblem: ";
					}
					int prefix = freeSum - Free[children[z]->getID()][postB[j]->getID()] + AllA[children[z]->getID()][postB[j]->getID()];
					int left = postB[j]->getRightmostForestNum() * (1 + subTreeSizeSum[z - 1]) + postB[j]->getSpecialForestNum() * (subTreeSizeSum[children.size() - 1] - subTreeSizeSum[z]);
					int right = postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					int sum = 0;
					if(left > right) sum = prefix + right;
					else sum = prefix + left;
					if(rightAMin > sum) {
						rightAMin = sum;
						rightAS->setKeyNode(children[z]->getID());
						if(left > right) rightAS->setDirection(0);
						else rightAS->setDirection(1);
					}
					if(DEBUG) {
						ou << to_string(sum) << " Direction: ";
						if(left > right) ou << "Right" << endl;
						else ou << "Left" << endl;
					}
				}
				int rightALeftmost = freeSum - Free[children[0]->getID()][postB[j]->getID()] + AllA[children[children.size() - 1]->getID()][postB[j]->getID()] + postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				if(rightAMin > rightALeftmost) {
					rightAS->setKeyNode(children[0]->getID());
					rightAS->setDirection(0);
				}
				if(DEBUG) {
					ou << "If select Tree A " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(rightALeftmost) << " Direction: Right" << endl;
				}
				RightA[postA[i]->getID()][postB[j]->getID()] = rightAMin;
				RightAStrategies[postA[i]->getID()][postB[j]->getID()] = rightAS;
				if(DEBUG) {
					ou << "RightAStrategies[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << rightAS->toString() << endl;
				}


				//Compute AllA inherited from the last left-right decomposition
				Strategy* allAS = new Strategy();
				allAS->setKeyNode(children[0]->getID());
				allAS->setTreeToDecompose(0);
				int allAmin = freeSum - Free[children[0]->getID()][postB[j]->getID()] + AllA[children[0]->getID()][postB[j]->getID()] + postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				if(DEBUG) {
					ou << "If select Tree A " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(allAmin) << endl;
				}
				for(int z = 1; z < children.size(); z++) {
					int sum = freeSum - Free[children[0]->getID()][postB[j]->getID()] + AllA[children[z]->getID()][postB[j]->getID()] + postB[j]->getSpecialForestNum() * (postA[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					if(DEBUG) {
					ou << "If select Tree A " << to_string(children[z]->getID()) << " #Subproblem: " << to_string(sum) << endl;
					}
					if(allAmin > sum) {
						allAmin = sum;
						allAS->setKeyNode(children[z]->getID());
					}
				}
				AllA[postA[i]->getID()][postB[j]->getID()] = allAmin;
				AllAStrategies[postA[i]->getID()][postB[j]->getID()] = allAS;
				if(DEBUG) {
					ou << "AllAStrategies[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << allAS->toString() << endl;
					ou << endl;
				}
			}
		}
	}

	for(int i = 0; i < postB.size(); i++) {
		vector<Node*> children = postB[i]->getChildren();
	
		for(int j = 0; j < postA.size(); j++) {
			if(DEBUG) {
				ou << "Calculate" << endl;
				ou << postA[j]->toString() << endl;
				ou << postB[i]->toString() << endl;
			}
			if(children.empty()) {
				LeftB[postA[j]->getID()][postB[i]->getID()] = postA[j]->getLeftmostForestNum();
				RightB[postA[j]->getID()][postB[i]->getID()] = postA[j]->getRightmostForestNum();
				AllB[postA[j]->getID()][postB[i]->getID()] = postA[j]->getSpecialForestNum();
				if(DEBUG) {
					ou << "LeftB[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << to_string(postA[j]->getLeftmostForestNum()) << endl;
					ou << "RightB[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << to_string(postA[j]->getRightmostForestNum()) << endl;
					ou << "AllB[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << to_string(postA[j]->getSpecialForestNum()) << endl;
				}
			} else {
				int freeSum = 0;
				vector<int> subTreeSizeSum;
				int prevSubTreeSizeSum = 0;
				
				for(int z = 0; z < children.size(); z++) {
					freeSum += Free[postA[j]->getID()][children[z]->getID()];
					prevSubTreeSizeSum += children[z]->getSubTreeSize();
					subTreeSizeSum.push_back(prevSubTreeSizeSum);
				}

				if(DEBUG) {
					ou << "Compute LeftB" << endl;
				}
				//Compute LeftB inherited from the last right decomposition.
				Strategy* leftBS = new Strategy();
				leftBS->setKeyNode(children[0]->getID());
				leftBS->setTreeToDecompose(1);
				leftBS->setDirection(0);
				if(DEBUG) {
					ou << "Free[" << postA[j]->getID() << "][" << children[0]->getID() << "] = " << Free[postA[j]->getID()][children[0]->getID()] << endl;
					ou << "LeftB[" << postA[j]->getID() << "][" << children[0]->getID() << "] = " << LeftB[postA[j]->getID()][children[0]->getID()] << endl;
				}
				int leftBLeftmost = freeSum - Free[postA[j]->getID()][children[0]->getID()] + LeftB[postA[j]->getID()][children[0]->getID()] + postA[j]->getLeftmostForestNum() * (postB[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				int leftBMin = leftBLeftmost;
				if(DEBUG) {
					ou << "If select Tree B " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(leftBLeftmost) << " Direction: Right" << endl;
				}
				for(int z = 1; z < children.size() - 1; z++) {
					if(DEBUG) {
						ou << "If select Tree B " << to_string(children[z]->getID()) << " #Subproblem: ";
					}
					int prefix = freeSum - Free[postA[j]->getID()][children[z]->getID()] + AllB[postA[j]->getID()][children[z]->getID()];
					int left = postA[j]->getSpecialForestNum() * (postB[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					int right = postA[j]->getLeftmostForestNum() * (1 + subTreeSizeSum[children.size() - 1] - subTreeSizeSum[z]) + postA[j]->getSpecialForestNum() * (subTreeSizeSum[z - 1]);
					int sum = 0;
					if(left > right) sum = prefix + right;
					else sum = prefix + left;
					if(DEBUG) {
						ou << to_string(sum) << " Direction: ";
						if(left > right) ou << "Right" << endl;
						else ou << "Left" << endl;
					}
					if(leftBMin > sum) {
						leftBMin = sum;
						leftBS->setKeyNode(children[z]->getID());
						if(left > right) leftBS->setDirection(0);
						else leftBS->setDirection(1);
					}
				}
				if(DEBUG) {
					ou << "Free[" << postA[j]->getID() << "][" << children[children.size() - 1]->getID() << "] = " << Free[postA[j]->getID()][children[0]->getID()] << endl;
					ou << "AllB[" << postA[j]->getID() << "][" << children[children.size() - 1]->getID() << "] = " << AllB[postA[j]->getID()][children[0]->getID()] << endl;
				}
				int leftBRightmost = freeSum - Free[postA[j]->getID()][children[children.size() - 1]->getID()] + AllB[postA[j]->getID()][children[children.size() - 1]->getID()] + postA[j]->getSpecialForestNum() * (postB[i]->getSubTreeSize() - children[children.size() - 1]->getSubTreeSize());
				if(DEBUG) {
					ou << "If select Tree B " << to_string(children[children.size() - 1]->getID()) << " #Subproblem: " << to_string(leftBRightmost) << " Direction: Left" << endl;
				}
				if(leftBMin > leftBRightmost) {
					leftBS->setKeyNode(children[children.size() - 1]->getID());
					leftBS->setDirection(1);
				}
				LeftB[postA[j]->getID()][postB[i]->getID()] = leftBMin;
				LeftBStrategies[postA[j]->getID()][postB[i]->getID()] = leftBS;
				if(DEBUG) {
					ou << "LeftBStrategies[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << leftBS->toString() << endl;
				}

				if(DEBUG) {
					ou << "Compute RightB" << endl;
				}


				//Compute RightB inherited from the last left decomposition.
				Strategy* rightBS = new Strategy();
				rightBS->setKeyNode(children[0]->getID());
				rightBS->setTreeToDecompose(1);
				rightBS->setDirection(1);
				int rightBRightmost = freeSum - Free[postA[j]->getID()][children[children.size() - 1]->getID()] + RightB[postA[j]->getID()][children[children.size() - 1]->getID()] + postA[j]->getRightmostForestNum() * (postB[i]->getSubTreeSize() - children[children.size() - 1]->getSubTreeSize());
				int rightBMin = rightBRightmost;
				if(DEBUG) {
					ou << "If select Tree B " << to_string(children[children.size() - 1]->getID()) << " #Subproblem: " << to_string(rightBRightmost) << " Direction: Left" << endl;
				}
				for(int z = 1; z < children.size() - 1; z++) {
					if(DEBUG) {
						ou << "If select Tree B " << to_string(children[z]->getID()) << " #Subproblem: ";
					}
					int prefix = freeSum - Free[postA[j]->getID()][children[z]->getID()] + AllB[postA[j]->getID()][children[z]->getID()];
					int left = postA[j]->getRightmostForestNum() * (1 + subTreeSizeSum[z - 1]) + postB[i]->getSpecialForestNum() * (subTreeSizeSum[children.size() - 1] - subTreeSizeSum[z]);
					int right = postA[j]->getSpecialForestNum() * (postB[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					int sum = 0;
					if(left > right) sum = prefix + right;
					else sum = prefix + left;
					if(rightBMin > sum) {
						rightBMin = sum;
						rightBS->setKeyNode(children[z]->getID());
						if(left > right) rightBS->setDirection(0);
						else rightBS->setDirection(1);
					}
					if(DEBUG) {
						ou << to_string(sum) << " Direction: ";
						if(left > right) ou << "Right" << endl;
						else ou << "Left" << endl;
					}
				}
				int rightBLeftmost = freeSum - Free[postA[j]->getID()][children[children.size() - 1]->getID()] + AllB[postA[j]->getID()][children[children.size() - 1]->getID()] + postA[j]->getLeftmostForestNum() * (postB[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				if(DEBUG) {
					ou << "If select Tree B " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(rightBLeftmost) << " Direction: Right" << endl;
				}
				if(rightBMin > rightBLeftmost) {
					rightBMin = rightBLeftmost;
					rightBS->setKeyNode(children[children.size() - 1]->getID());
					rightBS->setDirection(0);
				}
				RightB[postA[j]->getID()][postB[i]->getID()] = rightBMin;
				RightBStrategies[postA[j]->getID()][postB[i]->getID()] = rightBS;
				if(DEBUG) {
					ou << "rightBStrategies[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << rightBS->toString() << endl;
				}

				if(DEBUG) {
					ou << "Compute AllB" << endl;
				}

				//Compute AllA inherited from the last left-right decomposition
				Strategy* allBS = new Strategy();
				allBS->setKeyNode(children[0]->getID());
				allBS->setTreeToDecompose(1);
				int allBmin = freeSum - Free[postA[j]->getID()][children[0]->getID()] + AllA[postA[j]->getID()][children[0]->getID()] + postA[j]->getSpecialForestNum() * (postB[i]->getSubTreeSize() - children[0]->getSubTreeSize());
				if(DEBUG) {
					ou << "If select Tree B " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(allBmin) << endl;
				}
				for(int z = 1; z < children.size(); z++) {
					int sum = freeSum - Free[children[0]->getID()][postB[i]->getID()] + AllA[children[0]->getID()][postB[i]->getID()] + postA[j]->getSpecialForestNum() * (postB[i]->getSubTreeSize() - children[z]->getSubTreeSize());
					if(allBmin > sum) {
						allBmin = sum;
						allBS->setKeyNode(children[z]->getID());
					}
					if(DEBUG) {
					ou << "If select Tree B " << to_string(children[0]->getID()) << " #Subproblem: " << to_string(sum) << endl;
					}
				}
				AllB[postA[j]->getID()][postB[i]->getID()] = allBmin;
				AllBStrategies[postA[j]->getID()][postB[i]->getID()] = allBS;
				if(DEBUG) {
					ou << "AllBStrategies[" << to_string(postA[j]->getID()) << "][" << to_string(postB[i]->getID()) << "] = " << allBS->toString() << endl;
				}
			}
		}
	}

	if(DEBUG) {
		ou << "Compute Free" << endl;
	}


	for(int i = 0; i < postA.size(); i++) {
		for(int j = 0; j < postB.size(); j++) {
			
			int min = 0;
			int sum = 0;
			int prefix = 0;
			int left = 0;
			int right = 0;
			Strategy* freeS = new Strategy();

			vector<Node*> childrenA = postA[i]->getChildren();
			if(childrenA.empty()) {

				left = postB[j]->getRightmostForestNum();
				right = postB[j]->getLeftmostForestNum();
				if(left > right) {
					min = right;
					freeS->setDirection(0);
				} else {
					min = left;
					freeS->setDirection(1);
				}

				if(DEBUG) {
					ou << to_string(postA[i]->getID()) << " in Tree A has no child." << endl;
					ou << "If direction is left, then #Subproblem: " << to_string(left) << endl;
					ou << "If direction is right, then #Subproblem: " << to_string(right) << endl; 
				}
			}  else {
			int freeSumA = 0;
			vector<int> subTreeSizeSumA;
			int prevSubTreeSizeSumA = 0;
			for(int z = 0; z < childrenA.size(); z++) {
				freeSumA += Free[childrenA[z]->getID()][postB[j]->getID()];
				prevSubTreeSizeSumA += childrenA[z]->getSubTreeSize();
				subTreeSizeSumA.push_back(prevSubTreeSizeSumA);
			}
			

			prefix = freeSumA - Free[childrenA[0]->getID()][postB[j]->getID()] + LeftA[childrenA[0]->getID()][postB[j]->getID()] + postB[j]->getLeftmostForestNum() * (postA[i]->getSubTreeSize() - childrenA[0]->getSubTreeSize());
			right = 0;
			left = 0;
			sum = prefix + right;
			if(DEBUG) {
				ou << "If select Tree A " << to_string(childrenA[0]->getID()) << " #Subproblem: " << to_string(sum) << " Direction: Right";
			}
			if(min > sum) {
				freeS->setKeyNode(childrenA[0]->getID());
				freeS->setTreeToDecompose(0);
				freeS->setDirection(0);
				min = sum;
			}

			for(int z = 1; z < childrenA.size() - 1; z++) {
				if(DEBUG) {
					ou << "If select Tree A " << to_string(childrenA[z]->getID()) << " #Subproblem: ";
				}

				prefix = freeSumA - Free[childrenA[z]->getID()][postB[j]->getID()] + AllA[childrenA[z]->getID()][postB[j]->getID()];
				left = postB[j]->getRightmostForestNum() * (1 + subTreeSizeSumA[z - 1]) + postB[j]->getSpecialForestNum() * (subTreeSizeSumA[childrenA.size() - 1] - subTreeSizeSumA[z]);
				right = postB[j]->getLeftmostForestNum() * (1 + subTreeSizeSumA[childrenA.size() - 1] - subTreeSizeSumA[z]) + postB[j]->getSpecialForestNum() * (subTreeSizeSumA[z - 1]);
				if(left > right) sum = prefix + right;
				else sum = prefix + left;
				if(min > sum){
					min = sum;
					freeS->setKeyNode(childrenA[z]->getID());
					freeS->setTreeToDecompose(0);
					if(left > right)freeS->setDirection(0);
					else freeS->setDirection(1);
				}
				if(DEBUG) {
					ou << to_string(sum) << " Direction: ";
					if(left > right) ou << "Right" << endl;
					else ou << "Left" << endl;
				}

			}

			prefix = freeSumA - Free[childrenA[childrenA.size() - 1]->getID()][postB[j]->getID()] + RightA[childrenA[childrenA.size() - 1]->getID()][postB[j]->getID()] + postB[j]->getRightmostForestNum() * (postA[i]->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
			right = 0;
			left = 0;
			sum = prefix + left;
			if(DEBUG) {
				ou << "If select Tree A " << to_string(childrenA[0]->getID()) << " #Subproblem: " << to_string(sum) << " Direction: Left" << endl;
			}
			if(min > sum) {
				min = sum;
				freeS->setKeyNode(childrenA[childrenA.size() - 1]->getID());
				freeS->setTreeToDecompose(0);
				freeS->setDirection(1);
			}
			}
			
			vector<Node*> childrenB= postB[j]->getChildren();
			if(childrenB.empty()) {
				left = postA[i]->getRightmostForestNum();
				right = postA[i]->getLeftmostForestNum();
				if(min > right) {
					min = right;
					freeS->setDirection(0);
				}
				if(min > left) {
					min = left;
					freeS->setDirection(1);
				}
				if(DEBUG) {
					ou << to_string(postA[i]->getID()) << " in Tree B has no child." << endl;
					ou << "If direction is left, then #Subproblem: " << to_string(left) << endl;
					ou << "If direction is right, then #Subproblem: " << to_string(right) << endl; 
				}

			} else {
			int freeSumB = 0;
			vector<int> subTreeSizeSumB;
			int prevSubTreeSizeSumB = 0;
			for(int z = 0; z < childrenB.size(); z++) {
				freeSumB += Free[postA[i]->getID()][childrenB[z]->getID()];
				prevSubTreeSizeSumB += childrenB[z]->getSubTreeSize();
				subTreeSizeSumB.push_back(prevSubTreeSizeSumB);
			}


			prefix = freeSumB - Free[postA[i]->getID()][childrenB[0]->getID()] + LeftB[postA[i]->getID()][childrenB[0]->getID()] + postA[i]->getLeftmostForestNum() * (postB[j]->getSubTreeSize() - childrenB[0]->getSubTreeSize());
			left = 0;
			right = 0;
			sum = prefix + right;
			if(DEBUG) {
				ou << "If select Tree B " << to_string(childrenB[0]->getID()) << " #Subproblem: " << to_string(sum) << " Direction: Right";
			}
			if(min > sum) {
				freeS->setKeyNode(childrenB[0]->getID());
				freeS->setTreeToDecompose(1);
				freeS->setDirection(0);
			}
			for(int z = 1; z < childrenB.size() - 1; z++) {
				if(DEBUG) {
					ou << "If select Tree B " << to_string(childrenB[z]->getID()) << " #Subproblem: ";
				}

				prefix = freeSumB - Free[postA[i]->getID()][childrenB[z]->getID()] + AllB[postA[i]->getID()][childrenB[z]->getID()];
				left = postA[i]->getRightmostForestNum() * (1 + subTreeSizeSumB[z - 1]) + postA[i]->getSpecialForestNum() * (subTreeSizeSumB[childrenB.size() - 1] - subTreeSizeSumB[z]);
				right = postA[i]->getLeftmostForestNum() * (1 + subTreeSizeSumB[childrenB.size() - 1] - subTreeSizeSumB[z]) + postA[i]->getSpecialForestNum() * (subTreeSizeSumB[z - 1]);
				if(left > right) sum = prefix + right;
				else sum = prefix + left;
				if(min > sum) {
					min = sum;
					freeS->setKeyNode(childrenB[z]->getID());
					freeS->setTreeToDecompose(1);
					if(left > right) freeS->setDirection(0);
					else freeS->setDirection(1);
				}
				if(DEBUG) {
					ou << to_string(sum) << " Direction: ";
					if(left > right) ou << "Right" << endl;
					else ou << "Left" << endl;
				}
			}
			prefix = freeSumB - Free[postA[i]->getID()][childrenB[childrenB.size() - 1]->getID()] + RightB[postA[i]->getID()][childrenB[childrenB.size() - 1]->getID()] + postA[i]->getRightmostForestNum() * (postB[j]->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
			left = 0;
			right = 0;
			sum = prefix + left;
			if(DEBUG) {
				ou << "If select Tree B " << to_string(childrenB[childrenB.size() - 1]->getID()) << " #Subproblem: " << to_string(sum) << " Direction: Left";
			}
			if(min > sum) {
				min = sum;
				freeS->setKeyNode(childrenB[childrenB.size() - 1]->getID());
				freeS->setTreeToDecompose(1);
				freeS->setDirection(1);
			}
			}
			Free[postA[i]->getID()][postB[j]->getID()] = min;
			FreeStrategies[postA[i]->getID()][postB[j]->getID()] = freeS;
			if(DEBUG) {
				ou << "FreeStrategies[" << to_string(postA[i]->getID()) << "][" << to_string(postB[j]->getID()) << "] = " << freeS->toString() << endl;
			}
		}
	}
};*/

void TreeComparison::strategyComputation() {
	vector<Node*> preA = A_->getPreL();
	vector<Node*> preB = B_->getPreL();

	free(preA[0], preB[0]);
	if(DEBUG) {
		ou << "RESULT" << endl;
	/*	for(int i = 0; i < treeSizeA; i++) {
			for(int j = 0; j < treeSizeB; j++) {
				ou << Free[i][j] << " ";
			}
			ou << endl;
		}*/

		for(int i = 0; i < treeSizeA; i++) {
			for(int j = 0; j < treeSizeB; j++) {
				if(&FreeStrategies[i][j] != NULL) {
					ou << FreeStrategies[i][j].getLeaf() << " in ";
					if(FreeStrategies[i][j].getTreeToDecompose() == 0) ou << "A ";
					else ou << "B ";
				}
			}
			ou << endl;
		}
	}



};


int TreeComparison::free(Node* a, Node* b) {
	/*if(DEBUG) {
		ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
	}*/
/*	if(Free[a->getID()][b->getID()] != -1) {
		if(DEBUG) {
			ou << "Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") = " << to_string(Free[a->getID()][b->getID()]) << endl;
		}
		return Free[a->getID()][b->getID()];
	}*/
	
	/*if(DEBUG) {
		ou << a->toString() << endl;
		ou << b->toString() << endl;
	}*/
	if(Free[a->getID()][b->getID()] != -1) return Free[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	vector<Node*> childrenB = b->getChildren();
	Strategy freeS;
	int min = INT_MAX;
	int freeSumA = 0;
	vector<int> childrenSizeSumA;
	int freeSumB = 0;
	vector<int> childrenSizeSumB;

	if(childrenA.empty()) {
		int left = b->getRightmostForestNum();
		int right = b->getLeftmostForestNum();
		if(min > right) {
			min = right;
			freeS.setDirection(0);
			freeS.setTreeToDecompose(0);
			freeS.setKeyNode(a->getID());
			freeS.setLeaf(a->getID());
		} else if(min > left){
			min = left;
			freeS.setDirection(1);
			freeS.setTreeToDecompose(0);
			freeS.setKeyNode(a->getID());
			freeS.setLeaf(a->getID());
		}
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(a->getID()) << " in Tree A #Subproblem: " << to_string(min) << " Direction: ";
			if(left > right) ou << "Right" << endl;
			else ou << "Left" << endl;
		}
	}
	if(childrenB.empty()) {
		int left = a->getRightmostForestNum();
		int right = a->getLeftmostForestNum();
		if(min > right) {
			min = right;
			freeS.setDirection(0);
			freeS.setTreeToDecompose(1);
			freeS.setKeyNode(b->getID());
			freeS.setLeaf(b->getID());
		} else if(min > left){
			min = left;
			freeS.setDirection(1);
			freeS.setTreeToDecompose(1);
			freeS.setKeyNode(b->getID());
			freeS.setLeaf(b->getID());
		}

		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(b->getID()) << " in Tree B #Subproblem: " << to_string(min) << " Direction: ";
			if(left > right) ou << "Right" << endl;
			else ou << "Left" << endl;
		}
	} 

	int prev = 0;
	if(!childrenA.empty()) {
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
	}
	prev = 0;
	if(!childrenB.empty()) {
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
	}

	if(!childrenA.empty()) {
		int aleftmost = freeSumA - free(childrenA[0], b) + leftA(childrenA[0], b) + b->getLeftmostForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenA[0]->getID()) << "(leftmost) in Tree A #Subproblem: " << to_string(aleftmost) << " Direction: Right" << endl;
		}
		if(min > aleftmost) {
			min = aleftmost;
			freeS.setKeyNode(childrenA[0]->getID());
			freeS.setTreeToDecompose(0);
			freeS.setDirection(0);
			freeS.setLeaf(LeftAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
		}
		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getRightmostForestNum() * (1 + childrenSizeSumA[i - 1])  + b->getSpecialForestNum() * (childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]);
			int right = b->getLeftmostForestNum() * (1 + childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]) + b->getSpecialForestNum() * (childrenSizeSumA[i - 1]);
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				freeS.setKeyNode(childrenA[i]->getID());
				freeS.setTreeToDecompose(0);
				freeS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)freeS.setDirection(0);
				else freeS.setDirection(1);
			}
			if(DEBUG) {
				ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
				ou << "If select " << to_string(childrenA[i]->getID()) << " in Tree A #Subproblem: " << to_string(sum) << " Direction: ";
				if(left > right) ou << "Right" << endl;
				else ou << "Left" << endl;
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + rightA(childrenA[childrenA.size() - 1], b) + b->getRightmostForestNum() * (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenA[childrenA.size() - 1]->getID()) << "(rightmost) in Tree A #Subproblem: " << to_string(arightmost) << " Direction: Left" << endl;
		}
		if(min > arightmost) {
			min = arightmost;
			freeS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			freeS.setLeaf(RightAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			freeS.setTreeToDecompose(0);
			freeS.setDirection(1);
		}
	}

	if(!childrenB.empty()) {
		int bleftmost = freeSumB - free(a, childrenB[0]) + leftB(a, childrenB[0]) + a->getLeftmostForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenB[0]->getID()) << "(leftmost) in Tree B #Subproblem: " << to_string(bleftmost) << " Direction: Right" << endl;
		}
		if(min > bleftmost) {
			min = bleftmost;
			freeS.setKeyNode(childrenB[0]->getID());
			freeS.setLeaf(LeftBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			freeS.setTreeToDecompose(1);
			freeS.setDirection(0);
		}

		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getRightmostForestNum() * (1 + childrenSizeSumB[i - 1]) + a->getSpecialForestNum() * (childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]);
			int right = a->getLeftmostForestNum() * (1 + childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]) + a->getSpecialForestNum() * (childrenSizeSumB[i - 1]);
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				freeS.setKeyNode(childrenB[i]->getID());
				freeS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				freeS.setTreeToDecompose(1);
				if(left > right)freeS.setDirection(0);
				else freeS.setDirection(1);
			}
			if(DEBUG) {
				ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
				ou << "If select " << to_string(childrenB[i]->getID()) << " in Tree B #Subproblem: " << to_string(sum) << " Direction: ";
				if(left > right) ou << "Right" << endl;
				else ou << "Left" << endl;
			}
		}

		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + rightB(a, childrenB[childrenB.size() - 1]) + a->getRightmostForestNum() * (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		if(DEBUG) {
			ou << "Compute Free(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
			ou << "If select " << to_string(childrenB[childrenB.size() - 1]->getID()) << "(rightmost) in Tree B #Subproblem: " << to_string(brightmost) << " Direction: Left" << endl;
		}
		if(min > brightmost) {
			min = brightmost;
			freeS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			freeS.setLeaf(RightBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			freeS.setTreeToDecompose(1);
			freeS.setDirection(1);
		}
	}
	Free[a->getID()][b->getID()] = min;
	FreeStrategies[a->getID()][b->getID()] = freeS;
	if(DEBUG) {
		ou << "FreeS(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
		ou << FreeStrategies[a->getID()][b->getID()].toString() << endl;
	}
	return min;
};

int TreeComparison::leftA(Node* a, Node* b) {
	/*if(DEBUG) {
		ou << "Compute LeftA(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
	}*/
	/*if(LeftA[a->getID()][b->getID()] != -1) {
		if(DEBUG) {
			ou << "Compute LeftA(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") = " << to_string(LeftA[a->getID()][b->getID()]) << endl;
		}
		return LeftA[a->getID()][b->getID()];
	}*/
	if(LeftA[a->getID()][b->getID()] != -1) return LeftA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy leftAS;
	if(childrenA.empty()) {
		min = b->getLeftmostForestNum();
		leftAS.setKeyNode(a->getID());
		leftAS.setLeaf(a->getID());
		leftAS.setTreeToDecompose(0);
		leftAS.setDirection(0);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + leftA(childrenA[0], b) + b->getLeftmostForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(min >aleftmost) {
			min = aleftmost;
			leftAS.setKeyNode(childrenA[0]->getID());
			leftAS.setLeaf(LeftAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			leftAS.setTreeToDecompose(0);
			leftAS.setDirection(0);
		}

		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			int right = b->getLeftmostForestNum() * (1 + childrenSizeSumA[childrenSizeSumA.size() - 1] - childrenSizeSumA[i]) + b->getSpecialForestNum() * (childrenSizeSumA[i - 1]);
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				leftAS.setKeyNode(childrenA[i]->getID());
				leftAS.setTreeToDecompose(0);
				leftAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)leftAS.setDirection(0);
				else leftAS.setDirection(1);
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + allA(childrenA[childrenA.size() - 1], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		if(min > arightmost) {
			min = arightmost;
			leftAS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			leftAS.setLeaf(AllAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			leftAS.setTreeToDecompose(0);
			leftAS.setDirection(1);
		}
	}
	LeftA[a->getID()][b->getID()] = min;
	LeftAStrategies[a->getID()][b->getID()] = leftAS;
	return min;

};

int TreeComparison::rightA(Node* a, Node* b) {
	if(RightA[a->getID()][b->getID()] != -1) return RightA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy rightAS;
	if(childrenA.empty()) {
		min = b->getRightmostForestNum();
		rightAS.setKeyNode(a->getID());
		rightAS.setLeaf(a->getID());
		rightAS.setTreeToDecompose(0);
		rightAS.setDirection(1);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + allA(childrenA[0], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(min >aleftmost) {
			min = aleftmost;
			rightAS.setKeyNode(childrenA[0]->getID());
			rightAS.setLeaf(AllAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			rightAS.setTreeToDecompose(0);
			rightAS.setDirection(0);
		}

		for(int i = 1; i < childrenA.size() - 1; i++) {
			int prefix = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b);
			int left = b->getRightmostForestNum() * (1 + childrenSizeSumA[i - 1]) + b->getSpecialForestNum() * (childrenSizeSumA[childrenA.size() - 1] - childrenSizeSumA[i]);
			int right = b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				rightAS.setKeyNode(childrenA[i]->getID());
				rightAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				rightAS.setTreeToDecompose(0);
				if(left > right)rightAS.setDirection(0);
				else rightAS.setDirection(1);
			}
		}
		int arightmost = freeSumA - free(childrenA[childrenA.size() - 1], b) + rightA(childrenA[childrenA.size() - 1], b) + b->getLeftmostForestNum()* (a->getSubTreeSize() - childrenA[childrenA.size() - 1]->getSubTreeSize());
		if(min > arightmost) {
			min = arightmost;
			rightAS.setKeyNode(childrenA[childrenA.size() - 1]->getID());
			rightAS.setLeaf(AllAStrategies[childrenA[childrenA.size() - 1]->getID()][b->getID()].getLeaf());
			rightAS.setTreeToDecompose(0);
			rightAS.setDirection(1);
		}
	}
	RightA[a->getID()][b->getID()] = min;
	RightAStrategies[a->getID()][b->getID()] = rightAS;
	return min;
};

int TreeComparison::allA(Node* a, Node* b) {
	if(AllA[a->getID()][b->getID()] != -1) return AllA[a->getID()][b->getID()];
	vector<Node*> childrenA = a->getChildren();
	int min = INT_MAX;
	Strategy allAS;
	if(childrenA.empty()) {
		min = b->getSpecialForestNum();
		allAS.setKeyNode(a->getID());
		allAS.setLeaf(a->getID());
		allAS.setTreeToDecompose(0);
	} else {
		int freeSumA = 0;
		vector<int> childrenSizeSumA;
		int prev = 0;
		for(int i = 0; i < childrenA.size(); i++) {
			freeSumA += free(childrenA[i], b);
			prev += childrenA[i]->getSubTreeSize();
			childrenSizeSumA.push_back(prev);
		}
		int aleftmost = freeSumA - free(childrenA[0], b) + allA(childrenA[0], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[0]->getSubTreeSize());
		if(min >aleftmost) {
			min = aleftmost;
			allAS.setKeyNode(childrenA[0]->getID());
			allAS.setLeaf(AllAStrategies[childrenA[0]->getID()][b->getID()].getLeaf());
			allAS.setTreeToDecompose(0);
		}

		for(int i = 1; i < childrenA.size(); i++) {
			int sum = freeSumA - free(childrenA[i], b) + allA(childrenA[i], b) + b->getSpecialForestNum() * (a->getSubTreeSize() - childrenA[i]->getSubTreeSize());
			if(min > sum) {
				min = sum;
				allAS.setKeyNode(childrenA[i]->getID());
				allAS.setTreeToDecompose(0);
				allAS.setLeaf(AllAStrategies[childrenA[i]->getID()][b->getID()].getLeaf());
				if(left > right)allAS.setDirection(0);
				else allAS.setDirection(1);
			}
		}
	}
	AllA[a->getID()][b->getID()] = min;
	AllAStrategies[a->getID()][b->getID()] = allAS;
	return min;
};

int TreeComparison::leftB(Node* a, Node* b) {
	if(LeftB[a->getID()][b->getID()] != -1) return LeftB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy leftBS;
	if(childrenB.empty()) {
		min = a->getLeftmostForestNum();
		leftBS.setKeyNode(b->getID());
		leftBS.setLeaf(b->getID());
		leftBS.setTreeToDecompose(1);
		leftBS.setDirection(0);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + leftB(a, childrenB[0]) + a->getLeftmostForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(min >bleftmost) {
			min = bleftmost;
			leftBS.setKeyNode(childrenB[0]->getID());
			leftBS.setLeaf(LeftBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			leftBS.setTreeToDecompose(1);
			leftBS.setDirection(0);
		}

		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			int right = a->getLeftmostForestNum() * (1 + childrenSizeSumB[childrenSizeSumB.size() - 1] - childrenSizeSumB[i]) + a->getSpecialForestNum() * (childrenSizeSumB[i - 1]);
			int sum;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				leftBS.setKeyNode(childrenB[i]->getID());
				leftBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				leftBS.setTreeToDecompose(1);
				if(left > right)leftBS.setDirection(0);
				else leftBS.setDirection(1);
			}
		}
		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + allB(a, childrenB[childrenB.size() - 1]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		if(min > brightmost) {
			min = brightmost;
			leftBS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			leftBS.setLeaf(AllBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			leftBS.setTreeToDecompose(1);
			leftBS.setDirection(1);
		}
	}
	LeftB[a->getID()][b->getID()] = min;
	LeftBStrategies[a->getID()][b->getID()] = leftBS;
	return min;
};

int TreeComparison::rightB(Node* a, Node* b) {
	if(RightB[a->getID()][b->getID()] != -1) return RightB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy rightBS;
	if(childrenB.empty()) {
		min = a->getRightmostForestNum();
		rightBS.setKeyNode(b->getID());
		rightBS.setLeaf(b->getID());
		rightBS.setTreeToDecompose(1);
		rightBS.setDirection(1);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + allB(a, childrenB[0]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(min >bleftmost) {
			min = bleftmost;
			rightBS.setKeyNode(childrenB[0]->getID());
			rightBS.setLeaf(AllBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			rightBS.setTreeToDecompose(1);
			rightBS.setDirection(0);
		}

		for(int i = 1; i < childrenB.size() - 1; i++) {
			int prefix = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]);
			int left = a->getRightmostForestNum() * (1 + childrenSizeSumB[i - 1]) + a->getSpecialForestNum() * (childrenSizeSumB[childrenB.size() - 1] - childrenSizeSumB[i]);
			int right = a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			int sum = 0;
			if(left > right) sum = prefix + right;
			else sum = prefix + left;
			if(min > sum) {
				min = sum;
				rightBS.setKeyNode(childrenB[i]->getID());
				rightBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				rightBS.setTreeToDecompose(1);
				if(left > right)rightBS.setDirection(0);
				else rightBS.setDirection(1);
			}
		}
		int brightmost = freeSumB - free(a, childrenB[childrenB.size() - 1]) + rightB(a, childrenB[childrenB.size() - 1]) + a->getLeftmostForestNum()* (b->getSubTreeSize() - childrenB[childrenB.size() - 1]->getSubTreeSize());
		if(min > brightmost) {
			min = brightmost;
			rightBS.setKeyNode(childrenB[childrenB.size() - 1]->getID());
			rightBS.setLeaf(RightBStrategies[a->getID()][childrenB[childrenB.size() - 1]->getID()].getLeaf());
			rightBS.setTreeToDecompose(1);
			rightBS.setDirection(1);
		}
	}
	RightB[a->getID()][b->getID()] = min;
	RightBStrategies[a->getID()][b->getID()] = rightBS;
	return min;
};

int TreeComparison::allB(Node* a, Node* b) {
	if(AllB[a->getID()][b->getID()] != -1) return AllB[a->getID()][b->getID()];
	vector<Node*> childrenB = b->getChildren();
	int min = INT_MAX;
	Strategy allBS;
	if(childrenB.empty()) {
		min = a->getSpecialForestNum();
		allBS.setLeaf(b->getID());
		allBS.setKeyNode(b->getID());
		allBS.setTreeToDecompose(1);
	} else {
		int freeSumB = 0;
		vector<int> childrenSizeSumB;
		int prev = 0;
		for(int i = 0; i < childrenB.size(); i++) {
			freeSumB += free(a, childrenB[i]);
			prev += childrenB[i]->getSubTreeSize();
			childrenSizeSumB.push_back(prev);
		}
		int bleftmost = freeSumB - free(a, childrenB[0]) + allB(a, childrenB[0]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[0]->getSubTreeSize());
		if(min >bleftmost) {
			min = bleftmost;
			allBS.setKeyNode(childrenB[0]->getID());
			allBS.setLeaf(AllBStrategies[a->getID()][childrenB[0]->getID()].getLeaf());
			allBS.setTreeToDecompose(1);
		}

		for(int i = 1; i < childrenB.size(); i++) {
			int sum = freeSumB - free(a, childrenB[i]) + allB(a, childrenB[i]) + a->getSpecialForestNum() * (b->getSubTreeSize() - childrenB[i]->getSubTreeSize());
			if(min > sum) {
				min = sum;
				allBS.setKeyNode(childrenB[i]->getID());
				allBS.setLeaf(AllBStrategies[a->getID()][childrenB[i]->getID()].getLeaf());
				allBS.setTreeToDecompose(1);
				if(left > right)allBS.setDirection(0);
				else allBS.setDirection(1);
			}
		}
	}
	AllB[a->getID()][b->getID()] = min;
	AllBStrategies[a->getID()][b->getID()] = allBS;
	return min;
};
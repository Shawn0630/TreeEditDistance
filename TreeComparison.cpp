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

TreeComparison::TreeComparison(Tree* A, Tree* B, SimiMatrix costModel) {
	A_ = A;
	B_ = B;
	costModel_ = costModel;

	treeSizeA = A_->getTreeSize();
	treeSizeB = B_->getTreeSize();

	int maxSize = treeSizeA < treeSizeB? treeSizeB + 1 : treeSizeA + 1;
	fn = new int[maxSize + 1];
	ft = new int[maxSize + 1];

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
	delta = new float*[treeSizeA];
	hasVisited = new bool*[treeSizeA];
	d = new float*[treeSizeA];
	s = new float*[treeSizeA];
	t = new float*[treeSizeB];
	q = new float[treeSizeA];
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
		delta[i] = new float[treeSizeB];
		hasVisited[i] = new bool[treeSizeB];
		d[i] = new float[treeSizeB];
		s[i] = new float[treeSizeB];
	}

	for(int i = 0; i < treeSizeB; i++) {
		t[i] = new float[treeSizeB];
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
			delta[i][j] = 0.0f;
			hasVisited[i][j] = false;
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

void TreeComparison::deltaInit() {
	int treeSizeA = A_->getTreeSize();
	int treeSizeB = B_->getTreeSize();
	for(int i = 0; i < treeSizeA; i++) {
		Node* a = (*A_)[i];
		for(int j = 0; j < treeSizeB; j++) {
			Node* b = (*B_)[j];
			if(a->getSubTreeSize() == 1 && b->getSubTreeSize() == 1) {
				delta[i][j] = 0.0f;
			} else if(a->getSubTreeSize() == 1) {
				delta[i][j] = B_->preL_to_sumInsCost[j] - costModel_.ins(b->getLabel());
			} else if(b->getSubTreeSize() == 1) {
				delta[i][j] = A_->preL_to_sumDelCost[i] - costModel_.del(a->getLabel());
			}
		}
	}
};


void TreeComparison::computeSumInsAndDelCost(Tree* tree) {
	int treeSize = tree->getTreeSize();
	for(int i = 0; i < treeSize; i++) {
		int nodeForSum = treeSize - i - 1;// postOrder from bottom to up
		Node* node = (*tree)[nodeForSum];
		Node* parent = node->getParent();
	/*	cout << "Before" << endl;
		cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
		cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;
		cout << "ins " << node->getLabel() << " " << costModel_.ins(node->getLabel()) << endl;
		cout << "del " << node->getLabel() << " " << costModel_.del(node->getLabel()) << endl;
		cout << "After" << endl;*/
		tree->preL_to_sumInsCost[nodeForSum] += costModel_.ins(node->getLabel());
		tree->preL_to_sumDelCost[nodeForSum] += costModel_.del(node->getLabel());
/*		cout << "preL_to_sumInsCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumInsCost[nodeForSum]) << endl;
		cout << "preL_to_sumDelCost[" << to_string(nodeForSum) << "] = " << to_string(tree->preL_to_sumDelCost[nodeForSum]) << endl;*/
		if(parent != NULL) {
			/*cout << "Update Parent Before" << endl;
			cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
			cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;*/
			tree->preL_to_sumInsCost[parent->getID()] += tree->preL_to_sumInsCost[node->getID()];
			tree->preL_to_sumDelCost[parent->getID()] += tree->preL_to_sumDelCost[node->getID()];
/*			cout << "Update Parent After" << endl;
			cout << "preL_to_sumInsCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumInsCost[parent->getID()]) << endl;
			cout << "preL_to_sumDelCost[" << to_string(parent->getID()) << "] = " << to_string(tree->preL_to_sumDelCost[parent->getID()]) << endl;*/
		}
	}

};

float TreeComparison::gted(Node* a, Node* b) {
	
	if(hasVisited[a->getID()][b->getID()] == true) return 0.0f;
	hasVisited[a->getID()][b->getID()] = true;
	int treeSizeA = a->getSubTreeSize();
	int treeSizeB = b->getSubTreeSize();
	if(DEBUG) {
		ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
		ou << "treeSizeA = " << to_string(treeSizeA) << endl;
		ou << "treeSizeB = " << to_string(treeSizeB) << endl;
	}
	
	if ((treeSizeA == 1 || treeSizeB == 1)) {
      //return spf1(a, treeSizeA, b, treeSizeB);
		if(DEBUG) {
			ou << "return 0.0f" << endl;
		}
		return 0.0f;
    }


	int pathLeaf = FreeStrategies[a->getID()][b->getID()].getLeaf();
	int treeToDecompose = FreeStrategies[a->getID()][b->getID()].getTreeToDecompose();
	Node* currentPathNode = treeToDecompose == 0? (*A_)[pathLeaf] : (*B_)[pathLeaf];


	if(treeToDecompose == 0) { // decompose tree A
		Node* parent = currentPathNode->getParent();
		int pathType = getPathType(A_, a, pathLeaf);// 0 left 1 right 2 special
		while(parent != NULL && parent->getID() >= a->getID()) {
        	vector<Node*> children = parent->getChildren();
        	for(int i = 0; i < children.size(); i++) {
          		Node* child = children[i];
          		if(child->getID() != currentPathNode->getID()) {
          			if(DEBUG) {
          				ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
          				ou << "create problem in A " << "gted(" << to_string(child->getID()) << ", " << to_string(b->getID()) << ")" << endl;
          			}
            		gted(child, b);
          		}
        	}
        	parent = currentPathNode->getParent();
        	currentPathNode = parent;
        }

      	if (pathType == 0) {
        	//return spfL(a, b, false);
        	return 0.0f;
      	}
      	else if (pathType == 1) {
        	//return spfR(a, b, false);
        	return 0.0f;
      	}
      	return spfA(a, b, pathLeaf, pathType, false);
	} 

	else if(treeToDecompose == 1) {
		Node* parent = currentPathNode->getParent();
		int pathType = getPathType(B_, b, pathLeaf);
		while(parent != NULL && parent->getID() >= a->getID()) {
			vector<Node*> children = parent->getChildren();
			for(int i = 0; i < children.size(); i++) {
				Node* child = children[i];
				if(child->getID() != currentPathNode->getID()) {
					if(DEBUG) {
						ou << "gted(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ") ";
          				ou << "create problem in B " << "gted(" << to_string(a->getID()) << ", " << to_string(child->getID()) << ")" << endl;
          			}
					gted(a, child);
				}
			}
			parent = currentPathNode->getParent();
        	currentPathNode = parent;
		}

		if(pathType == 0) {
			//return spfL(b, a, true);
			return 0.0f;
		}
		else if(pathType == 1) {
			//return spfR(b, a, true);
			return 0.0f;
		}
		return spfA(b, a, pathLeaf, pathType, true);
	}
};

float TreeComparison::spfA(Node* a, Node* b, int leaf, int pathType, bool swap) {
	int endF = a->getID();
	int enF_in_preR = A_->preL_to_preR[endF];
	int endG = b->getID();
	int endG_in_preR = B_->preL_to_preR[endG];
	int sizeF = a->getSubTreeSize();
	int sizeG = b->getSubTreeSize();
	int endPathNode = leaf;
	int startPathNode = -1;
	int lFFirst, lFLast, lF;
	int rFFirst, rFLast, rF;
	int lGFirst, lGLast;
	int rGFirst, rGLast;

	if(DEBUG) {
		ou << "spfA(" << to_string(a->getID()) << ", " << to_string(b->getID()) << ")" << endl;
	}


	//loop A
	while(endPathNode >= endF) {
		int endPathNode_in_preR = A_->preL_to_preR[endPathNode];
		int startPathNode_in_preR = startPathNode == -1? 0x7fffffff : A_->preL_to_preR[startPathNode];

		int parent_of_endPathNode_preL = (*A_)[endPathNode]->getParent() == NULL? 0x7fffffff : (*A_)[endPathNode]->getParent()->getID();
		int parent_of_endPathNode_preR = (*A_)[endPathNode]->getParent() == NULL? 0x7fffffff : A_->preL_to_preR[parent_of_endPathNode_preL];

		bool hasLeftPart;
		bool hasRightPart;

		int lF_prev = endPathNode;

		if(startPathNode - endPathNode > 1) {
			hasLeftPart = true;
		} else {
			hasLeftPart = false;
		}
		if(startPathNode >= 0 && startPathNode_in_preR - endPathNode_in_preR > 1) {
			hasRightPart = true;
		} else {
			hasRightPart = false;
		}

		if(pathType == 1 || pathType == 2 && hasLeftPart) {
			// lFFirst and LFLast is important in this condition.
			// rGFirst and rGLast is to set in this stage.
			if(startPathNode == -1) {
				lFFirst = endPathNode;//the first node is the node on the path
				rFFirst = endPathNode_in_preR;// the first node is the node on the path
			} else {
				lFFirst = startPathNode - 1;//the first node is the node to the left of the path
				rFFirst = startPathNode_in_preR;//rFFirst set to the node on the path
			}
			if(!hasRightPart) {
				rFLast = endPathNode_in_preR;
			}
			lFLast = hasRightPart? endPathNode + 1 : endPathNode;

			rGLast = B_->preL_to_preR[endG];
			rGFirst = (rGLast + sizeG) - 1; // get the leftmost child in G

			fn[sizeG] = -1;
        	for (int i = endG; i < endG + sizeG; i++) {
            	fn[i] = -1;
            	ft[i] = -1;
        	}

			//loop B
			for(int rG = rGFirst; rG >= rGLast; rG--) {
				if(DEBUG) {
					ou << "new Round B" << endl;
				}
				int rG_in_preL = (B_)->preR_to_preL[rG];
				Node* parent = (*B_)[rG_in_preL]->getParent();
				int parent_of_rG_in_preL = parent == NULL? 0x7fffffff : parent->getID();
				int parent_of_rG_in_preR = parent == NULL? 0x7fffffff : B_->preL_to_preR[parent_of_rG_in_preL];
				lGFirst = B_->preR_to_preL[rG];// lGFirst is set to rGFirst;
				
				int rGminus1_in_preL = rG <= endG_in_preR? 0x7fffffff : B_->preR_to_preL[rG - 1];// rG should greater than endG_in_preR cause rG is the inner node of subtree enG
				int rGminus1_in_preR = rG <= endG_in_preR? 0x7fffffff : rG - 1;
				
				if (pathType == 1){
            		if (lGFirst == endG || rGminus1_in_preL != parent_of_rG_in_preL) {// parent not exist or not the rightmost child
              			lGLast = lGFirst;//lGLast is set to lGFirst
            		} else {
              			lGLast = parent_of_rG_in_preL + 1;//lGLast is set to the leftmost child of rG's parent
            		}
          		} else {
            		lGLast = lGFirst == endG? lGFirst : endG + 1;//lGLast is set to the leftmost child of the whole tree
          		}
			
          		updateFnArray(B_->preL_to_ln[lGFirst], lGFirst, endG); //stores the counter in D loop fn[ln] stores the start point
         		updateFtArray(B_->preL_to_ln[lGFirst], lGFirst); 
				//loop C
				for(int lF = lFFirst; lF >= lFLast; lF--) {
					if(DEBUG) {
						ou << "new round C" << endl;
					}
					int lG = lGFirst;
					int lF_in_preR = (A_)->preL_to_preR[lF];


					if(DEBUG) {
          				ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ")" << endl;
          				ou << "Save to S[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
					}

					int GcurrentForestSize = (*B_)[lG]->getSubTreeSize();
           		 	int GcurrentForestCost = (swap ? (B_)->preL_to_sumDelCost[lG] : (B_)->preL_to_sumInsCost[lG]); 


					bool fForestIsTree = lF_in_preR == rF;
					int lFSubtreeSize = (*A_)[lF]->getSubTreeSize();
					bool lFIsLeftSiblingOfCurrentPathNode = lF + lFSubtreeSize == startPathNode;
					bool lFIsConsecutiveNodeOfCurrentPathNode = startPathNode - lF == 1;

					int case1SLeftIndex, case1SRightIndex;//S[lF + 1, lG];
          			int case1TLeftIndex, case1TRightIndex;//T[lG, rG];
          			
          			int case2SLeftIndex, case2SRightIndex;//S[lF, lG];

          			int case3SLeftIndex, case3SRightIndex;
          			int case3TLeftIndex, case3TRightIndex;

          			case1SLeftIndex =  lF + 1;//fixed
          			case2SLeftIndex = lF;//fixed

          			case1TRightIndex = rG;//fixed

          			case3TRightIndex = rG;//fixed

          			float case1 = 0, case2 = 0, case3 = 0;
          			int case1_case, case2_case, case3_case;

          			float minCost = 0;

          			case1_case = 1;
          			case3_case = 1;

          			if (fForestIsTree) { // F_{lF,rF} is a tree.
              			if (lFSubtreeSize == 1) { // F_{lF,rF} is a single node.
              				// F_{lF, rF} - lF = null
                			case1_case = 3;
              			} else if (lFIsConsecutiveNodeOfCurrentPathNode) { // F_{lF,rF}-lF is the path node subtree.
              				// F_{lF, rF} - rF = tree
                			case1_case = 2;
              			}
              			case3 = 0;//F_{lF, rF} - F(lF) = null
              			if(DEBUG) {
          					ou << "case3 = 0" << endl; 
          				}
              			case3_case = 2;
              		} else {
              			if (lFIsConsecutiveNodeOfCurrentPathNode) {
                			// F_{lF, rF} - lF = tree
                			case1_case = 2;
              			}
              			//case3 = FcurrentForestCost - (swap ? (A_)->preL_to_sumInsCost[lF] : (A_)->preL_to_sumDelCost[lF]); // USE COST MODEL - Delete F_{lF,rF}-F_lF.
              			if(DEBUG) {
              				ou << "case3_case FcurrentForest - F(lF)" << endl;
              			}
              			
              			if (lFIsLeftSiblingOfCurrentPathNode) {
               				 case3_case = 3;
              			}
              		}

              		if(case3_case == 1) {
              			case3SLeftIndex = lF + lFSubtreeSize;
              		}

              		switch(case1_case) {
             			case 1: 
             				case1SRightIndex = lG;
             				//case1 = s[case1SLeftIndex][case1SRightIndex]; 
             				if(DEBUG) {
             					ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "]" << endl;
             				}
             				break;
              			case 2: 
              				case1TLeftIndex = lG;
              				//case1 = t[case1TLeftIndex][case1TRightIndex]; 
              				if(DEBUG) {
              					ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "]" << endl; 
              				}
              				break;
              			case 3: 
              				//case1 = GcurrentForestCost; 
              				if(DEBUG) {
              					ou << "case1_case3 " << to_string(GcurrentForestCost) << endl; 
              				}
              				break; // USE COST MODEL - Insert G_{lG,rG}.
            		}

            		case1 += (swap ? costModel_.ins((*A_)[lF]->getLabel()) : costModel_.del((*A_)[lF]->getLabel()));
            		//minCost = case1;

            		if (GcurrentForestSize == 1) { // G_{lG,rG} is a single node.
              			//case2 = FcurrentForestCost; // USE COST MODEL - Delete F_{lF,rF}.
              			if(DEBUG) {
              				ou << "case2_case1 FcurrentForestCost" << endl;
              			}
            		} else { // G_{lG,rG} is a tree.
              			//case2 = q[lF];
            			if(DEBUG) {
            				ou << "case2_case2 q[" << to_string(lF) << "]" << endl;
            			}
            		}

            		if (minCost < case3) {
            			//case3 += swap? d[lG][lF] : d[lF][lG];
            			if(DEBUG) {
            				if(swap) ou << "case3_case3 D[" << to_string(lG) << ", " << to_string(lF) << "]" << endl;
            				else ou << "case3_case3 D[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
            			}
            			if(minCost < case3) {
            				case3 += swap? costModel_.ren((*B_)[lG]->getLabel(), (*A_)[lF]->getLabel()) : costModel_.ren((*A_)[lF]->getLabel(), (*B_)[lG]->getLabel());
            			} 
            			if(minCost < case3) {
            				minCost = case3;
            			}
            		}
					lG = ft[lG];
					//loop D
					while (lG >= lGLast) {
						if(DEBUG) {
							ou << "Left (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ")" << endl;
							ou << "Save to S[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
						}

						GcurrentForestSize++;
						GcurrentForestCost += (swap ? costModel_.del((*B_)[lG]->getLabel()) : costModel_.ins((*B_)[lG]->getLabel()));
						minCost = 0;

						switch(case1_case) {
                			case 1:
                				case1SRightIndex = lG;
                				//case1 = s[case1SLeftIndex][case1SRightIndex] + (swap? costModel_.ins((*A_)[lF]->getLabel()) : costModel_.del((*A_)[lF]->getLabel())); 
                				if(DEBUG) {
                					ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "]" << endl;
                				}
                				break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
                			case 2: 
                				case1TLeftIndex = lG;
                				if(DEBUG) {
                					ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "]" << endl;
                				}
                				//case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*A_)[lF]->getLabel()) : costModel_.del((*A_)[lF]->getLabel())); 
                				break; // USE COST MODEL - Delete lF, leftmost root node in F_{lF,rF}.
                			case 3: 
                				//case1 = GcurrentForestCost + (swap? costModel_ins((*A_)[lF]->getLabel()) : costModel_.del((*A_)[lF]->getLabel()));
                				if(DEBUG) {
                					ou << "case1_case3 " << to_string(GcurrentForestCost) << endl;
                				}
                				break; // USE COST MODEL - Insert G_{lG,rG} and elete lF, leftmost root node in F_{lF,rF}.
              			}
              			//minCost = case1;

              			case2SRightIndex = fn[lG];
              			//case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*B_)[lG]->getLabel()) : costModel_.ins((*B_)[lG]->getLabel()));
              			if(DEBUG) {
              				ou << "case2 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "]" << endl;
              			}
              			if(case2 < minCost) {
              				//minCost = case2;
              			} 

              			//case3 = swap ? d[lG][lF] : d[lF][lG];
              			if(DEBUG) {
              				if(swap) {
              					ou << "case3 d[" << to_string(lG) << ", " << to_string(lF) << "]" << endl;
              				} else {
              					ou << "case3 d[" << to_string(lF) << ", " << to_string(lG) << "]" << endl;
              				}
              			}
              			//if (minCost < case3) {
                			switch(case3_case) {
                    			case 1: 
                    				case3SRightIndex = fn[lG] + (*B_)[lG]->getSubTreeSize();
                    				if(DEBUG) {
                    					ou << "case3 s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "]" << endl;
                    				}
                    				//case3 += s[case3SLeftIndex][case3SRightIndex]; 
                    				break;
                    			case 2: 
                    				//case3 += GcurrentForestCost - (swap ? (B_)->preL_to_sumDelCost[lG] : (B_)->preL_to_sumInsCost[lG]); 
                    				if(DEBUG) {
                    					ou << "case3 " << "GcurrentForestCost - G(lG)" << endl;
                    				}
                    				break; // USE COST MODEL - Insert G_{lG,rG}-G_lG.
                    			case 3: 
                    				case3TLeftIndex = fn[lG + (*B_)[lG]->getSubTreeSize() - 1];
                					if(DEBUG) {
                						ou << "case3 t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "]" << endl;
                					}
                    				//case3 += t[case3TLeftIndex][case3TRightIndex]; 
                    				break;
                			}
                			if (case3 < minCost) {
                  				case3 += (swap ? costModel_.ren((*B_)[lG]->getLabel(), (*A_)[lF]->getLabel()) : costModel_.ren((*A_)[lF]->getLabel(), (*B_)[lG]->getLabel())); // USE COST MODEL - Rename the leftmost root nodes in F_{lF,rF} and G_{lG,rG}.
                  				if (case3 < minCost) {
                    				minCost = case3;
                  				}
                			}
              			//}
						lG = ft[lG];
					}
					lF_prev = lF;
				}
				if(DEBUG) {
					ou << "rGminus1_in_preR = " << to_string(rGminus1_in_preR) << " rG = " << to_string(rG) << " parent_of_rG_in_preL = " << to_string(parent_of_rG_in_preL) << " parent_of_rG_in_preR = " << to_string(parent_of_rG_in_preR) << endl;
				}
				if(rGminus1_in_preR == parent_of_rG_in_preR && rGminus1_in_preR != 0x7fffffff) {
					if (!hasRightPart) {
              			if (hasLeftPart) {
              				if(swap) {
              					if(DEBUG) {
              						ou << "D[" << to_string(parent_of_rG_in_preL) << ", " << to_string(endPathNode) << "] = " << "S[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
              					}
              				} else {
              					if(DEBUG) {
              						ou << "D[" << to_string(endPathNode) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "S[" << to_string(lFLast + 1) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
              					}
              				}
              			}

              			if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {
                			if (swap) {
                				if(DEBUG) {
                					ou << "D[" << to_string(parent_of_rG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "S[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                				}
                			} else {
                				if(DEBUG) {
                					ou << "D[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_rG_in_preL) << "] = " << "S[" << to_string(lFLast) << ", " << to_string(rGminus1_in_preL + 1) << "]" << endl;
                				}
                			}
              			}
              		}

              		for (int lF = lFFirst; lF >= lFLast; lF--) {
              			if(DEBUG) {
              				ou << "Q[" << to_string(lF) << "] = " << "S[" << to_string(lF) << ", " << to_string(parent_of_rG_in_preL + 1) << "]" << endl;
              			}
            		}
				}

				for (int lG = lGFirst; lG >= lGLast; lG = ft[lG]) {
					if(DEBUG) {
						ou << "T[" << to_string(lG) << ", " << to_string(rG) << "] = " << "S[" << to_string(lFLast) << ", " << to_string(lG) << "]" << endl;
					}
          		}
			}
		}


		if (pathType == 0 || pathType == 2 && hasRightPart || pathType == 2 && !hasLeftPart && !hasRightPart) {
			if (startPathNode == -1) {
          		lFFirst = endPathNode;
          		rFFirst = A_->preL_to_preR[endPathNode];
        	} else {
          		rFFirst = A_->preL_to_preR[startPathNode] - 1;//the node right to the node on the path
          		lFFirst = endPathNode + 1;//lFirst is set to the node on the path
        	}

        	lFLast = endPathNode;
        	rFLast = A_->preL_to_preR[endPathNode];

        	lGLast = B_->preL_to_preR[endG];
        	lGFirst = (lGLast + sizeG) - 1;

        	fn[sizeG] = -1;
        	for (int i = endG; i < endG + sizeG + 2; i++) {
            	fn[i] = -1;
            	ft[i] = -1;
        	}

        	//loop B'
        	for (int lG = lGFirst; lG >= lGLast; lG--) {
        		if(DEBUG) {
					ou << "new Round B'" << endl;
				}
        		Node* parent = (*B_)[lG]->getParent();
        		int parent_of_lG_in_preL = parent == NULL? 0x7fffffff: parent->getID();
        		int parent_of_lG_in_preR = parent == NULL? 0x7fffffff : B_->preL_to_preR[parent->getID()];// not exist -1;
				rGFirst = B_->preL_to_preR[lG];
				int lG_in_preR = B_->preL_to_preR[lG];

				int lGminus1_in_preL = lG <= endG? 0x7fffffff : lG - 1;
				int lGminus1_in_preR = lG <= endG? 0x7fffffff : B_->preL_to_preR[lG - 1];

				 if (pathType == 0) {
           	 		if (lG == endG || lGminus1_in_preL != parent_of_lG_in_preL) {//parent not exists or not the leftmost child.
              			rGLast = rGFirst;
            		} else {
              			rGLast = parent_of_lG_in_preR + 1;
            		}
          		} else {// left and right
            		rGLast = rGFirst == endG_in_preR ? rGFirst : endG_in_preR;
          		}


				if(DEBUG) {
					ou << "updateFnArray(" << to_string(B_->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ", " << to_string(endG_in_preR) << ")" << endl;
				}
				updateFnArray(B_->preR_to_ln[rGFirst], rGFirst, endG_in_preR);
			
				if(DEBUG) {
					ou << "updateFtArray(" << to_string(B_->preR_to_ln[rGFirst]) << ", " << to_string(rGFirst) << ")" << endl;
				}
          		updateFtArray(B_->preR_to_ln[rGFirst], rGFirst);
          	
          		if(DEBUG) {
          			ou << "endG_in_preR = " << to_string(endG_in_preR) << endl;
          			ou << "start from rG = " << to_string(rGFirst) << endl;
          			ou << "FN" << endl;
          			for(int i = endG; i < endG + sizeG + 1; i++) {
          				ou << fn[i] << " ";
          			}
          			ou << endl;
          			ou << "FT" << endl;
          			for(int i = endG; i < endG + sizeG + 1; i++) {
          				ou << ft[i] << " ";
          			}
          			ou << endl;
          		}
          		lF = lF_prev;
          		// loop C'
          		for(int rF = rFFirst; rF >= rFLast; rF--) {
          			if(DEBUG) {
          				ou << "new Round C'" << endl;
          			}
          			int rG = rGFirst;
          			int rG_in_preL = (B_)->preR_to_preL[rG];
          			
          			if(rF == rFLast) lF = A_->preR_to_preL[rFLast]; 


          			if(DEBUG) {
          				ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ")" << endl;

          				ou << "Save to S[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
					}
          			

          			int rF_in_preL = (A_)->preR_to_preL[rF];
          			
          			bool FForestIsTree = lF == rF_in_preL;
          			int rFSubtreeSize = (*A_)[rF_in_preL]->getSubTreeSize();
          			
          			
          			int case1SLeftIndex, case1SRightIndex;//S[rF + 1, rG];
          			int case1TLeftIndex, case1TRightIndex;//T[lG, rG];

          			int case2SLeftIndex, case2SRightIndex;//S[rF, rG];

          			int case3SLeftIndex, case3SRightIndex;
          			int case3TLeftIndex, case3TRightIndex;
        

          			float case1 = 0, case2 = 0, case3 = 0;
          			int case1_case, case2_case, case3_case;

          			int GcurrentForestSize = (*B_)[lG]->getSubTreeSize();
          			float GcurrentForestCost = (swap ? (B_)->preL_to_sumDelCost[lG] : (B_)->preL_to_sumInsCost[lG]); // USE COST MODEL - reset to subtree insertion cost.

          			float minCost = 0;


          			case1SLeftIndex = rF + 1;//fixed
          			case1SRightIndex = rG;

          			case1TLeftIndex = lG;//fixed
          			case1TRightIndex = rG;

          			case2SLeftIndex = rF;//fixed


          			case3TLeftIndex = lG;//fixed


          			bool rFIsConsecutiveNodeOfCurrentPathNode;
          			bool rFIsRightSiblingOfCurrentPathNode;

          			case1_case = 1;
          			case3_case = 1;//otherwise

          			if (startPathNode > 0) {
              			rFIsConsecutiveNodeOfCurrentPathNode = startPathNode_in_preR - rF == 1;
              			rFIsRightSiblingOfCurrentPathNode = rF + rFSubtreeSize == startPathNode_in_preR;
            		} else {
              			rFIsConsecutiveNodeOfCurrentPathNode = false;//consecutiveNode use T;
              			rFIsRightSiblingOfCurrentPathNode = false;//
            		}

          			if(FForestIsTree) {
          				if(rFSubtreeSize == 1) {
          					// F_{lF, rF} - rF = null
          					//case1 = GcurrentForestCost;//sumG
          					case1_case = 3;
          				} else if(rFIsConsecutiveNodeOfCurrentPathNode) {
          					// F_{lF, rF} - rF = tree
          					//case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG];
          					case1_case = 2;
          				}
          				case3 = 0;
          				if(DEBUG) {
          					ou << "case3 = 0" << endl; 
          				}
          				case3_case = 2;// F_{lF, rF} - F(rF) = null
          			} else {
          				if (rFIsConsecutiveNodeOfCurrentPathNode) {// F_{lF, rF} - rF = the subforest to the left of the path
          					//case1 = t[case1TLeftIndex][case1TRightIndex];//T[lG, rG]
          					case1_case = 2;
          				} else {//otherwise
          					//case1 = s[case1SLeftIndex][case1SRightIndex];//S[rF + 1, rG];// have calculate
          					case1_case = 1;
          				}
          				//case3 = FcurrentForestCost - (swap ? (A_)->preL_to_sumInsCost[rF_in_preL] : (A_)->preL_to_sumDelCost[rF_in_preL]);// the first case in G should be G_{lG, rG} - l(rG) = null // F_{lF, rF} - F(rF), G_{lG, rG} - G(rG)
          				if(DEBUG) {
              				ou << "case3_case FcurrentForest - F(rF)" << endl;
              			}
          				if (rFIsRightSiblingOfCurrentPathNode) {
          					case3_case = 3; // use T
          				}
          			}

          			if (case3_case == 1) {
              			case3SLeftIndex = rF + rFSubtreeSize;//delete the whole rightmost tree//otherwise
            		}

            		switch(case1_case) {
            			case 1:
            				//case1 = s[case1SLeftIndex][case1SRightIndex];
            				if(DEBUG) {
            					ou << "case1_case1 = s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "]" << endl; 
            				}
            				break;
            			case 2:
            				//case1 = t[case1TLeftIndex][case1TRightIndex];
            				if(DEBUG) {
            					ou << "case1_case2 = t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "]" << endl;
            				}
            				break;
            			case 3:
            				//case1 = GcurrentForestCost;
            				if(DEBUG) {
            					ou << "case1_case3 = " << GcurrentForestCost << endl;
            				}
            				break;
            				
            		}
            		//case1 += (swap ? costModel_.ins((*A_)[rF]->getLabel()) : costModel_.del((*A_)[rF]->getLabel()));
            		minCost = case1;


          			if (GcurrentForestSize == 1) {// the first case in G should be a node or a tree
              			//case2 = FcurrentForestCost;
              			if(DEBUG) {
              				ou << "case2_case1 = FcurrentForestCost" << endl;
              			}
            		} else {
              			//case2 = q[rF];
              			if(DEBUG) {
              				ou << "case2_case2 = q[" << to_string(rF) << "]" << endl;
              			}
            		}

            		//case2 += (swap ? costModel_.del((*B_)[rG_in_preL]->getLabel()) : costModel_.ins((*B_)[rG_in_preL]->getLabel()));
            		if(case2 < minCost) {
            			minCost = case2;
            		}

            		//if(minCost < case3) { 
            			//case3 += swap ? d[rG_in_preL][rF_in_preL] : d[rF_in_preL][rG_in_preL];// F(rF) - rF
            			if(DEBUG) {
            				if(swap) ou << "case3_case3 D[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "]" << endl;
            				else ou << "case3_case3 D[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "]" << endl;
            			}
            			if(minCost < case3) {
            				 case3 += (swap ? costModel_.ren((*B_)[rG_in_preL]->getLabel(), (*A_)[rF_in_preL]->getLabel()) : costModel_.ren((*A_)[rF_in_preL]->getLabel(), (*B_)[rG_in_preL]->getLabel()));
            			}
            			if(case3 < minCost) {
            				minCost = case3;
            			}
            		//}

            		s[rF][rG] = minCost;

          			rG = ft[rG];
          			
          			// loop D'
          			while(rG >= rGLast) {// every G is a subforest not a subtree
          				rG_in_preL = (B_)->preR_to_preL[rG];

          				if(DEBUG) {
          					ou << "Right (" << to_string(lF) << ", " << to_string(rF) << ", " << to_string(lG) << ", " << to_string(rG) << ")" << endl;
          					ou << "Save to S[" << to_string(rF) << ", " << to_string(rG) << "]" << endl;
						}

						GcurrentForestSize++;
						GcurrentForestCost += (swap ? costModel_.del((*B_)[rG_in_preL]->getLabel()) : costModel_.ins((*B_)[rG_in_preL]->getLabel()));

						

						switch (case1_case) {
                			case 1:
                			//case1 = s[case1SLeftIndex][case1SRightIndex] + (swap ? costModel_.ins((*A_)[rF]->getLabel()) : costModel_.del((*A_)[rF]->getLabel())); 
                			case1SRightIndex = rG;
                			if(DEBUG) {
                				ou << "case1_case1 s[" << to_string(case1SLeftIndex) << ", " << to_string(case1SRightIndex) << "]" << endl; 
                			} 
                			break; 
                			case 2: 
                			case1TRightIndex = rG;
                			//case1 = t[case1TLeftIndex][case1TRightIndex] + (swap ? costModel_.ins((*A_)[rF]->getLabel()) : costModel_.del((*A_)[rF]->getLabel())); 
                			if(DEBUG) {
                				ou << "case1_case2 t[" << to_string(case1TLeftIndex) << ", " << to_string(case1TRightIndex) << "]" << endl;
                			}
                			break; 
                			case 3: 
                			//case1 = GcurrentForestCost + (swap ? costModel_.ins((*A_)[rF]->getLabel()) : costModel_.del((*A_)[rF]->getLabel())); 
                			if(DEBUG) {
                				ou << "case1_case3 " << GcurrentForestCost << endl;
                			}
                			break; 
              			}
              			minCost = case1;

              			case2SRightIndex = fn[rG];
              			if(DEBUG) {
              				ou << "case2_case3 s[" << to_string(case2SLeftIndex) << ", " << to_string(case2SRightIndex) << "]" << endl;
              			}
              			case2 = s[case2SLeftIndex][case2SRightIndex] + (swap ? costModel_.del((*B_)[rG_in_preL]->getLabel()) : costModel_.ins((*B_)[rG_in_preL]->getLabel()));//G is not a tree or a node for sure in D loop
              			if(case2 < minCost) {
              				minCost = case2;
              			}

              			//case3 = swap ? d[rG_in_preL][rF_in_preL] : d[rF_in_preL][rG_in_preL];//F_{rF} - rF, G_{rG} - rG
              			if(DEBUG) {
              				if(swap) {
              					ou << "case3_case d[" << to_string(rG_in_preL) << ", " << to_string(rF_in_preL) << "]" << endl;
              				} else {
              					ou << "case3_case d[" << to_string(rF_in_preL) << ", " << to_string(rG_in_preL) << "]" << endl;
              				}
              			}
              			//if(case3 < minCost) {
              				switch(case3_case) {
              					case 1: 
              					case3SRightIndex = fn[(rG + (*B_)[rG_in_preL]->getSubTreeSize()) - 1];
              					//case3 += s[case3SLeftIndex][case3SRightIndex];
              					if(DEBUG) {
              						ou << "case3_case1 s[" << to_string(case3SLeftIndex) << ", " << to_string(case3SRightIndex) << "]" << endl; 
              					}
              					break;
              					case 2: 
              					//case3 += GcurrentForestCost - (swap ? (B_)->preL_to_sumDelCost[rG_in_preL] : (B_)-?preL_to_sumInsCost[rG_in_preL]);
              					if(DEBUG) {
              						ou << "case3_case2 " << "GcurrentForestCost - G(rG) " << endl; 
              					}
              					break;
              					case 3: 
              					case3TRightIndex = fn[(rG + (*B_)[rG_in_preL]->getSubTreeSize()) - 1];
              					//case3 += t[case3TLeftIndex][case3TRightIndex];
              					if(DEBUG) {
              						ou << "case3_case3 t[" << to_string(case3TLeftIndex) << ", " << to_string(case3TRightIndex) << "]" << endl; 
              					}
              					break;
              				}
              				if(minCost < case3) {
              					case3 += (swap ? costModel_.ren((*B_)[rG_in_preL]->getLabel(), (*A_)[rF_in_preL]->getLabel()) : costModel_.ren((*A_)[rF_in_preL]->getLabel(), (*B_)[rG_in_preL]->getLabel()));
              					if(case3 < minCost) {
              						minCost = case3;
              					}
              				}
              			//}

              			s[rF][rG] = minCost;
          				rG = ft[rG];
          			}
          		}

          		if(lGminus1_in_preL == parent_of_lG_in_preL && lGminus1_in_preL != 0x7fffffff) { // lG is the leftmost child of its parent
    
          			if(hasRightPart) {
          				if(swap) {
          					if(DEBUG) {
          						ou << "D[" << to_string(parent_of_lG_in_preL) << ", " << to_string(endPathNode) << "] = " << "S[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "]" << endl; //rightmosts child of p(lG)
          					}
          				} else {
          					if(DEBUG) {
          						ou << "D[" << to_string(endPathNode) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "S[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "]" << endl; //rightmosts child of p(lG)
          					}
          				}
          			}
          			if (endPathNode > 0 && endPathNode == parent_of_endPathNode_preL + 1 && endPathNode_in_preR == parent_of_endPathNode_preR + 1) {//no left and right
          				if(swap) {
          					if(DEBUG) {
          						ou << "D[" << to_string(parent_of_lG_in_preL) << ", " << to_string(parent_of_endPathNode_preL) << "] = " << "S[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "]" << endl; //rightmosts child of p(lG)
          					}
          				} else {
          					if(DEBUG) {
          						ou << "D[" << to_string(parent_of_endPathNode_preL) << ", " << to_string(parent_of_lG_in_preL) << "] = " << "S[" << to_string(rFLast) << ", " << to_string(lGminus1_in_preR + 1) << "]" << endl; //rightmosts child of p(lG)
          					}
          				}
          			}

          			for (int rF = rFFirst; rF >= rFLast; rF--) {
          				if(DEBUG) {
          					ou << "Q[" << to_string(rF) << "] = " << "S[" << to_string(rF) << ", " << to_string(parent_of_lG_in_preR + 1) << "]" << endl;	
          				}
            		}
          		}
          		for (int rG = rGFirst; rG >= rGLast; rG = ft[rG]) {
            		if(DEBUG) {
            			ou << "T[" << to_string(lG) << ", " << to_string(rG) << "] = " << "S[" << to_string(rFLast) << ", " << to_string(rG) << "]" << endl;
            		}
          		}
        	}
		}
	rF = endPathNode_in_preR;//in D loop
	startPathNode = endPathNode;
    endPathNode = (*A_)[endPathNode] ->getParent() == NULL? -1 : (*A_)[endPathNode] ->getParent()->getID();	
	}
};

void TreeComparison::updateFnArray(int lnForNode, int node, int currentSubtreePreL) {
    if (lnForNode >= currentSubtreePreL) {
      fn[node] = fn[lnForNode];//the last leaf node whose next leaf is lnForNode
      fn[lnForNode] = node;// fn[lnfornode] points to the start point
      if(DEBUG) {
      	ou << "fn[" << to_string(node) << "] = fn[" << to_string(lnForNode) << "]" << endl; 
      	ou << "fn[" << to_string(lnForNode) << "] = " << to_string(node) << endl;
      }
    } else {
      int maxSize = treeSizeA < treeSizeB? treeSizeB + 1 : treeSizeA + 1;
      fn[node] = fn[maxSize];
      fn[maxSize] = node;
      if(DEBUG) {
      	ou << "fn[" << to_string(node) << "] = fn[" << to_string(maxSize) << "]" << endl; 
      	ou << "fn[" << to_string(maxSize) << "] = " << to_string(node) << endl;
      }
    }
}

void TreeComparison::updateFtArray(int lnForNode, int node) {
    ft[node] = lnForNode;
    if(DEBUG) {
    	ou << "ft[" << to_string(node) << "] = " << to_string(lnForNode) << endl;
    }
    if(fn[node] > -1) {
      ft[fn[node]] = node;
      if(DEBUG) {
      	ou << "ft[fn[" << to_string(node) << "]] = " << to_string(node) << endl;  
      }
    }
}


int TreeComparison::getPathType(Tree* tree, Node* node, int leaf) {
	if(tree->preL_to_lid[node->getID()] == leaf) return 0;
	else if(tree->preL_to_rid[node->getID()] == leaf) return 1;
	else return 2;
}

void TreeComparison::strategyComputation() {
	vector<Node*> preA = A_->getPreL();
	vector<Node*> preB = B_->getPreL();

	computeSumInsAndDelCost(A_);
	computeSumInsAndDelCost(B_);
	deltaInit();

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
	gted(preA[0], preB[0]);

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
#include "SimiMatrix.h"

#include <iostream>

using namespace std;

char& SimiMatrix::operator[] (int x){
	if(x >= 21) {
		cout << "Overflow" << endl;
		return base[0];
	}
	return base[i];

};

float& SimiMatrix::operator() (int x, int y) {
	if(x >= 50 || y >= 50) {
		cout << "Overflow" << endl;
		return simiMatrix[0][0];
	}
	return simiMatrix[x][y];

};
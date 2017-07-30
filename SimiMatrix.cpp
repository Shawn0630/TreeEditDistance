#include "SimiMatrix.h"

#include <iostream>

using namespace std;


SimiMatrix::SimiMatrix() {

};

char& SimiMatrix::operator[] (int x){
	if(x >= 21 || x < 0) {
		cout << to_string(x) << "Overflow" << endl;
		return base[0];
	}
	return base[x];
};

float& SimiMatrix::operator() (int x, int y) {
	if(x >= 50 || x < 0 || y >= 50 || y < 0) {
		cout << "Overflow" << endl;
		return simiMatrix[0][0];
	}
	return simiMatrix[x][y];
};

string SimiMatrix::toString(void) {
	string res = "";
	for(int i = 0; i < 21; i++) {
		res += base[i];// be careful of string + char[]
		res += " ";
	}
	res += "\n";
	for(int i = 0; i < 21; i++) {
		for(int j = 0; j < 21; j++) {
			res += to_string(simiMatrix[i][j]) + " ";
		}
		res += '\n';
	}
	return res;
};
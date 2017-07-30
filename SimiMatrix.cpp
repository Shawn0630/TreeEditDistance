#include "SimiMatrix.h"

#include <iostream>
#include <map>

using namespace std;


SimiMatrix::SimiMatrix() {

};

void SimiMatrix::setBase(char c, int i){
	base[c] = i;
};

int SimiMatrix::operator[] (char c) {
	if(base.find(c) == base.end()) return -1;
	return base[c];
}

float& SimiMatrix::operator() (int x, int y) {
	if(x >= 50 || x < 0 || y >= 50 || y < 0) {
		cout << "Overflow" << endl;
		return simiMatrix[0][0];
	}
	return simiMatrix[x][y];
};

string SimiMatrix::toString(void) {
	string res = "";
	/*for(int i = 0; i < 21; i++) {
		res += base[i];// be careful of string + char[]
		res += " ";
	}*/
	cout << base.size() << endl;
	map<char, int>::iterator it;
	for(it = base.begin(); it != base.end(); it++) {
		res += it->first;
		res += " ";
		res += to_string(it->second);
		res += "\n";
	}
	for(int i = 0; i < 21; i++) {
		for(int j = 0; j < 21; j++) {
			res += to_string(simiMatrix[i][j]) + " ";
		}
		res += '\n';
	}
	return res;
};
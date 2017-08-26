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

float SimiMatrix::del(char c) {
	if(base.find(c) == base.end()) return 0;
	int cid = base[c];
	return simiMatrix[cid][0];
};

float SimiMatrix::ins(char c) {
	if(base.find(c) == base.end()) return 0;
	int cid = base[c];
	return simiMatrix[cid][0];

};
float SimiMatrix::ren(char a, char b) {
	if(base.find(a) == base.end() || base.find(b) == base.end()) return 0;
	int aid = base[a];
	int bid = base[b];
	if(aid >= bid) return simiMatrix[aid][bid];
	else return simiMatrix[bid][aid];
};

string SimiMatrix::toString(void) {
	string res = "";
	/*for(int i = 0; i < 21; i++) {
		res += base[i];// be careful of string + char[]
		res += " ";
	}*/
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
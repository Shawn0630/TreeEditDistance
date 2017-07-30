#ifndef SIMIMATRIX_H
#define SIMIMATRIX_H

#include <string>
#include <map>

using namespace std;

class SimiMatrix {

private:
	map<char, int> base;
	float simiMatrix[50][50];

public:
	SimiMatrix();
	void setBase(char, int);
	int operator[] (char);
	float& operator() (int, int);
	float del(char);
	float ins(char);
	float ren(char, char);
	string toString(void);
};

#endif
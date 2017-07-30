#ifndef SIMIMATRIX_H
#define SIMIMATRIX_H

#include <string>

using namespace std;

class SimiMatrix {

private:
	char base[21];
	float simiMatrix[50][50];

public:
	SimiMatrix();
	char& operator[] (int);
	float& operator() (int, int);
	string toString(void);
};

#endif
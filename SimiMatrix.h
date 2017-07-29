#ifndef SIMIMATRIX_H
#define SIMIMATRIX_H


class SimiMatrix {

private:
	char base[21];
	float simiMatrix[50][50];

public:
	SimiMatrix();
	char& operator[] (int);
	float& operator() (int, int);

}

#endif
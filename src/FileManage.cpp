#include "FileManage.h"
#include "Errors.h"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

bool FileManage::instanceFlag_ = false;
FileManage* FileManage::instance_ = NULL;

FileManage::~FileManage() {
  instance_ = NULL;
  instanceFlag_ = false;
}

FileManage::FileManage() {
	rnaFileName_ = "";
	simiFileName_ = "";
}

FileManage* FileManage::getInstance(void) {
  if (!instanceFlag_) {
    instance_ = new FileManage();
    instanceFlag_ = true;
  }

  return instance_;
}

void FileManage::setRNAFileName(string rnaFileName) {
	rnaFileName_ = rnaFileName;
};


void FileManage::setSimiFileName(string simiFileName) {
	simiFileName_ = simiFileName;
};

vector<RNA> FileManage::readRNAsFromFile(void) {
	vector<RNA> res;
	if(rnaFileName_ == "") {
		readRNAFileError = true; 
		return res;
	}

	ifstream in;
	//cout << "rnaFileName_ = " << rnaFileName_ << endl;
	in.open(rnaFileName_);
	if(!in.is_open()) {
		readRNAFileError = true;
		return res;
	}
	string tmp;
	string RNA1Name, RNA2Name;
	RNA r1, r2;
	int RNA1Start, RNA2Start;
	int n1 = 2, n2 = 2;
	char c;
	
	getline(in, RNA1Name);
	r1.setRNAName(RNA1Name);
	in >> RNA1Start;
	// Read RNA1 Sequence
	do {
		in >> c;
		if(c >= 65 && c <= 90){
			r1[n1++] = c;   // the characters of the rna1
		}
	} while(c != '>');
	n1--;
	r1[1] = 'A';
	r1[++n1] = 'A';
	r1.setRNASize(n1);
	// Read RNA1 Secondary Structure
	do {
		in >> c;
		if(c == ')'){
			int start;
			int end;
			int length;
			in >> start >> end >> length;
			for(int i = start - RNA1Start + 1 + 1,j = end - RNA1Start + 1 + 1,k = 1; k <= length; i++, j--, k++){
			    r1(i) = j;   //pair of rna1
			    r1(j) = i;   //pair of rna1
		  	 }
		}
	} while(c != '>');
	r1(1) = n1;
	r1(n1) = 1; 
	// Read RNA1 Tertiary Structure
	do {
		in >> c;
	} while(c != '>');
	res.push_back(r1);

	getline(in, tmp);
	getline(in, tmp);

	getline(in, RNA2Name);
	r2.setRNAName(RNA2Name);
	in >> RNA2Start;
	// Read RNA2 Sequence
	do {
		in >> c;
		if(c >= 65 && c <= 90){
			r2[n2++] = c;   // the characters of the rna1
		}
	} while(c != '>');
	n2--;
	r2[1] = 'A';
	r2[++n2] = 'A';
	r2.setRNASize(n2);
	// Read RNA2 Secondary Structure
	do {
		in >> c;
		if(c == ')'){
			int start;
			int end;
			int length;
			in >> start >> end >> length;
			for(int i = start - RNA2Start + 1 + 1,j = end - RNA2Start + 1 + 1,k = 1; k <= length; i++, j--, k++){
			    r2(i) = j;   //pair of rna2
			    r2(j) = i;   //pair of rna2
		  	 }
		}
	} while(c != '>');
	r2(1) = n2;
	r2(n2) = 1; 
	// Read RNA2 Tertiary Structure
	do {
		in >> c;
	} while(c != '>');
	res.push_back(r2);
	in.close();
	return res;
};

SimiMatrix FileManage::readSimiFromFile(void) {
	SimiMatrix matrix;
	if(simiFileName_ == "") {
		readSimiFileError = true;
		return matrix;
	}

	ifstream in;
	string temp;
	char c;
	char index;
	int i = 0, j = 1;
	in.open(simiFileName_);
	//cout << "simiFileName_ = " << simiFileName_ << endl;
	if(!in.is_open()) {
		readSimiFileError = true;
		return matrix;
	}
	getline(in, temp);
	getline(in, temp);
	getline(in, temp);
	for(int k = 0; k < temp.length(); k++) {
		while(!(temp[k] >= 65 && temp[k] <= 90 || temp[k] == '-'))k++;
		if(k >= temp.length()) break;
		//matrix[i] = temp[k];
		index = temp[k];
		if(temp[k + 1] >= 65 && temp[k + 1] <= 90) {
			c = temp[k + 1];
			if(index == 'A' && c == 'A'){
				index = 'B';
			}
			else if(index == 'A' && c == 'G'){
				index = 'D';
			}
			else if(index == 'A' && c == 'C'){
				index = 'E';
			}
			else if(index == 'A' && c == 'U'){
				index = 'F';
			}
			else if(index == 'G' && c == 'A'){
				index = 'H';
			}
			else if(index == 'G' && c == 'G'){
				index = 'I';
			}
			else if(index == 'G' && c == 'C'){
				index = 'J';
			}
			else if(index == 'G' && c == 'U'){
				index = 'K';
			}
			else if(index == 'C' && c == 'A'){
				index = 'L';
			}
			else if(index == 'C' && c == 'G'){
				index = 'M';
			}
			else if(index == 'C' && c == 'C'){
				index = 'N';
			}
			else if(index == 'C' && c == 'U'){
				index = 'O';
			}
			else if(index == 'U' && c == 'A'){
				index = 'P';
			}
			else if(index == 'U' && c == 'G'){
				index = 'Q';
			}
			else if(index == 'U' && c == 'C'){
				index = 'R';
			}
			else if(index == 'U' && c == 'U'){
				index = 'S';
			}
			k++;
		}
		matrix.setBase(index, i);
		i++;
	}
	getline(in, temp);
	while(in >> c){
		if(j > 4) in >> c;
		for(int k = 0; k <= j; k++){
			float num = 0;
			in >> num;
			matrix(j, k) = num;
		}
		j++;
	}
	
	in.close();
	return matrix;

};

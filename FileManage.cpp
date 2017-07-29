#include "FileManage.h"

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
	if(rnaFileName_ == "") return res;
	vector<RNA> res;
	ifstream in;
	in.open(rnaFileName_);
	if(!in.is_open()) return res;
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
	ifstream in;
	in.open(simiFileName_);

	if(simiFileName == "") return matrix;
	if(!in.is_open()) return matrix;




};

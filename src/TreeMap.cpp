#include "TreeMap.h"
#include "Constant.h"

#include <iostream>

using namespace std;

TreeMap::TreeMap() {

}

TreeMap::TreeMap(Tree* F_, Tree* G_, SimiMatrix costModel_) {
	F = F_;
	G = G_;
	costModel = costModel_;
	
	FTreeSize = F->getTreeSize();
	GTreeSize = G->getTreeSize();

	FtoG = new int[FTreeSize];
	GtoF = new int[GTreeSize];

	fill_n(FtoG, FTreeSize, -1);
	fill_n(GtoF, GTreeSize, -1);

	counter = -1;
};

void TreeMap::setTreeF(Tree* F_) {
	F = F_;
};
void TreeMap::setTreeG(Tree* G_) {
	G = G_;
};
void TreeMap::setCostModel(SimiMatrix costModel_) {
	costModel = costModel_;
};

void TreeMap::init(void) {
	if(F == NULL || G == NULL) return;

	FTreeSize = F->getTreeSize();
	GTreeSize = G->getTreeSize();

	FtoG = new int[FTreeSize];
	GtoF = new int[GTreeSize];

	fill_n(FtoG, FTreeSize, -1);
	fill_n(GtoF, GTreeSize, -1);

	counter = -1;
};

void TreeMap::setMap(int Findex, int Gindex) {
	if(Findex > FTreeSize) {
		cout << "Findex = " << Findex << " FTreeSize = " << FTreeSize << endl; 
		cout << "TreeMap Overflow" << endl;
		return;
	}
	if(Gindex > GTreeSize) {
		cout << "Gindex = " << Gindex << " GTreeSzie = " << GTreeSize << endl;
		cout << "TreeMap Overflow" << endl;
		return;
	}

	if(Findex == -1) {
		GtoF[Gindex] = -1;
	} else if(Gindex == -1) {
		FtoG[Findex] = -1;
	} else {
		FtoG[Findex] = Gindex;
		GtoF[Gindex] = Findex;
	}

	return;
};

string TreeMap::toString(void) {
	string res = "";
	counter = 0;
	for(int i = 0; i < FTreeSize; i++) {
		res += to_string(FtoG[i]) + " ";
	}
	res += "\n";
	for(int i = 0; i < GTreeSize; i++) {
		res += to_string(GtoF[i]) + " ";
	}
	res += "\n";
	for(int i = 0; i < FTreeSize; i++) {
		if(FtoG[i] == -1) {
			res += to_string(i) + "(" + (*F)[i]->getLabel() + ") -> - \n";
			counter += costModel.del((*F)[i]->getLabel());
		} else {
			res += to_string(i) + "(" + (*F)[i]->getLabel() + ") -> " + to_string(FtoG[i]) + "(" + (*G)[FtoG[i]]->getLabel() + ") \n";
			counter += costModel.ren((*F)[i]->getLabel(), (*G)[FtoG[i]]->getLabel());
		}
	}
	for(int i = 0; i < GTreeSize; i++) {
		if(GtoF[i] != -1) continue;
		res += "- -> " + to_string(i) + "(" + (*G)[i]->getLabel() + ") \n";
		counter += costModel.ins((*G)[i]->getLabel());
	}

	return res;
}

int TreeMap::getCount() {
	//if(counter != -1) return counter;
	counter = 0;
	for(int i = 0; i < FTreeSize; i++) {
		if(FtoG[i] == -1) {
			counter += costModel.del((*F)[i]->getLabel());
			/*cout << costModel.del((*F)[i]->getLabel()) << endl;
			cout << "del " << (*F)[i]->getLabel() << " counter = " << counter << endl;*/
		} else {
			counter += costModel.ren((*F)[i]->getLabel(), (*G)[FtoG[i]]->getLabel());
			/*cout << costModel.ren((*F)[i]->getLabel(), (*G)[FtoG[i]]->getLabel()) << endl;
			cout << "ren " << (*F)[i]->getLabel()  << " -> " << (*G)[FtoG[i]]->getLabel() << " counter = " << counter << endl;*/
		}
	}
	for(int i = 0; i < GTreeSize; i++) {
		if(GtoF[i] != -1) continue;
		counter += costModel.ins((*G)[i]->getLabel());
		/*cout << costModel.ins((*G)[i]->getLabel()) << endl;
		cout << "ins " << (*G)[i]->getLabel() << " counter = " << counter << endl;*/
	}
	//cout << "Treemap counter = " << counter << endl;
	return counter;
};

//index F
int TreeMap::operator[](int Findex){
	if(Findex > FTreeSize) {
		cout << "Findex = " << Findex << " FTreeSize = " << FTreeSize << endl; 
		cout << "TreeMap Overflow" << endl;
		return 0;
	}

	return FtoG[Findex];
};
//index G
int TreeMap::operator()(int Gindex){
	if(Gindex > GTreeSize) {
		cout << "Gindex = " << Gindex << " GTreeSize = " << GTreeSize << endl; 
		cout << "TreeMap Overflow" << endl;
		return 0;
	}

	return GtoF[Gindex];
};

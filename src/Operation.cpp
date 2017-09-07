#include "Operation.h"


Operation::Operation() {

};

Operation::Operation(int from_ID, char from_label, int to_ID, char to_label) {
	from_ID_ = from_ID;
	from_label_ = from_label;
	to_ID_ = to_ID;
	to_label_ = to_label;
};


void Opeartion::setFromID(int from_ID) {
	from_ID_ = from_ID;
};
int Operation::getFromID(void) {
	return from_ID_;
};


void Opeartion::setFromLabel(char from_label) {
	from_label_ = from_label;
};
char Operation::getFromLabel(void) {
	return from_label_;
};


void Operation::setToID(int to_ID) {
	to_ID_ = to_ID;
};
int Operation::getToID(void) {
	return to_ID_;
};


void Operation::setToLabel(char to_label) {
	to_label_ = to_label;
};
char Operation::getToLabel(void) {
	return to_label_;
};


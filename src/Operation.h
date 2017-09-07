#ifndef OPERATION_H
#define OPERATION_H 1

class Operation {
private:
	int from_ID_;// -1 if insert to G
	int to_ID_;// -1 if delete from F
	char from_label_;// '-' if insert to G
	char to_label_;// '-' if delete from F
public:

	Operation();
	Operation(int, char, int, char);

	void setFromID(int);
	int getFromID(void);

	void setFromLabel(char);
	char getFromLabel(void);

	void setToID(int);
	int getToID(void);

	void setToLabel(char);
	char getToLabel(void);

};

#endif
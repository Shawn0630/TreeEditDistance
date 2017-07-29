#ifndef INCLUDED_FILEMANAGE_H
#define INCLUDED_FILEMANAGE_H 1

#include "Node.h"
#include "Tree.h"
#include "RNA.h"
#include "SimiMatrix.h"

#include <string>
#include <vector>
using namespace std;


class FileManage {

public:
  ~FileManage();
  static FileManage *getInstance(void);
  void setRNAFileName(string);
  void setSimiFileName(string);
  vector<RNA> readRNAsFromFile(void);
  SimiMatrix readSimiFromFile(void);

private:
  FileManage();
  string rnaFileName_;
  string simiFileName_;
  static FileManage *instance_;
  static bool instanceFlag_;
};

#endif /* INCLUDED_FILEMANAGE_H */
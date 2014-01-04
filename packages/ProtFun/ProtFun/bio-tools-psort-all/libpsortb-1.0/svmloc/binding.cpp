#include "binding.h"

SVM* createSVM(int st, int kt, int d, double g, double c0, double C, double nu, double e)
{
  return new SVM(st, kt, d, g, c0, C, nu, e);
}

void destroySVM(SVM* pSVM)
{
  delete pSVM;
}

int loadSVMModel(SVM* pSVM, char* filename)
{
  return pSVM->loadModel(filename);
}

int loadSVMFreqPattern(SVM* pSVM, char* filename)
{
  return pSVM->loadFreqPattern(filename);
}

double SVMClassify(SVM* pSVM, char* seq)
{
  return pSVM->classify(seq);
}


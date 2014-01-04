#ifndef __BINDING_H__
#define __BINDING_H__

#ifdef __cplusplus

#include "svmloc.h"

#else
  typedef 
    struct SVM
      SVM;
#endif

#ifdef __cplusplus
extern "C" {
#endif

  extern SVM* createSVM(int st, int kt, int d, double g, double c0, double C, double nu, double e);
  extern void destroySVM(SVM*);
  extern int loadSVMModel(SVM*, char*);
  extern int loadSVMFreqPattern(SVM*, char*);
  extern double SVMClassify(SVM*, char*);

#ifdef __cplusplus
}
#endif

#endif

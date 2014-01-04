#include <vector>
#include <map>

#include "bindings.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#ifdef __cplusplus
}
#endif

DataSet *_new_dataset(double l) {

  return new DataSet(l);
}

SVM *_new_svm(int st, int kt, int d, double g, double c0, double C,
	      double nu, double e) {

  return new SVM(st, kt, d, g, c0, C, nu, e);
}

MODULE = Bio::Tools::PSort::SVMLoc::DataSet	PACKAGE = Bio::Tools::PSort::SVMLoc::DataSet

DataSet *
_new_dataset(l)
     double l

double
DataSet::_getLabel()
    CODE:
      RETVAL = THIS->getLabel();
    OUTPUT:
      RETVAL

void
DataSet::_setLabel(l)
     double l
    CODE:
      THIS->setLabel(l);

double
DataSet::_getAttribute(k)
     int k
    CODE:
      RETVAL = THIS->getAttribute(k);
    OUTPUT:
      RETVAL

void
DataSet::_setAttribute(k,v)
     int k
     double v
    CODE:
      THIS->setAttribute(k,v);

void
DataSet::DESTROY()

MODULE = Bio::Tools::PSort::SVMLoc			PACKAGE = Bio::Tools::PSort::SVMLoc

SVM *
_new_svm(st,kt,d,g,c0,C,nu,e)
     int st
     int kt
     int d
     double g
     double c0
     double C
     double nu
     double e

double
SVM::_classify(seq)
     char *seq
    CODE:
      RETVAL = THIS->classify(seq);
    OUTPUT:
      RETVAL

int
SVM::_loadModel(filename)
     char *filename
    CODE:
      RETVAL = THIS->loadModel(filename);
    OUTPUT:
      RETVAL

int
SVM::_loadFreqPattern(filename)
     char *filename
    CODE:
      RETVAL = THIS->loadFreqPattern(filename);
    OUTPUT:
      RETVAL

double
SVM::_predict(ds)
     DataSet *ds
    CODE:
      RETVAL = THIS->predict(ds);
    OUTPUT:
      RETVAL

void
SVM::DESTROY()

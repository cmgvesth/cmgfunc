#ifdef __cplusplus
extern "C" {
#endif

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#ifdef __cplusplus
}
#endif

#include "hmm.h"
#include "globals.h"

HMM *hmm_new(char *filename, int df, int n2) {

  return new HMM(filename, df, n2);
}

HMMReport *hmm_report_new() {

  return new HMMReport();
}

HMMDomainHit *hmm_domain_hit_new(float s, double p, double e) {

  return new HMMDomainHit(s, p, e);
}

HMMGlobalHit *hmm_global_hit_new(float s, double p, double e) {

  return new HMMGlobalHit(s, p, e);
}


MODULE = Algorithm::HMM        PACKAGE = Algorithm::HMM::Hit::Domain

HMMDomainHit *
hmm_domain_hit_new(s, p, e)
     float s
     double p
     double e

double
HMMDomainHit::_getPValue()
    CODE:
        RETVAL = THIS->getPValue();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setPValue(double i)
    CODE:
        THIS->setPValue(i);

double
HMMDomainHit::_getEValue()
    CODE:
        RETVAL = THIS->getEValue();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setEValue(double i)
    CODE:
        THIS->setEValue(i);

float
HMMDomainHit::_getScore()
    CODE:
        RETVAL = THIS->getScore();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setScore(float i)
    CODE:
        THIS->setScore(i);

int
HMMDomainHit::_getSeqFrom()
    CODE:
        RETVAL = THIS->getSeqFrom();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setSeqFrom(i)
     int i
    CODE:
        THIS->setSeqFrom(i);

int
HMMDomainHit::_getSeqTo()
    CODE:
        RETVAL = THIS->getSeqTo();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setSeqTo(i)
     int i
    CODE:
        THIS->setSeqTo(i);

int
HMMDomainHit::_getSeqLength()
    CODE:
        RETVAL = THIS->getSeqLength();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setSeqLength(i)
     int i
    CODE:
        THIS->setSeqLength(i);

int
HMMDomainHit::_getHMMFrom()
    CODE:
        RETVAL = THIS->getHMMFrom();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setHMMFrom(i)
     int i
    CODE:
        THIS->setHMMFrom(i);

int
HMMDomainHit::_getHMMTo()
    CODE:
        RETVAL = THIS->getHMMTo();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setHMMTo(i)
     int i
    CODE:
        THIS->setHMMTo(i);

double
HMMDomainHit::_getParentPValue()
    CODE:
        RETVAL = THIS->getParentPValue();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setParentPValue(i)
     double i
    CODE:
        THIS->setParentPValue(i);

double
HMMDomainHit::_getParentEValue()
    CODE:
        RETVAL = THIS->getParentEValue();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setParentEValue(i)
     double i
    CODE:
        THIS->setParentEValue(i);

float
HMMDomainHit::_getParentScore()
    CODE:
        RETVAL = THIS->getParentScore();
    OUTPUT:
        RETVAL

void
HMMDomainHit::_setParentScore(i)
     float i
    CODE:
        THIS->setParentScore(i);


MODULE = Algorithm::HMM        PACKAGE = Algorithm::HMM::Hit::Global

HMMGlobalHit *
hmm_global_hit_new(s, p, e)
     float s
     double p
     double e

double
HMMGlobalHit::_getPValue()
    CODE:
        RETVAL = THIS->getPValue();
    OUTPUT:
        RETVAL

void
HMMGlobalHit::_setPValue(double i)
    CODE:
        THIS->setPValue(i);

double
HMMGlobalHit::_getEValue()
    CODE:
        RETVAL = THIS->getEValue();
    OUTPUT:
        RETVAL

void
HMMGlobalHit::_setEValue(double i)
    CODE:
        THIS->setEValue(i);

float
HMMGlobalHit::_getScore()
    CODE:
        RETVAL = THIS->getScore();
    OUTPUT:
        RETVAL

void
HMMGlobalHit::_setScore(float i)
    CODE:
        THIS->setScore(i);


MODULE = Algorithm::HMM        PACKAGE = Algorithm::HMM::Report

HMMReport *
hmm_report_new()

void
HMMReport::_addGlobalHit(h)
     HMMGlobalHit *h
    CODE:
        THIS->addGlobalHit(h);

int
HMMReport::_numGlobalHits()
    CODE:
        RETVAL = THIS->numGlobalHits();
    OUTPUT:
        RETVAL

HMMGlobalHit *
HMMReport::_getGlobalHit(i)
     int i
    CODE:
        RETVAL = THIS->getGlobalHit(i);
    OUTPUT:
        RETVAL

void
HMMReport::_addDomainHit(h)
     HMMDomainHit *h
    CODE:
        THIS->addDomainHit(h);

int
HMMReport::_numDomainHits()
    CODE:
        RETVAL = THIS->numDomainHits();
    OUTPUT:
        RETVAL

HMMDomainHit *
HMMReport::_getDomainHit(i)
     int i
    CODE:
        RETVAL = THIS->getDomainHit(i);
    OUTPUT:
        RETVAL

void
HMMReport::DESTROY()

MODULE = Algorithm::HMM        PACKAGE = Algorithm::HMM

HMM *
hmm_new(filename, df, n2)
     char *filename
     int df
     int n2

int
HMM::_save(filename, append, binary)
     char *filename
     int append
     int binary
    CODE:
        RETVAL = THIS->save(filename, append, binary);
    OUTPUT:
        RETVAL

int
HMM::_load(filename)
     char *filename
    CODE:
        RETVAL = THIS->load(filename);
    OUTPUT:
        RETVAL

int
HMM::_modelLoaded()
    CODE:
        RETVAL = THIS->modelLoaded();
    OUTPUT:
        RETVAL

int
HMM::_getDoForward()
    CODE:
        RETVAL = THIS->getDoForward();
    OUTPUT:
        RETVAL

void
HMM::_setDoForward(df)
     int df
    CODE:
        THIS->setDoForward(df);

int
HMM::_getDoNull2()
    CODE:
        RETVAL = THIS->getDoNull2();
    OUTPUT:
        RETVAL

void
HMM::_setDoNull2(n2)
     int n2
    CODE:
        THIS->setDoNull2(n2);

HMMReport*
HMM::_search(seq)
     char *seq
    CODE:
        RETVAL = THIS->search(seq);
    OUTPUT:
        RETVAL

int
HMM::_train(...)
    CODE:
	char **seqs = (char **)malloc(sizeof(char *) * items);
        if(seqs == NULL) {
	  RETVAL = 0;
        } else {
	  int i;
	  for(i = 1; i < items; i++) {
	    SV *sv = (SV *)ST(i);
	    seqs[i - 1] = (char *)SvPV_nolen((SV *)sv);
	  }
	  
	  RETVAL = THIS->train(seqs, items - 1);
	}
    OUTPUT:
        RETVAL

void
HMM::DESTROY()

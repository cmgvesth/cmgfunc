#ifndef __HMM_H__
#define __HMM_H__

extern "C" {
#include <float.h>
#include <structs.h>
#include <funcs.h>
#include <squid.h>
#include <msa.h>
#include <string.h>
}

#include <vector>

using namespace std;

class HMMHit {
 public:
  HMMHit(float s, double p, double e);
  double getPValue() { return pvalue; };
  void   setPValue(double p) { pvalue = p; };
  double getEValue() { return evalue; };
  void   setEValue(double e) { evalue = e; };
  float  getScore() { return score; };
  void   setScore(float s) { score = s; };
  ~HMMHit() { };
 private:
  double pvalue;
  double evalue;
  float  score;
};

class HMMGlobalHit : public HMMHit {
 public:
  HMMGlobalHit(float s, double p, double e) : HMMHit(s, p, e) { };
 private:
};

class HMMDomainHit : public HMMHit {
 public:
  HMMDomainHit(float s, double p, double e);
  void   setSeqFrom(int s) { seqfrom = s; };
  int    getSeqFrom() { return seqfrom; };
  void   setSeqTo(int s) { seqto = s; };
  int    getSeqTo() { return seqto; };
  void   setSeqLength(int s) { seqlen = s; };
  int    getSeqLength() { return seqlen; };
  void   setHMMFrom(int h) { hmmfrom = h; };
  int    getHMMFrom() { return hmmfrom; };
  void   setHMMTo(int h) { hmmto = h; };
  int    getHMMTo() { return hmmto; };
  double getParentPValue() { return ppvalue; };
  void   setParentPValue(double p) { ppvalue = p; };
  double getParentEValue() { return pevalue; };
  void   setParentEValue(double e) { pevalue = e; };
  float  getParentScore() { return pscore; };
  void   setParentScore(float s) { pscore = s; };
 private:
  int    seqfrom, seqto, seqlen, hmmfrom, hmmto;
  double ppvalue, pevalue;
  float  pscore;
};

class HMMReport {
 public:
  HMMReport() { };
  void addGlobalHit(HMMGlobalHit *h);
  HMMGlobalHit *getGlobalHit(int i);
  int numGlobalHits() { return global.size(); };
  void addDomainHit(HMMDomainHit *h);
  HMMDomainHit *getDomainHit(int i);
  int numDomainHits() { return domain.size(); };
  ~HMMReport();
 private:
  vector<HMMGlobalHit *> global;
  vector<HMMDomainHit *> domain;
};

class HMM {
 public:
  HMM() { init(NULL, 0, 1); };
  HMM(char *filename) { init(filename, 0, 1); };
  HMM(char *filename, int df, int n2) { init(filename, df, n2); };
  HMMReport *search(char *seq);
  int  train(char **seqs, int len);
  int  save(char *filename, int append, int binary);
  int  load(char *filename);
  int  modelLoaded() { return (model != NULL); };
  int  getDoForward() { return doForward; };
  void setDoForward(int df) { doForward = df; };
  int  getDoNull2() { return doNull2; };
  void setDoNull2(int n2) { doNull2 = n2; };
  ~HMM();
 private:
  void init(char *filename, int df, int n2);

  struct threshold_s thresh;
  struct plan7_s *model;
  int doForward;
  int doNull2;
};

#endif

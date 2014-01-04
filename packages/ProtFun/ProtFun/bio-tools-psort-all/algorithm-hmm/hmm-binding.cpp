#include "hmm.h"

// HMMHit constructors and methods 

HMMHit::HMMHit(float s, double p, double e) {

  pvalue = p;
  evalue = e;
  score  = s;
};

// HMMDomainHit constructors and methods  

HMMDomainHit::HMMDomainHit(float s, double p, double e) : HMMHit(s, p, e)  {

  ppvalue = pevalue = pscore = seqfrom = seqto = seqlen = hmmfrom = hmmto = 0;
}

// HMMReport constructors and methods 

void HMMReport::addGlobalHit(HMMGlobalHit *h) {

  if(h != NULL) global.push_back(h);
}

HMMGlobalHit *HMMReport::getGlobalHit(int i) {

  return ((i >= 0) && (i < global.size())) ? global[i] : NULL;
}

void HMMReport::addDomainHit(HMMDomainHit *h) {

  if(h != NULL) domain.push_back(h);
}

HMMDomainHit *HMMReport::getDomainHit(int i) {

  return ((i >= 0) && (i < domain.size())) ? domain[i] : NULL;
}

HMMReport::~HMMReport() {

  for(int i = 0; i < domain.size(); i++)
    if(domain[i] != NULL) delete domain[i];

  for(int i = 0; i < global.size(); i++)
    if(global[i] != NULL) delete global[i];
}

// HMM constructors and methods 

int HMM::save(char *filename, int append, int binary) {
  FILE *fp;

  if((model == NULL) || (filename == NULL)) return -1;

  fp = (append) ? fopen(filename, "a") : fopen(filename, "w");
  if(! fp) return -1;

  (binary) ? WriteBinHMM(fp, model) : WriteAscHMM(fp, model);

  fclose(fp);

  return 0;
}

int HMM::load(char *filename) {
  struct plan7_s *tmp_model = NULL;
  HMMFILE *hmmfp;

  if(strlen(filename) == 0) return -1;

  // Load the model from the specified file.
  if((hmmfp = HMMFileOpen(filename, "HMMERDB")) == NULL) return -1;

  if(! HMMFileRead(hmmfp, &tmp_model)) {
    HMMFileClose(hmmfp);
    return -1;
  }

  HMMFileClose(hmmfp);
  if(tmp_model == NULL) return -1;

  P7Logoddsify(tmp_model, 1);

  // If we had an old model already loaded, free it and make the new one
  // current.
  if(model != NULL) FreePlan7(model);
  model = tmp_model;

  SqdClean();

  return 0;
}

HMMReport *HMM::search(char *seq) {
  struct tophit_s *ghit, *dhit;
  struct p7trace_s *trace;
  double pvalue, evalue;
  int numhits;
  HMMReport *hmmrep;
  float score;
  int seqlen;
  char *dseq;

  // Ensure we have a valid model and sequence.
  if((model == NULL) ||(seq == NULL) || (!(seqlen = strlen(seq)))) return NULL;

  // Convert the sequence to the format the HMMER library expects.
  dseq = DigitizeSequence(seq, seqlen);

  // Calculate the raw scores for the sequence.
  if (P7ViterbiSize(seqlen, model->M) <= RAMLIMIT)
    score = P7Viterbi(dseq, seqlen, model, &trace);
  else
    score = P7SmallViterbi(dseq, seqlen, model, &trace);


  // Calculate the score using the forward algorithm if we need.
  if(doForward) {
    score = P7Forward(dseq, seqlen, model, NULL);
    if(doNull2) score -= TraceScoreCorrection(model, trace, dseq); 
  }

  // We're only looking at a single sequence, so the evalue will be the same
  // as the pvalue.
  evalue = pvalue = PValue(model, score);

  // If we met the threshold criteria, go and collect the global/domain
  // hit information we're going to return.
  if (score >= thresh.globT && evalue <= thresh.globE) {
    ghit = AllocTophits(200);
    dhit = AllocTophits(200);

    PostprocessSignificantHit(ghit, dhit, trace, model, dseq, seqlen, "SEQ",
			      NULL, NULL, doForward, score, doNull2, &thresh,
			      FALSE);
    numhits = ghit->num + dhit->num;


    // Allocate an array of hits that we'll return.
    hmmrep = new HMMReport();
    if(hmmrep == NULL) {
      P7FreeTrace(trace);
      FreeTophits(ghit);
      FreeTophits(dhit);
      free(dseq);
      SqdClean();
      return NULL;
    }


    // Sort the global hits by score, filtering out the ones that make the
    // cutoff and adding them to the array to be returned.
    FullSortTophits(ghit);
    for(int i = 0; i < ghit->num; i++) {
      double ev, pv;
      float sc;
      int ndom;

      GetRankedHit(ghit, i, &pv, &sc, NULL, NULL, NULL, NULL, NULL, NULL,
		   NULL, NULL, NULL, NULL, NULL, NULL, &ndom, NULL);
      ev = pv * (double)thresh.Z;

      if(ev <= thresh.globE && sc >= thresh.globT) {
	HMMGlobalHit *hit = new HMMGlobalHit(sc, pv, ev);
	hmmrep->addGlobalHit(hit);
      }
    }

    // Sort the domain hits by score, filtering out the ones that make the
    // cutoff and adding them to the array to be returned.
    FullSortTophits(dhit);
    for(int i = 0; i < dhit->num; i++) {
      int sqfrom, sqto, sqlen;
      double motherp, ev, pv;
      int hmmfrom, hmmto;
      float mothersc;
      int domidx;
      float sc;
      int ndom;

      GetRankedHit(dhit, i, 
		   &pv, &sc, &motherp, &mothersc,
		   NULL, NULL, NULL,
		   &sqfrom, &sqto, &sqlen,
		   &hmmfrom, &hmmto, NULL,
		   &domidx, &ndom,
		   NULL);
      ev = pv * (double)thresh.Z;

      if(motherp * (double)thresh.Z > thresh.globE || mothersc < thresh.globT)
	continue;
      else if (ev <= thresh.domE && sc >= thresh.domT) {
	HMMDomainHit *hit = new HMMDomainHit(sc, pv, ev);

	hit->setParentPValue(motherp);
	hit->setParentScore(mothersc);
	hit->setSeqFrom(sqfrom);
	hit->setSeqTo(sqto);
	hit->setSeqLength(sqlen);
	hit->setHMMFrom(hmmfrom);
	hit->setHMMTo(hmmto);

	hmmrep->addDomainHit(hit);
      }
    }
  }

  P7FreeTrace(trace);
  FreeTophits(ghit);
  FreeTophits(dhit);
  free(dseq);
  SqdClean();

  return hmmrep;
}

int HMM::train(char **seqs, int len) {
  MSA *msa;

  msa = (MSA *)malloc(sizeof(MSA));
  if(msa == NULL) return 0;

  for(int i = 0; i < len; i++) {
    printf("seq = %s\n", seqs[i]);
  }

  return 0;
}


HMM::~HMM() {

  if(model != NULL) FreePlan7(model);
}

void HMM::init(char *filename, int df, int n2) {

  thresh.globE   = 10.0;
  thresh.globT   = -FLT_MAX;
  thresh.domT    = -FLT_MAX;
  thresh.domE    = FLT_MAX;
  thresh.autocut = CUT_NONE;
  thresh.Z       = 1;

  doForward = df;
  doNull2   = n2;

  model = NULL;

  if(filename!= NULL) load(filename);
}

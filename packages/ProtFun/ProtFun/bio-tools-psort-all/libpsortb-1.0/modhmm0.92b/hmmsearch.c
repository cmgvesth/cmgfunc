#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "structs.h"
#include "funcs.h"
#include "cmdline_hmmsearch.h"


#define NORM_LOG_LIKE 0
#define LOGODDS 1

#define HMM 20
#define SEQ 21

#define LONGEST_SEQ -1 /* Note: this constant is also defined in read_seqs.c */
#define FIRST_SEQ 1

#define LEAD_SEQ 10
#define ALL_SEQS 11

#define MAX_LINE 500

#define HMMFILEERROR -10
#define REPFILEERROR -11
#define SEQERROR -12
#define PATHERROR -13

extern int verbose;


static struct hmm_multi_s hmm;
static struct hmm_multi_s retrain_hmm;
static struct msa_sequences_multi_s *msa_seq_infop;
static struct sequences_multi_s seq_info;
static struct replacement_letter_multi_s replacement_letters;
static struct null_model_multi_s null_model;
double *subst_mtxp;
double *subst_mtxp_2;
double *subst_mtxp_3;
double *subst_mtxp_4;
double *aa_freqs;
double *aa_freqs_2;
double *aa_freqs_3;
double *aa_freqs_4;


void get_null_model_multi(FILE *nullfile);
double get_nullmodel_score_multi(struct letter_s *seq, struct letter_s *seq_2, struct letter_s *seq_3,
				 struct letter_s *seq_4, int seq_len, int multi_scoring_method);
double get_nullmodel_score_msa_multi(int seq_len, int prf_mode, int use_prior,
				     int use_gap_shares, int normalize,
				     int scoring_method, int multi_scoring_method, double *aa_freqs,
				     double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
int get_scores(int use_labels, int scoring_method, int multi_scoring_method, helix_sites* hSites);
void get_post_prob_matrix(double **ret_post_prob_matrixp, double forw_score, struct forward_s *forw_mtx,
			  struct backward_s *backw_mtx, int seq_len);

helix_sites* get_helices(char* seq, char* hmmfilename, char* repfilename, char* path)
{

  /* command line options */
  helix_sites *hSites;
  FILE *hmmfile; /* file to read hmms from */
  FILE *seqnamefile; /* file for reading sequence names */
  FILE *replfile; /* file for reading special replacement letters */
  FILE *freqfile;
  int use_labels;
  int scoring_method;
  int multi_scoring_method;
  int helices;
  
  /* temporary variables */
  int i,j; /*standard loop indices */
 
  /*init some variables */
  hmmfile = NULL;
  seqnamefile = NULL;
  replfile = NULL;
  freqfile = NULL;
  use_labels = NO;
  scoring_method = DOT_PRODUCT;
  multi_scoring_method = JOINT_PROB;
  subst_mtxp = NULL;
  subst_mtxp_2 = NULL;
  subst_mtxp_3 = NULL;
  subst_mtxp_4 = NULL;
  aa_freqs = NULL;
  aa_freqs_2 = NULL;
  aa_freqs_3 = NULL;
  aa_freqs_4 = NULL;

  
  /* Create the structure to return the helix sites in */
  hSites = (helix_sites *) malloc_or_die(sizeof (helix_sites));
  hSites->helix_count = 0;
  hSites->helix = NULL;

  /* compulsory options */
  if(hmmfilename) {
    if((hmmfile = fopen(hmmfilename, "r")) == NULL) {
      hSites->helix_count = HMMFILEERROR;
      return hSites;
    }
  } else {
    hSites->helix_count = HMMFILEERROR;
    return hSites;
  }

  if(!path) {
    hSites->helix_count = PATHERROR;
    return hSites;
  }

  if(!seq) {
    hSites->helix_count = SEQERROR;
    return hSites;
  }
  
  if(repfilename) {
      if((replfile = fopen(repfilename, "r")) == NULL) {
	hSites->helix_count = REPFILEERROR;
	return hSites;
      }
  } else {
    hSites->helix_count = REPFILEERROR;
    return hSites;
  }

  
  /* get replacement letters */
  get_replacement_letters_multi(replfile, &replacement_letters);

  /* read null model - always null, don't bother */
  null_model.a_size = -1;
  null_model.nr_alphabets = -1;


  /* get hmm and score all seqs against them */
  /* get hmm from file */
  readhmm(hmmfile, &hmm, path);

  hmm.subst_mtx = subst_mtxp;
  hmm.subst_mtx_2 = subst_mtxp_2;
  hmm.subst_mtx_3 = subst_mtxp_3;
  hmm.subst_mtx_4 = subst_mtxp_4;
  hmm.replacement_letters = &replacement_letters;
  

  /* check sequence file for labels + check nolabel flag */
  use_labels = NO;
	
  /* allocate space for msa_seq_info structs */
  seq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
  seq_info.nr_alphabets = hmm.nr_alphabets;
  seq_info.nr_seqs = 1;
  seq_info.longest_seq = 0;
  seq_info.shortest_seq = INT_MAX;
  seq_info.avg_seq_len = 0;
	

  get_sequence_fasta_multi(seq, &seq_info, 0);
  hmm.alphabet_type = DISCRETE;

  seq_info.avg_seq_len = ((int)(seq_info.avg_seq_len / seq_info.nr_seqs));

	
  /* get score info for this seq  */

  helices = get_scores(use_labels, scoring_method, multi_scoring_method, hSites);
	
  /* deallocate seqinfo */
  free(((seq_info.seqs))->seq_1);

  free(msa_seq_infop);

      
  /* deallocate hmm_info */
  hmm_garbage_collection_multi(hmmfile, &hmm);


  /* deallocate replacement letters */
  if(replfile != NULL) {
    if(replacement_letters.nr_rl_1 > 0) {
      free(replacement_letters.letters_1);
      free(replacement_letters.probs_1);
    }
    fclose(replfile);

  }

  return hSites;

}


int find_helices(char* seq, char* hmmfilename, char* repfilename, char* path)
{
  /* command line options */
  FILE *hmmfile; /* file to read hmms from */
  FILE *seqnamefile; /* file for reading sequence names */
  FILE *replfile; /* file for reading special replacement letters */
  FILE *freqfile;
  int use_labels;
  int scoring_method;
  int multi_scoring_method;
  int helices;
  
  /* temporary variables */
  int i,j; /*standard loop indices */
 
  /*init some variables */
  hmmfile = NULL;
  seqnamefile = NULL;
  replfile = NULL;
  freqfile = NULL;
  use_labels = NO;
  scoring_method = DOT_PRODUCT;
  multi_scoring_method = JOINT_PROB;
  subst_mtxp = NULL;
  subst_mtxp_2 = NULL;
  subst_mtxp_3 = NULL;
  subst_mtxp_4 = NULL;
  aa_freqs = NULL;
  aa_freqs_2 = NULL;
  aa_freqs_3 = NULL;
  aa_freqs_4 = NULL;


  /* compulsory options */
  if(hmmfilename) {
    if((hmmfile = fopen(hmmfilename, "r")) == NULL) {
	return HMMFILEERROR;
    }
  } else {
      return HMMFILEERROR;
  }

  if(!path) {
    return PATHERROR;
  }

  if(!seq) {
      return SEQERROR;
  }
  
  if(repfilename) {
      if((replfile = fopen(repfilename, "r")) == NULL) {
	  return REPFILEERROR;
      }
  } else {
      return REPFILEERROR;
  }

  
  /* get replacement letters */
  get_replacement_letters_multi(replfile, &replacement_letters);

  /* read null model - always null, don't bother */
  null_model.a_size = -1;
  null_model.nr_alphabets = -1;


  /* get hmm and score all seqs against them */
  /* get hmm from file */
  readhmm(hmmfile, &hmm, path);

  hmm.subst_mtx = subst_mtxp;
  hmm.subst_mtx_2 = subst_mtxp_2;
  hmm.subst_mtx_3 = subst_mtxp_3;
  hmm.subst_mtx_4 = subst_mtxp_4;
  hmm.replacement_letters = &replacement_letters;
  

  /* check sequence file for labels + check nolabel flag */
  use_labels = NO;
	
  /* allocate space for msa_seq_info structs */
  seq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
  seq_info.nr_alphabets = hmm.nr_alphabets;
  seq_info.nr_seqs = 1;
  seq_info.longest_seq = 0;
  seq_info.shortest_seq = INT_MAX;
  seq_info.avg_seq_len = 0;
	

  get_sequence_fasta_multi(seq, &seq_info, 0);
  hmm.alphabet_type = DISCRETE;

  seq_info.avg_seq_len = ((int)(seq_info.avg_seq_len / seq_info.nr_seqs));

	
  /* get score info for this seq  */

  helices = get_scores(use_labels, scoring_method, multi_scoring_method, NULL);
	
  /* deallocate seqinfo */
  free(((seq_info.seqs))->seq_1);

  free(msa_seq_infop);

      
  /* deallocate hmm_info */
  hmm_garbage_collection_multi(hmmfile, &hmm);


  /* deallocate replacement letters */
  if(replfile != NULL) {
    if(replacement_letters.nr_rl_1 > 0) {
      free(replacement_letters.letters_1);
      free(replacement_letters.probs_1);      
    }
    fclose(replfile);
  }

  return helices;
}





/*************************************************************************************************/
void get_null_model_multi(FILE *nullfile)
{
  int i,j;
  char s[MAX_LINE];
  
  if(nullfile != NULL) {
    /* read nullfile */
    while(1) {
      if(fgets(s, MAX_LINE, nullfile) == NULL) {
	printf("Could not read null model file\n");
	exit(0);
      }
      if(s[0] == '#' || s[0] == '\n') {
	continue;
      }
      else {
	null_model.nr_alphabets = atoi(s);
	break;
      }
    }

    for(i = 0; i < null_model.nr_alphabets; i++) {
      while(1) {
	if(fgets(s, MAX_LINE, nullfile) == NULL) {
	  printf("Could not read null model file\n");
	  exit(0);
	}
	if(s[0] == '#' || s[0] == '\n') {
	  continue;
	}
	else {
	  switch(i) {
	  case 0:
	    null_model.a_size = atoi(s);
	    null_model.emissions = (double*)malloc_or_die(null_model.a_size * sizeof(double));
	    break;
	  case 1:
	    null_model.a_size_2 = atoi(s);
	    null_model.emissions_2 = (double*)malloc_or_die(null_model.a_size_2 * sizeof(double));
	    break;
	  case 2:
	    null_model.a_size_3 = atoi(s);
	    null_model.emissions_3 = (double*)malloc_or_die(null_model.a_size_3 * sizeof(double));
	    break;
	  case 3:
	    null_model.a_size_4 = atoi(s);
	    null_model.emissions_4 = (double*)malloc_or_die(null_model.a_size_4 * sizeof(double));
	    break;
	  }
	  break;
	}
      }
      j = 0;
      switch(i) {
      case 0:
	while(j < null_model.a_size) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      case 1:
	while(j < null_model.a_size_2) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_2 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;	
      case 2:
	while(j < null_model.a_size_3) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_3 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      case 3:
	while(j < null_model.a_size_4) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_4 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      }  
    }
    while(1) {
      if(fgets(s, MAX_LINE, nullfile) != NULL) {
	if(s[0] != '#' && s[0] != '\n') {
	  null_model.trans_prob = atof(s);
	  break;
	}
      }
      else {
	printf("Could not read null model file\n");
	exit(0);
      }
    }
  }
  else {
    null_model.a_size = -1;
    null_model.nr_alphabets = -1;
  }
}



/*********************score methods******************************************/
int get_scores(int use_labels, int scoring_method, int multi_scoring_method, helix_sites* hSites)
{
  struct forward_s *forw_mtx;
  struct one_best_s *one_best_mtx;
  struct backward_s *backw_mtx;
  struct forward_s *rev_forw_mtx;
  double *forw_scale, *rev_forw_scale, *scaling_factors;
  char *labels;
  struct letter_s *reverse_seq, *reverse_seq_2, *reverse_seq_3, *reverse_seq_4;
  struct helix_site *helixSite;

  int seq_len;
  double forward_score, backward_score, rev_forward_score, raw_forward_score;

  int a,b;
  int i,j,k;

  struct hmm_multi_s *hmmp;
  int IN_HELIX;
  int helices;

  hmmp = &hmm;


  /* get seq_length */
  seq_len = get_seq_length(seq_info.seqs->seq_1);


  /* if needed run forward */
  forward_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		&forw_mtx, &forw_scale, use_labels, multi_scoring_method);
  raw_forward_score = (forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
  forward_score = log10((forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
  for (j = seq_len; j > 0; j--) {
      forward_score = forward_score + log10(*(forw_scale + j));
  }

  /* if needed run backward */
  backward_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		 &backw_mtx, forw_scale, use_labels, multi_scoring_method);
  backward_score = log10((backw_mtx + get_mtx_index(0, 0, hmmp->nr_v))->prob);
  for (j = seq_len; j > 0; j--) {
      backward_score = backward_score + log10(*(forw_scale + j));
  }


  /* if needed run n-best */
  labels = (char*)malloc_or_die((seq_len + 1) * sizeof(char));
  one_best_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		 &one_best_mtx, &scaling_factors, use_labels, labels, multi_scoring_method);


  /* count helices + deallocate */
  j = 0;
  helices = 0;
  IN_HELIX = NO;
  while(1) {
      if(labels[j] == '\0') {
	  break;
      }
      else {
	  if((labels[j] == 'M') && (IN_HELIX == NO)) {
	      IN_HELIX = YES;
	      helices++;
	  } else if((labels[j] != 'M') && (IN_HELIX == YES)) {
	      IN_HELIX = NO;
	  }
	  j++;
      }
  }

  /* Check if we care where the helices are, if so, loop
     through and record their locations */

  if(hSites) {
    hSites->helix = (helix_sites *) malloc(sizeof (helix_sites) * helices);
    helixSite = hSites->helix;
    
    hSites->helix_count = helices;

    j = 0;
    helices = 0;
    IN_HELIX = NO;
    while(1) {
      if(labels[j] == '\0') {
	break;
      }
      else {
	if((labels[j] == 'M') && (IN_HELIX == NO)) {
	  IN_HELIX = YES;
	  helixSite->start = (j+1);
	  helices++;
	} else if((labels[j] != 'M') && (IN_HELIX == YES)) {
	  IN_HELIX = NO;
	  helixSite->end = (j);
	  helixSite++;
	}
	j++;
      }
    }
    
  }

  /* deallocate */
  free(labels);

  /* deallocate result matrix info */
    free(forw_mtx);
    free(forw_scale);
    free(backw_mtx);
    free(one_best_mtx);
    free(scaling_factors);

  return helices;

}


/*********************help functions***************************************/
double get_nullmodel_score_multi(struct letter_s *seq, struct letter_s *seq_2, struct letter_s *seq_3,
				 struct letter_s *seq_4, int seq_len, int multi_scoring_method)
{
  int avg_seq_len;
  double emiss_prob, emiss_prob_2, emiss_prob_3, emiss_prob_4;
  double e1, e2, e3, e4;
  double trans_prob;
  double null_model_score;
  int letter, letter_2, letter_3, letter_4;
  int l,i;
  double t_res;
  
  /* calculate null model score for seq */
  if(null_model.nr_alphabets < 0) {
    /* use default null model */
    emiss_prob = 1.0 / (double)hmm.a_size;
    if(hmm.nr_alphabets > 1) {
      emiss_prob_2 = 1.0 / (double)hmm.a_size_2;
    }
    if(hmm.nr_alphabets > 2) {
      emiss_prob_3 = 1.0 / (double)hmm.a_size_3;
    }
    if(hmm.nr_alphabets > 3) {
      emiss_prob_4 = 1.0 / (double)hmm.a_size_4;
    }
    trans_prob = (double)(seq_len)/(double)(seq_len + 1);
    if(multi_scoring_method == JOINT_PROB) {
      if(hmm.nr_alphabets == 1) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 2) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 3) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(emiss_prob_3) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 2) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(emiss_prob_3) + log10(emiss_prob_4)
				      + log10(trans_prob));
      }
    }
    else {
      /* use other multialpha scoring method, not implemented yet */
      printf("Error: only JOINT_PROB scoring is implemented\n");
      exit(0);
    }
  }
  else {
    /* use specified null model */
    null_model_score = 0.0;
    for(l = 0; l < seq_len; l++) {
      letter = get_alphabet_index(&seq[l], hmm.alphabet, hmm.a_size);
      if(hmm.nr_alphabets > 1) {
	letter_2 = get_alphabet_index(&seq_2[l], hmm.alphabet_2, hmm.a_size_2);
      }
      if(hmm.nr_alphabets > 2) {
	letter_3 = get_alphabet_index(&seq_3[l], hmm.alphabet_3, hmm.a_size_3);
      }
      if(hmm.nr_alphabets > 3) {
	letter_4 = get_alphabet_index(&seq_4[l], hmm.alphabet_4, hmm.a_size_4);
      }

      if(letter >= 0) {
	e1 = null_model.emissions[letter];
	//null_model_score += log10(null_model.emissions[letter]) + log10(null_model.trans_prob);
      }
      else {
	/* need replacement letters */
	letter = get_replacement_letter_index_multi(&seq[l], &replacement_letters, 1);
	if(letter >= 0) {
	  e1 = 0.0;
	  for(i = 0; i < hmm.a_size; i++) {
	    e1 += *(hmm.replacement_letters->probs_1 + get_mtx_index(letter, i, hmm.a_size)) * null_model.emissions[i];
	  }
	  //null_model_score += log10(t_res) + log10(null_model.trans_prob);
	}
	else {
	  printf("Could not find letter %s when scoring null model\n", &seq[l]);
	  return DEFAULT;
	}
      }
      if(hmm.nr_alphabets > 1) {
	if(letter_2 >= 0) {
	  e2 = null_model.emissions_2[letter_2];
	}
	else {
	  /* need replacement letters */
	  letter_2 = get_replacement_letter_index_multi(&seq_2[l], &replacement_letters, 2);
	  if(letter_2 >= 0) {
	    e2 = 0.0;
	    for(i = 0; i < hmm.a_size_2; i++) {
	      e2 += *(hmm.replacement_letters->probs_2 + get_mtx_index(letter_2, i, hmm.a_size_2)) * null_model.emissions_2[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq_2[l]);
	    return DEFAULT;
	  }
	}
      }
      if(hmm.nr_alphabets > 2) {
	if(letter_3 >= 0) {
	  e3 = null_model.emissions_3[letter_3];
	}
	else {
	  /* need replacement letters */
	  letter_3 = get_replacement_letter_index_multi(&seq_3[l], &replacement_letters, 3);
	  if(letter_3 >= 0) {
	    e3 = 0.0;
	    for(i = 0; i < hmm.a_size_3; i++) {
	      e3 += *(hmm.replacement_letters->probs_3 + get_mtx_index(letter_3, i, hmm.a_size_3)) * null_model.emissions_3[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq_3[l]);
	    return DEFAULT;
	  }
	}
      }
      if(hmm.nr_alphabets > 3) {
	if(letter_4 >= 0) {
	  e4 = null_model.emissions_4[letter_4];
	}
	else {
	  /* need replacement letters */
	  letter_4 = get_replacement_letter_index_multi(&seq_4[l], &replacement_letters, 4);
	  if(letter_4 >= 0) {
	    e4 = 0.0;
	    for(i = 0; i < hmm.a_size_4; i++) {
	      e4 += *(hmm.replacement_letters->probs_4 + get_mtx_index(letter_4, i, hmm.a_size_4)) * null_model.emissions_4[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq[l]);
	    return DEFAULT;
	  }
	}
      }

      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(e1) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(e2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(e3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(e4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }
    }
  }

  return null_model_score;
}






double get_nullmodel_score_msa_multi(int seq_len, int prf_mode, int use_prior,
				     int use_gap_shares, int normalize,
				     int scoring_method, int multi_scoring_method, double *aa_freqs, double *aa_freqs_2,
				     double *aa_freqs_3, double *aa_freqs_4)
{
  int avg_seq_len;
  double emiss_prob, emiss_prob_2, emiss_prob_3, emiss_prob_4;
  double trans_prob;
  double null_model_score;
  int letter;
  int k,l,p,m;
  double col_score, col_score_2, col_score_3, col_score_4;
  int using_default_null_model;
  double seq_normalizer, state_normalizer;

  using_default_null_model = NO;
  /* calculate null model score for seq */
  if(null_model.nr_alphabets < 0) {
    /* use default null model */
    emiss_prob = 1.0 / (double)hmm.a_size;
    if(hmm.nr_alphabets > 1) {
      emiss_prob_2 = 1.0 / (double)hmm.a_size_2;
    }
    if(hmm.nr_alphabets > 2) {
      emiss_prob_3 = 1.0 / (double)hmm.a_size_3;
    }
    if(hmm.nr_alphabets > 3) {
      emiss_prob_4 = 1.0 / (double)hmm.a_size_4;
    }
    trans_prob = (double)(seq_len)/(double)(seq_len + 1);
    null_model.a_size = hmm.a_size;
    null_model.emissions = (double*)malloc_or_die(null_model.a_size * sizeof(double));
    if(hmm.nr_alphabets > 1) {
      null_model.a_size_2 = hmm.a_size_2;
      null_model.emissions_2 = (double*)malloc_or_die(null_model.a_size_2 * sizeof(double));
    }
    if(hmm.nr_alphabets > 2) {
      null_model.a_size_3 = hmm.a_size_3;
      null_model.emissions_3 = (double*)malloc_or_die(null_model.a_size_3 * sizeof(double));
    }
    if(hmm.nr_alphabets > 3) {
      null_model.a_size_4 = hmm.a_size_4;
      null_model.emissions_4 = (double*)malloc_or_die(null_model.a_size_4 * sizeof(double));
    }
    for(k = 0; k < hmm.a_size; k++) {
      null_model.emissions[k] = emiss_prob;
      if(hmm.nr_alphabets > 1) {
	null_model.emissions_2[k] = emiss_prob_2;
      }
      if(hmm.nr_alphabets > 2) {
	null_model.emissions_3[k] = emiss_prob_3;
      }
      if(hmm.nr_alphabets > 3) {
	null_model.emissions_4[k] = emiss_prob_4;
      }
    }
    
    null_model.trans_prob = trans_prob;
    using_default_null_model = YES;
  }
  

  /* NOTE: must include other scoring methods here as well (copy from core_algorithms) */
  /* use specified null model */
  null_model_score = 0.0;
  if(scoring_method == DOT_PRODUCT) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_dp_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_dp_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_dp_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_dp_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == DOT_PRODUCT_PICASSO) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_dp_picasso_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_dp_picasso_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_dp_picasso_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_dp_picasso_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == PICASSO) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_picasso_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_picasso_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_picasso_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_picasso_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == PICASSO_SYM) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_picasso_sym_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_picasso_sym_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_picasso_sym_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_picasso_sym_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == SJOLANDER) {
    l = 0;
    while(1) {
       /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_sjolander_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
					   p, null_model.emissions,
					   0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_sjolander_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
						 p, null_model.emissions_2,
					       0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_sjolander_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					       p, null_model.emissions_3,
						 0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	  col_score_4 = get_sjolander_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
						 p, null_model.emissions_4,
						 0, normalize, msa_seq_infop->gap_shares);
	}
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	  if(hmm.nr_alphabets > 3) {
	    null_model_score += log10(col_score_4);
	  }
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == SJOLANDER_REVERSED) {
    l = 0;
    while(1) {
       /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_sjolander_reversed_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
						    p, null_model.emissions,
						    0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_sjolander_reversed_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
							p, null_model.emissions_2,
							0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_sjolander_reversed_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
							p, null_model.emissions_3,
							0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_sjolander_reversed_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
							p, null_model.emissions_4,
							0, normalize, msa_seq_infop->gap_shares);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	  if(hmm.nr_alphabets > 3) {
	    null_model_score += log10(col_score_4);
	  }
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  else if(scoring_method == SUBST_MTX_PRODUCT) {
    l = 0;
    while(1) {
      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_subst_mtx_product_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
						   p, null_model.emissions,
						   0, hmm.subst_mtx);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_subst_mtx_product_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
						       p, null_model.emissions_2,
						       0, hmm.subst_mtx_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_subst_mtx_product_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
						       p, null_model.emissions_3,
						       0, hmm.subst_mtx_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_subst_mtx_product_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
						       p, null_model.emissions_4,
						       0, hmm.subst_mtx_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }
      
      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }
    
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT) {
    printf("SMDP scoring not implemented yet\n");
    exit(0);
	   
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR) {
    printf("SMDPP scoring not implemented yet\n");
    exit(0);
  }
  
  if(using_default_null_model == YES) {
    free(null_model.emissions);
    null_model.a_size = -1;
  }
  return null_model_score;
}



void get_post_prob_matrix(double **ret_post_prob_matrixp, double forw_score, struct forward_s *forw_mtx,
			  struct backward_s *backw_mtx, int seq_len)
{
  int i,j;
  double *post_prob_matrix;
  double post_prob_score;
  
  /* must be freed by caller */
  post_prob_matrix = (double*)(malloc_or_die((seq_len+2) * hmm.nr_v * sizeof(double)));
  for(i = 0; i < seq_len + 2; i++) {
    for(j = 0; j < hmm.nr_v; j++) {
      //printf("forw_score = %f, forw_mtx = %f, backw_mtx = %f\n", forw_score,(forw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob,
      //     (backw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob);  
      post_prob_score = (forw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob * (backw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob /
	forw_score;
      *(post_prob_matrix + get_mtx_index(i,j,hmm.nr_v)) = post_prob_score;
    }
  }

  *ret_post_prob_matrixp = post_prob_matrix;
}

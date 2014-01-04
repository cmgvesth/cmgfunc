//#include "funcs.h"
#include "structs.h"

#include<stdio.h>

int main() {
  helix_sites *hSites;

  char *seq = "MMSSPHPMSSSRNTPLGVFYSLLACFYWGMVFVIPSMLGNFADLDIVLTRYSVFGICSLITILYKRSNIFKTVPFFLWKKGILWAFLINIAYYFGIAQAVRYSGSAVTVIIAGLAPIAILFYSNIKKKMLSYSFLLSMSGIIVVGIILSNVSEFQSESSSSLPLYLLGLGCVTAATSIWAGYIICNHDFLEQHSEISPDTWCHMLGISSLIICLPLIILGDTFGITHVTRNFLFHTPLSERCLFIVLCSAMGIFSSSRAIAAWNKASLHLSTALLGALLIFEPIFGWILSYLCKREMPSFQEGLGFFLMLGASLCLLLAQKKASEQETPSETLITTESDANJ";

  char *hmmfile = "./HMMS/S_TMHMM_0.92b.hmg";
  char *repfile = "./HMMS/replacement_letter_multi.rpl";
  char *path = "./HMMS/";
  int count = 0;
  int i = 0;
  int start, end;

  printf("Running get_helices\n");
  hSites = (helix_sites*) get_helices(seq, hmmfile, repfile, path);

  printf("Count of helicies\n");
  count = helix_count(hSites);
  printf("There are %d helices\n", count);

  for(i = 0; i < count; i++) {
    start = helix_start(hSites, i);
    end = helix_end(hSites, i);
    printf("Helix %d start: %d end: %d\n", i, start, end);
  }

  return 0;
}

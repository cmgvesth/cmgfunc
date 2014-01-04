#! /usr/freeware/bin/gawk -f

{
  if (FILENAME!=currfn) {
     i=0;
     j++;
  }
  line[j,++i] = $0;
  currfn = FILENAME;
}

END {	for (k=1; k<=i; k++) {
	    for (l=1; l<=j; l++)
	        printf(("\t" line[l,k]));
	    print "";
	}
}

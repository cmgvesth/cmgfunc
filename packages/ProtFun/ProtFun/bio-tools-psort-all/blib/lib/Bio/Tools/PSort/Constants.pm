package Bio::Tools::PSort::Constants;

require Exporter;

use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
use strict;

my @locs = qw(
CYTOPLASM
CYTOPLASMIC
CYTOPLASMICMEMBRANE 
PERIPLASMIC
OUTERMEMBRANE
EXTRACELLULAR
MITOCHONDRIAL
NUCLEAR
UNKNOWN
INNERMEMBRANE
MEMBRANE
CELLWALL
NONCYTOPLASMIC
OTHER);

my @secondlocs = qw(
HOSTASSOCIATED
HOSTCYTOPLASMIC
HOSTCYTOPLASMICMEMBRANE
FLAGELLAR
FIMBRIAL
T3SS
GASVESICLE
);

my @localizations = qw(
AllLocalizations
AllSecondarysLocalizations		       
);

my @common = qw(OK DONE ERROR);

%EXPORT_TAGS = (all           => [@common, @locs, @secondlocs, @localizations],
                codes         => \(@common),
                localizations => \@locs,
                secondlocalizations => \@secondlocs);

@EXPORT_OK   = (@common, @locs, @secondlocs, @localizations);

@ISA = qw(Exporter);

sub CYTOPLASM     { return 'Cytoplasmic'   }
sub CYTOPLASMIC     { return 'Cytoplasmic'   }
sub MEMBRANE { return 'CytoplasmicMembrane' }
sub INNERMEMBRANE { return 'CytoplasmicMembrane' }
sub EXTRACELLULAR { return 'Extracellular' }
sub CELLWALL { return 'Cellwall' }
sub UNKNOWN       { return 'Unknown'       }
sub OTHER {return 'Other' }
sub MITOCHONDRIAL { return 'Mitochondrial'} 
sub CYTOPLASMICMEMBRANE { return 'CytoplasmicMembrane' }
sub OUTERMEMBRANE { return 'OuterMembrane'}
sub PERIPLASMIC { return 'Periplasmic' }
sub NUCLEAR { return 'Nuclear' }
sub HOSTASSOCIATED { return 'HostAssociated' }
sub HOSTCYTOPLASMIC { return 'HostCytoplasmic' }
sub HOSTCYTOPLASMICMEMBRANE { return 'HostCytoplasmicMembrane' }
sub T3SS { return 'T3SS' }
sub FLAGELLAR { return 'Flagellar' }
sub FIMBRIAL { return 'Fimbrial' }
sub GASVESICLE { return 'GasVesicle' }
sub NONCYTOPLASMIC { return 'Non-Cytoplasmic' }

sub AllLocalizations { return @locs; }
sub AllSecondarysLocalizations { return @locs; }
1;

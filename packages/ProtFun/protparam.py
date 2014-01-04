from Bio.SeqUtils import ProtParam
from Bio.Data import IUPACData
from Bio import SeqIO
import yaml
import sys, os

#in vivo half-life dictionary
#Taken from http://web.expasy.org/protparam/protparam-doc.html
iv = {'A':(4.4*60, 20*60, 10*60), 'R':(60, 2, 2),
      'N':(1.4*60, 3, 10*60), 'D':(1.1*60,3,10*60),
      'C':(1.2*60, 20*60, 10*60), 'Q':(.8*60,10,10*60),
      'E':(60,30,10*60),'G':(30*60,20*60,10*60),
      'H':(3.5*60,10,10*60),'I':(20*60,30,10*60),
      'L':(5.5*60,3,2),'K':(1.3*60,3,2),
      'M':(30*60,20*60,10*60),'F':(1.1*60,3,2),
      'P':(20*60,20*60,0),'S':(1.9*60,20*60,10*60),
      'T':(7.2*60,20*60,10*60),'W':(2.8*60,3,2),
      'Y':(2.8*60,10,2),'V':(100*60,20*60,10*60)}

#Amino Acid Structures
aas = {'A':{'C':3,'H':7,'N':1,'O':2},  'R':{'C':6,'H':14,'N':4,'O':2},
       'N':{'C':4,'H':8,'N':2,'O':3},  'D':{'C':4,'H':7,'N':1,'O':4},
       'C':{'C':3,'H':7,'N':1,'O':2},  'E':{'C':5,'H':9,'N':1,'O':4},
       'Q':{'C':5,'H':10,'N':2,'O':3},  'G':{'C':2,'H':5,'N':1,'O':2},
       'H':{'C':6,'H':9,'N':3,'O':2},  'I':{'C':6,'H':13,'N':1,'O':2},
       'L':{'C':6,'H':13,'N':1,'O':2},  'K':{'C':6,'H':14,'N':2,'O':2},
       'M':{'C':5,'H':11,'N':1,'O':2,'S':1},  'F':{'C':9,'H':11,'N':1,'O':2},
       'P':{'C':5,'H':9,'N':1,'O':2},  'S':{'C':3,'H':7,'N':1,'O':3},
       'T':{'C':4,'H':9,'N':1,'O':3},  'W':{'C':11,'H':12,'N':2,'O':2},
       'Y':{'C':9,'H':11,'N':1,'O':3},  'V':{'C':5,'H':11,'N':1,'O':2}}
       
"""
Extends the ProteinAnalysis class to add some of the
missing ProtParam functions.
"""
class ProtFun(ProtParam.ProteinAnalysis):
    
    extinction_coefficient = None
    extinction_coefficient_no_cys = None
    aliphatic_index = None
    mol_weight = None
    atom_content = None
    atom_percent = None
    atom_length = None
    
    def molecular_weight(self):
        """
        Calculates the molecular weight. Override base method
        because it didn't have any caching, and calculated over
        the entire molecule one at a time...
        
        This method uses the cached molecule counts.
        """
        if( self.mol_weight ):
            return self.mol_weight
        
        counts = self.count_amino_acids()
        total_weight = 18.02 # add just one water molecule for the whole sequence
        for i in IUPACData.protein_weights:
            # remove a water molecule from the amino acid weight
            total_weight += counts[i] * (IUPACData.protein_weights[i] - 18.02)
        self.mol_weight = total_weight    
        return self.mol_weight
        
    def get_extinction_coefficient(self):
        """
        Calculates the extinction coefficient according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        
        (uses max cysteine pairs)
        """
        if ( self.extinction_coefficient ):
            return self.extinction_coefficient
        counts = self.count_amino_acids()
        tyr = 1490.0 * counts['Y']
        trp = 5500.0 * counts['W']
        cys = 125.0 * counts['C'] / 2.0 
        self.extinction_coefficient = tyr + trp + cys
        return self.extinction_coefficient
    
    def get_extinction_coefficient_no_cys(self):
        """
        Calculates the extinction coefficient according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        
        (uses no cysteine pairs)
        """
        if ( self.extinction_coefficient_no_cys ):
            return self.extinction_coefficient_no_cys
        counts = self.count_amino_acids()
        tyr = 1490.0 * counts['Y']
        trp = 5500.0 * counts['W']
        self.extinction_coefficient = tyr + trp
        return self.extinction_coefficient
        
    def get_absorbance(self):
        """
        Calculates the absorbance according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        
        (uses max cysteine pairs version of extinction)
        """
        return (self.get_extinction_coefficient() / self.molecular_weight())
        
    def get_absorbance_no_cys(self):
        """
        Calculates the absorbance according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        
        (uses no cysteine pairs version of extinction)
        """
        return (self.get_extinction_coefficient_no_cys() / self.molecular_weight())
        
    def get_in_vivo_half_life(self):
        """
        Calculates the in vivo half-life according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        
        Based on the 'N-end rule'
        """
        return iv.get( self.sequence[0], (0,0,0) )
        
    def get_aliphatic_index(self):
        """
        Calculates the aliphatic index according to:
        
        http://web.expasy.org/protparam/protparam-doc.html
        """
        if(self.aliphatic_index):
            return self.aliphatic_index
        mp = self.get_amino_acids_percent()
        a = 2.9
        b = 3.9
        self.aliphatic_index = mp['A'] * 100 + (a * mp['V'] * 100) + (b * (mp['I'] * 100 + mp['L'] * 100))
        return self.aliphatic_index
    
    def count_atoms(self):
        """
        Count atoms. Returns a dictionary.
        
        Counts the number of times each atom occurs in the protein
        sequence. Calculated once
        """
        if(self.atom_content):
            return self.atom_content
        atom_dict = dict([(k, 0) for k in IUPACData.atom_weights])
        counts = self.count_amino_acids()
        for aa in counts:
            for atom in aas[aa]:
                atom_dict[atom] += counts[aa] * aas[aa][atom]
        self.atom_content = atom_dict
        return self.atom_content
    
    def get_total_molecules(self):
        return self.length
    
    def get_total_atoms(self):
        """
        Get total number of atoms in the molecule (could be
        sped up with it's own translation table maybe -- just
        sum the amino * atomcount for each amino instead
        of counting the individual atoms and summing them)
        
        This is staying for now at least
        """
        if(self.atom_length):
            return self.atom_length
        counts = self.count_atoms()
        self.atom_length = sum(counts.values())
        return self.atom_length

    def get_number_negative(self):
        counts = self.count_amino_acids()
        return sum( [counts[k] for k in 'DE'] )
        
    def get_weight_negative(self):
        counts = self.count_amino_acids()
        return sum( [(counts[k]*IUPACData.protein_weights[k]) for k in 'DE'] )
    
    def get_number_positive(self):
        counts = self.count_amino_acids()
        return sum( [counts[k] for k in 'RHK'] )

    def get_weight_positive(self):
        counts = self.count_amino_acids()
        return sum( [(counts[k]*IUPACData.protein_weights[k]) for k in 'RHK'] )

    def get_number_aromatic(self):
        counts = self.count_amino_acids()
        return sum( [counts[k] for k in 'FYW'] )

    def get_weight_aromatic(self):
        counts = self.count_amino_acids()
        return sum( [(counts[k]*IUPACData.protein_weights[k]) for k in 'FYW'] )

    def get_number_nonpolar(self):
        counts = self.count_amino_acids()
        return sum( [counts[k] for k in 'GAVLIMFPCW'] )

    def get_weight_nonpolar(self):
        counts = self.count_amino_acids()
        return sum( [(counts[k]*IUPACData.protein_weights[k]) for k in 'GAVLIMFPCW'] )
        
    def get_number_uncharged(self):
        counts = self.count_amino_acids()
        return sum( [counts[k] for k in 'STNQY'] )

    def get_weight_uncharged(self):
        counts = self.count_amino_acids()
        return sum( [(counts[k]*IUPACData.protein_weights[k]) for k in 'STNQY'] )

def process_sequence( sequence ):
    pf = ProtFun( str(sequence.seq.upper()) )
    aa_dict = pf.count_amino_acids()
    aa_percnet = pf.get_amino_acids_percent()
    mw = pf.molecular_weight()
    ar = pf.aromaticity()
    ii = pf.instability_index()
    flex = pf.flexibility()
    average_flex = sum(flex)/len(flex)
    gravy = pf.gravy()
    iso_p = pf.isoelectric_point()
    seq_struc2 = pf.secondary_structure_fraction()
    seq_struc = {'helix':seq_struc2[0],'turn':seq_struc2[1],'sheet':seq_struc2[2]}
    ext_c_max = pf.get_extinction_coefficient()
    ext_c_min = pf.get_extinction_coefficient_no_cys()
    absorb_max = pf.get_absorbance()
    absorb_min = pf.get_absorbance_no_cys()
    half_life2 = pf.get_in_vivo_half_life()
    half_life = {'human':half_life2[0],'yeast':half_life2[1],'e.coli':half_life2[2]}
    aliph = pf.get_aliphatic_index()
    atoms = pf.count_atoms()
    total_molecules = pf.get_total_molecules()
    total_atoms = pf.get_total_atoms()
    number_neg = pf.get_number_negative()
    neg_weight = pf.get_weight_negative()
    number_pos = pf.get_number_positive()
    pos_weight = pf.get_weight_positive()
    num_arom = pf.get_number_aromatic()
    weight_arom = pf.get_weight_aromatic()
    num_nonpolar = pf.get_number_nonpolar()
    weight_nonpolar = pf.get_weight_nonpolar()
    num_uncharged = pf.get_number_uncharged()
    weight_uncharged = pf.get_weight_uncharged()
    seq_length = len(sequence)
    return {
            #'amino_acid_usage':aa_dict,
            'amino_acid_percentages':aa_percnet,
            'molecular_weight':mw,
            'aromaticity':ar,
            'instability_index':ii,
            #'flexibility_index':flex,
            'flexibility_average':average_flex,
            'gravy_coeff':gravy,
            'isoelectric_point':iso_p,
            'secondary_structure':seq_struc,
            'extinction_coeff_max':ext_c_max,
            'extinction_coeff_min':ext_c_min,
            'absorbtion_max':absorb_max,
            'absorbtion_min':absorb_min,
            'in_vivo_halflife':half_life,
            'aliphatic_index':aliph,
            #'atom_usage':atoms,
            'total_molecules':total_molecules,
            'total_atoms':total_atoms,
            'number_negative':number_neg,
            'weight_negative':neg_weight,
            'number_positive':number_pos,
            'weight_positive':pos_weight,
            'number_aromatic':num_arom,
            'weight_aromatic':weight_arom,
            'number_nonpolar':num_nonpolar,
            'weight_nonpolar':weight_nonpolar,
            'number_uncharged':num_uncharged,
            'weight_uncharged':weight_uncharged,
            'sequence_length':seq_length,
           }
           
def format_results( results_dict ):
    return yaml.dump(results_dict, width=sys.maxint)
        
def usage():
    print "usage:"
    print "\tpython protparam.py <input.fsa> [output]"
    print "Description:"
    print "  Calculates Expasy ProtParam calculations on input data"
    print "  as well as Amino Acid Comp (replaces AACOMP tool and "
    print "  Expasy ProtParam [which is web only])"
    print "Arguments:"
    print "  <input.fsa> - the input file to run calculations on"
    print "  [output.fsa] - the output file to save results to, if"
    print "                 not specified, results are saved to:"
    print "                 <input.fsa>.protparam"

def main(args):
    if( len(args) != 2 and len(args) != 3 ):
        usage()
        exit( 1 )
    (root, ext) = os.path.splitext( args[1] )
    #save_as = args[1] + ".protparam"
    #if len(args) == 3:
    #    save_as = args[2]
    
    print "# Reading input file: '%s'" % args[1]   
    if( args[1] is "-" ):
        input_handle = sys.stdin
    else:
        input_handle = open( args[1], 'r' )

    fasta = SeqIO.parse( input_handle, 'fasta' )
    print "# Success"
    
    #print "# Will save results in: '%s'" % save_as
    #if os.path.exists( save_as ):
    #    print("# Output file '%s' already exists, will overwrite." % save_as)

    #if( args[1] is "-" ):
    #    save_file = open( "tmp_protparam", 'w')
    #else:
    #    save_file = open( save_as, 'w')

    count = 0    
    dict_results = dict()

    for sequence in fasta:
        count = count + 1
	
        dict_results[sequence.id]=process_sequence( sequence )
        if (count % 500 == 0):
            print "# Finished processing count sequences %d" % count 
            
        result_string = format_results(dict_results)
        
        if( args[1] is "-" ):
            print result_string
        else:
            #print >> save_file, result_string
            print result_string
        dict_results = dict()
    print "# Finish"
    #save_file.close()
    return 0

if __name__ == '__main__':
    try:
        main( sys.argv )
    except KeyboardInterrupt, e:
        pass


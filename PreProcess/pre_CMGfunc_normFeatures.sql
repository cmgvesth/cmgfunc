use proteinfunction;
update feat set nweight_aromatic=weight_aromatic/molecular_weight;
update feat set nweight_negative=weight_negative/molecular_weight;
update feat set nweight_nonpolar=weight_nonpolar/molecular_weight;
update feat set nweight_positive=weight_positive/molecular_weight; 
update feat set nweight_uncharged=weight_uncharged/molecular_weight;

update feat set nin_vivo_halflife_ecoli=in_vivo_halflife_ecoli/6000;
update feat set nisoelectric_point=isoelectric_point/14;
update feat set nmolecular_weight=molecular_weight/4713080;
update feat set nsequence_length=sequence_length/40000;
update feat set nExpAA=ExpAA/5000;
update feat set nFirst60=First60/60;
update feat set nPredHel=PredHel/300;

update feat set nCytoplasmic_Score_n=Cytoplasmic_Score_n/10;
update feat set nCytoplasmicMembrane_Score_n=CytoplasmicMembrane_Score_n/10;
update feat set nPeriplasmic_Score_n=Periplasmic_Score_n /10;
update feat set nOuterMembrane_Score_n=OuterMembrane_Score_n/10;
update feat set nExtracellular_Score_n= Extracellular_Score_n/10;
update feat set nFinal_Score_n=Final_Score_n/10;
update feat set nCytoplasmic_Score_p=Cytoplasmic_Score_p/10;
update feat set nCytoplasmicMembrane_Score_p=CytoplasmicMembrane_Score_p/10;
update feat set nCellwall_Score_p=Cellwall_Score_p/10;
update feat set nExtracellular_Score_p=Extracellular_Score_p/10;
update feat set nFinal_Score_p=Final_Score_p/10;


update feat set nseq_low_total=seq_low_total/sequence_length;
update feat set nseq_low_avg=seq_low_avg/sequence_length;
update feat set nseq_low_coverage=seq_low_coverage/sequence_length;
update feat set nseq_high_total=seq_high_total/sequence_length;
update feat set nseq_high_avg=seq_high_avg/sequence_length;
update feat set nseq_high_coverage=seq_high_coverage/100;
update feat set naromaticity=aromaticity/sequence_length;
update feat set ntotal_atoms=total_atoms/sequence_length;
update feat set ninstability_index=instability_index/250;

update feat set nflexibility_average=flexibility_average/1.1;
update feat set naliphatic_index=aliphatic_index/400;
update feat set nextinction_coeff_min=(extinction_coeff_min)/1800;
update feat set nextinction_coeff_max=(extinction_coeff_max)/1800;

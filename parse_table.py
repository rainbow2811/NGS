#!/usr/bin/python

import os
import sys
import getopt
import re

correct_input_header = '\t'.join(\
['CHROM','POS','ID','REF','ALT','GEN[0].GT','GEN[1].GT','GEN[2].GT','GEN[0].AD','GEN[1].AD','GEN[2].AD',\
'ANN[*].GENE','ANN[*].FEATUREID','ANN[*].FEATURE','ANN[*].BIOTYPE','ANN[*].EFFECT','ANN[*].IMPACT','ANN[*].HGVS_C','ANN[*].HGVS_P','ANN[*].ERRORS',\
'CLNSIG','dbNSFP_SIFT_pred'])

output_header = '\t'.join(\
['#Chrom','Pos','rs_ID','Ref','Alt','NA12878_genotype','NA12891_genotype','NA12892_genotype','NA12878_refDP','NA12878_altDP','NA12891_refDP','NA12891_altDP','NA12892_refDP','NA12892_altDP',\
'Gene','NM_id','Effect','Impact','HGVS_c','HGVS_p','Clinvar_significance','SIFT_prediction'])


def usage(err):
	title = '''
 usage: %s [ table file from SnpSift-ExtractFields ]
	'''%sys.argv[0]
	sys.stderr.write(title)
	sys.stderr.write(str(err)+'\n')

if __name__=="__main__":
	if len(sys.argv) < 2:
		usage('\n');sys.exit(2)

	input_table = arg = sys.argv[1]

	if not os.path.isfile(input_table):
		usage('\nfile "%s" not found'%input_table);sys.exit(2)

	finp = open(input_table)
	header = finp.readline().strip()
	if not header == correct_input_header:
		err = '''
 wrong header found!
 <correct header>

%s'''%('\n'.join(correct_input_header.split('\t')))
		finp.close()
		usage(err)
		sys.exit(2)


        print output_header
	for line in finp:
		token = line.strip().split('\t')
		chrom,pos,id,ref,alt,na12878_gt,na12891_gt,na12892_gt,na12878_ad,na12891_ad,na12892_ad,gene_sum,featureid_sum,feature_sum,biotype_sum,effect_sum,impact_sum,hgvsc_sum,hgvsp_sum,error_sum,clnsig,sift_pred = token
		snpeff_set = dict()
		#--tri-allele filter --
		if ',' in alt:
			continue
		#--no-call include filter --

		if './.' in (na12878_gt,na12891_gt,na12892_gt):
			continue
		#----------------------
		try:
			na12878_refdp,na12878_altdp = na12878_ad.split(',')
		except:
			continue
		try:
			na12891_refdp,na12891_altdp = na12891_ad.split(',')
		except:
			continue
		try:
			na12892_refdp,na12892_altdp = na12892_ad.split(',')
		except:
			continue

		#--snpeff parsing by transcript id --
		gene_set = gene_sum.split(',')
		featureid_set = featureid_sum.split(',')
		feature_set = feature_sum.split(',')
		biotype_set = biotype_sum.split(',')
		effect_set = effect_sum.split(',')
		impact_set = impact_sum.split(',')
		hgvsc_set = hgvsc_sum.split(',')
		hgvsp_set = hgvsp_sum.split(',')
		error_set = error_sum.split(',')

                print '#'

		for i in range(len(gene_set)):
                        print i
			gene = gene_set[i]
			featureid = featureid_set[i]
			feature = feature_set[i]
			biotype = biotype_set[i]
			effect = effect_set[i]
			impact = impact_set[i]
			hgvsc = hgvsc_set[i]
			hgvsp = hgvsp_set[i]
			error = error_set[i]
			#----filter--------
			if re.search('^WARNING',error):  # abnormal transcript
				continue
			if not biotype == 'protein_coding': # intergenic
				continue
			if not feature == 'transcript': # only transcript
				continue
			#-------------------
			snpeff_set[featureid] = (gene,effect,impact,hgvsc,hgvsp)

                print '#'
		#---filter ----
		if snpeff_set == {}: # 0 proper transcript
			continue
		selected_featureid = sorted(snpeff_set.keys(),key=lambda x:float(x.replace('NM_','').split('.')[0]))[0]
		selected_gene,selected_effect,selected_impact,selected_hgvsc,selected_hgvsp = snpeff_set[selected_featureid]

		text_set = (chrom,pos,id,ref,alt,na12878_gt,na12891_gt,na12892_gt,na12878_refdp,na12878_altdp,na12891_refdp,na12891_altdp,na12892_refdp,na12892_altdp,\
selected_gene,selected_featureid,selected_effect,selected_impact,selected_hgvsc,selected_hgvsp,clnsig,sift_pred)

	#	print '\t'.join(['%s']*len(text_set))%text_set

	finp.close()
	


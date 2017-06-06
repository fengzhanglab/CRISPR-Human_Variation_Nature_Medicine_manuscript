#David Scott, MIT, 2016
import sys
import csv
import gzip
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import vcf_obj_v0_6 as vcf
import geno_gen_obj_v0_6_1 as gg

import vmb

def revcomp(seq):
	seq = seq.lower()
	seq = seq.replace('a','T')
	seq = seq.replace('t','A')
	seq = seq.replace('c','G')
	seq = seq.replace('g','C')
	seq = seq.upper()
	return seq[::-1]

def nt_query_gg_all(gg_obj_file,out_file):

	with open(gg_obj_file,'rb') as gg_handle:
		gg_obj = pickle.load(gg_handle)
		
		rf_split = gg_obj_file.split('_')
		chrid = rf_split[len(rf_split)-2]
		offset = int(rf_split[len(rf_split)-1].replace('.gg',''))
		print offset

		with open(out_file,'wb') as csvfile:
			mywriter = csv.writer(csvfile, delimiter=',')	

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'

			#skip overlap if offset > 0
			if offset == 0:
				ref_lb = 0
			else:
				ref_lb = 1000
			ref_ub = len(gg_obj.ref)
							
			out = gg_obj.get_ntvar_counts(ref_lb,ref_ub)

			mywriter.writerow([chrid,offset]+out['A']+out['C']+out['G']+out['T'])						

if __name__ == '__main__':
	gg_obj_file = sys.argv[1]
	out_file = sys.argv[2]
	# threshold = sys.argv[5]
	# search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),threshold)
	search_seqs = nt_query_gg_all(gg_obj_file,out_file)	



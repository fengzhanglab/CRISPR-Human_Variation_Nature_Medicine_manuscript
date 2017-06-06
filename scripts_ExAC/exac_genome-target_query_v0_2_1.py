#David Scott, MIT, 2016
import sys
import csv
import gzip
import numpy
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import vcf_obj_v0_6 as vcf
import gff3_obj_v0_6 as gff3
import exac_gen_obj_v0_6_1 as gg

import vmb

def revcomp(seq):
	seq = seq.lower()
	seq = seq.replace('a','T')
	seq = seq.replace('t','A')
	seq = seq.replace('c','G')
	seq = seq.replace('g','C')
	seq = seq.upper()
	return seq[::-1]

#def PAM_compile_gg(gg_obj_file,targ_len,PAM_orient,PAM_list,threshold):
def target_query_gg_all(gg_obj_file,gff3_file,PAM_orient,PAM_list,target_len,out_file):
	PAM_len = len(PAM_list[0])

	PAM_dict = {}
	for pam in PAM_list:
		PAM_dict[str(pam)] = pam
	PAM_dict_rc = {}
	for pam in PAM_list:
		PAM_dict_rc[revcomp(str(pam))] = pam

	gff3_pct = {}
	gff3_file_list = gff3_file.replace('.gff3','').split('_')
	offset = int(gff3_file_list[len(gff3_file_list)-1])
	with open(gff3_file,'rb') as fin:
		for line in fin:
			if line[0] == '#':
				continue			
			row_obj = gff3.GFF3_row(line)
			# print str(row_obj)
			#protein coding transcript (PCT) and exon or UTR (EU)
			if (row_obj.PCT and row_obj.EU):
				if row_obj.INFO['transcript_id'] not in gff3_pct:
					gff3_pct[row_obj.INFO['transcript_id']] = gff3.PCT(row_obj)
				else:
					gff3_pct[row_obj.INFO['transcript_id']].add_element(row_obj)					

	gff3_lines = []
	for pct_id in gff3_pct:
		gff3_lines+=gff3_pct[pct_id].get_coding()
		# print str(gff3_pct[pct_id])

	# print gff3_lines

	with open(gg_obj_file,'rb') as gg_handle:
		gcr_obj = gg.Geno_CRISPR(pickle.load(gg_handle))		
		with open(out_file,'wb') as csvfile:
			mywriter = csv.writer(csvfile, delimiter=',')	

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'

			for line_els in gff3_lines:
				# if line[0] == '#':
				# 	continue

				# line_els = line.strip().split()
				if line_els[2] == 'exon':
					ref_lb = int(line_els[3])-offset-1
					ref_ub = int(line_els[4])-offset-1
					if ref_lb < 0:
						ref_lb = 0
					if ref_ub > gcr_obj.len:
						ref_ub = gcr_obj.len

					af_bins = ([[0.00001,[0]*(target_len+PAM_len)],[0.0001,[0]*(target_len+PAM_len)],[0.001,[0]*(target_len+PAM_len)],
						[0.01,[0]*(target_len+PAM_len)],[0.1,[0]*(target_len+PAM_len)],[1,[0]*(target_len+PAM_len)]])	
					hetf_bins = ([[0.00001,[0]*(target_len+PAM_len)],[0.0001,[0]*(target_len+PAM_len)],[0.001,[0]*(target_len+PAM_len)],
						[0.01,[0]*(target_len+PAM_len)],[0.1,[0]*(target_len+PAM_len)],[1,[0]*(target_len+PAM_len)]])										
					[target_inds, target_vars] = gcr_obj.get_var_targets_del_pams(PAM_orient,PAM_list,target_len,ref_lb,ref_ub)
					n_targ = len(target_inds)
					n_var = len(target_vars)
					for loc in target_inds:
						max_var = 0
						max_var_els = []							

						#change 161006 taking max only
						if str(loc) in target_vars:
							# print [[el[0],el[1].get_af_adj()] for el in target_vars[str(loc)]]
							for el in target_vars[str(loc)]:
								if el[1].get_af_adj() > max_var:
									max_var = el[1].get_af_adj()
									max_var_els = [el]
								elif el[1].get_af_adj() == max_var:
									max_var_els+=[el]
							
							#pam start is index 0 at this point
							#for equal af, take closest to PAM

							# print '##########'
							# print [[el[0],el[1].get_af_adj()] for el in max_var_els]

							max_var_els.sort(key=lambda x: x[0])

							# print '###'
							# print [[el[0],el[1].get_af_adj()] for el in max_var_els]
							# print '##########'

							el = max_var_els[0]
							if PAM_orient == 'R':
								#orient PAM at right hand side
								ind = (target_len+PAM_len)-(el[0]+1)
							elif PAM_orient == 'L':
								ind = el[0]

							af = el[1].get_af_adj()
							hetf = el[1].get_hetf()

							#change 161006 taking max only
							varct = 1
							if af < af_bins[0][0]:
								af_bins[0][1][ind]+=1
								hetf_bins[0][1][ind]+=hetf/varct
							elif af < af_bins[1][0]:
								af_bins[1][1][ind]+=1
								hetf_bins[1][1][ind]+=hetf/varct
							elif af < af_bins[2][0]:
								af_bins[2][1][ind]+=1
								hetf_bins[2][1][ind]+=hetf/varct
							elif af < af_bins[3][0]:
								af_bins[3][1][ind]+=1
								hetf_bins[3][1][ind]+=hetf/varct
							elif af < af_bins[4][0]:
								af_bins[4][1][ind]+=1
								hetf_bins[4][1][ind]+=hetf/varct
							elif af <= af_bins[5][0]:
								af_bins[5][1][ind]+=1
								hetf_bins[5][1][ind]+=hetf/varct
						
					mywriter.writerow(line_els+[el[1] for el in af_bins]+
						[el[1] for el in hetf_bins]+[n_var,n_targ,'var'])																							

if __name__ == '__main__':
	gg_obj_file = sys.argv[1]
	gff3_file = sys.argv[2]
	target_len = int(sys.argv[3])
	PAM_orient = sys.argv[4]
	PAM_list = sys.argv[5]
	out_file = sys.argv[6]
	# threshold = sys.argv[5]
	# search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),threshold)
	search_seqs = target_query_gg_all(gg_obj_file,gff3_file,PAM_orient,PAM_list.split(','),target_len,out_file)	



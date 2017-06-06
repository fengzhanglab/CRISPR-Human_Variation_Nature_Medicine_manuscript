#David Scott, MIT, 2016
import sys
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

#def PAM_compile_gg(gg_obj_file,targ_len,PAM_orient,PAM_list,threshold):
def PAM_compile_gg(gg_obj_file,targ_len,PAM_orient,PAM_list,out_file):
	PAM_len = len(PAM_list[0])

	PAM_dict = {}
	for pam in PAM_list:
		PAM_dict[str(pam)] = pam
	PAM_dict_rc = {}
	for pam in PAM_list:
		PAM_dict_rc[revcomp(str(pam))] = pam

	with open(gg_obj_file,'rb') as gg_handle:
		with open(out_file,'wb') as fout:

			ggs_obj = gg.Geno_slice(pickle.load(gg_handle))	
			
			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'
			
			print ggs_obj.ref
			# print ggs_obj.var_track

			len_slice = targ_len+(2*PAM_len)

			ggs_obj.slice(0, len_slice)

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'

			ct = 0
			ct2 = 0
			search_seqs = []
			for i in range(0,(len(ggs_obj.ref)-len_slice)):
				
				gg_seqs = ggs_obj.compile_geno_seqs()
				gg_seqs += [gg.Geno_seq(ggs_obj.ref[i:(i+len_slice)],ggs_obj.ind-ggs_obj.len,-1,['0'])]
				for gg_seq in gg_seqs:

					if PAM_orient == 'R':
						if gg_seq.seq[0:PAM_len] in PAM_dict_rc:
							# if (levenshtein(revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)]),
							# 	target) <= threshold):
							# fout.write('>'+gg_seq.seq[0:PAM_len]+
							# 	'|'+PAM_orient+'|'+'BS'+'|'+str(gg_seq)+'\n')
							# fout.write(revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)])+'\n')
							search_seqs.append([revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)]),
								gg_seq.seq[0:PAM_len],PAM_orient,'BS',gg_seq.ind,gg_seq.allele,gg_seq.names])
							ct2+=1

						if gg_seq.seq[(len_slice-PAM_len):len_slice] in PAM_dict:
							# if (levenshtein(gg_seq.seq[PAM_len:(PAM_len+targ_len)],
							# 	target) <= threshold):
							# fout.write('>'+gg_seq.seq[(len_slice-PAM_len):len_slice]+
							# 	'|'+PAM_orient+'|'+'TS'+'|'+str(gg_seq)+'\n')
							# fout.write(gg_seq.seq[PAM_len:(PAM_len+targ_len)]+'\n')
							search_seqs.append([gg_seq.seq[PAM_len:(PAM_len+targ_len)],
								gg_seq.seq[(len_slice-PAM_len):len_slice],PAM_orient,
								'TS',gg_seq.ind,gg_seq.allele,gg_seq.names])							
							ct2+=1

					elif PAM_orient == 'L':
						if gg_seq.seq[0:PAM_len] in PAM_dict:
							# if (levenshtein(gg_seq.seq[PAM_len:(PAM_len+targ_len)],
							# 	target) <= threshold):							
							# fout.write('>'+gg_seq.seq[0:PAM_len]+
							# 	'|'+PAM_orient+'|'+'TS'+'|'+str(gg_seq)+'\n')
							# fout.write(gg_seq.seq[PAM_len:(PAM_len+targ_len)]+'\n')
							search_seqs.append([gg_seq.seq[PAM_len:(PAM_len+targ_len)],
								gg_seq.seq[0:PAM_len],PAM_orient,'TS',gg_seq.ind,gg_seq.allele,gg_seq.names])							
							ct2+=1

						if gg_seq.seq[(len_slice-PAM_len):len_slice] in PAM_dict_rc:
							# if (levenshtein(revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)]),
							# 	target) <= threshold):								
							# fout.write('>'+gg_seq.seq[(len_slice-PAM_len):len_slice]+
							# 	'|'+PAM_orient+'|'+'BS'+'|'+str(gg_seq)+'\n')
							# fout.write(revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)])+'\n')
							search_seqs.append([revcomp(gg_seq.seq[PAM_len:(PAM_len+targ_len)]),
								gg_seq.seq[(len_slice-PAM_len):len_slice],PAM_orient,
								'BS',gg_seq.ind,gg_seq.allele,gg_seq.names])							
							ct2+=1
					ct+=1

				ggs_obj.ref_base_shift() 

				if (i%100000) == 0:
					print '####################'
					print 'vmb memory:' + str(vmb.memory())
					print 'vmb resident: ' + str(vmb.resident())
					print '####################'
					print str(float(i)/(float(len(ggs_obj.ref))-float(len_slice)))
			print ct2

			pickle.dump(search_seqs,fout)

if __name__ == '__main__':
	gg_obj_file = sys.argv[1]
	target_len = int(sys.argv[2])
	PAM_orient = sys.argv[3]
	PAM_list = sys.argv[4]
	out_file = sys.argv[5]
	# threshold = sys.argv[5]
	# search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),threshold)
	search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),out_file)	



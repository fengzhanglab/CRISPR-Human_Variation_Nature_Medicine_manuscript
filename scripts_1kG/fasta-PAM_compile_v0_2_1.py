#David Scott, MIT, 2016
import sys
import gzip
import pickle
from Bio.Seq import Seq
from Bio import SeqIO

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
def PAM_compile_fasta(in_file,targ_len,PAM_orient,PAM_list,out_file):
	PAM_len = len(PAM_list[0])

	PAM_dict = {}
	for pam in PAM_list:
		PAM_dict[str(pam)] = pam
	PAM_dict_rc = {}
	for pam in PAM_list:
		PAM_dict_rc[revcomp(str(pam))] = pam

	len_slice = targ_len+(2*PAM_len)

	ct = 0
	ct2 = 0
	search_seqs = []
	with open(in_file,'r') as fin:
		with open(out_file,'wb') as fout:	
			for record in SeqIO.parse(fin, "fasta"):
				seq = str(record.seq)
				name = str(record.name)

				for i in range(0,(len(seq)-len_slice)):
					seq_temp = seq[i:(i+len_slice)]

					if PAM_orient == 'R':
						if seq_temp[0:PAM_len] in PAM_dict_rc:
							fout.write('>'+name+'|'+str(i)+'|'+revcomp(seq_temp[0:PAM_len])+
								'|'+PAM_orient+'|'+'BS'+'\n')
							fout.write(revcomp(seq_temp[PAM_len:(PAM_len+targ_len)])+'\n')
							# search_seqs.append([revcomp(seq_temp[PAM_len:(PAM_len+targ_len)]),
							# 	seq_temp[0:PAM_len],PAM_orient,'BS'])
							ct2+=1

						if seq_temp[(len_slice-PAM_len):len_slice] in PAM_dict:
							fout.write('>'+name+'|'+str(i)+'|'+seq_temp[(len_slice-PAM_len):len_slice]+
								'|'+PAM_orient+'|'+'TS'+'\n')
							fout.write(seq_temp[PAM_len:(PAM_len+targ_len)]+'\n')
							# search_seqs.append([seq_temp[PAM_len:(PAM_len+targ_len)],
							# 	seq_temp[(len_slice-PAM_len):len_slice],PAM_orient,'TS'])							
							ct2+=1

					elif PAM_orient == 'L':
						if seq_temp[0:PAM_len] in PAM_dict:						
							fout.write('>'+name+'|'+str(i)+'|'+seq_temp[0:PAM_len]+
								'|'+PAM_orient+'|'+'TS'+'\n')
							fout.write(seq_temp[PAM_len:(PAM_len+targ_len)]+'\n')
							# search_seqs.append([seq_temp[PAM_len:(PAM_len+targ_len)],
							# 	seq_temp[0:PAM_len],PAM_orient,'TS'])							
							ct2+=1

						if seq_temp[(len_slice-PAM_len):len_slice] in PAM_dict_rc:								
							fout.write('>'+name+'|'+str(i)+'|'+revcomp(seq_temp[(len_slice-PAM_len):len_slice])+
								'|'+PAM_orient+'|'+'BS'+'\n')
							fout.write(revcomp(seq_temp[PAM_len:(PAM_len+targ_len)])+'\n')
							# search_seqs.append([revcomp(seq_temp[PAM_len:(PAM_len+targ_len)]),
							# 	seq_temp[(len_slice-PAM_len):len_slice],PAM_orient,'BS'])							
							ct2+=1
						ct+=1


					if (i%100000) == 0:
						print '####################'
						print 'vmb memory:' + str(vmb.memory())
						print 'vmb resident: ' + str(vmb.resident())
						print '####################'
						print str(float(i)/(float(len(seq))-float(len_slice)))
				print ct2

if __name__ == '__main__':
	in_file = sys.argv[1]
	target_len = int(sys.argv[2])
	PAM_orient = sys.argv[3]
	PAM_list = sys.argv[4]
	out_file = sys.argv[5]
	# threshold = sys.argv[5]
	# search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),threshold)
	search_seqs = PAM_compile_fasta(in_file,target_len,PAM_orient,PAM_list.split(','),out_file)	



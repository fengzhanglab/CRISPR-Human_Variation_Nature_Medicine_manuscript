#David Scott, MIT, 2016
import sys
import csv
import gzip
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import vcf_obj_v0_6 as vcf
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
def nt_query_gg_all(gg_obj_file,gff3_file,out_file):

	with open(gg_obj_file,'rb') as gg_handle:
		gg_obj = pickle.load(gg_handle)		
		with open(out_file,'wb') as csvfile:
			mywriter = csv.writer(csvfile, delimiter=',')	

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'

			gff3_file_list = gff3_file.replace('.gff3','').split('_')
			offset = int(gff3_file_list[len(gff3_file_list)-1])
			with open(gff3_file,'rb') as fin:
				
				for line in fin:
					if line[0] == '#':
						continue

					line_els = line.strip().split()
					if line_els[2] == 'exon':
						ref_lb = int(line_els[3])-offset-1
						ref_ub = int(line_els[4])-offset-1
						if ref_lb < 0:
							ref_lb = 0
						if ref_ub > len(gg_obj.ref):
							ref_ub = len(gg_obj.ref)
										
						out = gg_obj.get_ntvar_counts(ref_lb,ref_ub)

						mywriter.writerow(line_els+out['A']+out['T']+out['C']+out['G']+out['Cnot']+out['Gnot'])						

if __name__ == '__main__':
	gg_obj_file = sys.argv[1]
	gff3_file = sys.argv[2]
	out_file = sys.argv[3]
	# threshold = sys.argv[5]
	# search_seqs = PAM_compile_gg(gg_obj_file,target_len,PAM_orient,PAM_list.split(','),threshold)
	search_seqs = nt_query_gg_all(gg_obj_file,gff3_file,out_file)	



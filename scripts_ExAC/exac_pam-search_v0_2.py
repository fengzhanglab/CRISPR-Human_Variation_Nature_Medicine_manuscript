#David Scott, MIT, 2016
import sys
import csv
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

def levenshtein(s1, s2):
	# adapted from:
	# en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
	if len(s1) < len(s2):
		return levenshtein(s2, s1)

	# len(s1) >= len(s2)
	if len(s2) == 0:
		return len(s1)

	previous_row = range(len(s2) + 1)
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			# j+1 instead of j since previous_row and current_row are one character longer
			insertions = previous_row[j + 1] + 1 
			# than s2
			deletions = current_row[j] + 1       
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	
	return previous_row[-1]

#def PAM_compile_gg(gg_obj_file,targ_len,PAM_orient,PAM_list,threshold):
def leven_search_pam(pam_obj_file,target,threshold,out_file):
	target = str.upper(target)
	targ_len = len(target)

	d = 0
	hit = []
	with open(pam_obj_file,'rb') as pam_handle:
		with open(out_file,'wb') as csvfile:
			pam_obj = pickle.load(pam_handle)	

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'

			ct = 0
			for el in pam_obj:
				d = levenshtein(el[0].upper(),target)
				if d <= threshold:
					hit.append(el+[d])
				
				ct+=1
				if (ct%1000) == 0:
					print '####################'
					print 'vmb memory:' + str(vmb.memory())
					print 'vmb resident: ' + str(vmb.resident())
					print '####################'
					print str(float(ct)/float(len(pam_obj)))

			mywriter = csv.writer(csvfile, delimiter=',')
			for el in hit:
				mywriter.writerow(el)

if __name__ == '__main__':
	pam_obj_file = sys.argv[1]
	target = sys.argv[2]
	threshold = int(sys.argv[3])
	out_file = sys.argv[4]
	search_seqs = leven_search_pam(pam_obj_file,target,threshold,out_file)	



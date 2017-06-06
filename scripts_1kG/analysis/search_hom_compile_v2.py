import os
import csv
import sys
import numpy
from Bio.Seq import Seq
from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

#allele number [female, male, allele count]
#allele count is updated below
hg_ploidy ={'1':[2,2,0],
'2':[2,2,0],
'3':[2,2,0],
'4':[2,2,0],
'5':[2,2,0],
'6':[2,2,0],
'7':[2,2,0],
'8':[2,2,0],
'9':[2,2,0],
'10':[2,2,0],
'11':[2,2,0],
'12':[2,2,0],
'13':[2,2,0],
'14':[2,2,0],
'15':[2,2,0],
'16':[2,2,0],
'17':[2,2,0],
'18':[2,2,0],
'19':[2,2,0],
'20':[2,2,0],
'21':[2,2,0],
'22':[2,2,0],
'X':[2,1,0],
'Y':[0,1,0]}

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

def linear(s1, s2):
	dist=0
	for i in range(0,len(s1)):
		if not (s1[i] == s2[i]):
			dist+=1
	return dist

def main(fin,fin_names,fin_males,fout):

	print 'importing data'

	#read names of 1kG individuals
	namect = 0
	names = []
	name_gend = {}
	with open(fin_names, 'rb') as txtin:
		names = txtin.read().splitlines()
		namect = len(names)

	#read names of 1kG males
	malect=0
	males = []
	females = []
	male_dict = {}
	with open(fin_males, 'rb') as txtin:
		males = txtin.read().splitlines()
		malect=len(males)
		for name in males:
			male_dict[name] = name

	for name in names:
		if name in male_dict:
			name_gend[name] = 1
		else:
			females+=[name]
			name_gend[name] = 0            

	if not ((len(males)+len(females)) == namect):
		print 'error: malect+femalect != namect'
		print 'exiting'
		print len(females)
		print len(males)
		print namect
		exit(1)

	#update population allele counts for each chrm
	for i in range(1,23):
		hg_ploidy[str(i)][2] = 2*namect
	hg_ploidy['X'][2] = (2*namect)-malect
	hg_ploidy['Y'][2] = malect

	########################################
	#time to process!
	########################################

	PAM_len = 0
	PAM_orient = ''
	PAM_str = ''
	pam_temp = ''
	pam_flag = 0

	targ_len = 0
	targ_temp = ''
	targ_specs = []

	d = 0
	len_slice = 0
	seq_temp = ''
	ref_temp = ''
	nmm = 0
	
	ct = 0
	rawdata_loc = 0
	rawdata_info = ''
	rawdata_chrid = ''
	rawdata_offset = 0
	rawdata_allele = 0
	rawdata_allelect = 0
	rawdata_names = ''

	alt_name_dict = {}

	row_out = []
	AC = 0
	AC_Het = 0
	AC_Hom = 0
	AN = 0
	names_out = []	
	with open(fin, 'rb') as csvin:
		csvreader = csv.reader(csvin, delimiter=',')
		with open(fout,'wb') as csvout:
			mywriter = csv.writer(csvout, delimiter=',')
			########################################
			# translate single alleles to pop alleles
			########################################
			for row in csvreader: 
				ct+=1
				if (ct % 1000) == 0:
					print ct 

				#pull ref targets and build targ_specs
				PAM_orient = row[11]
				PAM_temp = row[10]
				PAM_str = row[12]
				targ_temp = row[9]

				rawdata_chrid = row[7]
				rawdata_loc = int(row[13])

				#get allele
				rawdata_allele = int(row[14])
				rawdata_allelect = hg_ploidy[rawdata_chrid][2]
				rawdata_names = row[15].replace('[','').replace(']','')
				rawdata_names = rawdata_names.replace("'",'').replace(' ','').split(',')

				idstr = (str(rawdata_chrid)+'|'+str(rawdata_loc)+'|'+
					PAM_orient+'|'+PAM_str+'|'+PAM_temp+'|'+targ_temp)

				if idstr not in alt_name_dict:
					alt_name_dict[idstr] = [[],[],{}]
					alt_name_dict[idstr][rawdata_allele] = row
					for name in rawdata_names:
						if name not in alt_name_dict[idstr][2]:
							alt_name_dict[idstr][2][name] = 1
						else:
							alt_name_dict[idstr][2][name]+=1
				else:
					alt_name_dict[idstr][rawdata_allele] = row
					for name in rawdata_names:
						if name not in alt_name_dict[idstr][2]:
							alt_name_dict[idstr][2][name] = 1
						else:
							alt_name_dict[idstr][2][name]+=1

			for idstr in alt_name_dict:
				if ((len(alt_name_dict[idstr][0]) > 0) and 
					(len(alt_name_dict[idstr][1]) > 0)):
					row = alt_name_dict[idstr][0]
					if not ((str(alt_name_dict[idstr][0][0:14]) == 
						str(alt_name_dict[idstr][1][0:14])) and
						(str(alt_name_dict[idstr][0][16:26]) == 
						str(alt_name_dict[idstr][1][16:26]))):
						print 'error: non-equivalent paired alleles'
						print alt_name_dict[idstr][0][0:14]
						print alt_name_dict[idstr][1][0:14]
						print alt_name_dict[idstr][0][16:26]
						print alt_name_dict[idstr][1][16:26]						
						print 'exiting'
						exit(1)
				elif len(alt_name_dict[idstr][0]) > 0:
					row = alt_name_dict[idstr][0]
				elif len(alt_name_dict[idstr][1]) > 0:
					row = alt_name_dict[idstr][1]
	
				rawdata_chrid = row[7]

				row_out = []
				row_out+=row[0:14]
				targ_specs = row[16:22]+row[23:26]
			
				AC = 0
				AC_Het = 0
				AC_Hom = 0
				AN = int(row[22])
				names_out = []
				##############################
				if rawdata_chrid == 'Y':
					for name in alt_name_dict[idstr][2]:
						if name not in male_dict:
							print 'error: non male name in Y allele'
							# print 'FEMALES:'
							# print females
							# print 'MALES:'
							# print males
							# print 'MALE DICT:'
							# print male_dict
							print 'IDSTR:'
							print idstr
							print 'ALT NAMES:'
							print alt_name_dict[idstr]
							print 'NAME:'
							print name
							print 'exiting'
							exit(1) 
						else:
							AC+=alt_name_dict[idstr][2][name]
							#No Het/Hom for Y chr
							# if alt_name_dict[idstr][2][name] == 2:
							# 	AC_Hom+=alt_name_dict[idstr][2][name]
							# elif alt_name_dict[idstr][2][name] == 1:
							# 	AC_Het+=alt_name_dict[idstr][2][name]
							names_out+=[name]												

				##############################
				elif (rawdata_chrid == 'X'):
					for name in alt_name_dict[idstr][2]:
						if name in male_dict:
							AC+=alt_name_dict[idstr][2][name]
							#No Het/Hom for male X chr
							# if alt_name_dict[idstr][2][name] == 2:
							# 	AC_Hom+=alt_name_dict[idstr][2][name]
							# elif alt_name_dict[idstr][2][name] == 1:
							# 	AC_Het+=alt_name_dict[idstr][2][name]
							names_out+=[name]
						else:
							AC+=alt_name_dict[idstr][2][name]
							if alt_name_dict[idstr][2][name] == 2:
								AC_Hom+=alt_name_dict[idstr][2][name]
							elif alt_name_dict[idstr][2][name] == 1:
								AC_Het+=alt_name_dict[idstr][2][name]
							names_out+=[name]							

				##############################
				else:
					for name in alt_name_dict[idstr][2]:
						AC+=alt_name_dict[idstr][2][name]
						if alt_name_dict[idstr][2][name] == 2:
							AC_Hom+=alt_name_dict[idstr][2][name]
						elif alt_name_dict[idstr][2][name] == 1:
							AC_Het+=alt_name_dict[idstr][2][name]
						names_out+=[name]							

				row_out+=[AC,AC_Het,AC_Hom,AN,names_out]
				row_out+=targ_specs
				mywriter.writerow(row_out)

fin = sys.argv[1]
fin_names = sys.argv[2]
fin_males = sys.argv[3]
fout = sys.argv[4]

main(fin,fin_names,fin_males,fout)



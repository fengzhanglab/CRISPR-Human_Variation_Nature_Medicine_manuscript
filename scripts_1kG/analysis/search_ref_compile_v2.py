import os
import csv
import sys
import numpy
from Bio.Seq import Seq
from Bio import SeqIO

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

def main(fin,fin_names,fin_males,targct,refin,pamin,fout):

	print 'importing data'

	ids = (['transcript','targchr','targloc','PAM','PAMori','str','targ','searchID',
		'sMM','sTarg','sPAM','sPAMori','sstr','sloc','snames'])
	
	PAM_offset_cut = 3

	#prepare PAM dicts
	PAM_list = pamin.strip().split(',')
	PAM_dict = {}
	for pam in PAM_list:
		PAM_dict[str(pam)] = pam
	PAM_dict_rc = {}
	for pam in PAM_list:
		PAM_dict_rc[revcomp(str(pam))] = pam

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

	#prepare reference dict
	ref_dict = {}
	with open(refin,'r') as fasta:
		for record in SeqIO.parse(fasta, "fasta"):
			print record.id
			ref_dict[record.id] = str(record.seq) 

	########################################
	#time to process!
	########################################
	gRNA_loc = 0
	gRNA_chrid = ''
	gRNA_orient = ''
	gRNA_str = ''

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

	ref_lines = []
	alt_name_dict = {}
	with open(fin, 'rb') as csvin:
		csvreader = csv.reader(csvin, delimiter=',')
		with open(fout,'wb') as csvout:
			mywriter = csv.writer(csvout, delimiter=',')
			for row in csvreader:
				ct+=1
				if (ct % 1000) == 0:
					print ct 
				
				nmm = int(row[8]) 
				if nmm > 3:
					continue

				#pull ref targets and build targ_specs
				PAM_len = len(row[10])
				PAM_orient = row[11]
				PAM_str = row[12]
				targ_temp = row[6]
				targ_len = len(row[6])

				rawdata_info = row[7]
				rawdata_chrid = rawdata_info.split('_')[1]
				rawdata_offset = int(rawdata_info.split('_')[2].split('.')[0])
				rawdata_loc = int(row[13])+rawdata_offset

				#get allele
				rawdata_allele = int(row[14])
				rawdata_allelect = hg_ploidy[rawdata_chrid][2]
				rawdata_names = row[15].replace('[','').replace(']','')
				rawdata_names = rawdata_names.replace("'",'').replace(' ','').split(',')

				#pull ref seq from genomic fasta
				len_slice = targ_len+(2*PAM_len)
				seq_temp = ref_dict[rawdata_chrid][rawdata_loc:(rawdata_loc+len_slice)]

				pam_temp = ''
				ref_temp = ''
				pam_flag = 0
				if (PAM_orient == 'R') and (PAM_str == 'BS'):
					rawdata_loc+=(PAM_offset_cut+PAM_len)
					pam_temp = seq_temp[0:PAM_len]
					ref_temp = revcomp(seq_temp[PAM_len:(PAM_len+targ_len)])                    
					if seq_temp[0:PAM_len] in PAM_dict_rc:
						pam_flag = 1
				elif (PAM_orient == 'R') and (PAM_str == 'TS'):
					rawdata_loc+=((len_slice-(PAM_offset_cut+PAM_len))+1)
					pam_temp = seq_temp[(len_slice-PAM_len):len_slice]
					ref_temp = seq_temp[PAM_len:(PAM_len+targ_len)]                    
					if seq_temp[(len_slice-PAM_len):len_slice] in PAM_dict:
						pam_flag = 1
				elif (PAM_orient == 'L') and (PAM_str == 'TS'):
					rawdata_loc+=PAM_offset_cut+PAM_len
					pam_temp = seq_temp[0:PAM_len]
					ref_temp = seq_temp[PAM_len:(PAM_len+targ_len)]                    
					if seq_temp[0:PAM_len] in PAM_dict:                     
						pam_flag = 1
				elif (PAM_orient == 'L') and (PAM_str == 'BS'):
					rawdata_loc+=((len_slice-(PAM_offset_cut+PAM_len))+1)
					pam_temp = seq_temp[(len_slice-PAM_len):len_slice]
					ref_temp = revcomp(seq_temp[PAM_len:(PAM_len+targ_len)])                    
					if seq_temp[(len_slice-PAM_len):len_slice] in PAM_dict_rc:                              
						pam_flag = 1

				#update fields to hold chr and cut site
				row[7] = rawdata_chrid
				row[13] = rawdata_loc

				#index 1 is PAM present 1/0
				d = linear(targ_temp,ref_temp)
				targ_specs = [d, pam_flag, ref_temp, pam_temp, PAM_orient, PAM_str]

				#annotate on vs off-targets
				not_ontarg = 1
				gRNA_loc = str(row[2])
				gRNA_chrid = str(row[1])
				gRNA_orient = str(row[4])
				gRNA_str = str(row[5])

				# print row
				# if nmm == 0:
				# 	print '#####'
				# 	print 'gRNA:'+str(gRNA_chrid)+'|'+str(gRNA_loc)+'|'+gRNA_orient+'|'+gRNA_str
				# 	print 'targ:'+str(rawdata_chrid)+'|'+str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str
				# 	print '#####'

				if ((str(gRNA_chrid)+'|'+str(gRNA_loc)+'|'+gRNA_orient+'|'+gRNA_str) == 
					(str(rawdata_chrid)+'|'+str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)): 
					print 'gRNA target found'
					print str(gRNA_chrid)+'|'+str(gRNA_loc)+'|'+gRNA_orient+'|'+gRNA_str
					print str(rawdata_chrid)+'|'+str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str
					print row
					not_ontarg = 0

				#Do not write reference lines on first pass
				if str(row[len(row)-1]) == '[\'0\']':
					# print 'reference allele'
					if not (d == int(row[8])):
						print 'error: reference allele targ distance mismatch'
						print '####################'
						print 'row: '+str(row)
						print 'targ_specs: '+str(targ_specs)
						print '####################'
						print 'exiting'
						exit(1)
					else:
						targ_specs+=[rawdata_allelect,targct,not_ontarg,'ref']
						ref_lines+=[row+targ_specs]
						continue

				alt_allele = row[9]+row[10]
				ref_allele = ref_temp+pam_temp
				if ((alt_allele == ref_allele) and
					(not (str(row[len(row)-1]) == '[\'0\']'))):
					continue
					# print 'identical ref and alt alleles: not appending'
					# print 'alt variant in non-pam pam search flank'
					# print '####################'
					# print 'row: '+str(row)
					# print 'targ_specs: '+str(targ_specs)
					# print '####################'  
				else:
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
					str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr not in alt_name_dict:
						alt_name_dict[idstr] = {}
						for name in rawdata_names:
							alt_name_dict[idstr][name] = name
					else:
						for name in rawdata_names:
							alt_name_dict[idstr][name] = name

					targ_specs+=[rawdata_allelect,targct,not_ontarg,'alt']
					mywriter.writerow(row+targ_specs)                                      

			########################################
			# translate reference candidates to 1kG alleles
			########################################
			ct = 0
			for row in ref_lines: 
				ct+=1
				if (ct % 1000) == 0:
					print ct 

				#pull ref targets and build targ_specs
				PAM_orient = row[11]
				PAM_str = row[12]

				rawdata_chrid = row[7]
				rawdata_loc = int(row[13])

				##############################
				if rawdata_chrid == 'Y':
					rawdata_allele = 0
					row[14] = rawdata_allele
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
						str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr not in alt_name_dict:
						row[15] = males
					else:
						for name in alt_name_dict[idstr]:
							if name not in male_dict:
								print 'error: non male name in Y allele'
								# print 'FEMALES:'
								# print females
								# print 'MALES:'
								# print males
								# print 'MALE DICT:'
								# print male_dict
								# print 'IDSTR:'
								# print idstr
								# print 'ALT NAMES:'
								# print alt_name_dict[idstr]
								print 'NAME:'
								print name
								print 'exiting'
								exit(1) 
						#add names not already allocated for allele
						#these are presumed then to be reference                    
						row[15] = []
						for name in males:
							if name not in alt_name_dict[idstr]:
								row[15]+=[name]
					
					# print '#####Y'
					# print rawdata_chrid
					# print len(row[15])

					if len(row[15]) > 0:
						mywriter.writerow(row)

				##############################
				elif (rawdata_chrid == 'X'):
					#compile all names in male alternate alleles for targ
					male_names_temp = {}
					rawdata_allele = 0
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
						str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr in alt_name_dict:
						for name in alt_name_dict[idstr]:
							if name in male_dict:
								if name in male_names_temp:
									print 'error: male name present on > 1 X allele'
									# print 'MALE DICT:'
									# print male_dict
									# print 'IDSTR:'
									# print idstr
									# print 'ALT NAMES:'
									# print alt_name_dict[idstr]
									print 'NAME:'
									print name
									# print 'exiting'
									# exit(1) 	
								else:								
									male_names_temp[name] = name
					rawdata_allele = 1
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
						str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr in alt_name_dict:
						for name in alt_name_dict[idstr]:
							if name in male_dict:
								if name in male_names_temp:
									print 'error: male name present on > 1 X allele'
									# print 'MALE DICT:'
									# print male_dict
									# print 'IDSTR:'
									# print idstr
									# print 'ALT NAMES:'
									# print alt_name_dict[idstr]
									print 'NAME:'
									print name
									# print 'exiting'
									# exit(1) 	
								else:								
									male_names_temp[name] = name

					#add remaining unincluded males to first allele
					rawdata_allele = 0
					row[14] = rawdata_allele
					if len(male_names_temp) == 0:
						#copy list not reference
						row[15] = males[:]
					else:
						#add names not already allocated for allele
						#these are presumed then to be reference                    
						row[15] = []
						for name in males:
							if name not in male_names_temp:
								row[15]+=[name]

					# print '#####X0m'
					# print len(males)					
					# print rawdata_chrid
					# print len(row[15])

					#add remaining unincluded females to first allele
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
						str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr not in alt_name_dict:
						row[15]+=females
					else:
						#add names not already allocated for allele
						#these are presumed then to be reference                    
						for name in females:
							if name not in alt_name_dict[idstr]:
								row[15]+=[name]

					# print '#####X0mf'
					# print len(males)					
					# print rawdata_chrid
					# print len(row[15])

					if len(row[15]) > 0:
						mywriter.writerow(row)

					#add remaining unincluded females to second allele
					rawdata_allele = 1
					row[14] = rawdata_allele
					idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
						str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
					if idstr not in alt_name_dict:
						row[15] = females
					else:
						#add names not already allocated for allele
						#these are presumed then to be reference                    
						row[15] = []
						for name in females:
							if name not in alt_name_dict[idstr]:
								row[15]+=[name]	
					
					# print '#####X1'					
					# print rawdata_chrid
					# print len(row[15])

					if len(row[15]) > 0:
						mywriter.writerow(row)

				##############################
				else:
					for i in range(0,2):
						rawdata_allele = i
						row[14] = rawdata_allele
						idstr = (str(rawdata_allele)+'|'+str(rawdata_chrid)+'|'+
							str(rawdata_loc)+'|'+PAM_orient+'|'+PAM_str)
						if idstr not in alt_name_dict:
							row[15] = names
						else:
							#add names not already allocated for allele
							#these are presumed then to be reference                    
							row[15] = []
							for name in names:
								if name not in alt_name_dict[idstr]:
									row[15]+=[name]
					
						# print '#####else'
						# print rawdata_chrid
						# print len(row[15])

						if len(row[15]) > 0:
							mywriter.writerow(row)

fin = sys.argv[1]
fin_names = sys.argv[2]
fin_males = sys.argv[3]
targct = sys.argv[4]
refin = sys.argv[5]
pamin = sys.argv[6]
fout = sys.argv[7]

main(fin,fin_names,fin_males,targct,refin,pamin,fout)



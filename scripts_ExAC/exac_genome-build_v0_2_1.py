#David Scott, MIT, 2016
import sys
import gzip
import math
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import vcf_obj_v0_6 as vcf
import exac_gen_obj_v0_6_1 as gg

import vmb

def split_chrom_ref_vcf_gff3(ref_filename,vcf_filename,gff3_filename):
	
	with open(ref_filename,'r') as fin:
		for record in SeqIO.parse(fin, "fasta"):
			print record.id

			print '####################'
			print 'vmb memory:' + str(vmb.memory())
			print 'vmb resident: ' + str(vmb.resident())
			print '####################'
			
			with open(ref_filename.replace('.fa',
				'_'+str(record.id)+".fa"),'wb') as fout:
				SeqIO.write(record, fout, "fasta")			

	# with gzip.open(vcf_filename,'r') as fin:
	# 	head_flag = True
	# 	vcf_head = []
	# 	CHROM0 = ''
	# 	CHROM = ''
	# 	ct = 0		
	# 	for line in fin:
	# 		if head_flag:
	# 			if line[0] == '#':
	# 				print line
	# 				vcf_head.append(line)
	# 				continue
	# 			else:
	# 				head_flag = False

	# 		CHROM = str(line).strip().split()[0]
	# 		tPOS = str(line).strip().split()[1]

	# 		if ct == 0:
	# 			CHROM0 = CHROM
	# 			fout = open(vcf_filename.replace('.vep.vcf.gz',
	# 				'_'+str(CHROM)+".vep.vcf"),'wb')
	# 			print 'opened: '+vcf_filename.replace('.vep.vcf.gz',
	# 				'_'+str(CHROM)+".vep.vcf")
	# 			print 'writing head'				
	# 			for line_tmp in vcf_head:
	# 				fout.write(line_tmp)				

	# 		print 'CHROM0: '+str(CHROM0)
	# 		print 'CHROM: '+str(CHROM)
	# 		print 'tPOS: '+str(tPOS)

	# 		if not (CHROM == CHROM0):
	# 			fout.close()
	# 			print 'closed: '+vcf_filename.replace('.vep.vcf.gz',
	# 				'_'+str(CHROM0)+".vep.vcf")
	# 			CHROM0 = CHROM
	# 			fout = open(vcf_filename.replace('.vep.vcf.gz',
	# 				'_'+str(CHROM)+".vep.vcf"),'wb')
	# 			print 'opened: '+vcf_filename.replace('.vep.vcf.gz',
	# 				'_'+str(CHROM)+".vep.vcf")
	# 			print 'writing head'				
	# 			for line_tmp in vcf_head:
	# 				fout.write(line_tmp)

	# 		fout.write(line)
	# 		ct+=1
	# 	fout.close()

	# with gzip.open(gff3_filename,'r') as fin:
	# 	head_flag = True
	# 	gff3_head = []
	# 	lines_out = []
	# 	CHROM0 = ''
	# 	CHROM = ''
	# 	ct = 0	

	# 	for line in fin:
	# 		if line[0] == '#':
	# 			print line
	# 			continue

	# 		CHROM = str(line).strip().split()[0].replace('chr','')
	# 		tPOS = int(str(line).strip().split()[3])

	# 		if ct == 0:
	# 			CHROM0 = CHROM
	# 			fout = open(gff3_filename.replace('.gff3.gz',
	# 				'_'+str(CHROM)+".gff3"),'wb')
	# 			print 'opened: '+gff3_filename.replace('.gff3.gz',
	# 				'_'+str(CHROM)+".gff3")
	# 			print 'writing head'				
	# 			for line_tmp in gff3_head:
	# 				fout.write(line_tmp)				

	# 		print 'CHROM0: '+str(CHROM0)
	# 		print 'CHROM: '+str(CHROM)
	# 		print 'tPOS: '+str(tPOS)

	# 		if not (CHROM == CHROM0):
	# 			lines_out.sort(key=lambda x: x[0])
	# 			for el in lines_out:
	# 				fout.write(el[1])
	# 			lines_out = []
	# 			fout.close()
	# 			print 'closed: '+gff3_filename.replace('.gff3.gz',
	# 				'_'+str(CHROM0)+".gff3")
	# 			CHROM0 = CHROM
	# 			fout = open(gff3_filename.replace('.gff3.gz',
	# 				'_'+str(CHROM)+".gff3"),'wb')
	# 			print 'opened: '+gff3_filename.replace('.gff3.gz',
	# 				'_'+str(CHROM)+".gff3")
	# 			print 'writing head'				
	# 			for line_tmp in gff3_head:
	# 				fout.write(line_tmp)

	# 		lines_out.append([tPOS,line])
	# 		ct+=1

	# 	lines_out.sort(key=lambda x: x[0])
	# 	for el in lines_out:
	# 		fout.write(el[1])
	# 	lines_out = []
	# 	fout.close()			

def split_gg_ref_vcf_gff3(ref_filename,vcf_filename,gff3_filename,cov_filename,bin_size,bin_overlap):
	
	# with open(ref_filename,'r') as fin:
	# 	for record in SeqIO.parse(fin, "fasta"):
	# 		print record.id
	# 		# with open(record.id+".fa",'w') as fo:
	# 		# 	SeqIO.write(record, fo, "fasta")			
	# 		print len(record)
	# 		print (bin_size-bin_overlap)
	# 		print math.ceil(float(len(record))/float(bin_size-bin_overlap))
	# 		print int(math.ceil(float(len(record))/float(bin_size-bin_overlap)))
	# 		print range(int(math.ceil(float(len(record))/float(bin_size-bin_overlap))))

	# 		bin_start = 0
	# 		write_bin = int(math.ceil(float(len(record))/float(bin_size-bin_overlap))) - 1
	# 		for i in range(int(math.ceil(float(len(record))/float(bin_size-bin_overlap)))):
	# 			print i
	# 			print bin_start
	# 			if (i == write_bin):
	# 				print 'writing'
	# 				fout = open(ref_filename.replace(
	# 					'.fa','_'+str(bin_start)+'.fa'),'wb')
	# 				SeqIO.write(record[bin_start:(bin_start+bin_size)],
	# 					fout,"fasta")
	# 				fout.close()
	# 			bin_start = bin_start+(bin_size-bin_overlap)

	# tPOS = 0
	# bin_end = 0
	# bin_start = 0
	# head_flag = True
	# vcf_head = []
	# vcf_vars = []
	# vcf_vars_overlap = []
	# with open(vcf_filename,'r') as fin:
	# 	for line in fin:
	# 		if head_flag:
	# 			if line[0] == '#':
	# 				print line
	# 				vcf_head.append(line)
	# 				continue
	# 			else:
	# 				head_flag = False

	# 		tPOS = int(str(line).strip().split()[1])
	# 		print tPOS
	# 		if (tPOS-bin_end) > bin_size:
	# 			if len(vcf_vars) > 0:
	# 				fout = open(vcf_filename.replace(
	# 					'.vep.vcf','_'+str(bin_start)+'.vep.vcf'),'wb')
	# 				for line_tmp in vcf_head:
	# 					fout.write(line_tmp)
	# 				for line_tmp in vcf_vars:
	# 					fout.write(line_tmp)
	# 				fout.close()
	# 			while bin_end < tPOS:
	# 				bin_start = bin_start+(bin_size-bin_overlap)
	# 				bin_end = bin_start+bin_size
	# 			vcf_vars = []
	# 			vcf_vars_overlap = []
				
	# 		elif (tPOS-bin_start) > bin_size:
	# 			if len(vcf_vars) > 0:
	# 				fout = open(vcf_filename.replace(
	# 					'.vep.vcf','_'+str(bin_start)+'.vep.vcf'),'wb')
	# 				for line_tmp in vcf_head:
	# 					fout.write(line_tmp)
	# 				for line_tmp in vcf_vars:
	# 					fout.write(line_tmp)
	# 				fout.close()
	# 			bin_start = bin_start+(bin_size-bin_overlap)
	# 			bin_end = bin_start+bin_size
	# 			vcf_vars = vcf_vars_overlap
	# 			vcf_vars_overlap = []

	# 		if (tPOS-bin_start) > (bin_size-bin_overlap):
	# 			vcf_vars.append(line)
	# 			vcf_vars_overlap.append(line)
	# 		else:
	# 			vcf_vars.append(line)

	# 	if len(vcf_vars) > 0:
	# 		fout = open(vcf_filename.replace(
	# 			'.vep.vcf','')+'_'+str(bin_start)+'.vep.vcf','wb')
	# 		for line_tmp in vcf_head:
	# 			fout.write(line_tmp)
	# 		for line_tmp in vcf_vars:
	# 			fout.write(line_tmp)
	# 		fout.close()

	# tPOS = 0
	# bin_end = 0
	# bin_start = 0
	# head_flag = True
	# gff3_head = []
	# gff3_vars = []
	# gff3_vars_overlap = []
	# with open(gff3_filename,'r') as fin:
	# 	for line in fin:
	# 		if head_flag:
	# 			if line[0] == '#':
	# 				print line
	# 				gff3_head.append(line)
	# 				continue
	# 			else:
	# 				head_flag = False

	# 		tPOS = int(str(line).strip().split()[3])
	# 		print tPOS
	# 		if (tPOS-bin_end) > bin_size:
	# 			if len(gff3_vars) > 0:
	# 				fout = open(gff3_filename.replace(
	# 					'.gff3','_'+str(bin_start)+'.gff3'),'wb')
	# 				for line_tmp in gff3_head:
	# 					fout.write(line_tmp)
	# 				for line_tmp in gff3_vars:
	# 					fout.write(line_tmp)
	# 				fout.close()
	# 			while bin_end < tPOS:
	# 				bin_start = bin_start+(bin_size-bin_overlap)
	# 				bin_end = bin_start+bin_size
	# 			gff3_vars = []
	# 			gff3_vars_overlap = []
				
	# 		elif (tPOS-bin_start) > bin_size:
	# 			if len(gff3_vars) > 0:
	# 				fout = open(gff3_filename.replace(
	# 					'.gff3','_'+str(bin_start)+'.gff3'),'wb')
	# 				for line_tmp in gff3_head:
	# 					fout.write(line_tmp)
	# 				for line_tmp in gff3_vars:
	# 					fout.write(line_tmp)
	# 				fout.close()
	# 			bin_start = bin_start+(bin_size-bin_overlap)
	# 			bin_end = bin_start+bin_size
	# 			gff3_vars = gff3_vars_overlap
	# 			gff3_vars_overlap = []

	# 		if (tPOS-bin_start) > (bin_size-bin_overlap):
	# 			gff3_vars.append(line)
	# 			gff3_vars_overlap.append(line)
	# 		else:
	# 			gff3_vars.append(line)

	# 	if len(gff3_vars) > 0:
	# 		fout = open(gff3_filename.replace(
	# 			'.gff3','')+'_'+str(bin_start)+'.gff3','wb')
	# 		for line_tmp in gff3_head:
	# 			fout.write(line_tmp)
	# 		for line_tmp in gff3_vars:
	# 			fout.write(line_tmp)
	# 		fout.close()

	tPOS = 0
	bin_end = 0
	bin_start = 0
	head_flag = True
	cov_head = []
	cov_vars = []
	cov_vars_overlap = []
	with gzip.open(cov_filename,'r') as fin:
		for line in fin:
			if head_flag:
				if line[0] == '#':
					print line
					cov_head.append(line)
					continue
				else:
					head_flag = False
					
			tPOS = int(str(line).strip().split()[1])
			print tPOS
			if (tPOS-bin_end) > bin_size:
				if len(cov_vars) > 0:
					fout = open(cov_filename.replace(
						'.txt.gz','_'+str(bin_start)+'.txt'),'wb')
					for line_tmp in cov_head:
						fout.write(line_tmp)
					for line_tmp in cov_vars:
						fout.write(line_tmp)
					fout.close()
				while bin_end < tPOS:
					bin_start = bin_start+(bin_size-bin_overlap)
					bin_end = bin_start+bin_size
				cov_vars = []
				cov_vars_overlap = []
				
			elif (tPOS-bin_start) > bin_size:
				if len(cov_vars) > 0:
					fout = open(cov_filename.replace(
						'.txt.gz','_'+str(bin_start)+'.txt'),'wb')
					for line_tmp in cov_head:
						fout.write(line_tmp)
					for line_tmp in cov_vars:
						fout.write(line_tmp)
					fout.close()
				bin_start = bin_start+(bin_size-bin_overlap)
				bin_end = bin_start+bin_size
				cov_vars = cov_vars_overlap
				cov_vars_overlap = []

			if (tPOS-bin_start) > (bin_size-bin_overlap):
				cov_vars.append(line)
				cov_vars_overlap.append(line)
			else:
				cov_vars.append(line)

		if len(cov_vars) > 0:
			fout = open(cov_filename.replace(
				'.txt.gz','')+'_'+str(bin_start)+'.txt','wb')
			for line_tmp in cov_head:
				fout.write(line_tmp)
			for line_tmp in cov_vars:
				fout.write(line_tmp)
			fout.close()

def build_geno_genome(ref_filename, vcf_filename, cov_filename):
	rf_split = ref_filename.split('_')
	offset = int(rf_split[len(rf_split)-1].replace('.fa',''))
	print offset

	gg_obj = object()
	with open(ref_filename,'r') as fin:
		for record in SeqIO.parse(fin, "fasta"):
			print record.id
			# with open(record.id+".fa",'w') as fo:
			# 	SeqIO.write(record, fo, "fasta")			
			gg_obj = gg.Geno_genome(str(record.seq))

	with open(ref_filename.replace('.fa','')+'.gg','wb') as fout:
		if len(vcf_filename) > 0:

			vcf_head = []
			vcf_people = []
			head_flag = True
			with open(vcf_filename,'r') as fin:
				last_line = ''
				for line in fin:
					if head_flag:
						if line[0] == '#':
							print line
							vcf_head.append(line)
							last_line = line
							continue
						else:
							head_flag = False

					#test block
					# print '########################################'
					# print str(line.strip().split())
					# print str(last_line.strip().split())
					# print '########################################'

					variant_row = vcf.Vcf_variant_row(str(line).strip().split(),offset)
					if variant_row.tFILTER == 'PASS':
						variant_row.compile_var_specs()
						# variant_row.print_var_specs()
						variants = variant_row.get_variants()
						gg_obj.add_variants(variants)

			loc = 0
			cov = 0
			q20 = 0
			row = []
			cov_head = []
			head_flag = True
			with open(cov_filename,'r') as fin:
				for line in fin:
					if head_flag:
						if line[0] == '#':
							print line
							print str(line).strip().split()[1]
							print str(line).strip().split()[2]
							print str(line).strip().split()[8]
							cov_head.append(line)
							continue
						else:
							head_flag = False

					#test block
					# print '########################################'
					# print str(line.strip().split())
					# print str(last_line.strip().split())
					# print '########################################'

					row = str(line).strip().split()
					loc = int(row[1])
					cov = float(row[2])
					q20 = float(row[8])
					if (cov*q20) > 10:
						print loc
						gg_obj.add_mask(loc,offset)

		pickle.dump(gg_obj,fout)

if __name__ == '__main__':

	method = sys.argv[1]
	print "geno_genome-build_v0_1.py"
	print method
	
	if method == 'split_chrom':
		fasta_filename = sys.argv[2]
		vcf_filename = sys.argv[3]
		gff3_filename = sys.argv[4]
		split_chrom_ref_vcf_gff3(fasta_filename,vcf_filename,gff3_filename)	
	elif method == 'split_vcf':
		fasta_filename = sys.argv[2]
		vcf_filename = sys.argv[3]
		gff3_filename = sys.argv[4]
		cov_filename = sys.argv[5]
		bin_size = int(sys.argv[6])
		bin_overlap = int(sys.argv[7])
		split_gg_ref_vcf_gff3(fasta_filename,vcf_filename,gff3_filename,cov_filename,bin_size,bin_overlap)
	elif method == 'build_gg':
		fasta_filename = sys.argv[2]
		if len(sys.argv) > 3:
			vcf_filename = sys.argv[3]
			cov_filename = sys.argv[4]
		else:
			vcf_filename = ''
			cov_filename = ''
		build_geno_genome(fasta_filename,vcf_filename,cov_filename)



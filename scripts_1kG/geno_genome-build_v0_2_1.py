#David Scott, MIT, 2016
import sys
import gzip
import math
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import vcf_obj_v0_6 as vcf
import geno_gen_obj_v0_6_1 as gg

def split_ref_genome(ref_filename):
	
	with open(ref_filename,'r') as fin:
		for record in SeqIO.parse(fin, "fasta"):
			print record.id
			with open(ref_filename.replace('.fa',
				'_'+str(record.id)+".fa"),'wb') as fout:
				SeqIO.write(record, fout, "fasta")			

def split_gg_ref_vcf(ref_filename,vcf_filename,bin_size,bin_overlap):
	
	with open(ref_filename,'r') as fin:
		for record in SeqIO.parse(fin, "fasta"):
			print record.id
			# with open(record.id+".fa",'w') as fo:
			# 	SeqIO.write(record, fo, "fasta")			
			print len(record)
			print (bin_size-bin_overlap)
			print math.ceil(float(len(record))/float(bin_size-bin_overlap))
			print int(math.ceil(float(len(record))/float(bin_size-bin_overlap)))
			print range(int(math.ceil(float(len(record))/float(bin_size-bin_overlap))))

			bin_start = 0
			write_bin = int(math.ceil(float(len(record))/float(bin_size-bin_overlap))) - 1
			for i in range(int(math.ceil(float(len(record))/float(bin_size-bin_overlap)))):
				print i
				print bin_start
				if (i == write_bin):
					print 'writing'
					fout = open(ref_filename.replace(
						'.fa','_'+str(bin_start)+'.fa'),'wb')
					SeqIO.write(record[bin_start:(bin_start+bin_size)],
						fout,"fasta")
					fout.close()
				bin_start = bin_start+(bin_size-bin_overlap)

	# tPOS = 0
	# bin_end = 0
	# bin_start = 0
	# head_flag = True
	# vcf_head = []
	# vcf_vars = []
	# vcf_vars_overlap = []
	# with gzip.open(vcf_filename,'r') as fin:
	# 	for line in fin:
	# 		if head_flag:
	# 			if line[0] == '#':
	# 				print line
	# 				vcf_head.append(line)
	# 				continue
	# 			else:
	# 				head_flag = False

	# 		tPOS = int(str(line).strip().split()[1])
	# 		if (tPOS-bin_end) > bin_size:
	# 			if len(vcf_vars) > 0:
	# 				fout = open(vcf_filename.replace(
	# 					'.vcf.gz','_'+str(bin_start)+'.vcf'),'wb')
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
	# 					'.vcf.gz','_'+str(bin_start)+'.vcf'),'wb')
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
	# 			'.vcf.gz','')+'_'+str(bin_start)+'.vcf','wb')
	# 		for line_tmp in vcf_head:
	# 			fout.write(line_tmp)
	# 		for line_tmp in vcf_vars:
	# 			fout.write(line_tmp)
	# 		fout.close()

def build_geno_genome(ref_filename, vcf_filename):
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

	head_flag = True
	vcf_head = []
	vcf_people = []
	with open(ref_filename.replace('.fa','')+'.gg','wb') as fout:
		if len(vcf_filename) > 0:
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
							vcf_people = vcf.Vcf_row(str(last_line).strip().split())
							head_flag = False

					#test block
					# print '########################################'
					# print str(line.strip().split())
					# print str(last_line.strip().split())
					# print '########################################'

					variant_row = vcf.Vcf_variant_row(str(line).strip().split(),
						offset,vcf_people.data)
					if variant_row.tFILTER == 'PASS':
						variant_row.compile_var_specs()
						# variant_row.print_var_specs()
						variants = variant_row.get_variants()
						gg_obj.add_variants(variants)

		pickle.dump(gg_obj,fout)

if __name__ == '__main__':

	method = sys.argv[1]
	print "geno_genome-build_v0_1.py"
	print method
	
	if method == 'split_ref':
		fasta_filename = sys.argv[2]
		split_ref_genome(fasta_filename)	
	elif method == 'split_vcf':
		fasta_filename = sys.argv[2]
		vcf_filename = sys.argv[3]
		bin_size = int(sys.argv[4])
		bin_overlap = int(sys.argv[5])
		split_gg_ref_vcf(fasta_filename,vcf_filename,bin_size,bin_overlap)
	elif method == 'build_gg':
		fasta_filename = sys.argv[2]
		if len(sys.argv) > 3:
			vcf_filename = sys.argv[3]
		else:
			vcf_filename = ''
		build_geno_genome(fasta_filename,vcf_filename)



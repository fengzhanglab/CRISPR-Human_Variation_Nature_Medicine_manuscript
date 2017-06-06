#David Scott, MIT, 2016
import sys
import gzip


class Variant(object):
	"""docstring for Variant"""
	def __init__(self, location, var_type, var_seq, var_allele):
		super(Variant, self).__init__()
		self.loc = int(location)
		self.seq = str(var_seq)
		self.type = str(var_type)
		self.allele = int(var_allele)
		self.idstr = (str(self.allele)+':'+
			str(self.loc)+self.seq+self.type)
		self.names = []
		
	def add_name(self,name):
		self.names.append(name)

	def add_names(self,name_list):
		self.names+=name_list	

class Vcf_row(object):
	"""docstring for Vcf_row"""
	def __init__(self, row_list):
		super(Vcf_row, self).__init__()
		self.CHROM = str(row_list[0])
		self.tPOS = str(row_list[1])
		self.tID = str(row_list[2])
		self.tREF = str(row_list[3])
		self.tALT = str(row_list[4])
		self.tQUAL = str(row_list[5])
		self.tFILTER = str(row_list[6])
		self.tINFO = str(row_list[7])
		self.tFORMAT = str(row_list[8])
		self.data = row_list[9:len(row_list)]
		self.len = len(self.data)

class Vcf_variant_row(Vcf_row):
	"""docstring for Vcf_variant_row"""
	def __init__(self, row_list, offset, data_names_list):
		super(Vcf_variant_row, self).__init__(row_list)
		self.tPOS = int(self.tPOS)-int(offset)
		self.data_names = data_names_list
		self.var_dict = {}

	def var_dict_empty(self):
		if len(self.var_dict) > 0:
			return False
		else:
			return True

	def get_snp(self,tREF_list,tALT_list,allele_id):
		#snp

		########################################
		#checks for structure of insertion
		########################################
		#check that if len REF is greater than 1
		#the 3p extension of REF beyond 1
		#should exist at the 3p of the snp ALT seq
		#case of complex locus consisting of mult snp/ins/del alleles
		if len(tREF_list[0]) > 1: 
			ind1 = (len(tALT_list[allele_id])-
				(len(tREF_list[0])-1))
			ind2 = len(tALT_list[allele_id])
			if not (tALT_list[allele_id][ind1:ind2] == 
				tREF_list[0][1:len(tREF_list[0])]):
				print "COMPLEX SNP reference 3p flank mismatch"

				#test block
				print '########################################'
				print 'CHROM: '+str(self.CHROM)
				print 'location: '+str(self.tPOS)
				print 'REF allele: '+str(tREF_list[0])
				print 'ALT allele: '+str(tALT_list[allele_id])
				print '########################################'

				#160714: report, but do not exit
				# sys.exit(1)
			else:
				print '########## SNP COMPLEX CASE ##########'

		#return [location, type, seq]
		return [self.tPOS, 'SNP', tALT_list[allele_id][0]]

	def get_indel_del(self,tREF_list,tALT_list,allele_id):
		#deletion

		########################################
		#checks for structure of deletion
		#
		#check if first base of string always matches
		#checks if in complex cases illustrated below
		#that indel always right justifies
		#	if deletion insertion, case where first base may not
		#	REF: TCAGAT
		#	ALT: TC-TAT (del of 2, ins 1 OR del 1, snp 1)
		#	-> 2  CA  C
		#	+  3  G   T
		#	REF: TC-AGAT
		#	ALT: TCTTTAT (del of 2, ins 3 OR ins 1, snp 2)
		#	-> 2  C   CT
		#	+  3  A   T
		#	+  4  G   T
		########################################
		if not (tREF_list[0][0] == tALT_list[allele_id][0]):
			print "INDEL:DEL reference 5p flank mismatch"

			#test block
			print '########################################'
			print 'CHROM: '+str(self.CHROM)
			print 'location: '+str(self.tPOS)
			print 'REF allele: '+str(tREF_list[0])
			print 'ALT allele: '+str(tALT_list[allele_id])
			print '########################################'

			#160714: report, but do not exit
			# sys.exit(1)

		seq = tREF_list[0][1:(
				len(tREF_list[0])-(
					len(tALT_list[allele_id])-1))]

		#return [location, type, seq]
		return [self.tPOS+1, 'INDEL:DEL', seq]		

	def get_indel_ins(self,tREF_list,tALT_list,allele_id):
		#insertion

		########################################
		#checks for structure of insertion
		########################################

		#check that if len REF is greater than 1
		#the 3p extension of REF beyond 1
		#should exist at the 3p of the insertion ALT seq
		#case of complex locus consisting of mult snp/ins/del alleles
		if len(tREF_list[0]) > 1: 
			ind1 = (len(tALT_list[allele_id])-
				(len(tREF_list[0])-1))
			ind2 = len(tALT_list[allele_id])
			if not (tALT_list[allele_id][ind1:ind2] == 
				tREF_list[0][1:len(tREF_list[0])]):
				print "COMPLEX INDEL:INS reference 3p flank mismatch"

				#test block
				print '########################################'
				print 'CHROM: '+str(self.CHROM)
				print 'location: '+str(self.tPOS)
				print 'REF allele: '+str(tREF_list[0])
				print 'ALT allele: '+str(tALT_list[allele_id])
				print '########################################'

				#160714: report, but do not exit
				# sys.exit(1)														
			else:
				print '########## INS COMPLEX CASE ##########'

		if not (tREF_list[0][0] == tALT_list[allele_id][0]):
			print "INDEL:INS reference 5p flank mismatch"

			#test block
			print '########################################'
			print 'CHROM: '+str(self.CHROM)
			print 'location: '+str(self.tPOS)
			print 'REF allele: '+str(tREF_list[0])
			print 'ALT allele: '+str(tALT_list[allele_id])
			print '########################################'
			
			seq = tALT_list[allele_id]

		else:
			seq = tALT_list[allele_id][1:(
					len(tALT_list[allele_id])-(
						len(tREF_list[0])-1))]			

		#return [location, type, seq]
		return [self.tPOS, 'INDEL:INS', seq]

	def compile_var_specs(self):
		########################################
		#currently only handling SNP, INDEL
		#not currently handling long INS/DEL
		#future versions should include all ##ALT variants
		########################################
		tINFO_list = self.tINFO.split(';')
		tINFO_dict = {}
		for el in tINFO_list:
			els = el.split('=')
			if len(els) > 1:
				tINFO_dict[els[0]] = els[1]

		if not (('SNP' in tINFO_dict['VT']) or 
			('INDEL' in tINFO_dict['VT'])):
			return []

		tREF_list = self.tREF.split(',')
		tALT_list = self.tALT.split(',')
		
		#tREF list should be of length 1
		if len(tREF_list) > 1:
			print "tREF_list length greater than 1 for variant ID:"
			print self.tID
			print tREF_list
			print tALT_list			

			#160714: report, but do not exit
			# sys.exit(1)

		tALT_list_lens = [len(el) for el in tALT_list]		

		#INDEL:DEL structure check
		#ONLY INDEL case where len(tREF_list[0]) > 1
		#AND min ALT len should reflect full DEL

		#160714 turned off
		#see case: X	16001576	.	ATT	ATTT
		
		# if (len(tREF_list[0]) > 1) and (min(tALT_list_lens) > 1):
		# 	print "INDEL:DEL not rooted at tPOS"
		# 	print self.tID
		# 	print tREF_list
		# 	print tALT_list
		# 	#raise exception
		# 	sys.exit(1)

		var_loc0 = 0
		var_seq0 = ''
		var_loc = 0
		var_seq = ''
		var_type = ''	
		var_name = ''
		var_idstr = ''
		allele_ct = 0

		el_list = []
		ct = 0

		#test block
		print '##self.tPOS: '+str(self.tPOS) 
		# print '########################################'
		# print self.data
		# print self.data_names
		# print '########################################'

		for el in self.data:
			#skip locations that have no call
			if '.' in el:

				#test block
				# print '########################################'
				# print 'allele string: '+str(el)
				# print 'location: '+str(self.tPOS)
				# print 'data_name: '+str(self.data_names[ct])		
				# print '########################################'

				ct+=1
				continue

			el_list = el.split('|')

			allele_ct = 0
			for allele_id in el_list:
				allele_id = int(allele_id)-1
				
				if allele_id == -1:
					allele_ct+=1
					continue

				var_name = self.data_names[ct]	

				#triage INDEL, SNP
				if len(tREF_list[0]) == len(tALT_list[allele_id]): 
					OUT = self.get_snp(tREF_list,tALT_list,allele_id)
					if len(OUT)==3:
						[var_loc, var_type, var_seq] = OUT

						var_idstr = (str(var_loc)+str(var_type)+str(var_seq)+str(allele_ct))
						if var_idstr not in self.var_dict:
							self.var_dict[var_idstr] = ([var_loc,var_type,var_seq,allele_ct,[]])
							self.var_dict[var_idstr][4].append(var_name)
						else:
							self.var_dict[var_idstr][4].append(var_name)

				elif len(tREF_list[0]) > len(tALT_list[allele_id]):
					OUT = self.get_indel_del(tREF_list,tALT_list,allele_id)
					if len(OUT)==3:
						[var_loc0, var_type, var_seq0] = OUT

						#test block
						# print '########################################'
						# print 'var_loc0: '+str(var_loc0)
						# print 'var_type: '+str(var_type)
						# print 'var_seq0: '+str(var_seq0)
						# print '########################################'

						for i in range(0,len(var_seq0)):
							var_loc = var_loc0+i
							var_seq = var_seq0[i]
	
							#test block
							# print '########################################'
							# print 'var_ind: '+str(i)
							# print 'var_loc: '+str(var_loc)
							# print 'var_seq: '+str(var_seq)
							# print '########################################'
	
							var_idstr = (str(var_loc)+str(var_type)+str(var_seq)+str(allele_ct))
							if var_idstr not in self.var_dict:
								self.var_dict[var_idstr] = ([var_loc,var_type,var_seq,allele_ct,[]])
								self.var_dict[var_idstr][4].append(var_name)
							else:
								self.var_dict[var_idstr][4].append(var_name)	

				elif len(tREF_list[0]) < len(tALT_list[allele_id]):
					OUT = self.get_indel_ins(tREF_list,tALT_list,allele_id)
					if len(OUT)==3:
						[var_loc, var_type, var_seq] = OUT

						var_idstr = (str(var_loc)+str(var_type)+str(var_seq)+str(allele_ct))
						if var_idstr not in self.var_dict:
							self.var_dict[var_idstr] = ([var_loc,var_type,var_seq,allele_ct,[]])
							self.var_dict[var_idstr][4].append(var_name)
						else:
							self.var_dict[var_idstr][4].append(var_name)											

				#test block
				# print '########################################'
				# print 'allele string: '+str(el)
				# print 'location: '+str(self.tPOS)
				# print 'REF allele: '+str(tREF_list[0])
				# print 'ALT allele: '+str(tALT_list[allele_id])
				# print 'count: '+str(ct)
				# print 'var_loc: '+str(var_loc)
				# print 'var_type: '+str(var_type)
				# print 'var_seq: '+str(var_seq)
				# print 'var_name: '+str(var_name)
				# print '########################################'						

				allele_ct+=1

			ct+=1

	def print_var_specs(self):
		if not self.var_dict_empty():
			for idstr in self.var_dict:
				var = self.var_dict[idstr]
				print 'location: '+str(var[0])
				print 'type: '+str(var[1])
				print 'sequence: '+str(var[2])
				print 'allele: '+str(var[3])
				print 'names: '+str(var[4])

	def get_variants(self):
		variants = []
		if not self.var_dict_empty():
			for idstr in self.var_dict:
				var = self.var_dict[idstr]
				variant = Variant(var[0], var[1], var[2], var[3])
				variant.add_names(var[4])
				variants.append(variant)
		return variants



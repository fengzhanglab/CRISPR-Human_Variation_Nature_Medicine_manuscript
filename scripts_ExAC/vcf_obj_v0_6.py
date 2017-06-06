#David Scott, MIT, 2016
import sys
import gzip

class Vep(object):
	def __init__(self,vep_list):
		super(Vep, self).__init__()	
		self.Consequence = vep_list[0].split('&')
		self.IMPACT = vep_list[1]
		self.SYMBOL = vep_list[2]
		self.Gene = vep_list[3]
		self.Feature_type = vep_list[4]
		self.Feature = vep_list[5]
		self.BIOTYPE = vep_list[6]
		self.EXON = vep_list[7]
		self.INTRON = vep_list[8]
		self.HGVSc = vep_list[9]
		self.HGVSp = vep_list[10]
		self.cDNA_position = vep_list[11]
		self.CDS_position = vep_list[12]
		self.Protein_position = vep_list[13]
		self.Amino_acids = vep_list[14]
		self.Codons = vep_list[15]
		self.Existing_variation = vep_list[16]
		self.ALLELE_NUM = vep_list[17]
		self.DISTANCE = vep_list[18]
		self.STRAND = vep_list[19]
		self.VARIANT_CLASS = vep_list[20]
		self.MINIMISED = vep_list[21]
		self.SYMBOL_SOURCE = vep_list[22]
		self.HGNC_ID = vep_list[23]
		self.CANONICAL = vep_list[24]
		self.TSL = vep_list[25]
		self.CCDS = vep_list[26]
		self.ENSP = vep_list[27]
		self.SWISSPROT = vep_list[28]
		self.TREMBL = vep_list[29]
		self.UNIPARC = vep_list[30]
		self.SIFT = vep_list[31]
		self.PolyPhen = vep_list[32]
		self.DOMAINS = vep_list[33]
		self.HGVS_OFFSET = vep_list[34]
		self.GMAF = vep_list[35]
		self.AFR_MAF = vep_list[36]
		self.AMR_MAF = vep_list[37]
		self.ASN_MAF = vep_list[38]
		self.EAS_MAF = vep_list[39]
		self.EUR_MAF = vep_list[40]
		self.SAS_MAF = vep_list[41]
		self.AA_MAF = vep_list[42]
		self.EA_MAF = vep_list[43]
		self.CLIN_SIG = vep_list[44]
		self.SOMATIC = vep_list[45]
		self.PHENO = vep_list[46]
		self.PUBMED = vep_list[47]
		self.MOTIF_NAME = vep_list[48]
		self.MOTIF_POS = vep_list[49]
		self.HIGH_INF_POS = vep_list[50]
		self.MOTIF_SCORE_CHANGE = vep_list[51]
		self.LoF_info = vep_list[52]
		self.LoF_flags = vep_list[53]
		self.LoF_filter = vep_list[54]
		self.LoF = vep_list[55]
		self.context = vep_list[56]
		self.ancestral = vep_list[57]

class Variant(object):
	"""docstring for Variant"""
	def __init__(self, location, var_type, var_ref, var_seq, var_allele):
		super(Variant, self).__init__()
		self.loc = int(location)
		self.ref = str(var_ref)
		self.seq = str(var_seq)
		self.type = str(var_type)
		self.allele = int(var_allele)
		self.idstr = (str(self.allele)+':'+
			str(self.loc)+self.seq+self.type)
		self.names = []
		self.veps = {}		
		self.AC=0
		self.AC_AFR=0
		self.AC_AMR=0
		self.AC_Adj=0
		self.AC_EAS=0
		self.AC_FIN=0
		self.AC_Het=0
		self.AC_Hom=0
		self.AC_NFE=0
		self.AC_OTH=0
		self.AC_SAS=0
		self.AN=0
		self.AN_AFR=0
		self.AN_AMR=0
		self.AN_Adj=0
		self.AN_EAS=0
		self.AN_FIN=0
		self.AN_NFE=0
		self.AN_OTH=0
		self.AN_SAS=0
		self.Het_AFR=0
		self.Het_AMR=0
		self.Het_EAS=0
		self.Het_FIN=0
		self.Het_NFE=0
		self.Het_OTH=0
		self.Het_SAS=0
		self.Hom_AFR=0
		self.Hom_AMR=0
		self.Hom_EAS=0
		self.Hom_FIN=0
		self.Hom_NFE=0
		self.Hom_OTH=0
		self.Hom_SAS=0
		
	def add_name(self,name):
		self.names.append(name)

	def add_names(self,name_list):
		self.names+=name_list

	def update_info(self,var_info_dict):
		self.AC=var_info_dict['AC']
		self.AC_AFR=var_info_dict['AC_AFR']
		self.AC_AMR=var_info_dict['AC_AMR']
		self.AC_Adj=var_info_dict['AC_Adj']
		self.AC_EAS=var_info_dict['AC_EAS']
		self.AC_FIN=var_info_dict['AC_FIN']
		self.AC_Het=var_info_dict['AC_Het']
		self.AC_Hom=var_info_dict['AC_Hom']
		self.AC_NFE=var_info_dict['AC_NFE']
		self.AC_OTH=var_info_dict['AC_OTH']
		self.AC_SAS=var_info_dict['AC_SAS']
		self.AN=var_info_dict['AN']
		self.AN_AFR=var_info_dict['AN_AFR']
		self.AN_AMR=var_info_dict['AN_AMR']
		self.AN_Adj=var_info_dict['AN_Adj']
		self.AN_EAS=var_info_dict['AN_EAS']
		self.AN_FIN=var_info_dict['AN_FIN']
		self.AN_NFE=var_info_dict['AN_NFE']
		self.AN_OTH=var_info_dict['AN_OTH']
		self.AN_SAS=var_info_dict['AN_SAS']
		self.Het_AFR=var_info_dict['Het_AFR']
		self.Het_AMR=var_info_dict['Het_AMR']
		self.Het_EAS=var_info_dict['Het_EAS']
		self.Het_FIN=var_info_dict['Het_FIN']
		self.Het_NFE=var_info_dict['Het_NFE']
		self.Het_OTH=var_info_dict['Het_OTH']
		self.Het_SAS=var_info_dict['Het_SAS']
		self.Hom_AFR=var_info_dict['Hom_AFR']
		self.Hom_AMR=var_info_dict['Hom_AMR']
		self.Hom_EAS=var_info_dict['Hom_EAS']
		self.Hom_FIN=var_info_dict['Hom_FIN']
		self.Hom_NFE=var_info_dict['Hom_NFE']
		self.Hom_OTH=var_info_dict['Hom_OTH']
		self.Hom_SAS=var_info_dict['Hom_SAS']

	def get_af_adj(self):
		return float(self.AC_Adj)/float(self.AN_Adj)

	def get_af_hom_adj(self):
		return float(self.AC_Hom)/float(self.AC_Hom+self.AC_Het)

	def get_af_het_adj(self):
		return float(self.AC_Het)/float(self.AC_Hom+self.AC_Het)

	def get_hetf(self):
		if int(self.AC_Hom) == 0:
			return 1
		elif int(self.AC_Het) == 0:
			return 0
		else:
			return float(self.AC_Het)/float(self.AC_Hom+self.AC_Het)	

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

class Vcf_variant_row(Vcf_row):
	"""docstring for Vcf_variant_row"""
	def __init__(self, row_list, offset):
		super(Vcf_variant_row, self).__init__(row_list)
		self.tPOS = int(self.tPOS)-int(offset)
		self.var_dict = {}

	def var_dict_empty(self):
		if len(self.var_dict) > 0:
			return False
		else:
			return True

	def compile_info_dict(self,allele_id):
		var_info_list = self.tINFO.strip().split(';')
		var_info_dict = {}
		el_name = ''
		el_list = []

		for el in var_info_list:
			if '=' not in el:
				continue

			#test_block
			# print '########################################'
			# print 'el: '+str(el)
			# print 'el split: '+str(el.split('=')[1].strip().split(','))
			# print '########################################'

			el_name = el.split('=')[0]
			el_list = el.split('=')[1].split(',')
			if el_name[0:2] == 'AC': 
				var_info_dict[el_name] = el_list[allele_id]	
			elif el_name[0:2] == 'AN': 
				var_info_dict[el_name] = el_list[0]	
			elif el_name[0:3] == 'Het': 
				var_info_dict[el_name] = el_list[allele_id]	
			elif el_name[0:3] == 'Hom': 
				var_info_dict[el_name] = el_list[allele_id]	

			if el_name == 'Hom_SAS':
				break												
		return var_info_dict

	def compile_vep_obj(self,allele_id):
		var_info_list = self.tINFO.strip().split(';')
		vep_allele_dict = {}
		vep_info_ind = 0
		vep_el_len = 58
		allele_ind = 17

		#register vep data
		for i in range(0,len(var_info_list)):
			if var_info_list[i][0:4] == 'CSQ=': 
				vep_info_ind = i
				break

		list_temp = []
		vep_temp = object
		vep_info_list = var_info_list[vep_info_ind].strip().split('|')
		for i in range(1,len(vep_info_list)):			
			list_temp+=[vep_info_list[i]]
			if (i%vep_el_len) == 0:

				#test block
				# print '########################################'
				# print 'list_temp: '+str(list_temp)
				# print 'allele_id: '+str(allele_id)
				# print 'list_temp[allele_ind+1]: '+str(list_temp[allele_ind])
				# print '########################################'

				if int(list_temp[allele_ind]) == (allele_id+1):
					vep_temp = Vep(list_temp)

					#test block
					# print '########################################'
					# print 'vep_temp: '+str(vep_temp)					
					# print 'vep_temp.Consequence: '+str(vep_temp.Consequence)
					for el in vep_temp.Consequence:

						# print 'el: '+str(el)

						if el not in vep_allele_dict:
							vep_allele_dict[el] = [vep_temp]
						else:
							vep_allele_dict[el]+=[vep_temp]
					# print '########################################'
				
				list_temp = []	
		return vep_allele_dict

	def get_snp(self,tREF_list,tALT_list,allele_id):
		#snp

		########################################
		#checks for structure of snp
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
		return [self.tPOS, 'SNP', tREF_list[0][0], tALT_list[allele_id][0]]

	def compile_var_specs(self):
		########################################
		#currently only handling SNP, INDEL
		#not currently handling long INS/DEL
		#future versions should include all ##ALT variants
		########################################

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

		#test block
		print '##self.tPOS: '+str(self.tPOS) 

		for allele_id in range(0,len(tALT_list)):	

			#triage INDEL, SNP
			if len(tREF_list[0]) == len(tALT_list[allele_id]): 
				OUT = self.get_snp(tREF_list,tALT_list,allele_id)
				if len(OUT)==4:
					[var_loc, var_type, var_ref, var_seq] = OUT

					var_idstr = (str(var_loc)+str(var_type)+str(var_ref)+str(var_seq)+str(allele_id))
					self.var_dict[var_idstr] = ([var_loc,var_type,var_ref,var_seq,allele_id,{},[]])
					self.var_dict[var_idstr][5] = self.compile_info_dict(allele_id)
					self.var_dict[var_idstr][6] = self.compile_vep_obj(allele_id)										

	def print_var_specs(self):
		if not self.var_dict_empty():
			for idstr in self.var_dict:
				var = self.var_dict[idstr]
				print 'location: '+str(var[0])
				print 'type: '+str(var[1])
				print 'reference: '+str(var[2])
				print 'sequence: '+str(var[3])
				print 'allele: '+str(var[4])
				print 'info: '+str(var[5])

	def get_variants(self):
		variants = []
		if not self.var_dict_empty():
			for idstr in self.var_dict:
				var = self.var_dict[idstr]
				variant = Variant(var[0], var[1], var[2], var[3], var[4])				

				variant.update_info(var[5])
				variant.veps = var[6]
				variants.append(variant)
		return variants



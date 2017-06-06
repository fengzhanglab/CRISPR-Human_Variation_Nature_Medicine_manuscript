#David Scott, MIT, 2016
import sys
import gzip

class GFF3_row(object):
	"""docstring for GFF3_row"""
	def __init__(self, GFF3_line):
		super(GFF3_row, self).__init__()
		row_list = GFF3_line.strip().split()
		self.CHROM = str(row_list[0])
		self.SOURCE = str(row_list[1])
		self.TYPE = str(row_list[2])
		self.LB = int(row_list[3])
		self.UB = int(row_list[4])
		self.DOT1 = str(row_list[5])
		self.STR = str(row_list[6])
		self.DOT2 = str(row_list[7])
		self.INFO = self.parse_info(row_list[8])
		self.PCT = False
		if 'transcript_type' in self.INFO:
			if self.INFO['transcript_type'] == 'protein_coding':
				self.PCT = True

		self.EU = False
		if ((self.TYPE == 'exon') or (self.TYPE == 'UTR')):
			self.EU = True

	def parse_info(self,info_str):
		info_list = info_str.strip().split(';')
		info_dict = {}
		el_list = []
		for el in info_list:
			el_list = el.split('=')
			info_dict[el_list[0]] = el_list[1]
		return info_dict

	def get_exon_list(self):
		if not (self.TYPE == 'exon'):
			print 'error: get_exon_list called for non exon'
			print self.TYPE
			print 'exiting'
			exit(1)

		out = [self.CHROM]
		out+=[self.INFO['gene_id']]
		out+=[self.TYPE]
		out+=[self.LB]
		out+=[self.UB]			
		out+=[self.STR]
		out+=[self.INFO['transcript_id']]
		out+=[self.INFO['transcript_name']]
		out+=[self.INFO['exon_id']]
		return out

	def __str__(self):
		out = [self.CHROM]
		out+=[self.SOURCE]
		out+=[self.TYPE]
		out+=[self.LB]
		out+=[self.UB]			
		out+=[self.STR]
		return str(out)

class PCT(object):
	"""docstring for PCT"""
	def __init__(self, GFF3_row):
		super(PCT, self).__init__()
		self.ID = GFF3_row.INFO['transcript_id']
		self.CHROM = GFF3_row.CHROM
		self.LB = GFF3_row.LB
		self.UB = GFF3_row.UB
		self.STR = GFF3_row.STR
		self.utrs = []
		self.exons = []
		self.coding = []		
		if GFF3_row.TYPE == 'UTR':
			self.utrs+=[GFF3_row]
		elif GFF3_row.TYPE == 'exon':
			self.exons+=[GFF3_row]

	def add_element(self,GFF3_row):
		if not (GFF3_row.LB <= GFF3_row.UB):
			print 'error: GFF3_row.LB greater than GFF3_row.UB'
			print GFF3_row
			print 'exiting'
			exit(1)

		if self.LB > GFF3_row.LB:
			self.LB = GFF3_row.LB
		if self.UB < GFF3_row.UB:
			self.UB = GFF3_row.UB

		if GFF3_row.TYPE == 'UTR':
			self.utrs+=[GFF3_row]
		elif GFF3_row.TYPE == 'exon':
			self.exons+=[GFF3_row]		

	def compile_coding(self):
		self.exons.sort(key=lambda x: x.LB)
		for i in range(len(self.exons)):
			exon = self.exons[i]
			for j in range(len(self.utrs)):
				if exon.LB == self.utrs[j].LB:
					exon.LB = self.utrs[j].UB+1
				if exon.UB == self.utrs[j].UB:
					exon.UB = self.utrs[j].LB-1			

			if (exon.UB - exon.LB) > 0:
				self.coding+=[exon]

	def get_coding(self):
		self.compile_coding()

		out = []
		for i in range(len(self.coding)):
			out+=[self.coding[i].get_exon_list()] 
		return out

	def __str__(self):
		out = [self.ID]
		out+=[self.CHROM]
		out+=[self.LB]
		out+=[self.UB]			
		out+=[self.STR]
		out+=[[str(el) for el in self.utrs]]
		out+=[[str(el) for el in self.exons]]
		out+=[[str(el) for el in self.coding]]	
		return str(out)


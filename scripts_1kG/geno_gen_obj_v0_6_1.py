#David Scott, MIT, 2016
import sys
import gzip
import vmb

class Geno_seq(object):
	"""docstring for Geno_seq"""
	def __init__(self, sequence, index, allele, names):
		super(Geno_seq, self).__init__()
		self.seq = str(sequence)
		self.ind = int(index)
		self.allele = int(allele)
		self.names = names

	def __str__(self):
		string = str(self.ind)
		for name in self.names:
			string+='|'+str(name)
		return string

class Geno_genome(object):
	"""docstring for Geno_genome"""
	def __init__(self, reference):
		super(Geno_genome, self).__init__()
		self.ref = reference
		self.var_track = [[] for _ in range(len(reference))]

	def add_variant(self,variant):
		self.var_track[varaiant.loc-1].append(variant)

	def add_variants(self,variant_list):
		for variant in variant_list:
			#var_loc must be -1 (1-justified loc not index)
			if (variant.loc-1) < len(self.var_track):
				try:
					self.var_track[variant.loc-1].append(variant)
				except Exception, err:
					#exc block
					print '########################################'				
					print 'len self.var_track: '+str(len(self.var_track))
					print 'variant.loc: '+str(variant.loc)
					print '########################################'
					print Exception, err
					sys.exit(1)

	def get_ntvar_counts(self,ref_lb,ref_ub):
		out = {'A':[0,0],'C':[0,0],'G':[0,0],'T':[0,0]}
		for i in range(ref_lb,ref_ub):
			if self.ref[i] in out:
				#total base count
				out[self.ref[i]][1]+=1
				if len(self.var_track[i]) > 0:
					#base has variation
					out[self.ref[i]][0]+=1
		return out

	def get_nt_counts(self,ref_lb,ref_ub):
		out = {'A':0,'C':0,'G':0,'T':0}
		for i in range(ref_lb,ref_ub):
			if self.ref[i] in out:
				out[self.ref[i]]+=1
		return out

class Geno_slice(object):
	"""docstring for Geno_slice"""
	def __init__(self, geno_genome):
		super(Geno_slice, self).__init__()
		self.ref = geno_genome.ref
		self.var_track = geno_genome.var_track

		print '####################'
		print 'vmb memory:' + str(vmb.memory())
		print 'vmb resident: ' + str(vmb.resident())
		print '####################'
		
		self.ind = 0
		self.len = 0
		self.name_dict = {}
		self.geno_dict = {}

	def slice(self, index, len_slice):
		self.ind = index+len_slice
		self.len = len_slice
		self.name_dict = {}
		self.geno_dict = {}

		#add variants in initial var track to name_dict
		for var_list in self.var_track[(self.ind-self.len):self.ind]: 
			if len(var_list) > 0:
				for el in var_list:
					self.append_name_dict(el)

	def ref_base_shift(self):
		#only goes forward, not back
		var_list_pop = self.var_track[self.ind-self.len]

		if len(var_list_pop) > 0:
			for el in var_list_pop:
				self.pop_name_dict(el)

		if len(self.var_track[self.ind]) > 0:
			for el in self.var_track[self.ind]:
				self.append_name_dict(el)

		self.ind+=1

	def pop_name_dict(self,variant):

		#test block
		# print '########################################'
		# print self.ind
		# print self.name_dict
		# print variant

		for name in variant.names:
			ind = 0

			# print name
			# print variant.idstr

			for el in self.name_dict[name]:

				# print el.idstr
				# print '########################################'				
				if el.idstr == variant.idstr:
					self.name_dict[name].pop(ind)
					break
				ind+=1
			
			if len(self.name_dict[name]) == 0:
				del self.name_dict[name] 

	def append_name_dict(self,variant):
		for name in variant.names:
			if name not in self.name_dict:
				self.name_dict[name] = [variant]
			else:
			 	self.name_dict[name].append(variant)

		#test block
		# print '########################################'
		# print 'name_dict: '
		# print self.name_dict
		# print '########################################'

	def compile_geno_dict(self):
		self.geno_dict = {}
		for name in self.name_dict:
			self.name_dict[name].sort(key=lambda x: x.idstr)
			allele = self.name_dict[name][0].allele
			idstr = ''
			ct = [0,0]
			for el in self.name_dict[name]:
				if el.allele == allele:
					idstr+=el.idstr
					ct[1]+=1
				else:
					if idstr not in self.geno_dict:
						self.geno_dict[idstr] = (
							[self.name_dict[name][ct[0]:ct[1]],[]])
						self.geno_dict[idstr][1].append(name)
					else:
						self.geno_dict[idstr][1].append(name)
					allele = el.allele
					idstr = el.idstr
					ct[0] = ct[1]
					ct[1]+=1

			if idstr not in self.geno_dict:
				self.geno_dict[idstr] = (
					[self.name_dict[name][ct[0]:ct[1]],[]])
				self.geno_dict[idstr][1].append(name)
			else:
				self.geno_dict[idstr][1].append(name)

		#test block
		# print '########################################'
		# print 'geno_dict: '
		# print self.geno_dict
		# print '########################################'

	def compile_geno_seqs(self):
		geno_seqs = []
		self.compile_geno_dict()
		for idstr in self.geno_dict:
			self.geno_dict[idstr][0].sort(key=lambda x: x.loc)
			seq = list(self.ref[(self.ind-self.len):(self.ind)])
			seq_vars = self.geno_dict[idstr][0]
			ind_append = self.ind
			offset = 0			

			for var in seq_vars:
				#var_loc must be -1 (1-justified loc not index)
				ind = (var.loc-1) - (self.ind-self.len) + offset			

				if var.type == 'SNP':

					seq[ind] = var.seq

					if not (''.join(seq[ind:(ind+len(var.seq))]) == var.seq):
						print 'error: post seq mod snp mismatch! Exiting!'
						print '########################################'
						print 'idstr: '+idstr
						print 'reference: '+''.join(seq)
						print 'index: '+str(ind)
						print 'variant type: '+var.type
						print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
						print 'ALT allele: '+var.seq
						print '########################################'
						exit(1)

				elif (var.type == 'INDEL:INS') or (var.type == 'INDEL:DEL'):		
					seq = []
					break

			if len(seq) > 0:
				for ind in range(0,(len(seq)-self.len)+1):
					geno_seqs.append(Geno_seq(''.join(seq[ind:ind+self.len]),self.ind-self.len,
						self.geno_dict[idstr][0][0].allele,self.geno_dict[idstr][1]))

		return geno_seqs



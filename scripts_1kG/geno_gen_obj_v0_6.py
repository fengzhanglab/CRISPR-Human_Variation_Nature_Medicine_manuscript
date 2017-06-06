#David Scott, MIT, 2016
import sys
import gzip

class Geno_seq(object):
	"""docstring for Geno_seq"""
	def __init__(self, sequence, index, names):
		super(Geno_seq, self).__init__()
		self.seq = str(sequence)
		self.ind = int(index)
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

class Geno_slice(object):
	"""docstring for Geno_slice"""
	def __init__(self, geno_genome):
		super(Geno_slice, self).__init__()
		self.ref = geno_genome.ref
		self.var_track = geno_genome.var_track
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
		self.var_track.append(self.var_track[self.ind])

		if len(var_list_pop) > 0:
			for el in var_list_pop:
				self.pop_name_dict(el)

		if len(self.var_track[self.ind]) > 0:
			for el in self.var_track[self.ind]:
				self.append_name_dict(el)

		#test block
		print self.ind
		if len(self.var_track[self.ind]) > 0:
			print '########################################'
			# print 'var_track: '+str(self.var_track)
			print 'appended var_list: '+str(self.var_track[self.ind])
			print 'popped var_list: '+str(var_list_pop)		
			print '########################################'

		self.ind+=1

	def pop_name_dict(self,variant):
		for name in variant.names:
			ind = 0
			for el in self.name_dict[name]:
				if el.idstr == variant.idstr:
					self.name_dict[name].pop(ind)
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
			seqs = [list(self.ref[(self.ind-self.len):(self.ind)])]
			offset = 0			
			ind_append = self.ind

			#test block
			print '########################################'
			print 'self.ref len: '+str(len(self.ref))
			print 'self.ind: '+str(self.ind)
			print 'idstr: '+idstr
			print 'idstr_names: '+str(self.geno_dict[idstr][1])
			if len(self.geno_dict[idstr][0]) > 1:
				print '########## MULT VAR ##########'
			print '########################################'

			seq_vars = self.geno_dict[idstr][0]
			for var in seq_vars:
				#var_loc must be -1 (1-justified loc not index)
				ind = (var.loc-1) - (self.ind-self.len) + offset

				#test block
				if (ind >= len(seq)) or (ind < 0):
					print 'error: index out of bounds!'
					print '########################################'
					print 'idstr: '+idstr
					print 'reference: '+''.join(seq)
					print 'index: '+str(ind)
					print 'self.ind: '+str(self.ind)
					print 'offset: '+str(offset)
					print 'len(seq): '+str(len(seq))
					print 'variant loc: '+str(var.loc)
					print 'variant type: '+var.type
					print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
					print 'ALT allele: '+var.seq
					print '########################################'

				#test block
				print 'index calculation'
				print '########################################'
				print '(var.loc-1): '+str(var.loc-1)
				print '(self.ind-self.len): '+str(self.ind-self.len)
				print 'offset: '+str(offset)
				print 'index: '+str(ind)
				print '########################################'				

				if var.type == 'SNP':

					for i in range(0,len(seqs)):

						#test block
						# print 'before seq mod'
						# print '########################################'
						# print 'idstr: '+idstr
						# print 'reference: '+''.join(seqs[i])
						# print 'index: '+str(ind)
						# print 'variant type: '+var.type
						# print 'REF allele: '+''.join(seqs[i][ind:(ind+len(var.seq))])
						# print 'ALT allele: '+var.seq
						# print '########################################'
	 
						seq[i][ind] = var.seq

						if not (''.join(seq[ind:(ind+len(var.seq))]) == var.seq):
							print 'error: post seq mod SNP mismatch! Exiting!'
							print '########################################'
							print 'idstr: '+idstr
							print 'reference: '+''.join(seq)
							print 'index: '+str(ind)
							print 'variant type: '+var.type
							print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
							print 'ALT allele: '+var.seq
							print '########################################'
							exit(1)						

						#test block
						# print 'after seq mod'
						# print '########################################'
						# print 'idstr: '+idstr
						# print 'reference: '+''.join(seqs[i])
						# print 'index: '+str(ind)
						# print 'variant type: '+var.type
						# print 'REF allele: '+''.join(seqs[i][ind:(ind+len(var.seq))])
						# print 'ALT allele: '+var.seq
						# print '########################################'

				elif var.type == 'INDEL:DEL':
					#v0_6 forward, all deletions should be len 1

					#skip sequences starting with a deletion
					#	any ensuing sequence will be sampled after passing del
					if ind == 0:
						seqs = []
						break

					#rules for deletion handling
					#v0_6 forward, all deletions should be len 1
					#1. sequence should be maintained at length self.len
					#2. append ensuing bases after del to end of seq
					#3. if no variant for allele is encountered in ensuing
					#	then just append next reference base
					#4. keep track of ensuing base index for cases of mult del
					#5. if variant is encountered for allele in ensuing bases:
					#	that matches all names associated with current allele
					#	then append variable to seq_vars for allele										
					#STILL NEED TO DO:
					#if variant is encountered for in ensuing bases:
					#	that does not match all names associated with current allele
					#	branch alleles for new allelic variant combinations created
					#NOTE: currently these sequences are thrown out
					#	prevents introduction of false positive
					#	possibly creating false negative
					seq_append = ''
					while len(seq_append) < 1:
						if len(self.var_track[ind_append]) > 0:
							for var_temp in self.var_track[ind_append]:

								#test block
								print '########################################'
								print 'variant type: '+var_temp.type
								print 'variant loc: '+str(var_temp.loc)
								print '########################################'

								if not (var_temp.allele == var.allele):
									seq_append = self.ref[ind_append]
									continue

								#greedy implementation, need to streamline
								ct_match = 0
								allele_match = True
								var_temp_names_dict = {}
								for name_temp in var_temp.names:
									var_temp_names_dict[name_temp] = ''

								#check how many name match between the variant
								#	and names in the current allele
								for name_temp in self.geno_dict[idstr][1]:
									if name_temp in var_temp_names_dict:
										ct_match+=1

								#if name matches between var and curr allele equal
								#	then the variant is part of current allele
								if ct_match == len(self.geno_dict[idstr][1]):
									if not var_temp.type == 'INDEL:DEL':
										seq_append = self.ref[ind_append]
										seq_vars+=[var_temp]
										break
								#if name matches between var and curr allele > 0
								#	but not equal, then the variant branches allele										
								elif ct_match > 0:
									#variant is encountered for in ensuing bases:
									#	that does not match all curr allele names
									#	alleles need to be branched, delete curr allele
									seq = ''
									break
								else:
									seq_append = self.ref[ind_append]																	
						else:
							seq_append = self.ref[ind_append]
						ind_append+=1

						if len(seq) == 0:
								#variant is encountered for in ensuing bases:
								#	that does not match all curr allele names
								#	alleles need to be branched, delete curr allele
								break						

					if len(seq) == 0:
						print '########################################'
						print 'variant is encountered for in ensuing bases:'
						print '    that does not match all curr allele names'
						print '    alleles need to be branched, delete curr allele'
						print '########################################'
						break

					seq+=list(seq_append)
					
					#test block
					print 'before seq mod'
					print '########################################'
					print 'idstr: '+idstr
					print 'reference: '+''.join(seq)
					print 'index: '+str(ind)
					print 'variant type: '+var.type
					print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
					print 'ALT allele: '+var.seq
					print 'self.ind: '+str(self.ind)
					print 'append_seq: '+seq_append
					print '########################################'

					if not (''.join(seq[ind:(ind+len(var.seq))]) == var.seq):
						print 'error: pre seq mod INDEL:DEL mismatch! Exiting!'
						print '########################################'
						print 'idstr: '+idstr
						print 'reference: '+''.join(seq)
						print 'index: '+str(ind)
						print 'variant type: '+var.type
						print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
						print 'ALT allele: '+var.seq
						print '########################################'
						exit(1)							

					del seq[ind:(ind+len(var.seq))]
					offset-=len(var.seq)

					#test block
					print 'after seq mod'
					print '########################################'
					print 'idstr: '+idstr
					print 'reference: '+''.join(seq)
					print 'index: '+str(ind)
					print 'variant type: '+var.type
					print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
					print 'ALT allele: '+var.seq
					print '########################################'

				elif var.type == 'INDEL:INS':

					for i in range(0,len(seqs)):
						#test block
						# print 'before seq mod'
						# print '########################################'
						# print 'idstr: '+idstr
						# print 'reference: '+''.join(seqs[i])
						# print 'index: '+str(ind)
						# print 'variant type: '+var.type
						# print 'REF allele: '+''.join(seqs[i][ind:(ind+len(var.seq))])
						# print 'ALT allele: '+var.seq
						# print '########################################'

						seqs[i][(ind+1):(ind+1)] = list(var.seq)

						if not (''.join(seq[ind:(ind+len(var.seq))]) == var.seq):
							print 'error: post seq mod INDEL:INS mismatch! Exiting!'
							print '########################################'
							print 'idstr: '+idstr
							print 'reference: '+''.join(seq)
							print 'index: '+str(ind)
							print 'variant type: '+var.type
							print 'REF allele: '+''.join(seq[ind:(ind+len(var.seq))])
							print 'ALT allele: '+var.seq
							print '########################################'
							exit(1)						

						#test block
						# print 'after seq mod'
						# print '########################################'
						# print 'idstr: '+idstr
						# print 'reference: '+''.join(seqs[i])
						# print 'index: '+str(ind)
						# print 'variant type: '+var.type
						# print 'REF allele: '+''.join(seqs[i][(ind+1):(ind+len(var.seq)+1)])
						# print 'ALT allele: '+var.seq
						# print '########################################'				
			
					offset+=len(var.seq)

			if len(seqs) > 0:
				for seq in seqs:
					#given insertion, can have sequence of length
					#	greater than self.len, and this seq is
					#	broken down into subseqs of length self.len for return
					for ind in range(0,(len(seq)-self.len)+1):
						geno_seqs.append(Geno_seq(''.join(seq[ind:ind+self.len]),
							self.ind-self.len,self.geno_dict[idstr][1]))

		return geno_seqs



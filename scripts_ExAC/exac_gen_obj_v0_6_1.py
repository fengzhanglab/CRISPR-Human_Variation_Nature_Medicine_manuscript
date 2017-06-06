#David Scott, MIT, 2016
import sys
import gzip
import vmb

def revcomp(seq):
	seq = seq.lower()
	seq = seq.replace('a','T')
	seq = seq.replace('t','A')
	seq = seq.replace('c','G')
	seq = seq.replace('g','C')
	seq = seq.upper()
	return seq[::-1]

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
		self.mask_track = [0 for _ in range(len(reference))] 

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

	def add_mask(self,loc,offset):
		self.mask_track[(loc-offset)-1] = 1

	def get_ntvar_counts(self,ref_lb,ref_ub):
		out = {'A':[0,0],'C':[0,0],'G':[0,0],'T':[0,0],'Gnot':[0,0],'Cnot':[0,0]}
		for i in range(ref_lb,ref_ub):
			if self.mask_track[i]:
				#total base count
				out[self.ref[i]][1]+=1
				if len(self.var_track[i]) > 0:
					#base has variation
					out[self.ref[i]][0]+=1
				if i>ref_lb:
					if ((self.ref[i] == 'G') and 
						(not (self.ref[i-1] == 'C'))):
						out['Gnot'][1]+=1
						if len(self.var_track[i]) > 0:
							out['Gnot'][0]+=1
				if i<(ref_ub-1):
					if ((self.ref[i] == 'C') and 
						(not (self.ref[i+1] == 'G'))):
						out['Cnot'][1]+=1
						if len(self.var_track[i]) > 0:
							out['Cnot'][0]+=1						
		return out

	def get_nt_counts(self,ref_lb,ref_ub):
		out = {'A':0,'C':0,'G':0,'T':0}
		for i in range(ref_lb,ref_ub):
			if self.mask_track[i]:
				out[self.ref[i]]+=1
		return out

class Geno_CRISPR(object):
	"""docstring for Geno_genome"""
	def __init__(self, geno_genome):
		super(Geno_CRISPR, self).__init__()
		self.ref = geno_genome.ref
		self.len = len(self.ref)
		self.var_track = geno_genome.var_track
		self.mask_track = geno_genome.mask_track

	def get_del_PAMS(self,PAM_list,ref_lb,ref_ub):
		#have to separate TS and BS PAMS bc,
		#	if a TS is converted to BS, the target it not conserved
		PAM_len = len(PAM_list[0])

		PAM_dict = {}
		for pam in PAM_list:
			PAM_dict[str(pam)] = pam
		PAM_dict_rc = {}
		for pam in PAM_list:
			PAM_dict_rc[revcomp(str(pam))] = pam

		pam_inds = []
		del_vars = {}
		motif_temp = ''
		motif_list = []
		for i in range(ref_lb,(ref_ub-PAM_len)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+PAM_len)]) == PAM_len):
				continue

			motif = self.ref[i:(i+PAM_len)]

			#test block
			# print 'get_del_PAMS'
			# print '########################################'
			# print 'motif_ref: '+motif
			# print 'i: '+str(i)

			if (motif in PAM_dict):
				pam_inds+=[i]
				for j in range(0,PAM_len):

					# print 'j: '+str(j)

					for var in self.var_track[i+j]:
						if var.type == 'SNP':
							if (('synonymous_variant' in var.veps) or
								('missense_variant' in var.veps)):
								motif_list = list(motif)
								
								# print 'motif_list[j]: '+motif_list[j]
								# print 'var.ref: '+var.ref
								# print 'var.seq: '+var.seq

								if not (motif_list[j] == var.ref):
									print 'error: post seq mod snp mismatch! Exiting!'
									print '########################################'
									print 'var.idstr: '+var.idstr
									print 'reference: '+''.join(motif_list)
									print 'index: '+str(j)
									print 'variant type: '+var.type
									print 'REF allele: '+var.ref
									print 'ALT allele: '+var.seq
									print '########################################'
									exit(1)

								motif_list[j] = var.seq

								# print 'motif_alt: '+''.join(motif_list)

								if not (''.join(motif_list) in PAM_dict):

									# print 'append del_vars: '+str(var)

									if i not in del_vars:
										del_vars[i] = [var]
									else:
										del_vars[i]+=[var]

						# print '########################################'

		for i in range(ref_lb,(ref_ub-PAM_len)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+PAM_len)]) == PAM_len):
				continue			

			motif = self.ref[i:(i+PAM_len)]

			#test block
			# print 'get_del_PAMS-RC'
			# print '########################################'
			# print 'motif_ref: '+motif
			# print 'i: '+str(i)

			if (motif in PAM_dict_rc):
				pam_inds+=[i]
				for j in range(0,PAM_len):

					# print 'j: '+str(j)

					for var in self.var_track[i+j]:
						if var.type == 'SNP':
							if (('synonymous_variant' in var.veps) or
								('missense_variant' in var.veps)):
								motif_list = list(motif)
								
								# print 'motif_list[j]: '+motif_list[j]
								# print 'var.ref: '+var.ref
								# print 'var.seq: '+var.seq

								if not (motif_list[j] == var.ref):
									print 'error: post seq mod snp mismatch! Exiting!'
									print '########################################'
									print 'var.idstr: '+var.idstr
									print 'reference: '+''.join(motif_list)
									print 'index: '+str(j)
									print 'variant type: '+var.type
									print 'REF allele: '+var.ref
									print 'ALT allele: '+var.seq
									print '########################################'
									exit(1)

								motif_list[j] = var.seq

								# print 'motif_alt: '+''.join(motif_list)

								if not (''.join(motif_list) in PAM_dict_rc):

									# print 'append del_vars: '+str(var)

									if i not in del_vars:
										del_vars[i] = [var]
									else:
										del_vars[i]+=[var]

						# print '########################################'

		return [pam_inds, del_vars]		

	def get_add_PAMS(self,PAM_list,ref_lb,ref_ub):
		#have to separate TS and BS PAMS bc,
		#	if a TS is converted to BS, the target it not conserved
		PAM_len = len(PAM_list[0])

		PAM_dict = {}
		for pam in PAM_list:
			PAM_dict[str(pam)] = pam
		PAM_dict_rc = {}
		for pam in PAM_list:
			PAM_dict_rc[revcomp(str(pam))] = pam

		no_pam_inds = []
		add_vars = {}
		motif_temp = ''
		motif_list = []
		for i in range(ref_lb,(ref_ub-PAM_len)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+PAM_len)]) == PAM_len):
				continue

			motif = self.ref[i:(i+PAM_len)]

			#test block
			# print 'get_add_PAMS'
			# print '########################################'
			# print 'motif_ref: '+motif
			# print 'i: '+str(i)

			if not (motif in PAM_dict):
				no_pam_inds+=[i]
				for j in range(0,PAM_len):

					# print 'j: '+str(j)

					for var in self.var_track[i+j]:
						if var.type == 'SNP':
							if (('synonymous_variant' in var.veps) or
								('missense_variant' in var.veps)):
								motif_list = list(motif)
								
								# print 'motif_list[j]: '+motif_list[j]
								# print 'var.ref: '+var.ref
								# print 'var.seq: '+var.seq

								if not (motif_list[j] == var.ref):
									print 'error: post seq mod snp mismatch! Exiting!'
									print '########################################'
									print 'var.idstr: '+var.idstr
									print 'reference: '+''.join(motif_list)
									print 'index: '+str(j)
									print 'variant type: '+var.type
									print 'REF allele: '+var.ref
									print 'ALT allele: '+var.seq
									print '########################################'
									exit(1)

								motif_list[j] = var.seq

								# print 'motif_alt: '+''.join(motif_list)

								if (''.join(motif_list) in PAM_dict):

									# print 'append add_vars: '+str(var)

									if i not in add_vars:
										add_vars[i] = [var]
									else:
										add_vars[i]+=[var]

						# print '########################################'

		for i in range(ref_lb,(ref_ub-PAM_len)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+PAM_len)]) == PAM_len):
				continue

			motif = self.ref[i:(i+PAM_len)]

			#test block
			# print 'get_add_PAMS-RC'
			# print '########################################'
			# print 'motif_ref: '+motif
			# print 'i: '+str(i)

			if not (motif in PAM_dict_rc):
				no_pam_inds+=[i]
				for j in range(0,PAM_len):

					# print 'j: '+str(j)

					for var in self.var_track[i+j]:
						if var.type == 'SNP':
							if (('synonymous_variant' in var.veps) or
								('missense_variant' in var.veps)):
								motif_list = list(motif)
								
								# print 'motif_list[j]: '+motif_list[j]
								# print 'var.ref: '+var.ref
								# print 'var.seq: '+var.seq

								if not (motif_list[j] == var.ref):
									print 'error: post seq mod snp mismatch! Exiting!'
									print '########################################'
									print 'var.idstr: '+var.idstr
									print 'reference: '+''.join(motif_list)
									print 'index: '+str(j)
									print 'variant type: '+var.type
									print 'REF allele: '+var.ref
									print 'ALT allele: '+var.seq
									print '########################################'
									exit(1)

								motif_list[j] = var.seq

								# print 'motif_alt: '+''.join(motif_list)

								if (''.join(motif_list) in PAM_dict_rc):

									# print 'append add_vars: '+str(var)

									if i not in add_vars:
										add_vars[i] = [var]
									else:
										add_vars[i]+=[var]

						# print '########################################'

		return [no_pam_inds, add_vars]	

	def get_var_targets(self,PAM_orient,PAM_list,target_len,ref_lb,ref_ub):
		PAM_len = len(PAM_list[0])

		PAM_dict = {}
		for pam in PAM_list:
			PAM_dict[str(pam)] = pam
		PAM_dict_rc = {}
		for pam in PAM_list:
			PAM_dict_rc[revcomp(str(pam))] = pam

		target_inds = []
		target_vars = {}
		motif_temp = ''
		motif_list = []
		len_slice = target_len+(2*PAM_len)
		for i in range(ref_lb,(ref_ub-len_slice)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+len_slice)]) == len_slice):
				continue

			########################################
			########################################
			if PAM_orient == 'R':
				########################################
				if self.ref[i:(i+PAM_len)] in PAM_dict_rc:

					#test block
					print 'get_var_targets-RPAM-RC'
					print '########################################'
					print 'motif_ref: '+self.ref[i:(i+PAM_len)]
					print 'i: '+str(i)

					target_inds+=[i]
					for j in range(PAM_len,PAM_len+target_len):

						print 'j: '+str(j)

						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):
									
									print 'self.ref[i+j]: '+self.ref[i+j]
									print 'var.ref: '+var.ref
									print 'var.seq: '+var.seq

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = j-(PAM_len)

									print 'var_targ_ind: '+str(var_targ_ind)
									
									if i not in target_vars:
										target_vars[i] = [[var_targ_ind,var]]
									else:
										target_vars[i]+=[[var_targ_ind,var]]

							print '########################################'

				########################################
				if self.ref[(i+len_slice-PAM_len):(i+len_slice)] in PAM_dict:

					#test block
					print 'get_var_targets-RPAM'
					print '########################################'
					print 'motif_ref: '+self.ref[(i+len_slice-PAM_len):(i+len_slice)]
					print 'i: '+str(i)

					target_inds+=[i]
					for j in range(PAM_len,PAM_len+target_len):

						print 'j: '+str(j)

						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):
									
									print 'self.ref[i+j]: '+self.ref[i+j]
									print 'var.ref: '+var.ref
									print 'var.seq: '+var.seq

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = target_len-(j-(PAM_len))-1

									print 'var_targ_ind: '+str(var_targ_ind)

									if i not in target_vars:
										target_vars[i] = [[var_targ_ind,var]]
									else:
										target_vars[i]+=[[var_targ_ind,var]]

							print '########################################'

			########################################
			########################################
			elif PAM_orient == 'L':
				########################################
				if self.ref[i:(i+PAM_len)] in PAM_dict:

					target_inds+=[i]
					for j in range(PAM_len,PAM_len+target_len):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = j-(PAM_len)
									if i not in target_vars:
										target_vars[i] = [[var_targ_ind,var]]
									else:
										target_vars[i]+=[[var_targ_ind,var]]

				########################################
				if self.ref[(i+len_slice-PAM_len):(i+len_slice)] in PAM_dict_rc:

					target_inds+=[i]
					for j in range(PAM_len,PAM_len+target_len):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = target_len-(j-(PAM_len))-1
									if i not in target_vars:
										target_vars[i] = [[var_targ_ind,var]]
									else:
										target_vars[i]+=[[var_targ_ind,var]]

		return [target_inds, target_vars]

	def get_var_targets_del_pams(self,PAM_orient,PAM_list,target_len,ref_lb,ref_ub):
		PAM_len = len(PAM_list[0])

		PAM_dict = {}
		for pam in PAM_list:
			PAM_dict[str(pam)] = pam
		PAM_dict_rc = {}
		for pam in PAM_list:
			PAM_dict_rc[revcomp(str(pam))] = pam

		target_inds = []
		target_vars = {}
		motif_temp = ''
		motif_list = []
		len_slice = target_len+(2*PAM_len)
		for i in range(ref_lb,(ref_ub-len_slice)):
			#confirm that all nucleotides in PAM motif contain coverage
			if not (sum(self.mask_track[i:(i+len_slice)]) == len_slice):
				continue

			########################################
			########################################
			if PAM_orient == 'R':
				########################################
				if self.ref[i:(i+PAM_len)] in PAM_dict_rc:
					ref_seq = revcomp(self.ref[i:(i+PAM_len+target_len)])

					target_inds+=[[i,'BS',ref_seq]]
					for j in range(0,PAM_len+target_len):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = j

									if j < PAM_len:
										motif_list = list(self.ref[i:(i+PAM_len)])
										motif_list[j] = var.seq
										if not (''.join(motif_list) in PAM_dict_rc):
											if str([i,'BS',ref_seq]) not in target_vars:
												target_vars[str([i,'BS',ref_seq])] = [[var_targ_ind,var]]
											else:
												target_vars[str([i,'BS',ref_seq])]+=[[var_targ_ind,var]]
									else:
										if str([i,'BS',ref_seq]) not in target_vars:
											target_vars[str([i,'BS',ref_seq])] = [[var_targ_ind,var]]
										else:
											target_vars[str([i,'BS',ref_seq])]+=[[var_targ_ind,var]]											

				########################################
				if self.ref[(i+len_slice-PAM_len):(i+len_slice)] in PAM_dict:
					ref_seq = self.ref[(i+PAM_len):(i+len_slice)]

					target_inds+=[[i+PAM_len,'TS',ref_seq]]
					for j in range(PAM_len,len_slice):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = (target_len+PAM_len)-(j-(PAM_len))-1

									if j >= (len_slice-PAM_len):
										motif_list = list(self.ref[(i+len_slice-PAM_len):(i+len_slice)])
										motif_list[j-(len_slice-PAM_len)] = var.seq
										if not (''.join(motif_list) in PAM_dict):
											if str([i+PAM_len,'TS',ref_seq]) not in target_vars:
												target_vars[str([i+PAM_len,'TS',ref_seq])] = [[var_targ_ind,var]]
											else:
												target_vars[str([i+PAM_len,'TS',ref_seq])]+=[[var_targ_ind,var]]
									else:
										if str([i+PAM_len,'TS',ref_seq]) not in target_vars:
											target_vars[str([i+PAM_len,'TS',ref_seq])] = [[var_targ_ind,var]]
										else:
											target_vars[str([i+PAM_len,'TS',ref_seq])]+=[[var_targ_ind,var]]

			########################################
			########################################
			elif PAM_orient == 'L':
				########################################
				if self.ref[i:(i+PAM_len)] in PAM_dict:
					ref_seq = self.ref[i:(i+PAM_len+target_len)]
					
					target_inds+=[[i,'TS',ref_seq]]
					for j in range(0,PAM_len+target_len):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = j

									if j < PAM_len:
										motif_list = list(self.ref[i:(i+PAM_len)])
										motif_list[j] = var.seq
										if not (''.join(motif_list) in PAM_dict):
											if str([i,'TS',ref_seq]) not in target_vars:
												target_vars[str([i,'TS',ref_seq])] = [[var_targ_ind,var]]
											else:
												target_vars[str([i,'TS',ref_seq])]+=[[var_targ_ind,var]]
									else:
										if str([i,'TS',ref_seq]) not in target_vars:
											target_vars[str([i,'TS',ref_seq])] = [[var_targ_ind,var]]
										else:
											target_vars[str([i,'TS',ref_seq])]+=[[var_targ_ind,var]]

				########################################
				if self.ref[(i+len_slice-PAM_len):(i+len_slice)] in PAM_dict_rc:
					ref_seq = revcomp(self.ref[(i+PAM_len):(i+len_slice)])
					
					target_inds+=[[i+PAM_len,'BS',ref_seq]]
					for j in range(PAM_len,len_slice):
						for var in self.var_track[i+j]:
							if var.type == 'SNP':
								if (('synonymous_variant' in var.veps) or
									('missense_variant' in var.veps)):

									if not (self.ref[i+j] == var.ref):
										print 'error: post seq mod snp mismatch! Exiting!'
										print '########################################'
										print 'var.idstr: '+var.idstr
										print 'index: '+str(j)
										print 'variant type: '+var.type
										print 'REF allele: '+var.ref
										print 'ALT allele: '+var.seq
										print '########################################'
										exit(1)

									var_targ_ind = (target_len+PAM_len)-(j-(PAM_len))-1

									if j >= (len_slice-PAM_len):
										motif_list = list(self.ref[(i+len_slice-PAM_len):(i+len_slice)])
										motif_list[j-(len_slice-PAM_len)] = var.seq
										if not (''.join(motif_list) in PAM_dict_rc):
											if str([i+PAM_len,'BS',ref_seq]) not in target_vars:
												target_vars[str([i+PAM_len,'BS',ref_seq])] = [[var_targ_ind,var]]
											else:
												target_vars[str([i+PAM_len,'BS',ref_seq])]+=[[var_targ_ind,var]]
									else:
										if str([i+PAM_len,'BS',ref_seq]) not in target_vars:
											target_vars[str([i+PAM_len,'BS',ref_seq])] = [[var_targ_ind,var]]
										else:
											target_vars[str([i+PAM_len,'BS',ref_seq])]+=[[var_targ_ind,var]]

		return [target_inds, target_vars]

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
					geno_seqs.append(Geno_seq(''.join(seq[ind:ind+self.len]),
						self.ind-self.len,self.geno_dict[idstr][1]))

		return geno_seqs



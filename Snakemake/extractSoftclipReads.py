import sys
import pysam 
import re 

bamfile='/home/xchen2/FOXR2/hgg.RNA.l1hs/pipeline.v4.nopos97/tmp.bam'
sclipfile='/home/xchen2/FOXR2/hgg.RNA.l1hs/pipeline.v4.nopos97/tmp.fq'

def extractSoftclipReads(bamfile,sclipfile):
	fq = ''
	file=open(bamfile,'r')
	file=pysam.AlignmentFile(file,'rb')
	sclipOut=open(sclipfile,'w')

	for read in file.fetch():
		seq = read.query_sequence
		rname = file.getrname(read.rname)
		qual = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))
		qname = read.qname
		flag = read.flag
		pos = read.pos
		read = str(read).split()
		cigar = read[5]

		num_M = re.findall("\d+[IDSM]",cigar)
		num_int = re.split("[SDMI]", cigar)
		num_int = num_int[0:-1]
		num_int = list(map(int,num_int))

		pos_ref = [1 for i in range(len(num_M))]
		pos_seq = [1 for i in range(len(num_M))]

		for i in range(len(num_M)):
			if "I" in num_M[i] or 'S' in num_M[i] :
				pos_ref[i] = 0
			if "D" in num_M[i] :
				pos_seq[i] = 0

		pre_sum_ref = [pos_ref[i]*num_int[i] for i in range(len(num_int))]
		cum_sum_ref = [sum(pre_sum_ref[:i])+pos for i in range(1, len(num_int)+1)]
		cum_sum_ref = [pos] + cum_sum_ref 

		pre_sum_seq = [pos_seq[i]*num_int[i] for i in range(len(num_int))]
		cum_sum_seq = [sum(pre_sum_seq[:i]) for i in range(1, len(num_int)+1)]
		cum_sum_seq = [0] + cum_sum_seq 

		split_seq = [seq[i:j] for i, j in zip(cum_sum_seq[:-1], cum_sum_seq[1:])]
		split_qual = [qual[i:j] for i, j in zip(cum_sum_seq[:-1], cum_sum_seq[1:])]

		for i in range(len(num_M)):
			if len(split_seq[i])>= 15 and 'S' in num_M[i] :
				fq += '@%s chr:%s_pos:%s_cigar:%s_flag:%s\n' %(qname,rname,str(cum_sum_ref[i]+1),cigar,str(flag))
				fq += '%s\n' %split_seq[i]
				fq += '+\n'
				fq += '%s\n' %split_qual[i]

	sclipOut.write(fq)
	sclipOut.close()
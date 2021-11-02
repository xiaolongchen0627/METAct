#!/hpcf/authorized_apps/rhel7_apps/python/install/3.7.0/bin/python
import sys
import pysam 
import re 
import os 
import pandas as pd

def extractSoftclipReads(bamfile,sclipfile,keepall=False,sclip_len=20):
	fq = ''
	file=open(bamfile,'r')
	file=pysam.AlignmentFile(file,'rb')
	sclipOut=open(sclipfile,'w')
	sclip_len = int(sclip_len)
	for read in file.fetch():
		seq = read.query_sequence
		rname = file.getrname(read.rname)
		qual = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))
		qname = read.qname
		flag = read.flag
		pos = read.pos
		read = str(read).split()
		cigar = read[5]
		num_M = re.findall("\d+[IDNSM]",cigar)
		num_IDS = re.findall("\d+[IDS]",cigar)
		if (len(num_IDS) >= 2) and (keepall == False ):   ######only keep reads have softclipped reads on one end and without [ID] 
			continue
		else : 
			num_int = re.split("[SDNMI]", cigar)
			num_int = num_int[0:-1]
			num_int = list(map(int,num_int))

			pos_ref = [1 for i in range(len(num_M))]
			pos_seq = [1 for i in range(len(num_M))]

			for i in range(len(num_M)):
				if "I" in num_M[i] or 'S' in num_M[i] :
					pos_ref[i] = 0
				if "D" in num_M[i] or 'N' in num_M[i] :
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
				if len(split_seq[i])>= sclip_len and 'S' in num_M[i] :
					fq += '@%s_chr:%s_pos:%s_cigar:%s_flag:%s\n' %(qname,rname,str(cum_sum_ref[i]+1),cigar,str(flag))
					fq += '%s\n' %split_seq[i]
					fq += '+\n'
					fq += '%s\n' %split_qual[i]

	sclipOut.write(fq)
	sclipOut.close()

def filter_hg38Sam(sam_file,psl_file,output_file):
    psl = open(psl_file, 'r')
    #header_list=["match","mis_Match","rep_Match","Ns","Q_Gap_count","Q_Gap_bases",
    #"T_Gap_count","T_Gap_bases","strand","Q_Name","Q_Size","Q_Start","Q_End",
    #"T_Name","T_Size","T_Start","T_End","block_Count","blockSizes","qStarts","tStarts"]
    psl = pd.read_csv(psl,sep="\t",low_memory=False)
    Q_Name_info = psl["Q_Name"].str.split("\\.",expand=True)
    psl["Q_Name_Chr"] = Q_Name_info[0]
    psl["Q_Name_Start"] = [int(i) for i in Q_Name_info[1]]
    psl["Q_Name_End"] = [int(i) for i in Q_Name_info[2]]
    psl["Q_Name_ID"] = Q_Name_info[3]
    psl["Q_Name_Strand"] = Q_Name_info[4] 
    #print(["Q_Name_Strand"])
    samfile=open(str(sam_file),'r')
    samfile=pysam.AlignmentFile(samfile,'r')
    output = pysam.AlignmentFile(output_file,'wb',template=samfile)
    for read in samfile.fetch():        
        if read.has_tag("NM") and ((read.get_tag("NM") == 1) or (read.get_tag("NM") == 0)):
            NM = read.get_tag("NM")
            rname = samfile.getrname(read.rname)
            qname = read.qname 
            pos = int(read.pos)
            read_str = str(read).split()
            cigar = read_str[5]
            read_l1hs_pos = int(re.findall("pos:(\d*)_cigar",read_str[0])[0])
            subset_psl = psl[(psl["Q_Name_Chr"] == rname) & (psl["Q_Name_Start"] >= pos - 20000) & (psl["Q_Name_End"] <= pos + 20000) &
            ( read_l1hs_pos >= psl["T_Start"] ) & (read_l1hs_pos <= psl["T_End"])]
            read.tags=read.tags+[("L1",len(subset_psl))]            
            output.write(read)
    samfile.close()
    output.close()

# ref_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/SJHGG030242_D1.bam'
# query_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030242_D1.hg38.bam'
# output_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030242_D1.hg38.tmp.bam' 

def scoring_alignment(ref_file,query_file,output_file):
ref = open(str(ref_file),'r')
query = open(str(query_file),'r')
output = open(str(output_file), 'w')

ref = pysam.AlignmentFile(ref,'rb')
query = pysam.AlignmentFile(query,'rb')
i=0
for queryRead in query.fetch():
	queryRname = query.getrname(queryRead.rname)
	queryQname = queryRead.qname
	queryPos = int(queryRead.pos)
	#print(queryRead.flag)
	
	for refRead in ref.fetch(queryRname,queryPos-50,queryPos+50):
		refQname = refRead.qname
		refPos = int(refRead.pos)
		print([i,queryRname,refQname,queryQname,queryRead.pos,refPos])
	i = i + 1
	if i==2:
		break 




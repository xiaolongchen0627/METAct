#!/hpcf/authorized_apps/rhel7_apps/python/install/3.7.0/bin/python
import sys
import pysam 
import re 
import os 
import pandas as pd

def extractSoftclipReads(bamfile,sclipfile,sclip_len=20):
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
		pos = int(read.pos)
		read = str(read).split()
		cigar = read[5]
		num_M = re.findall("\d+[IDNSM]",cigar)
		num_IDS = re.findall("\d+[IDS]",cigar)
		num_S = [int(i) for i in re.findall("(\d+)[S]",cigar)]
		if len(num_S)>=1: 
			if ((len(num_IDS) >= 2) or (min(num_S)>5 and len(num_S)>=2) ) :   ######only keep reads have softclipped reads on one end and without [ID] , allow maximum 5 reads if both end have that 
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

def filter_hg38Sam(sam_file,psl_file,output_file,focus_site = 97,target = 'L1HS' ):
    psl = open(psl_file, 'r')
    #header_list=["match","mis_Match","rep_Match","Ns","Q_Gap_count","Q_Gap_bases",
    #"T_Gap_count","T_Gap_bases","strand","Q_Name","Q_Size","Q_Start","Q_End",
    #"T_Name","T_Size","T_Start","T_End","block_Count","blockSizes","qStarts","tStarts"]
    try : 
    	focus_site = int(focus_site)
    except:
    	focus_site = 'all'
    if target.upper() == 'L1HS':
    	psl = pd.read_csv(psl,sep="\t",low_memory=False)
    	Q_Name_info = psl["Q_Name"].str.split("\\.",expand=True)
    	psl["Q_Name_Chr"] = Q_Name_info[0]
    	psl["Q_Name_Start"] = [int(i) for i in Q_Name_info[1]]
    	psl["Q_Name_End"] = [int(i) for i in Q_Name_info[2]]
    	psl["Q_Name_ID"] = Q_Name_info[3]
    	psl["Q_Name_Strand"] = Q_Name_info[4] 
    
    samfile=open(str(sam_file),'r')
    samfile=pysam.AlignmentFile(samfile,'r')
    output = pysam.AlignmentFile(output_file,'wb',template=samfile)
    for read in samfile.fetch():        
        if read.has_tag("NM") and ((read.get_tag("NM") == 1) or (read.get_tag("NM") == 0)) : 
            read_str = str(read).split()
            qual=[int(i) for i in read.query_qualities]
			qual = round(sum(qual)/len(qual),2)
            cigar = read_str[5]
            sclip = [int(i) for i in re.findall('(\d+)S',cigar)]
            Mapped_reads_len = [int(i) for i in re.findall('(\d+)M',cigar)][0]
            Qname = read.qname
            Pos = int(read.pos)
            flag = read.flag
            BP = Pos if (flag&16 == 0 ) else Pos + Mapped_reads_len ### BreakPoints from the alignments.
            if (not re.findall("[ID]",cigar)):	
            	if (len(sclip)<=1):
            		sclip_len = 0 if len(sclip) == 0 else sclip[0]
            		if sclip_len <=3 :
            			rname = samfile.getrname(read.rname)
            			pos = int(read.pos)	
            			read_l1hs_pos = int(re.findall("pos:(\d*)_cigar",read_str[0])[0])
            			if focus_site != 'all':
            				if abs(read_l1hs_pos - focus_site) > 5 :  ###only extract 97 SD for now . and allow a few ambiguity 
            					continue
            				else :
            					if target.upper() == 'L1HS':
	            					subset_psl = psl[(psl["Q_Name_Chr"] == rname) & (psl["Q_Name_Start"] >= pos - 20000) & (psl["Q_Name_End"] <= pos + 20000) & ( read_l1hs_pos >= psl["T_Start"] ) & (read_l1hs_pos <= psl["T_End"])]
	            					L1_within = psl[(psl["Q_Name_Chr"] == rname) & (pos >= psl["Q_Name_Start"]) & (pos <=psl["Q_Name_End"])]
	            					read.tags = read.tags+[("OD",":".join(subset_psl["Q_Name_ID"]))]  ######## older L1s 
	            					read.tags = read.tags+[("WI",":".join(L1_within["Q_Name_ID"]))] ###### if there is anything within one of the older L1s
	            				read.tags = read.tags+[('BP',BP)] ######breakpoints of sclip reads         
	            				read.tags = read.tags+[("AQ",qual)] ########### Average Quality of sclip reads 	            			
	            				output.write(read)
	            		else:
	            			if target.upper() == 'L1HS':
	            				subset_psl = psl[(psl["Q_Name_Chr"] == rname) & (psl["Q_Name_Start"] >= pos - 20000) & (psl["Q_Name_End"] <= pos + 20000) & ( read_l1hs_pos >= psl["T_Start"] ) & (read_l1hs_pos <= psl["T_End"])]
	            				L1_within = psl[(psl["Q_Name_Chr"] == rname) & (pos >= psl["Q_Name_Start"]) & (pos <=psl["Q_Name_End"])]
	            				read.tags = read.tags+[("OD",":".join(subset_psl["Q_Name_ID"]))]  ######## older L1s 
	            				read.tags = read.tags+[("WI",":".join(L1_within["Q_Name_ID"]))] ###### if there is anything within one of the older L1s
	            			read.tags = read.tags+[('BP',BP)] ######breakpoints of sclip reads         
	            			read.tags = read.tags+[("AQ",qual)] ########### Average Quality of sclip reads 	            			
	            			output.write(read)

    samfile.close()
    output.close()

#ref_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/SJHGG030242_D1.bam'
#query_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030242_D1.hg38.bam'
#filter_sam='/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030242_D1.hg38.tmp.sam' 
#

def scoring_alignment(ref_file,query_file,filter_sam):
	summary_file = filter_sam + '.summary.txt'
	fastq_output = re.sub('.sam$','.fq.name.txt',filter_sam) 
	
	#summary_file.write('SeqName\tL1HS.pos\treference_breakpoints\treads_support\treads_from_different_pair\thg38_reads\tsclip_Seqs\tref_Seqnames\tref_Seqs\n')


	ref = open(str(ref_file),'r')
	query = open(str(query_file),'r')

	filter_sam = open(str(filter_sam),'w')

	ref = pysam.AlignmentFile(ref,'rb')
	query = pysam.AlignmentFile(query,'rb')
	filter_sam = pysam.AlignmentFile(filter_sam,'w',template=ref)

	readSupportQueryBP = {}
	readSupportRefBP = {} 
	rcSupportQueryBP = {}
	rcSupportRefBP = {} 
	seqSupportQueryBP = {}
	seqSupportRefBP = {} 

	processedRef=[]
	chrs=["chr"+str(i+1) for i in list(range(22))]
	chrs.append('chrX')
	chrs.append('chrY')
	#chrs=['chrX']
	for queryRead in query.fetch():	
		queryRname = query.getrname(queryRead.rname)
		if queryRname in chrs:
			queryQname = queryRead.qname
			queryPos = int(queryRead.pos)
			#print(queryQname)
			queryL1pos = [int(i) for i in re.findall('pos:(\d+)_',queryQname)][0]
			queryL1cigar = [i for i in re.findall('_cigar:(.*)_flag:',queryQname)][0]			
			queryL1cigar_num = re.findall("\d+[IDNSM]",queryL1cigar)
			if len(queryL1cigar_num ) <= 4 : ##simple but useful also 
				queryCigar = str(queryRead).split()[5]
				queryMapped_reads_len = [int(i) for i in re.findall('(\d+)M',queryCigar)][0]
				
				querySeq = queryRead.query_sequence
				queryBP = int(queryRead.get_tag('BP'))
				queryPSL = str(queryRead.get_tag('OD'))
				rcSupportQueryBP[''.join([queryRname,':',str(queryBP)])]=rcSupportQueryBP.setdefault(''.join([queryRname,':',str(queryBP)]),0) + 1
				seqSupportQueryBP.setdefault(''.join([str(queryRname),':',str(queryBP)]),[]).append(str(querySeq))
				readSupportQueryBP.setdefault(''.join([str(queryRname),':',str(queryBP)]),[]).append(str(queryQname))
				
				for refRead in ref.fetch(queryRname,queryBP-100,queryBP+100): ##search for reads that witin 200 bp region of the breakpoints.
					refQname = refRead.qname
					same_read = 0
			
					if not refQname in processedRef : 
						refPos = int(refRead.pos)
						refSeq = refRead.query_sequence
						refCigar = str(refRead).split()[5]
						refMapped_reads_len = sum([int(i) for i in re.findall('(\d+)[MD]',refCigar)])
						same_read = 0  # if the reads are the same one .  
						if (queryBP > refPos) and (queryBP < refPos + refMapped_reads_len) : 
						# reads contain the Breakpoints. treat as prefiltered WT reads
							if re.match(refQname,queryQname):
								# if the reads can be perfectly mapped in the original RNA-seq data or 
								# containes few softclipped reads 
								same_read = 1 
							else :
								rcSupportRefBP[''.join([queryRname,':',str(queryBP)])]=rcSupportRefBP.setdefault(''.join([queryRname,':',str(queryBP)]),0) + 1
								seqSupportRefBP.setdefault(''.join([str(queryRname),':',str(queryBP)]),[]).append(str(refSeq))
								readSupportRefBP.setdefault(''.join([str(queryRname),':',str(queryBP)]),[]).append(str(refQname))
							num_M = re.findall("\d+[IDNSM]",refCigar)
							refRead.tags=refRead.tags+[("MH",same_read)]+[("PT",len(num_M))]+[('BP',queryBP)]+[("OD",queryPSL)]		
							filter_sam.write(refRead)
							#fastq_output.write(re)
						processedRef.append(refQname)


	filter_sam.close()
	ref.close()
	query.close()
	#summary_file = 
	#summary_file='/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030242_D1.hg38.tmp.summary.txt' 
	summary_file = open(str(summary_file),'w')
	summary_file.write('SeqName\tL1HS.pos\treference_breakpoints\treads_support\treads_from_different_pair\thg38_reads\tsclip_Seqs\tref_Seqnames\tref_Seqs\n')

	for key in rcSupportQueryBP.keys():
		if(rcSupportQueryBP[key])>=2:
			reads = [re.sub('_chr.*','',i) for i in readSupportQueryBP[key]]
			cigar = [re.findall('_cigar:(.*)_flag:',i)[0] for i in readSupportQueryBP[key]]
			pos = [int(re.findall('pos:(\d+)_',i)[0]) for i in readSupportQueryBP[key]]
			realSupport = len(list(set([re.sub('/.*','',i) for i in readSupportQueryBP[key]])))
			
			if key in readSupportRefBP:
				hg38reads = rcSupportRefBP[key]
				hg38seq = ';;'.join(list(set([re.sub('/.*','',i) for i in seqSupportRefBP[key]])))
				hg38seqname = ';;'.join(list(set([re.sub('/.*','',i) for i in readSupportRefBP[key]])))
			else:
				hg38reads = 0
				hg38seq = ''
				hg38seqname = ''
			for i in range(len(reads)) :
				summary_file.write('%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n' %(reads[i],pos[i],key,rcSupportQueryBP[key],realSupport,hg38reads,seqSupportQueryBP[key][i],hg38seqname,hg38seq))

	summary_file.close()

	fastq_output = open(str(fastq_output),'w')

	processedRef = list(set([re.sub('\/%d+','',i) for i in processedRef])) ###remove /1 and /2 of the read name 
	if len(processedRef)>= 1 :
		for qname in processedRef :
			fastq_output.write('%s\n' %qname)

	fastq_output.close()

					#print(len(processedRef))


				#if same_read ==1 :
					#summary_file.write('SeqName\tL1HS.pos\treference_chr\tbreakpoints\tref_alignment_cigar\tsclip_seq\ttotal_seq\n')
					#summary_file.write('%s\t%d\t%s\t%s\t%s\t%s\t%s\n' %(queryQname,queryL1pos,queryRname,queryBP,refCigar,querySeq,refSeq) )
					#i = i + 1
					#print(refRead.tags)
	#if i > 10:
		#break 



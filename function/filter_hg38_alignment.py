import pysam
import os 
import re
import pandas as pd

#psl="/home/xchen2/FOXR2/rmsk.L1HS.consensus.psl"
#sam="/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030851_D1.hg38.sam"
#output="/home/xchen2/FOXR2/hgg.RNA.l1hs/hg38/SJHGG030851_D1.hg38.output1.bam"

def filter_hg38Sam(samfile,psl,output):
    psl = open(psl, 'r')
    psl = pd.read_csv(psl,sep="\t",skiprows=5,names=["match","mis_Match","rep_Match","Ns","Q_Gap_count","Q_Gap_bases",
    "T_Gap_count","T_Gap_bases","strand","Q_Name","Q_Size","Q_Start","Q_End",
    "T_Name","T_Size","T_Start","T_End","block_Count","blockSizes","qStarts","tStarts"])
    Q_Name_info = psl["Q_Name"].str.split("\\.",expand=True)
    psl["Q_Name_Chr"] = Q_Name_info[0]
    psl["Q_Name_Start"] = [int(i) for i in Q_Name_info[1]]
    psl["Q_Name_End"] = [int(i) for i in Q_Name_info[2]]
    psl["Q_Name_ID"] = Q_Name_info[3]
    psl["Q_Name_Strand"] = Q_Name_info[4] 
    samfile=open(samfile,'r')
    samfile=pysam.AlignmentFile(samfile,'r')
    output = pysam.AlignmentFile(output,'wb',template=samfile)
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

filter_hg38Sam(snakemake.input,snakemake.params.psl,snakemake.output)


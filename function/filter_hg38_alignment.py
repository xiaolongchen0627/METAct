import pysam
import os 
import re
import pandas as pd

psl="/home/xchen2/FOXR2/rmsk.L1HS.consensus.psl"
sam="/home/xchen2/FOXR2/hgg.RNA.l1hs/pipeline.v4.nopos97/SJHGG030242_D1.hg38.sam"

def filter_hg38Sam(psl,samfile,output):
    file=open(samfile,'r')
    file=pysam.AlignmentFile(file,'r')

    psl = open(psl, 'r')
    psl = pd.read_csv(psl,sep="\t",skiprows=5,names=["match","mis_Match","rep_Match","Ns","Q_Gap_count","Q_Gap_bases",
        "T_Gap_count","T_Gap_bases","strand","Q_Name","Q_Size","Q_Start","Q_End",
        "T_Name","T_Size","T_Start","T_End","block_Count","blockSizes","qStarts","tStarts"])

    Q_Name_info = psl["Q_Name"].str.split("\\.",expand=True)
    
    psl["Q_Name_Chr"] = Q_Name_info[0]
    psl["Q_Name_Start"] = Q_Name_info[1]
    psl["Q_Name_End"] = Q_Name_info[2]
    psl["Q_Name_ID"] = Q_Name_info[3]
    psl["Q_Name_Strand"] = Q_Name_info[4]

    for read in file.fetch():



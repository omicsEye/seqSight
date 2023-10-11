import pandas as pd
import csv
import os

# file_path = "seqSight_BGC_All_Results/BGC_Frankel/SRR5930493-sam-report.tsv"
# samp1 = pd.read_csv(file_path,delimiter='\t',header=0,skiprows=1)
# temp1 = samp1[["Genome","Final Guess"]]
#
# # renamelst =  ["Genome_ID","SampleName/BestHitNumber"]
# renamelst =  ["Genome_ID","SampleName/BestHitNumber"]
# temp1.columns = renamelst
# # print(temp1)
def merge():
    # Get all file names from the result folder
    file_path = "/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/seqSight_BGC_All_Results/BGC_Frankel/"
    path_list = os.listdir(file_path)
    path_name = []
    for i in path_list:
      path_name.append(i.split("-sam")[0])
    print(sorted(path_name))

    # Set the first one as base
    samp2 = pd.read_csv(sorted(path_name)[0], delimiter='\t', header=0, skiprows=1)
    s2 = samp2[["Genome", "Final Guess"]]
    # lst = ["Genome_ID", "SRR5930493"]
    # s2.columns = lst

    # Apply for the rest others
    for i in sorted(path_name)[1:]:
        if i != ".DS_Store":
            print(i)
            samp = pd.read_csv("seqSight_BGC_All_Results/BGC_Frankel/"+i+"-sam-report.tsv", delimiter='\t', header=0,skiprows=1)
            temp = samp[["Genome","Final Guess"]]
            renamelst =  ["Genome_ID",i]
            temp.columns = renamelst
            s2 = pd.merge(s2,temp,how='outer',on=['Genome_ID'])
    # pd.set_option('display.max_rows',None)
    print(s2)
    s2.set_index('Genome_ID',inplace=True)
    # s2['count'] = s2.apply(lambda x: x.count(), axis=1)
    # s2.loc["count1"] = s2.count()

    # put the path you wish to save the output
    s2.to_csv("seqSight_BGC_All_Results/BGC_Frankel_Merged.tsv")
    return

import pandas as pd
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import ast

Pseudo_seq_df = ref_df = pd.read_csv('ref_OBPs_2.csv')
Refrence_seq_df = pd.read_csv('seq_OBP.csv')

def homology_matched_sequences(Pseudo_seq_df,Refrence_seq_df):

    data_df = Pseudo_seq_df
    data_df = Pseudo_seq_df.T
    ref_df = Refrence_seq_df
    ref_df['Pseudo_sites'] = ref_df['Pseudo_sites'].apply(ast.literal_eval)
    
    pseudo_df = pd.DataFrame(columns = ["PDB_sites","Pseudo_seq",'OBP_sites'])
    
    matrix = matlist.blosum62
    
    for OBP,row_1 in data_df.iterrows():
    
        max_score = 0
        name_OBP = OBP
        print(name_OBP)
        seq_OBP = row_1[0]
        print(seq_OBP + '\n')
        max_seq_OBPs = []
        max_ref_OBPs = []
        max_pseudo_indexs = []
        max_pseudo_OBPs = []
        PDB_IDs = []
    
        for index,row in ref_df.iterrows():
            ref_OBP = ref_df.iloc[index,1]
            alignments = pairwise2.align.globaldx(seq_OBP, ref_OBP, matrix)[0]
            score = alignments[2]
            
            if score > max_score:
                max_score = score
                max_seq_OBPs = [alignments[0]]
                max_ref_OBPs = [alignments[1]]
                max_pseudo_indexs = [ref_df.iloc[index,2]]
                max_pseudo_OBPs = [ref_df.iloc[index,3]]
                PDB_IDs = [ref_df.iloc[index,0]]
                
            elif score == max_score:
                max_seq_OBPs += [alignments[0]]
                max_ref_OBPs += [alignments[1]]
                max_pseudo_indexs += [ref_df.iloc[index,2]]
                max_pseudo_OBPs += [ref_df.iloc[index,3]]
                PDB_IDs += [ref_df.iloc[index,0]]
                
        for i in range(len(PDB_IDs)):
            
            max_ref_index = []
            max_seq_index = []
            max_index_int = 1
            max_seq_OBP = max_seq_OBPs[i]
            max_ref_OBP = max_ref_OBPs[i]
            max_pseudo_index = max_pseudo_indexs[i]
            max_pseudo_OBP = max_pseudo_OBPs[i]
            PDB_ID = PDB_IDs[i]
        
            for char in max_seq_OBP:
                        
                if char == '-':
                    
                    max_seq_index += ['-'] 
                
                else:
                    
                    max_seq_index += [max_index_int] 
                    max_index_int += 1
            
            max_index_int = 1
            
            for char in max_ref_OBP:
                        
                if char == '-':
                    
                    max_ref_index += [0] 
                
                else:
                    
                    max_ref_index += [max_index_int] 
                    max_index_int += 1
                        
            pseduo_pos = -1
            pseduo_new_seq = ''
            OBP_sites = []
            
            for p in max_pseudo_index:
            
                pos_p = max_ref_index.index(p)
                OBP_sites += [max_seq_index[pos_p]]
                
                if max_seq_OBP[pos_p] == '-':
                    pseduo_new_seq += max_ref_OBP[pos_p]
                else:
                    pseduo_new_seq += max_seq_OBP[pos_p]
                    
            print(pseduo_new_seq)
            print(max_pseudo_OBP + '\n')
            
            pseudo_df.at[name_OBP + '_' + PDB_ID, "PDB_sites"] = max_pseudo_index
            pseudo_df.at[name_OBP + '_' + PDB_ID, "Pseudo_seq"] = pseduo_new_seq
            pseudo_df.at[name_OBP + '_' + PDB_ID, "OBP_sites"] = OBP_sites
    
    



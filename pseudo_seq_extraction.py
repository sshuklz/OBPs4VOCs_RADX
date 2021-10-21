import pandas as pd
import urllib.request
import numpy as np
import os
import time

t = time.time()

PDB_list = ['1KX9','3L4L','3QME','3L4A','3D73','3D74','3FE6','3FE8','3FE9','3CYZ',
            '3CZ0','3BFA','3BFB','2H8V','3CAB','3CDN','3BJH','3CZ1','3BFH',
            '3D75','3D76','3D77','3D78','4Z39','6QQ4','3N7H','5EL2','4FQT',
            '3K1E','1N8U','1N8V','4IJ7','5DIC','3B6X','6HHE','3OGN','6P2E',
            '6OPB','6OGH','6OMW','3R72','2WC5','3S0D','3S0E','3RZS','3S0B',
            '3B88','3B87','3B86','2QDI','2GTE','1OOF','6OTL','6OII','4PT1',
            '1OOG','3Q8I','3VB1','4F7F','3V2L','6VQ5','1OOH','4INX','4INW',
            '2P71','2P70','1DQE','1ORG','2WCH','2WCJ','2WCM','2WCL','2WC6',
            '3R1O','3R1P','1P28','1OW4']

PDB_list_test = [
            '2P71','2P70','1DQE','1ORG','2WCH','2WCJ','2WCM','2WCL','2WC6',
            '3R1O','3R1P','1P28','1OW4']

def Pseudo_seq_extraction(PDB_ID, 
                          Method = 'Distance', 
                          Ions = False, 
                          Sorting = 'TSP',
                          Padding = 'X',
                          contact_sites_cat = []):
    
    print('\nextracting pseudo sequences for: ' + PDB_ID)
    
    IONS = ['CA', 'CO', 'MG', 'MN', 'FE', 'CU', 'ZN', 'FE2', 'FE3', 'FE4', 
            'LI', 'NA', 'K', 'RB', 'SR', 'CS', 'BA', 'CR', 'NI', 'FE1',
            'NI', 'RU', 'RU1', 'RH', 'RH1', 'PD', 'AG', 'CD', 'LA', 'W',
            'W1', 'OS', 'IR', 'PT', 'PT1', 'AU', 'HG', 'CE', 'PR', 'SM',
            'EU', 'GD', 'TB', 'YB', 'LU', 'AL', 'GA', 'IN', 'SB', 'TL', 'PB',
            'CL', 'IOD', 'BR','MO', 'RE', 'HO', 'NH4', 'MOH','EDO']
    
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
             'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), 
                (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), 
                (60, 66), (76, 78),(78, 80)]
    
    names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid',
             'resseq', 'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor',
             'element', 'charge']
    
    if not (Method == 'PLIP') | (Method == 'Distance') | (Method == 'Test'):
        
        raise ValueError('Method must be: "PLIP", or "Distance"') 
        
    if not (Ions == True) | (Ions == False):
        
        raise ValueError('Ions must be: True, or False') 
        
    if not (Sorting == 'N2C') | (Sorting == 'TSP'):
        
        raise ValueError('Sorting must be: N2C or TSP') 
    
    file = os.path.join('./PDB/', PDB_ID + '.pdb')
    
    if os.path.isdir('./PDB') == False:
        
        print('Creating PDB folder')
        os.mkdir('./PDB')
    
    if not os.path.isfile('./PDB/' + PDB_ID + '.pdb'):
        
        print('Downloading PDB for ' + PDB_ID)
        url = 'https://files.rcsb.org/download/' + PDB_ID + '.pdb'
        urllib.request.urlretrieve(url, file)
    
    pdb_df = pd.read_fwf(file, names=names, colspecs=colspecs) 
    seq_df = pdb_df.loc[pdb_df['ATOM'] == 'ATOM']
    seq_df = seq_df.astype({'resseq': int, 'x': float,'y': float,'z': float})
    lig_df = pdb_df.loc[pdb_df['ATOM'] == 'HETATM']
    lig_df = lig_df[lig_df.resname != 'HOH']
    lig_df = lig_df.astype({'resseq': int, 'x': float,'y': float,'z': float})
    main_chain = lig_df.iloc[0,5]
    seq_df = seq_df[seq_df.chainid == main_chain]
    lig_df = lig_df[lig_df.chainid == main_chain]
    
    res_num = 0 ; AA_num = 0
    seq_index = [] ; Full_seq = ''
    
    for index, row in seq_df.iterrows():
        
        if res_num != row['resseq']:
            
            if row['chainid'] == main_chain:
            
                res_num = row['resseq']
                Full_seq += d3to1[row['resname']]
                AA_num += 1
                seq_index.append([AA_num, int(res_num)]) 
                
    seq_index_df = pd.DataFrame(seq_index, columns=['AA_num','res_num'])
    
    def pseudo_sequence_generator(contact_sites,pseudo_seq_sites):
        
        pseudo_sequence = ''
        pseudo_sequence_padded = ''
        print('\nfull_sequence: ' + str(Full_seq))
        
        if Method != 'Test':
            
            contact_sites = list(set(contact_sites))
        
            for index, row in seq_index_df.iterrows(): 
               
                if row[1] in contact_sites:
                    
                    pseudo_seq_sites += [row[0]]
                    
        else:
            
            for index, row in seq_index_df.iterrows(): 
               
                if row[0] in pseudo_seq_sites:
                    
                    contact_sites += [row[1]]
        
        if Sorting == 'N2C':
        
            pseudo_seq_sites.sort()
            print('\npseudo_sequence_sites: ' + str(pseudo_seq_sites))
            
        if Sorting == 'TSP':
        
            from functools import lru_cache
            from typing import Dict, Optional, Tuple 
            import numpy as np
        
            def nearest_neighbor(d):

                n = d.shape[0]
                idx = np.arange(n)
                path = np.empty(n, dtype=int)
                mask = np.ones(n, dtype=bool)
            
                last_idx = 0
                path[0] = last_idx
                mask[last_idx] = False
                
                for k in range(1, n):
                    
                    last_idx = idx[mask][np.argmin(d[last_idx, mask])]
                    path[k] = last_idx
                    mask[last_idx] = False
                    
                return path
            
            def two_opt(d, path, verbose=False):
                """https://en.wikipedia.org/wiki/2-opt
                """
                path = np.array(path)
            
                edges = np.stack([path[:-1], path[1:]])
                min_path_cost = np.sum(d[tuple(edges)])
                n = d.shape[0]
               
                while True:
                    found_new = False
                    for i in range(n - 1):
                        for k in range(i + 2, n + 1):
                            new_path = np.concatenate([path[:i], path[i:k][::-1], path[k:]])
                            edges = np.stack([new_path[:-1], new_path[1:]])
                            path_cost = np.sum(d[tuple(edges)])
                            if path_cost < min_path_cost:
                                if verbose:
                                    print(
                                        "Found better path ({} > {})".format(
                                            min_path_cost, path_cost
                                        )
                                    )
                                path = new_path
                                min_path_cost = path_cost
                                # Go back to outmost loop
                                found_new = True
                                break
                        if found_new:
                            break
                    if not found_new:
                        break
                return path
        
            def Held_Karp_TSP(dist_mat, maxsize: Optional[int] = None):

                set_N = frozenset(range(1, dist_mat.shape[0]))
                stored_values: Dict[Tuple, int] = {}
                @lru_cache(maxsize=maxsize)
                
                def path(node_initial: int, set_N: frozenset):
                    
                    if not set_N:
                        
                        return dist_mat[node_initial, 0]
            
                    cost = [(node_J, dist_mat[node_initial, node_J] + 
                             path(node_J, set_N.difference({node_J}))) 
                             for node_J in set_N]
                    
                    nmin, min_cost = min(cost, key=lambda x: x[1])
                    stored_values[(node_initial, set_N)] = nmin
                    return min_cost
            
                smallest_path = path(0, set_N)
                node_initial = 0 ; arrangement = [0]
                
                while set_N:
                    
                    node_initial = stored_values[(node_initial, set_N)]
                    arrangement.append(node_initial)
                    set_N = set_N.difference({node_initial})
            
                print('\nshortest path: ' + str(smallest_path) + str(arrangement))
                return arrangement
        
            def Drop_largest_edge(TSP_FC):
        
                LE = 0; BI = 0
        
                for p in range(len(TSP_FC)-1):
                    
                    edge = pseudo_dist_df.iat[TSP_FC[p],TSP_FC[p+1]]
                    print('edge in path: ' + str(edge))
                    
                    if edge > LE:
                        
                        LE = edge; BI = p
                    
                edge = pseudo_dist_df.iat[TSP_FC[-1],TSP_FC[0]]
                print('edge in path: ' + str(edge))

                if edge > LE:
                        
                    return TSP_FC
                    
                else:
                    
                    TSP_FC = TSP_FC[BI+1:] + TSP_FC[:BI+1]
                    print('drop largest edge: ' +str(TSP_FC))
                    return TSP_FC

            pseudo_seq_sites_TSP = []
            pseudo_seq_sites.sort()
            pseudo_seq_df = seq_df.loc[pdb_df['name'] == 'CA']
            pseudo_seq_df = pseudo_seq_df.iloc[:,[6,8,9,10]]
            
            pseudo_seq_df = pseudo_seq_df.loc[
                pseudo_seq_df['resseq'].isin(map(int, contact_sites))]
            
            for index1, row1 in pseudo_seq_df.iterrows(): 
           
                for index2, row2 in seq_index_df.iterrows(): 
           
                    if row2[1] == row1[0]:
                        
                        pseudo_seq_df.at[index1, 'resseq'] = row2[0]
            
            pseudo_dist_df = pd.DataFrame(index = pseudo_seq_sites,
                                          columns = pseudo_seq_sites)
            
            for site_1 in range(len(pseudo_seq_sites)):
                
                C1_site = np.array([pseudo_seq_df.iat[site_1,1],
                                   pseudo_seq_df.iat[site_1,2],
                                   pseudo_seq_df.iat[site_1,3]])
            
                for site_2 in range(len(pseudo_seq_sites)):

                    C2_site = np.array([pseudo_seq_df.iat[site_2,1],
                                        pseudo_seq_df.iat[site_2,2],
                                        pseudo_seq_df.iat[site_2,3]])
                    
                    dist = np.sqrt(np.sum((C1_site-C2_site)**2, axis=0))
                    pseudo_dist_df.iat[site_1,site_2] = dist
            
            print(pseudo_dist_df)
            pseudo_dist_arr = pseudo_dist_df.to_numpy()
            
            if len(pseudo_seq_sites) < 25:
            
                TSP_sort = Held_Karp_TSP(pseudo_dist_arr)
            
            else:
                
                TSP_sort = two_opt(pseudo_dist_arr, nearest_neighbor(pseudo_dist_arr), verbose=False)
            
            TSP_sort = Drop_largest_edge(TSP_sort)
            
            for i in TSP_sort:
                
                pseudo_seq_sites_TSP += [pseudo_seq_sites[i]]
            
            pseudo_seq_sites = pseudo_seq_sites_TSP
            print('\npseudo_sequence_sites_NN: ' + str(pseudo_seq_sites))
        
        for n in pseudo_seq_sites:
                
            pseudo_sequence += Full_seq[n - 1] 
        
        if Padding != None:
            
            pad = []
            
            for j in range(len(pseudo_seq_sites)-1):
                
                pad_count = np.floor(pseudo_dist_df.iat[j,j+1] / 3.5)
                pad += [int(pad_count) * Padding]
            
            pad += ['']; p = 0
            
            for n in pseudo_seq_sites:
                
                pseudo_sequence_padded += Full_seq[n - 1]
                pseudo_sequence_padded += pad[p]; p += 1   
            
            print('pseudo_sequence: ' + str(pseudo_sequence))
            print('pseudo_sequence_padded: ' + str(pseudo_sequence_padded))
            return [Full_seq, pseudo_seq_sites, pseudo_sequence, pseudo_sequence_padded]

        print('pseudo_sequence: ' + str(pseudo_sequence))
        return [Full_seq, pseudo_seq_sites, pseudo_sequence]
    
    def Hard_compute_method(Dmin = 4.0, max_sites = 12):
        
        print('Distance method'); sites = max_sites + 1
        
        while sites > max_sites:
        
            contact_sites = [] ; pseudo_seq_sites = []
            
            def Distance(p_site,l_site,p_res,contact_sites,Dmin):
                
                dist = np.sqrt(np.sum((p_site-l_site)**2, axis=0))
                
                if dist < Dmin:
                    
                    contact_sites += [p_res]
                
            for index1, row1 in seq_df.iterrows():
                
                p_site = np.array([row1[8], row1[9], row1[10]])
            
                for index2, row2 in lig_df.iterrows():
                
                    if Ions == False:
                        
                        if row2[4] not in IONS: 
                            
                            l_site = np.array([row2[8], row2[9], row2[10]])
                            Distance(p_site,l_site,row1[6],contact_sites,Dmin)
                        
                    else:
                        
                        l_site = np.array([row2[8], row2[9], row2[10]])
                        Distance(p_site,l_site,row1[6],contact_sites,Dmin)
             
            sites = len(list(set(contact_sites)))
            print('Max_dist: ' +  str(round(Dmin, 3)) + ' Sites found: ' + str(sites))
            
            if sites - max_sites >= 4:
                Dmin -= 0.05
            
            elif sites - max_sites >= 2:
                Dmin -= 0.02
                
            else:
                Dmin -= 0.005
                
        return pseudo_sequence_generator(contact_sites,pseudo_seq_sites)
    
    def PLIP_method():
        
        from plip.structure.preparation import PDBComplex
        
        print('PLIP method')
        mol = PDBComplex() ; mol.load_pdb(file) # Load the PDB file into PLIP class
        lig = []; contact_res = []; contact_sites = []; pseudo_seq_sites = []
        
        def binding_interaction(lig,contact_res,contact_sites):
            
            mol.analyze()
            lig = ':'.join([ligand.hetid,ligand.chain,str(ligand.position)])
            bonds = mol.interaction_sets[lig] # Contains all interactions
            
            hp = [hp.resnr for hp in bonds.hydrophobic_contacts]
            
            wb = [wb.resnr for wb in bonds.water_bridges]
            
            ha = [h.resnr for h in bonds.halogen_bonds]
            
            mc = [metalc.resnr for metalc in bonds.metal_complexes]
            
            pi_s = [pistack.resnr for pistack in bonds.pistacking]
            
            pi_c = [picat.resnr for picat in bonds.pication_paro +
                    bonds.pication_laro]
            
            hb = [hbond.resnr for hbond in bonds.hbonds_pdon +
                  bonds.hbonds_ldon]
            
            sb = [sbridge.resnr for sbridge in bonds.saltbridge_lneg +
                  bonds.saltbridge_pneg]
            
            contact_res = pi_s + pi_c + sb + hp + hb + wb + ha + mc
            contact_sites += contact_res
        
        for ligand in mol.ligands:
            
            if Ions == False:
            
                if ligand.hetid not in IONS:
                    
                    binding_interaction(lig,contact_res,contact_sites)
                    
            else:
                
                binding_interaction(lig,contact_res,contact_sites)
                
        return pseudo_sequence_generator(contact_sites,pseudo_seq_sites)
            
    if Method == 'PLIP':
    
        return PLIP_method()

    if Method == 'Distance':
    
        return Hard_compute_method()
    
    if Method == 'Test':
    
        return pseudo_sequence_generator(contact_sites = [], 
                                         pseudo_seq_sites = contact_sites_cat)
    
PDB_df = pd.DataFrame(columns=['PDB_ID','Full_seq','Pseudo_sites','Pseudo_seq','Pseudo_seq_padded'])
PDB_df_sorted = pd.DataFrame(columns=['PDB_ID','Full_seq','Pseudo_sites','Pseudo_seq','Pseudo_seq_padded'])

PDB_df_length = 0
res_match = [] 

for PDB in PDB_list:
    Result = [PDB] + Pseudo_seq_extraction(PDB_ID = PDB)
    print ('elapsed time: ' + str(time.time() - t))
    PDB_df.loc[PDB_df_length] = Result
    PDB_df_length += 1
    
PDB_df = PDB_df.sort_values(by=['Full_seq'])
PDB_df = PDB_df.reset_index()
PDB_df = PDB_df.drop(['index'], axis=1)
PDB_df_length = 0
PDB_mult = []
seq_last = ' '
    
print('making combined dataframe')

for index,row in PDB_df.iterrows():

    if index > 1:
    
        if row[1] != seq_last:
    
            res_contact = [i for n, i in enumerate(res_match) if i not in res_match[:n]]
            res_match = []
            print(res_contact)
            Result = [PDB] + Pseudo_seq_extraction(
                Method = 'Test', PDB_ID = PDB, contact_sites_cat = res_contact)
            Result[0] = PDB_mult
            PDB_df_sorted.loc[PDB_df_length] = Result
            PDB_mult = []
            res_contact = []
            PDB_df_length += 1
    
        res_match += row[2] 
        PDB_mult += [row[0]] 
        
        if index + 1 == len(PDB_df):
            
            res_contact = [i for n, i in enumerate(res_match) if i not in res_match[:n]]
            res_match = []
            print(res_contact)
            Result = [PDB] + Pseudo_seq_extraction(
                Method = 'Test', PDB_ID = PDB, contact_sites_cat = res_contact)
            Result[0] = PDB_mult
            PDB_df_sorted.loc[PDB_df_length] = Result
            PDB_mult = []
            res_contact = []
            PDB_df_length += 1
    
    else:
        
        res_match += row[2] 
        PDB_mult += [row[0]]
            
    PDB = row[0]
    seq_last = row[1]
    

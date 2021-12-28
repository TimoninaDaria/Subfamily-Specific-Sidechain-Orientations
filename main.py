import numpy as np
import os
import module
from Bio import pairwise2
import ali_count
from multiprocessing import Pool
import rmsdm
import math
import knee
import copy
from sklearn.cluster import OPTICS
from sklearn.cluster import DBSCAN
import hdbscan
import ntpath
from sklearn.metrics import silhouette_score
from scipy.spatial import distance
from tqdm import tqdm 
from tabulate import tabulate
import sys
import multiprocessing
import pickle

class Amino_acid:
    def __init__(self, name, coords, protein_name, number):
        self.coords = coords
        self.name = name
        self.protein_name = protein_name
        self.number = number
        self.calculate_start_end()
    
    def print_warning_lack_of_atom(self, *atoms):
        atoms_lack = [atom for atom in atoms if atom not in self.coords]
        if len(atoms_lack) != 0:
            print("Warning: Lack of atom(s) {} in file {} for amino acid {}".format(atoms_lack, self.protein_name, self.number)) 
        
    def calculate_start_end(self):
        if self.name == 'ALA':
            self.start = self.coords.get('CA', None)
            self.end = self.coords.get('CB', None)
            self.print_warning_lack_of_atom('CA', 'CB')
        elif self.name == 'ARG':
            self.start = self.coords.get('CD', None)
            if 'NH1' in self.coords and 'NH2' in self.coords:
                self.end = (self.coords['NH1'] + self.coords['NH2']) / 2
            else: 
                self.end = None
            self.print_warning_lack_of_atom('CD', 'NH1', 'NH2')
        elif self.name == 'ASN':
            self.start = self.coords.get('CB', None)
            if 'OD1' in self.coords and 'ND2' in self.coords:
                self.end = (self.coords['OD1'] + self.coords['ND2']) / 2
            else:
                self.end = None
            self.print_warning_lack_of_atom('CB', 'OD1', 'ND2')
        elif self.name == 'ASP':
            self.start = self.coords.get('CB', None)
            if 'OD1' in self.coords and 'OD2' in self.coords:
                self.end = (self.coords['OD1'] + self.coords['OD2']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CB', 'OD1', 'OD2')
        elif self.name == 'CYS':
            self.start = self.coords.get('CA', None)
            self.end = self.coords.get('SG', None)
            self.print_warning_lack_of_atom('CA', 'SG')
        elif self.name == 'GLN':
            self.start = self.coords.get('CG', None)
            if 'NE2' in self.coords and 'OE1' in self.coords:
                self.end = (self.coords['NE2'] + self.coords['OE1']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CG', 'NE2', 'OE1')
        elif self.name == 'GLU':
            self.start = self.coords.get('CG', None)
            if 'OE1' in self.coords and 'OE2' in self.coords:
                self.end = (self.coords['OE1'] + self.coords['OE2']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CG', 'OE2', 'OE2')
        elif self.name == 'GLY':
            self.start = self.coords.get('CA', None)
            if 'N' in self.coords and 'C' in self.coords:
                self.end = (self.coords['N'] + self.coords['C']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CA', 'N', 'C')
        elif self.name == 'HIS':
            self.start =  self.coords.get('CG', None)
            if 'CE1' in self.coords and 'NE2' in self.coords:
                self.end = (self.coords['CE1'] + self.coords['NE2']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CG', 'CE1', 'NE2')
        elif self.name == 'ILE':
            self.start = self.coords.get('CB', None)
            self.end = self.coords.get('CD1', None)
            self.print_warning_lack_of_atom('CB', 'CD1')
        elif self.name == 'LEU':
            self.start = self.coords.get('CB', None)
            if 'CD1' in self.coords and 'CD2' in self.coords:
                self.end = (self.coords['CD1'] + self.coords['CD2']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CB', 'CD1', 'CD2')
        elif self.name == 'LYS':
            self.start =  self.coords.get('CG', None)
            self.end =  self.coords.get('NZ', None)
            self.print_warning_lack_of_atom('CG', 'NZ')
        elif self.name == 'MET':
            self.start = self.coords.get('CA', None)
            self.end = self.coords.get('SD', None)
            self.print_warning_lack_of_atom('CA', 'SD')
        elif self.name == 'PHE':
            self.start = self.coords.get('CG', None)
            self.end = self.coords.get('CZ', None)
            self.print_warning_lack_of_atom('CG', 'CZ')
        elif self.name == 'PRO':
            self.start = self.coords.get('CA', None) 
            if 'CG' in self.coords and 'CD' in self.coords:
                self.end = (self.coords['CG'] + self.coords['CD']) / 2 
            else:
                self.start = None
            self.print_warning_lack_of_atom('CG', 'CA', 'CD')
        elif self.name == 'SER':
            self.start = self.coords.get('CA', None) 
            self.end = self.coords.get('OG', None)
            self.print_warning_lack_of_atom('CA', 'OG')
        elif self.name == 'THR':
            self.start = self.coords.get('CA', None) 
            if 'CG2' in self.coords and 'OG1' in self.coords:
                self.end = (self.coords['CG2'] + self.coords['OG1']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CA', 'CG2', 'OG1')
        elif self.name == 'TRP':
            self.start = self.coords.get('CD1', None)
            if 'CZ2' in self.coords and 'CZ3' in self.coords:
                self.end = (self.coords['CZ2'] + self.coords['CZ3']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CD1', 'CZ2', 'CZ3')
        elif self.name == 'TYR':
            self.start = self.coords.get('CG', None)
            self.end = self.coords.get('CZ', None)
            self.print_warning_lack_of_atom('CG', 'CZ')
        elif self.name == 'VAL': 
            self.start = self.coords.get('CA', None)
            if 'CG1' in self.coords and 'CG2' in self.coords:
                self.end = (self.coords['CG1'] + self.coords['CG2']) / 2 
            else:
                self.end = None
            self.print_warning_lack_of_atom('CA', 'CG1', 'CG2')
        else:
            self.start = None
            self.end = None

class PDB:
    def __init__(self):
        self.amino_acid_num_to_name = {}
        self.amino_acid_num_to_atom_name_to_coords = {}
        self.seq = []
        self.fasta_name = None 
    
    def glue(self, seq):
        glue_seq = []
        for i, val in enumerate(seq):
            if i == 0: glue_seq.append(val)
            else:
                if val == seq[i-1]: pass
                else: glue_seq.append(val)
        return glue_seq
    
    def name_of_at(self, line):
        return(line[12:16].replace(' ',''))

    def num_of_aa(self, line):
        return(int(line[22:26]))  

    def name_of_ch(self, line):
        return line[21:22] 
    
    def name_of_aa(self, line):
        return line[16:20] .replace(' ','')       

    def num_of_at(self, line):
        return line[6:11].replace(' ','') 
    
    def reading(self, filename):
        self.filename = os.path.abspath(filename)
        self.base_filename = ntpath.basename(os.path.abspath(filename))
        chain = 'No chain'
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("ATOM"):
                    if chain != self.name_of_ch(line) and chain != 'No chain':
                        raise Exception('Exception: One pdb-file should contain one chain. File {} contains more, than one chain'.format(filename))
                    
                    chain = self.name_of_ch(line)
                    self.amino_acid_num_to_name[self.num_of_aa(line)] = self.name_of_aa(line)
                    self.seq.append(self.num_of_aa(line))
                    
                    if self.num_of_aa(line) in self.amino_acid_num_to_atom_name_to_coords:
                        self.amino_acid_num_to_atom_name_to_coords[self.num_of_aa(line)][self.name_of_at(line)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    else:
                        self.amino_acid_num_to_atom_name_to_coords[self.num_of_aa(line)] = {}
                        self.amino_acid_num_to_atom_name_to_coords[self.num_of_aa(line)][self.name_of_at(line)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                else:
                    pass
        self.seq = self.glue(self.seq)    
        
        self.amino_acids = {}
        for amino_acid_num in self.amino_acid_num_to_name:
            self.amino_acids[amino_acid_num] = Amino_acid(self.amino_acid_num_to_name[amino_acid_num], self.amino_acid_num_to_atom_name_to_coords[amino_acid_num], self.base_filename, amino_acid_num)

class PDBS_alignment:
    def __init__(self, dir_with_aligned_pdbs, filename_fasta):
        self.dir_with_aligned_pdbs = dir_with_aligned_pdbs
        self.filename_fasta = filename_fasta
        self.molecules = []         
        self.alignment_fasta = []
        self.alignment = []
        self.dirs = []   
        
        self.get_fasta_alignment()
        self.check_if_alignment_is_ok()
    
    def check_if_alignment_is_ok(self):
        list_of_bad_backbone_aa = {}
        for i in range(len(self.alignment)):
            for col in range(len(self.alignment[0])):
                if self.alignment[i][col] != '-':
                    if all(k in self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][col]] for k in ('N', 'CA', 'C', 'O')):
                        continue
                    else:
                        if i in list_of_bad_backbone_aa:
                            list_of_bad_backbone_aa[i].append(self.alignment[i][col])
                        else:
                            list_of_bad_backbone_aa[i] = [self.alignment[i][col]]
        for i in list_of_bad_backbone_aa:
            print('Warning: PDB entry {} contains incomplete set of backbone atoms for amino acid(s)'.format(self.dirs[i]), list_of_bad_backbone_aa[i])
        
        for i in range(len(self.alignment)):
            for j in range(i, len(self.alignment)):
                rmsds = []
                for col in range(len(self.alignment[0])):
                    if self.alignment[i][col] == '-' or self.alignment[j][col] == '-':
                        continue
                   
                    if 'CA' in self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][col]] and 'CA' in self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][col]]:
                        P1 = []
                        P2 = []
                        P1.append(self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][col]]['CA'])
                        P2.append(self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][col]]['CA'])
                        rmsds.append(rmsdm.rmsd(P1, P2))
                if len(rmsds) != 0 :
                    mean_rmsd = sum(rmsds) / len(rmsds)
                    if mean_rmsd > 9:
                         print('Warning: Average RMSD between {} and {} is {} what may indicate poor alignment quality of the two proteins'.format(self.dirs[i], self.dirs[j],str(mean_rmsd)))
    
    def get_conserv_positions(self, perc_of_gaps=0.05, max_content_of_mismatch=0.05, mismatch_threshold=None):
        
        def align_analysis_rec(heap, align, cut, RMSD_MATR):
            i_max = np.unravel_index(RMSD_MATR.argmax(), RMSD_MATR.shape)[0]
            j_max = np.unravel_index(RMSD_MATR.argmax(), RMSD_MATR.shape)[1]
            max = np.max(RMSD_MATR)
            if max < cut: 
                return
            else: 
                max1 = RMSD_MATR[i_max].sum()
                max2 = RMSD_MATR[j_max].sum()
            if max1 > max2:
                align[i_max][heap] = '-'
                RMSD_MATR[i_max][:] = 0
                RMSD_MATR[:][i_max] = 0
            else: 
                align[j_max][heap] = '-' 
                RMSD_MATR[j_max][:] = 0
                RMSD_MATR[:][j_max] = 0
            align_analysis_rec(heap, align, cut, RMSD_MATR)
        
        def RMSD_MATR_one_column(heap):
            RMSD_MATR = np.zeros((len(self.alignment),len(self.alignment))) 
            for i in range(len(self.alignment)):
                for j in range(i, len(self.alignment)): 
                    if self.alignment[i][heap] == '-' or self.alignment[j][heap] == '-': 
                        RMSD_MATR[i][j] = 0
                        RMSD_MATR[j][i] = 0
                    else:
                        if all(k in self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][heap]] for k in ('N', 'CA', 'C', 'O')) and all(k in self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][heap]] for k in ('N', 'CA', 'C', 'O')):
                            P1 = []
                            P2 = []
                            P1.append(self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][heap]]['N'])
                            P1.append(self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][heap]]['CA'])
                            P1.append(self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][heap]]['C'])
                            P1.append(self.molecules[i].amino_acid_num_to_atom_name_to_coords[self.alignment[i][heap]]['O'])
                            P2.append(self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][heap]]['N'])
                            P2.append(self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][heap]]['CA'])
                            P2.append(self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][heap]]['C'])
                            P2.append(self.molecules[j].amino_acid_num_to_atom_name_to_coords[self.alignment[j][heap]]['O'])
                            RMSD_MATR[i][j] = rmsdm.rmsd(P1, P2)
                            RMSD_MATR[j][i] = RMSD_MATR[i][j]
                        else:
                            RMSD_MATR[i][j] = 0
                            RMSD_MATR[j][i] = 0   
            return RMSD_MATR
        
        def align_analysis(heap, align, cut):
            RMSD_MATR = RMSD_MATR_one_column(heap)
            align_analysis_rec(heap, align, cut, RMSD_MATR) 
        
        def get_max_rmsd_in_column(heap):
            RMSD_MATR = RMSD_MATR_one_column(heap)   
            max = np.max(RMSD_MATR)
            return max
        
        def num_of_let(alignment, col): 
            pos_with_letters = []
            for j in range(len(alignment)):
                if alignment[j][col] != '-': 
                    pos_with_letters.append(j)
            return len(pos_with_letters)
        
        def get_elbow(rmsds):
            y_axis = sorted(rmsds)
            Perc_5 = math.floor(len(y_axis) * 0.05)
            x_axis = [i for i in range(len(rmsds))]
            data = np.asarray(list(zip(x_axis[Perc_5:len(x_axis) - Perc_5],y_axis[Perc_5:len(x_axis) - Perc_5])))
            elbow_index = knee.find_elbow(data, knee.get_data_radiant(data))
            cut = data[elbow_index][1] 
            if cut == 0 or cut == y_axis[Perc_5]:
                cut = y_axis[-1]
            print('Info: The automatically selected cut-off value to discriminate spatially aligned from misaligned residues is {} A'.format(round(cut,3)))
            return cut
        
        if mismatch_threshold is None:
            heaps = []
            for col in range(len(self.alignment[0])):
                if num_of_let(self.alignment, col) >= len(self.alignment) * (1 - perc_of_gaps):
                    heaps.append(col)
            assert not len(heaps) == 0, "It seems that the number of common core residues is not sufficient for cluster analysis. You should increase the 'max_content_of_gaps' parameter or reconsider your input alignment"
            rmsds = [get_max_rmsd_in_column(heap) for heap in heaps]
            mismatch_threshold = get_elbow(rmsds)
        align_tmp = copy.deepcopy(self.alignment)
        heaps = []
        for col in range(len(self.alignment[0])):
            if num_of_let(self.alignment, col) > len(self.alignment) * (1 - max_content_of_mismatch):
                heaps.append(col)
        for heap in heaps:
            align_analysis(heap, align_tmp, mismatch_threshold)
        heaps = []
        for col in range(len(align_tmp[0])):
            if num_of_let(align_tmp, col) > len(align_tmp) * (1 - max_content_of_mismatch):
                heaps.append(col)
        assert len(heaps) > 1, "It seems that the number of common core residues is not sufficient for cluster analysis. You should increase the 'max_content_of_mismatch' parameter or reconsider your input alignment"
        return heaps
    
    def get_fasta_alignment(self):
        def convertion(a):    
            aminoac = {
            "G" : 'GLY',
            "L" : 'LEU',
            "Y" : 'TYR',
            "S" : 'SER',
            "E" : 'GLU',
            "Q" : 'GLN',
            "D" : 'ASP',
            "N" : 'ASN',
            "F" : 'PHE',
            "A" : 'ALA',
            "K" : 'LYS',
            "R" : 'ARG',
            "H" : 'HIS',
            "C" : 'CYS',
            "V" : 'VAL',
            "P" : 'PRO',
            "W" : 'TRP',
            "I" : 'ILE',
            "M" : 'MET',
            "T" : 'THR',
            "GLY" : 'G',
            "LEU" : 'L',
            "TYR" : 'Y',
            "SER" : 'S',
            "GLU" : 'E',
            "GLN" : 'Q',
            "ASP" : 'D',
            "ASN" : 'N',
            "PHE" : 'F',
            "ALA" : 'A',
            "LYS" : 'K',
            "ARG" : 'R',
            "HIS" : 'H',
            "CYS" : 'C',
            "VAL" : 'V',
            "PRO" : 'P',
            "HYP" : 'hP',
            "TRP" : 'W',
            "ILE" : 'I',
            "MET" : 'M',
            "THR" : 'T',
            }
            if a in aminoac:
                b = aminoac[a]
            else:
                b = 'X'
            return b
    
        with open(self.filename_fasta, 'r') as fasta:
            names_of_prot_in_fasta = module.names_of_prot(fasta)
            seq_from_fasta = module.seq_from_fasta(fasta)
        dirs_in_pdb_folder = os.listdir(self.dir_with_aligned_pdbs)
        assert len(names_of_prot_in_fasta) == len(dirs_in_pdb_folder), 'Error: number of proteins in fasta-file is not equal to number of proteins in pdb-files folder'
        for i, name_in_fasta in enumerate(names_of_prot_in_fasta):
            if name_in_fasta + '.pdb' in dirs_in_pdb_folder:
                dir = name_in_fasta + '.pdb'
            elif name_in_fasta.replace(':','_')+'.pdb' in dirs_in_pdb_folder:
                dir = name_in_fasta.replace(':', '_') + '.pdb'
            elif name_in_fasta.split(':')[0]+'.pdb' in dirs_in_pdb_folder:
                dir = name_in_fasta.split(':')[0]+'.pdb'
            else:
                raise Exception('Warning: Lack of file ' + names_of_prot_in_fasta[i] + '.pdb in input folder')
            self.alignment.append([])
            self.dirs.append(dir)
            self.alignment_fasta.append(seq_from_fasta[i])
            mol = PDB()
            mol.reading(os.path.join(self.dir_with_aligned_pdbs, dir))
            mol.fasta_name = name_in_fasta
            self.molecules.append(mol) 
            seq_from_pdb =''
            for j in mol.seq:
                seq_from_pdb += convertion(mol.amino_acid_num_to_name[j])
            alignment_tmp = pairwise2.align.globalxx(seq_from_fasta[i].replace("-", ""), seq_from_pdb)
            c = [i for i in range(len(seq_from_fasta[i].replace("-", "")))]
            c_res = ali_count.fasta_bind_site(c, alignment_tmp[0][0], alignment_tmp[0][1])
            k = 0
            for j in seq_from_fasta[i]:
                if j == '-':
                    self.alignment[-1].append('-')
                    continue
                elif j != '-':
                    if c_res[k] != -1:
                        self.alignment[-1].append(mol.seq[c_res[k]])
                    elif c_res[k] == -1:
                        self.alignment[-1].append('-')
                    k += 1

class One_column_clustering:
    def __init__(self, col, alignment, method='hdbscan'):
        self.col = col
        self.alignment = alignment
        self.method = method
        
    def clustering(self, input):
        clust_matr = []
        self.mols_and_aa = []
        self.names_of_aa = []
        self.ref =  self.alignment.alignment[input.ref_num][self.col]
        for i in range(len(self.alignment.alignment)):
            num_of_aa = self.alignment.alignment[i][self.col]
            if num_of_aa != '-' and self.alignment.molecules[i].amino_acids[num_of_aa].start is not None and self.alignment.molecules[i].amino_acids[num_of_aa].end is not None:
                self.names_of_aa.append(self.alignment.molecules[i].amino_acids[num_of_aa].name)
                self.mols_and_aa.append((self.alignment.molecules[i], num_of_aa)) 
                clust_matr.append(np.hstack((self.alignment.molecules[i].amino_acids[num_of_aa].start, self.alignment.molecules[i].amino_acids[num_of_aa].end)))
        clust_matr = np.array(clust_matr)
        if self.method == 'optics':
            clusterer = OPTICS(metric='euclidean', n_jobs=input.cpu_threads, min_samples=input.min_samples)
        elif self.method == 'hdbscan':
            clusterer = hdbscan.HDBSCAN(metric='euclidean', min_cluster_size=input.min_cluster_size, min_samples=input.min_samples) 
        elif self.method == 'dbscan':
            clusterer = DBSCAN(metric='euclidean', n_jobs=input.cpu_threads, eps=input.eps, min_samples=input.min_samples)
        db = clusterer.fit(clust_matr)
        self.lab = db.labels_
        self.number_of_clusters = len(set(self.lab)) - 1 if -1 in set(self.lab) else len(set(self.lab))
        self.is_ok = len(set(self.lab)) > 1
        if -1 not in set(self.lab) and len(set(self.lab)) == 2 or  -1 in set(self.lab) and len(set(self.lab)) == 3:
            self.sil = silhouette_score(clust_matr[self.lab != -1], self.lab[self.lab != -1], metric='euclidean')
            dist_matr = np.array([[distance.euclidean(clust_matr[i], clust_matr[j]) for i in range(len(clust_matr))] for j in range(len(clust_matr))])
            mean_diams_clusters = [dist_matr[self.lab == i].T[self.lab == j].mean() for i in set(self.lab) for j in set(self.lab) if i != j and i != -1 and j != -1] 
            self.diam = max(mean_diams_clusters) 
            self.pre_score = self.sil * self.diam 
        else: 
            self.sil = None
            self.diam = None
            self.pre_score = None
        
    def printing_py_result(self, input):
        for i in range(len(self.lab)):
            if i == input.ref_num:
                pref = 'Ref_'
            else:
                pref = ''
            f.write('cmd.show("sticks", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            if self.lab[i] == 0:
                f.write('cmd.color("lightpink", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 1:
                f.write('cmd.color("palegreen", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 2:
                f.write('cmd.color("lightblue", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 3:
                f.write('cmd.color("paleyellow", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 4:
                f.write('cmd.color("wheat", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 5:
                f.write('cmd.color("palecyan", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 6:
                f.write('cmd.color("lightorange", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == 7:
                f.write('cmd.color("bluewhite", "{} and resi {}")\n'.format(pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            elif self.lab[i] == -1:
                f.write('cmd.color("{}", "{} and resi {}")\n'.format(self.lab[i] + 3, pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))
            else:
                f.write('cmd.color("{}", "{} and resi {}")\n'.format(self.lab[i] + 3, pref + self.mols_and_aa[i][0].base_filename, self.mols_and_aa[i][1]))

class Input:
    def __init__(self, sys_argv):
        self.dict_of_sys_argv = dict([arg.split('=') for arg in sys_argv[1:]])
        if len(self.dict_of_sys_argv) == 0:
            print('''
Usage:   python3 main.py aligned_pdbs=</path/to/folder> aligned_fasta=</path/to/file> output=</path/to/folder> [options]
Example: python3 main.py aligned_pdbs=./input_pdbs aligned_fasta=./input.fasta output=./output

Mandatory input parameters:
===========================
aligned_pdbs=<string>  # Path to folder with aligned protein 3D-structures as separate files in the PDB format (each file should represent one chain)
aligned_fasta=<string> # Path to the corresponding sequence representation of the alignment in the FASTA format
output=<string>        # Path to folder to store results

Utilization of computing resources: 
===================================
cpu_threads=<int> # Number of parallel CPU threads to utilize (the default is "all" physically available)

Cluster analysis methods:
=========================
method=hdbscan              # Use HDBSCAN automatic method (default) 
method=optics               # Use OPTICS automatic method
method=dbscan eps=<float>   # Use DBSCAN method for manual fine-tuning of the results by specifying the ‘eps’ value (eps > 0)

Cluster size:
=========================
min_samples=<int>           # The 'min_samples' parameter of HDBSCAN, OPTICS, and DBSCAN (the number of points in a neighborhood 
                              for a point to be considered as the cluster core)
min_cluster_size=<int>      # The HDBSCAN 'min_cluster_size' parameter to regulate the minimal size of a cluster

Selection of common core positions:
===================================
max_content_of_gaps=<int>      # Define the allowed gap content in alignment column, in % (at most 5%, by default)
max_content_of_mismatch=<int>  # Define the allowed 3D-mismatch content in alignment column, in % (at most 5%, by default)
mismatch_threshold=<float>     # Define the cut-off value to discriminate spatially aligned from misaligned residues, in angstroms (selected automatically)

Printing result:
===================================
number_of_result_resids=<int>     # Number of residues with the greatest score to print in RESULT-files (10, by default) 
ref=<string>                      # Set the reference protein by its name in the FASTA alignment (default is the firs protein)

PyMol parameters:
=================
compile_pymol_pse=false # Disable compilation of PyMol sessions with 3D-annotation of results
              ''')
            quit()
        
        self.path_to_program = os.path.dirname(os.path.abspath(sys_argv[0]))
        self.pdbflex_extract()
        self.input_preparation()

    def pdbflex_extract(self):
        self.pdbflex_sils = []
        self.pdbflex_diams = []
        for i in range(0,100):
            if 'PDBFLEX_{}_sil.txt'.format(i) in os.listdir(os.path.join(os.path.join(self.path_to_program, 'PDBFlex'))):
                with open(os.path.join(os.path.join(self.path_to_program, 'PDBFlex'), 'PDBFLEX_{}_sil.txt'.format(i)), 'rb') as f:
                    self.pdbflex_sils.append(np.loadtxt(f))
                with open(os.path.join(os.path.join(self.path_to_program, 'PDBFlex'), 'PDBFLEX_{}_diam.txt'.format(i)), 'rb') as f:
                    self.pdbflex_diams.append(np.loadtxt(f))
            
    def input_preparation(self):
        assert 'aligned_fasta' in self.dict_of_sys_argv and 'aligned_pdbs' in self.dict_of_sys_argv, 'Not enough input arguments to run the program'
        self.aligned_pdbs = os.path.abspath(self.dict_of_sys_argv['aligned_pdbs'])
        self.aligned_fasta = os.path.abspath(self.dict_of_sys_argv['aligned_fasta'])
        
        assert 'output' in self.dict_of_sys_argv, 'Choose folder for output'
        self.output = os.path.abspath(self.dict_of_sys_argv['output'])

        if 'compile_pymol_pse' in self.dict_of_sys_argv:
            self.create_pse = self.dict_of_sys_argv['compile_pymol_pse'] == 'true'
        else:
            self.create_pse = True

        self.cpu_threads = int(self.dict_of_sys_argv.get('cpu_threads', multiprocessing.cpu_count()))
        
        self.method_of_clustering = self.dict_of_sys_argv.get('method', 'hdbscan')
        assert  self.method_of_clustering  in ['optics', 'dbscan', 'hdbscan'], 'Method must be optics, dbscan or hdbscan'

        
        if 'eps' in self.dict_of_sys_argv:
            self.eps_for_clustering = float(self.dict_of_sys_argv['eps'])
        else:
            self.eps_for_clustering = None
            if self.method_of_clustering == 'dbscan':
                raise Exception('If you use dbscan-method, you should choose eps')
        
        self.min_cluster_size = int(self.dict_of_sys_argv.get('min_cluster_size', 5))
        
        if self.method_of_clustering == 'hdbscan':
            self.min_samples = self.dict_of_sys_argv.get('min_samples', None)
        else:
            self.min_samples = self.dict_of_sys_argv.get('min_samples', 5)
        if self.min_samples is not None:
            self.min_samples = int(self.min_samples)

        self.perc_of_gaps = float(self.dict_of_sys_argv.get('max_content_of_gaps', 5))/100.

        self.max_content_of_mismatch = float(self.dict_of_sys_argv.get('max_content_of_mismatch', 5))/100.
   
        self.number_of_result_resids = int(self.dict_of_sys_argv.get('number_of_result_resids', 10))
        
        self.ref = self.dict_of_sys_argv.get('ref', None)

        if 'mismatch_threshold' in self.dict_of_sys_argv:
            self.mismatch_threshold = float(self.dict_of_sys_argv['mismatch_threshold'])
        else:
            self.mismatch_threshold = None 

        if os.path.exists(self.output) and os.path.isdir(self.output):
            assert not os.listdir(self.output), 'Folder {} exists and is not empty'.format(self.output)
        else:
            os.makedirs(self.output)

        assert os.path.exists(self.aligned_fasta), 'No fasta file {} in here'.format(self.aligned_fasta)
        assert os.path.exists(self.aligned_pdbs) and os.path.isdir(self.aligned_pdbs), 'No folder {} in here'.format(self.aligned_pdbs)
        
def statistic(cls, input):
    
    all_diams = np.array([cl.diam for cl in cls if cl.diam is not None])
    all_sils = np.array([cl.sil for cl in cls if cl.sil is not None])
    
    if len(all_diams) == 0:
        for cl in cls:
            cl.score = None
        return
    
    min_sil = min(min(map(np.min, input.pdbflex_sils)), np.min(all_sils))
    max_sil = max(max(map(np.max, input.pdbflex_sils)), np.max(all_sils))
    min_diam = min(min(map(np.min, input.pdbflex_diams)), np.min(all_diams))
    max_diam = max(max(map(np.max, input.pdbflex_diams)), np.max(all_diams))

    input.scaler_sils_pdbflex = list(map(lambda k: (k - min_sil)/(max_sil - min_sil),  input.pdbflex_sils))
    input.scaler_diams_pdbflex = list(map(lambda k: (k - min_diam)/(max_diam - min_diam),  input.pdbflex_diams))
    
    for cl in cls:
        if cl.diam:
            cl.diam_scaler = (cl.diam - min_diam)/(max_diam - min_diam)
            cl.sil_scaler = (cl.sil - min_sil)/(max_sil - min_sil)
        else:
            cl.diam_scaler = None
            cl.sil_scaler = None
            
    mean_sil_diam = np.mean([np.median(input.scaler_sils_pdbflex[i]*input.scaler_diams_pdbflex[i]) for i in range(len(input.scaler_sils_pdbflex))])
    std_sil_diam = np.std([np.median(input.scaler_sils_pdbflex[i]*input.scaler_diams_pdbflex[i]) for i in range(len(input.scaler_sils_pdbflex))])

    for cl in cls:
        if cl.diam:
            cl.spec = cl.sil_scaler*cl.diam_scaler
            cl.score = (cl.sil_scaler*cl.diam_scaler - mean_sil_diam)/std_sil_diam
        else:
            cl.score = None
            cl.spec = None

if __name__ == "__main__":
    
    v = 1.0
        
    input = Input(sys.argv)

    alignment = PDBS_alignment(input.aligned_pdbs, input.aligned_fasta)
    
    if input.ref is None:
        input.ref_num = 0
        input.ref = alignment.molecules[0].fasta_name
    else:
        fasta_names = [alignment.molecules[i].fasta_name for i in range(len(alignment.molecules))]
        assert input.ref in fasta_names, 'No reference protein in fasta file'
        input.ref_num = fasta_names.index(input.ref)
    
    heaps = alignment.get_conserv_positions(perc_of_gaps=input.perc_of_gaps, max_content_of_mismatch=input.max_content_of_mismatch, mismatch_threshold=input.mismatch_threshold)
        
    cls = []
    for heap in tqdm(heaps):
        cl = One_column_clustering(heap, alignment, input.method_of_clustering)
        cl.clustering(input)
        cls.append(cl)

    assert len([cl.is_ok for cl in cls if cl.is_ok]) > 0, 'No specific positions found' 

    statistic(cls, input)

    f = open(os.path.join(input.output, 'RESULT.py'), 'w')
    f.write('cmd.bg_color("white")\n')
    for i in range(len(alignment.molecules)):
        if i == input.ref_num:
            pref = 'Ref_'
        else:
            pref = ''

        f.write('cmd.load(r"{}", "{}")\n'.format(alignment.molecules[i].filename, pref + alignment.molecules[i].base_filename))
        f.write('cmd.color("silver", "{}")\n'.format(pref + alignment.molecules[i].base_filename))        
        
    cls_sorted = sorted(cls, key=lambda x: x.score if x.score is not None else -9999, reverse=True)
    for cl in cls_sorted[0:min(input.number_of_result_resids, len(cls_sorted))]:
        if cl.is_ok:
            cl.printing_py_result(input) 

    f.write('cmd.show("cartoon")\n')
    f.write('cmd.hide("lines")\n')
    f.write('cmd.hide("nonbonded")\n')
    f.write('cmd.set("stick_radius", "0.3")\n')
    f.write('cmd.color("chartreuse", "hetatm and element C")\n')
    f.write('cmd.show("sticks", "hetatm")\n')
    f.write('cmd.color("orange", "hetatm and not element C+H+N+O+S+P")\n')
    f.write('cmd.show("spheres", "polymer and hetatm and not element C+H+N+O+S+P")\n')
    f.write('cmd.set("sphere_scale", "0.5", "hetatm")\n')
    f.write('cmd.set("cartoon_gap_cutoff", "0")\n')
    f.write('cmd.save(r"'+os.path.join(input.output, 'RESULT.pse")'))
    f.close()
    if input.create_pse:
        os.system('pymol -qc '+os.path.join(input.output, 'RESULT.py'))
    table = [["Rank", "Specificity", "Z-Score", "Number of clusters", "Residue number in reference protein PDB {}".format(input.ref)]]
    i = 1
    for _, cl in enumerate(cls_sorted[0:min(input.number_of_result_resids, len(cls_sorted))]):
        if cl.is_ok:
            if cl.score is not None:
                table.append([i, round(cl.spec, 3), round(cl.score, 3), cl.number_of_clusters, cl.ref])
                i += 1
            else:
                table.append([i, 'N/A', 'N/A',  cl.number_of_clusters, cl.ref])
                i += 1
    with open(os.path.join(input.output, 'RESULT.txt'), 'w') as f:
        f.write(tabulate(table)+"\n")
    
    with open(os.path.join(input.output, 'RESULT_extended_info.txt'), 'w') as f:
        for i, cl in enumerate(cls_sorted[0:min(input.number_of_result_resids, len(cls_sorted))]):
            if cl.is_ok:
                f.write('cluster #{}:\n'.format(i + 1))
                for label in set(cl.lab):
                    for j in range(len(cl.lab)):
                        if cl.lab[j] == label:
                            f.write('label : {}, name of prot: {}\n'.format(label, cl.mols_and_aa[j][0].base_filename))
                f.write('\n')


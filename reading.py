import numpy as np

def glue(seq):
    glue_seq = []
    for i, val in enumerate(seq):
        if i == 0: glue_seq.append(val)
        else:
            if val == seq[i-1]: pass
            else: glue_seq.append(val)
    return glue_seq
    
def name_of_at(line):
    return(line[12:16].replace(' ',''))

def num_aa(line):
    return(int(line[22:26]))  

def name_of_ch(line):
    return line[21:22] 
    
def name_of_aa(line):
    return line[16:20] .replace(' ','')       

def num_of_at(line):
    return line[6:11].replace(' ','') 

def readprot(filename):
    dict1 = {}
    dict2 = {}
    seq = []
    chain = 'No chain'
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if chain != name_of_ch(line) and chain != 'No chain':
                    print('Warning: One pdb-file should contain one chain. File '+filename+' contains more, than one chain')
                    raise Exception
                else: 
                    chain =  name_of_ch(line)
                aminoac = name_of_aa(line)
                dict1[num_aa(line)] = name_of_aa(line)
                seq.append(num_aa(line))
                try:
                    dict2[num_aa(line)][name_of_at(line)] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                except Exception:
                    dict2[num_aa(line)] = {}
                    dict2[num_aa(line)][name_of_at(line)] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            else:
                pass
    
    return dict1, dict2, glue(seq)

def seq_rename(seq):
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
    "hP" : 'HYP',
    "W" : 'TRP',
    "I" : 'ILE',
    "M" : 'MET',
    "T" : 'THR',
    "hK" : 'HYL', 
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
    "HYL" : 'hK', 
    }
    a = []
    for i in seq:
        try: a.append(aminoac[i])
        except Exception: a.append('X')
    return a




  

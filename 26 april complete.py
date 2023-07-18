#!/usr/bin/env python
# coding: utf-8

# In[2]:


#generate ORF
import pickle
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

Entrez.email='francisrobitaille07@gmail.com'

cytstop = ['TAA','TAG','TGA']
cytstart = ['ATG']
mitstop = ['TAA','TAG', 'AGA', 'AGG'] # ribosome stalling patterns? # UGA for humanin?
#mitnouvstop = ['AGA','AGG']
mitstart = ['ATG','ATT','ATA'] #is there molecular evidence, or could it point to some other initiation mechanism?

#replace Z with real AA value
#should I combine mitstop and mitnouvstop?? I probably will not find

#what does AUA code , I or M. Same with AUU, it seems like AUU needs to code for I to get the same refPROT (ND1)

#I'll have to have 2 codes, one for translation and one to create proteins, B can't pop up

#MAJOR PROBLEM: now mit_alt is used which assigns R for AGA and AGG. Is that ok (as per Kienzle et al)
# Or should I have a mit with AGA and AGG with stop and a cytoplasmic code (export). Or a third? Resolve this urgent matter

mit = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "B", "AUC": "I", "AUA": "M", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "W", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

#AUU is directly I instead of provisional B (this code is for double checking presence of refORF)
mit_for_translation_of_ref = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "M", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "W", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

mit_alt = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "B", "AUC": "I", "AUA": "M", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "W", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

cyt = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

class Sequence:
    def __init__(self, name, dna_seq, start, end):
    # Start and end are defined for each strand. L coord is not equal H coord. See H to L functions
        self.name = name
        self.start = start
        self.end = end
        self.dna_seq = str(dna_seq)
        self.rna_seq = self.dna_seq.replace("T", "U")
        
class ORF:
    def __init__(self, sequence, mother, start, end, reading_frame):
        
        #a relative indice to differentiate the reading frame of a protein (is prot X on same RF as Y)
        self.reading_frame = reading_frame
        self.sequence = sequence
        
        #mother is the strand (tRNA punctuated or not) from where the protein comes from, ex: H_1 or H_full (no tRNA punctuation)
        self.mother = mother
        
        #The longest protein in a sequence between * codons gets assigned a list of smaller sequences
        #it comprises. Ex: MAAMF*: MAAMF* is the main orf and MF* is in the list of sister proteins. They are
        #of course all in the same reading frame
        
        self.list_sisters = []
    
        
        #start and end relative to the frame of reference of mother
        self.start = start
        self.end = end
        self.dna = None
        self.rna = None
        #self.hasstop = None

#defining functions

def seperate_codons(s):
    char = ','
    n = 3
    return char.join(s[i:i+n] for i in range(0, len(s), n))
    
    #ex: string = ND6.rna_seq

    #new_string = seperate_codons(string)
    #print(new_string)


def transcribe(dna_seq):
    #T,U normally
    rna_seq = dna_seq.replace("T", "U")
    return rna_seq

def transcribe_2(dna_seq):
    #T,U normally
    rna_seq = dna_seq.replace("U", "T")
    return rna_seq
        
def translate_all_3frames (mother, code):
    #seq is RNA
    seq = mother.rna_seq
    translated_3_frames = []
    for frame in range (3):
        reading_frame = ''
        for i in range(frame, len(seq),3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                if codon in code:
                    reading_frame += code[codon]
        translated_3_frames.append(reading_frame)
    return translated_3_frames



#split the AA sequence if '*' and reconstruct the stop codon if it was there in the first place
def split (orfs):
    split = orfs.split('*')
    split_wise = [seq + '*' for seq in split[:-1]]
    if orfs[-1] == '*':
        split_wise.append(split[-1] + '*')
    else:
        split_wise.append(split[-1])
    return split_wise

def create_proteins (orfs, mother):
    
    #minimum len in AA
    min_len = 10
    final_prots = []
    
    #for each reading frame (3)
    
    compteur_RF = 0
    for i in range(3):
        
        
        # j is each split sequence of AA
        RF_pos = 0
        
        #split separates MYK*KLKK into MYK and KLKK
        for j in split(orfs[i]):
            longest_prot = 0
            longest_protein = None
            
            # k is each amino acid in j
            k = 0
            while k < len(j):
                
                protein = ''
                
                #if it starts with AUA or AUG, just replace downstream AUU's with 'I'
                
                if j[k] == 'M' or j[k] == 'B':
                    
                    protein = j[k:len(j)]
                    
                    #if it starts with AUU, replace it with 'M' and replace downstream AUU's with 'I':
                    if protein[0] == 'B':
                        protein = 'M' + protein[1:]
                        
                    if len(protein) > min_len:
                        protein = protein.replace('B','I')
                        
                        #makes sur the protein is the longest. If not, assigns itself as sister
                        #protein to the longest protein
                        
                        if longest_prot == 0:
                            longest_protein = ORF(protein, mother.name, mother.start + compteur_RF + (RF_pos+ k )*3, mother.start + compteur_RF + (RF_pos+ k +len(protein))*3, compteur_RF)
                            final_prots.append(longest_protein)
                            
                        if longest_prot != 0:
                            longest_protein.list_sisters.append(protein)
                        longest_prot += 1
                    
                k += 1
            RF_pos += len(j)
        compteur_RF += 1
        
    return final_prots

def find_reforf(protein, ref_H):
    # ref_orf is a Sequence object
    # prot is a ORF object
    for prot in protein:
        
        is_ref = 0
    
    #sometimes the actual refPROT is shorter because it's actually in the sister list of our ORF, 
    #hence the 'in' and not '=='
        for ref_orf in ref_H:
            ref_prot = translate2(ref_orf.rna_seq, mitstart, mitstop, mit_for_translation_of_ref)
            
            if 'L' in prot.mother:
                # [1:] to account for 2 translation tables (first codon sometimes fucks up (I instead of M in ref_prot))
                if ref_prot[1:] in prot.sequence:
                    
                    
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print( 'ref protein: ', ref_orf.name, ', location: ', prot.mother,'Reading_frame: ', prot.reading_frame)
                    print(' start:         ',L_to_H(prot.start),', end:         ',  L_to_H(prot.end)) 
                    print(' genbank start: ', ref_orf.start,    ', genbank end: ', ref_orf.end) 
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
            
            else:
                if ref_prot[1:] in prot.sequence:
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print('ref protein: ', ref_orf.name, ', location: ', prot.mother, 'Reading_frame: ', prot.reading_frame)
                    print(' start:         ',prot.start,', end:         ',  prot.end) 
                    print(' genbank start: ', ref_orf.start,    ', genbank end: ', ref_orf.end)  
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
                    
        if is_ref == 0:
            print(prot.mother, ', start: ', prot.start, ', end: ', prot.end, ', Reading frame: ', prot.reading_frame)
            print('sequence: ', prot.sequence)
            if len(prot.list_sisters) == 0:
                print('No sister proteins')
            else:
                print('sisters: ', prot.list_sisters)
            print('')
            
def codon_count_3frames (mother, ref_prot_RF):
    
    total_count = {}
    codons = {'AGA' : 0, 'AGG' : 0, 'UAG' : 0, 'UAA' : 0, 'AUG' : 0, 'AUU' : 0, 'AUA' : 0, 'UGA': 0}
    #seq is RNA
    seq = mother.rna_seq
    translated_3_frames = []
    for frame in range (3):
        reading_frame = ''
        for i in range(frame, len(seq),3):
            codon = seq[i:i+3]
            
            if codon == 'AGA':
                codons[codon]+=1
            if codon == 'AGG':
                codons[codon]+=1
            if codon == 'UAA':
                codons[codon]+=1
            if codon == 'UAG':
                codons[codon]+=1
            if codon == 'AUG':
                codons[codon]+=1
            if codon == 'AUA':
                codons[codon]+=1
            if codon == 'AUU':
                codons[codon]+=1
            
        if ref_prot_RF == 0:
            if frame == 0:
                total_count['RF'] = codons
            if frame == 1:
                total_count['+1'] = codons
            if frame == 2:
                total_count['+2'] = codons
            
        if ref_prot_RF == 1:
            if frame == 0:
                total_count['+2'] = codons
            if frame == 1:
                total_count['+RF'] = codons
            if frame == 2:
                total_count['+1'] = codons
                
        if ref_prot_RF == 2:
            if frame == 0:
                total_count['+1'] = codons
            if frame == 1:
                total_count['+2'] = codons
            if frame == 2:
                total_count['RF'] = codons
                
        codons = {'AGA' : 0, 'AGG' : 0, 'UAG' : 0, 'UAA' : 0, 'AUG' : 0, 'AUU' : 0, 'AUA' : 0}
    return total_count

def other_count_3frames (mother, ref_prot_RF):
    
    total_count = {}
    codons = {'GUG' : 0, 'CUG':0, 'UGA': 0}
    #seq is RNA
    seq = mother.rna_seq
    translated_3_frames = []
    for frame in range (3):
        reading_frame = ''
        for i in range(frame, len(seq),3):
            codon = seq[i:i+3]
            
            if codon == 'GUG':
                codons[codon]+=1
            if codon == 'CUG':
                codons[codon]+=1
            if codon == 'UGA':
                codons[codon]+=1
            
                    
        #concevoir tous les scénarios RF = 0, 1, 2
                                    #Frame = 0, 1, 2
        if ref_prot_RF == 0:
            if frame == 0:
                total_count['RF'] = codons
            if frame == 1:
                total_count['+1'] = codons
            if frame == 2:
                total_count['+2'] = codons
            
        if ref_prot_RF == 1:
            if frame == 0:
                total_count['+2'] = codons
            if frame == 1:
                total_count['+RF'] = codons
            if frame == 2:
                total_count['+1'] = codons
                
        if ref_prot_RF == 2:
            if frame == 0:
                total_count['+1'] = codons
            if frame == 1:
                total_count['+2'] = codons
            if frame == 2:
                total_count['RF'] = codons
                
        codons = {'GUG' : 0, 'CUG':0, 'UGA': 0}
        
    return total_count

#both of these functions are identical bc reciprocal (H to L and L to H)
def L_to_H(pos_L):
    
    if 0<= pos_L <=410:
        pos_H = 410 - pos_L
    elif 410 <= pos_L <= 16569: #13749
        pos_H = 16979 - pos_L
    else:
        pos_H ='out a luck'
    return pos_H

def H_to_L(pos_H):
    
    if 0<= pos_H <=410:
        pos_L = 410 - pos_H
    elif 410 <= pos_H <= 16569:
        pos_L = 16979 - pos_H
    else:
        pos_L ='out a luck'
    return pos_L

#need another translate function for mitochondria (maybe have the codon table as input)

def translate2(sequence, start_codons, stop_codons, codon_table):
    protein_sequence = ''
    i = 0
    while i <= len(sequence):
        codon = sequence[i:i+3]
        #print (len(sequence))
        #print (codon)
        if codon in start_codons and len(protein_sequence) == 0:
            protein_sequence += 'M'
        elif codon in codon_table:
            protein_sequence += codon_table[codon]
        elif codon in stop_codons:
            protein_sequence += '*'
            break
        i += 3
    return protein_sequence


#returns the complementary strand in the correct order for translation (5' -> 3')
def reverse_complement(sequence):
    # Define a dictionary for complementing nucleotides
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
    # Reverse the sequence
    reversed_sequence = sequence[::-1]
    # Complement each nucleotide and join them to form the reverse complement
    reverse_complement = ''.join([complement[nucleotide] for nucleotide in reversed_sequence])
    return reverse_complement

#http://scikit-bio.org/docs/0.4.2/generated/skbio.alignment.global_pairwise_align_protein.html#skbio.alignment.global_pairwise_align_protein

#NW_score = global_pairwise_align_protein(seq1, seq2, gap_open_penalty=10, gap_extend_penalty=0.5, substitution_matrix=BLOSUM62, penalize_terminal_gaps=False)
#pseudo_parent_prot_pairs[pseudo_alt].update({'sortest_path_len':sp, 'NW_score':NW_score[1]})
#not used as of now

def protein_aligner (query, list_of_prots):
    
    # Initialize a PairwiseAligner object
    aligner = PairwiseAligner()

    # Loop through the list of sequences to find the best match
    best_score = 0
    best_seq = None
    
    for seq in list_of_prots:
    # Use the aligner object to get the optimal global alignment between the two sequences
        alignments = aligner.align(query, seq)
    # Get the alignment score of the first alignment
        alignment_score = alignments[0].score
    # If the current alignment score is better than the previous best score, update the best score and best sequence
        if alignment_score > best_score:
            best_score = alignment_score
            best_seq = seq
            
    print('best match:   ',best_seq)
    print('query peptide:',query)
    print('match score:', best_score)
    
    return best_score


# In[ ]:





# In[3]:


# generate the H strand and L strand and also enforce tRNA punctuation

handle = Entrez.efetch(db='nucleotide',id='NC_012920.1', rettype='fasta')
record = SeqIO.read(handle, 'fasta')

#refer to : https://docs.google.com/document/d/1RO4yXNNFRrOhQN3_erMYJJWUxRD6nK-XdxhKZbS8blU/edit
#for all litterature backed choices

#generate H-strand, from the upstream HSP (around 561) -> entire length -> 535 (last hance at a stop codon)).
H_strand_1 = record.seq[561:16569] 
H_strand_2 = record.seq[0:535]

H_full =  Sequence('H_full', Seq(str(H_strand_1) + str(H_strand_2)), 561, 535)

#generate one continuous L-strand from the LSP (410) -> entire strand -> 409.
#First generate two strands as such: 410 -> 0 and 16569 -> 409, then join them together

#flips the sequence and corrects the Us to Ts
L_strand_1 = transcribe_2(reverse_complement(record.seq[0:410]))
L_strand_2 = transcribe_2(reverse_complement(record.seq[409:16569]))

#same as above
L_full = Sequence('L_full', Seq(str(L_strand_1) + str(L_strand_2)), H_to_L(410), H_to_L(409))
   
#Heavy Strand tRNA punctuated: note that individual tRNAs are not considered as potential mRNAs post cleavage

# upstream HSP to Phe (577 to 647). Too short (560 to 577): will scrap H1
#H_1 = Sequence('H_1', record.seq[560:576] , 560, 576)

# Phe (577 to 647) to Val (1602 to 1670)
H_2 = Sequence('H_2', record.seq[648:1601], 648, 1601)

# Val (1602 to 1670) to Leu 1 (3230 to 3304)
H_3 = Sequence('H_3', record.seq[1671:3229] , 1617, 3229)

# Leu 1 (3230 to 3304) to Isoleu (4263 to 4331)
H_4 = Sequence('H_4', record.seq[3305:4262] , 3305, 4262)

# Isoleu (4263 to 4331) to Met(4402 to 4469) - complement is tRNA Glu
H_5 = Sequence('H_5', record.seq[4332:4401] , 4332, 4401)

# Met (4402 to 4469) to trypto (5512 to 5579)
H_6 = Sequence('H_6', record.seq[4469:5511] , 4470,5511)

# Trypto (5512 to 5579) to AspAcid (7518 to 7585)
H_7 = Sequence('H_7', record.seq[5580:7517] , 5580,7517)

# AspAcid (7518 to 7585) to Lys (8295 to 8364)
H_8 = Sequence('H_8', record.seq[7585:8294] , 7586, 8294)

# Lys (8295 to 8364) to Gly (9991 to 10058)
H_9 = Sequence('H_9', record.seq[8365:9990] , 8365, 9990)

# Gly (9991 to 10058) to Arg (10405 to 10469)
H_10 = Sequence('H_10', record.seq[10058:10404] , 10059, 10404)

# Arg (10405 to 10469) to his (12138 to 12206) until Leucine 2 (.. to 12236)
H_11 = Sequence('H_11', record.seq[10469:12137] , 10470, 12137)

# Leucine 2 (.. to 12236) to Threonine (15888 to 15954)
H_12 = Sequence('H_12', record.seq[12237:15887] , 12237, 15887)

# Threonine (15888 to 15954) NOT to end of H-strand transcription (TAS @ 16166)
# Because there is still significant transcription beyond TAS
H_13_1 = record.seq[15955:16569]
H_13_2 = record.seq[0:559]
H_13 = Sequence('H_13', H_13_1 + H_13_2 , 15955, 559)

H_tRNA_enforced = [H_2,H_3,H_4,H_5,H_6,H_7,H_8,H_9,H_10,H_11,H_12,H_13]

#Light Strand:

# in the reverse order
# prendre juste ce qui est pas dans le tRNA 

L_1_1 = transcribe_2(reverse_complement(record.seq[0:407]))
L_1_2 = transcribe_2(reverse_complement(record.seq[16024:16569]))

# LSP (400) to Proline (16023 - 15956) 
L_1 = Seq(str(L_1_1) + str(L_1_2))
L_1 = Sequence('L_1', L_1 , H_to_L(407), H_to_L(16024))

# Pro (16023 - 15956) to Glutamic acid (14742 - 14674)
L_2 = Sequence('L_2', Seq(transcribe_2(reverse_complement(record.seq[14743:15955]))), H_to_L(15955), H_to_L(14743))

# Glutamic acid (14742 - 14674) to Serine 1 (7514 - 7446)
L_3 = Sequence('L_3', Seq(transcribe_2(reverse_complement(record.seq[7515:14673]))), H_to_L(14673), H_to_L(7515))

# Serine 1 (7514 - 7446) to ANCY (5903 - 5580)
L_4 = Sequence('L_4', Seq(transcribe_2(reverse_complement(record.seq[5904:7445]))), H_to_L(7445), H_to_L(5904))

# ANCY (5903 - 5580) to Glutamine (4402 - 4331)
#check_reforf(L_5,mitstart,mitstop,mit)
L_5 = Sequence('L_5', Seq(transcribe_2(reverse_complement(record.seq[4403:5578]))), H_to_L(5578), H_to_L(4403))


# Glutamine (4402 - 4331) to LSP (see doc for choices)
L_6 = Sequence('L_6', Seq(transcribe_2(reverse_complement(record.seq[408:4330]))), H_to_L(4330), H_to_L(408))

L_tRNA_enforced = [L_1,L_2,L_3,L_4,L_5,L_6]


# In[4]:


# ALL KNOWN REFORFS AND ALTORFS

#generate sequences and check that my mitochondrial code 'mit' actually reconstitutes all refprots
def check_reforf(seq, start, stop, code):
    #print(seq)
    prot_seq = translate(transcribe(seq), start, stop, code)
    #print(prot_seq)
    return prot_seq

# ALL 3 are broken up by AGA's
HUMANIN = 'MAPRGFSCLLLLTSEIDLPVKRRA'
#Sequence('HUMANIN', record.seq[2632:2707] , 2632, 2707)

ALTND4 = 'MRHNYNKLHLPTTNRPKIAHCMLFNQPHSPRSNSHSHPNPLKLHRRSHSHNRPRAYILITILPSKLKLRTHSQSHHNPLS'
#Sequence('ALTND4', record.seq[11556:11856] , 11556,11856)
#M*HNYNKLHLPTTN*PKIAHCMLFNQPHSPRSNSHSHPNPLKLHRRSHSHNRPRAYILITILPSKLKLRTHSQSHHNPLS*

ALTCO1 = 'MCNNLLHSNTHHNRRLWQLTSSPNNRCPRYGVSPHKQHKLLTLTSLSPTPARICYSGGRSRNRLNSLPSLSRELLPPWSLRRPNHLLLTPSRCLLYLRGHQFHHNNYQYKTPCHNPMPNAPLRLIRPNHSSPTSPISPSPSCWHHYTTNRPQPQHHLLRPRRRRRPHSMPTPILIFRSPWSLYSYPTRLRNNLPYCNLLLRKKRTIWMHRYGLSYDINWLPRVYRVSTPYIYSRNRRRHTSMFHLRYHNHRYPHRRQSI'
#Sequence('ALTCO1', record.seq[6088:6868] , 6088,6868)
#ICNNLLHSNTHHNR*LWQLTSSPNNRCPRYGVSPHKQHKLLTLTSLSPTPARICYSGGRS*N*LNSLPSLS*ELLPPWSLR*PNHLLLTPS*CLLYL*GHQFHHNNYQYKTPCHNPMPNAPLRLIRPNHSSPTSPISPSPSCWHHYTTN*PQPQHHLLRPRR***PHSMPTPILIFRSPWSLYSYPT*LRNNLPYCNLLLRKK*TIWMH*YGLSYDINWLP*VYRVSTPYIYS*N*R*HTSMFHLRYHNHRYPHRRQSI*

SHLP1 = 'MCHWAGGASNTGDARGDVFGKQAG'
#Sequence('SHLP1', Seq(transcribe_2(reverse_complement(record.seq[2484:2559]))), H_to_L(2559), H_to_L(2484))
# MCHWAGGASNTGDA*GDVFGKQA

SHLP2 = 'MGVKFFTLSTRFFPSVQRAVPLWTNS'
#Sequence('SHLP2', Seq(transcribe_2(reverse_complement(record.seq[2087:2168]))), H_to_L(2168), H_to_L(2087))
# MGVKFFTLST*FFPSVQ*AVPLWTN

SHLP3 = 'MLGYNFSSFPCGTISIAPGFNFYRLYFIWVNGLAKVVW'
#Sequence('SHLP3', Seq(transcribe_2(reverse_complement(record.seq[1702:1819]))), H_to_L(1819), H_to_L(1702))
#'MLGYNFSSFPCGTMSIAPGFNFYRLYFIWVNGLAKVVW*'

SHLP4 = 'MLEVMFLVNRRGKICRVPFTFFNLSL'
#Sequence('SHLP4', Seq(transcribe_2(reverse_complement(record.seq[2441:2522]))), H_to_L(2522), H_to_L(2441))
# MLEVMFLVN*RGKICRVPFTFFNLS

SHLP5 = 'MYCSEVGFCSEVAPTEIFNAGLVV'
#Sequence('SHLP5', Seq(transcribe_2(reverse_complement(record.seq[2779:2854]))), H_to_L(2854), H_to_L(2779))
# 'MYCSEVGFCSEVAPTEIFNAGLVV*'

SHLP6 = 'MLDQDIPMVQPLLKVRLFND'
#Sequence('SHLP6', record.seq[2989:3052], 2992, 3049)
# 'MLDQDIPMVQPLLKVRLFND*'
    
MOTSC = 'MRWQEMGYIFYPRKLR'
#Sequence('SHLP6', record.seq[1342:1393], 1342, 1392)
# M*WQEMGYIFYP*KLR*

SHMOOSE = 'MPPCLTTWLSQLLKDNSYPLVLGPKNFGATPNKSNNHAHYYNHPNPDFPNSPHPYHPR'
#Sequence('SHMOOSE', record.seq[12233:12410] , 12230,12407)

ND1 = 'MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT'
#Sequence('ND1', record.seq[3306:4262] , 3306, 4262)

# 1,2. ('AUU' = 'I')
ND2 = 'MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLNVSLLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL'
#Sequence('ND2', record.seq[4469:5511] , 4469, 5511)
# 1,2 ('first AUU' = 'M'). the translate function needs to account for that. 

CO1 = 'MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQPGNLLGNDHIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSLLLLLASAMVEAGAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKWSAAVLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTTWNILSSVGSFISLTAVMLMIFMIWEAFASKRKVLMVEEPSMNLEWLYGCPPPYHTFEEPVYMKS'
#Sequence('CO1', record.seq[5903:7446] , 5903, 7446)

CO2 = 'MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFLTLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGHQWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVLHSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKIFEMGPVFTL'
#Sequence('CO2', record.seq[7585:8269] , 7585, 8269)

ATP8 = 'MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKMKNYNKPWEPKWTKICSLHSLPPQS'
#Sequence('ATP8', record.seq[8365:8572] , 8365, 8572)

ATP6 = 'MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQWLIKLTSKQMMTMHNTKGRTWSLMLVSLIIFIATTNLLGLLPHSFTPTTQLSMNLAMAIPLWAGTVIMGFRSKIKNALAHFLPQGTPTPLIPMLVIIETISLLIQPMALAVRLTANITAGHLLMHLIGSATLAMSTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT'
#Sequence('ATP6', record.seq[8526:9207] , 8526, 9207)

CO3 = 'MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLMLGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFYHSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMIQALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTICFIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWGS'
#Sequence('CO3', record.seq[9206:9990] , 9206, 9990)

ND3 = 'MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPMSPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSLAYEWLQKGLDWTE'
#Sequence('ND3', record.seq[10058:10404] , 10058, 10404)

ND4L = 'MPLIYMNIMLAFTISLLGMLVYRSHLMSSLLCLEGMMLSLFIMATLMTLNTHSLLANIVPIAMLVFAACEAAVGLALLVSISNTYGLDYVHNLNLLQC'
#Sequence('ND4L', record.seq[10469:10766] , 10469, 10766)

ND4 = 'MLKLIVPTIMLLPLTWLSKKHMIWINTTTHSLIISIIPLLFFNQINNNLFSCSPTFSSDPLTTPLLMLTTWLLPLTIMASQRHLSSEPLSRKKLYLSMLISLQISLIMTFTATELIMFYIFFETTLIPTLAIITRWGNQPERLNAGTYFLFYTLVGSLPLLIALIYTHNTLGSLNILLLTLTAQELSNSWANNLMWLAYTMAFMVKMPLYGLHLWLPKAHVEAPIAGSMVLAAVLLKLGGYGMMRLTLILNPLTKHMAYPFLVLSLWGMIMTSSICLRQTDLKSLIAYSSISHMALVVTAILIQTPWSFTGAVILMIAHGLTSSLLFCLANSNYERTHSRIMILSQGLQTLLPLMAFWWLLASLANLALPPTINLLGELSVLVTTFSWSNITLLLTGLNMLVTALYSLYMFTTTQWGSLTHHINNMKPSFTRENTLMFMHLSPILLLSLNPDIITGFSS'
#Sequence('ND4', record.seq[10759:12137] , 10759, 12137)

ND5 ='MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIMEFSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYARADANTAAIQAILYNRIGDIGFILALAWFILHSNSWDPQQMALLNANPSLTPLLGLLLAAAGKSAQLGLHPWLPSAMEGPTPVSALLHSSTMVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLMMVTIGINQPHLAFLHICTHAFFKAMLFMCSGSIIHNLNNEQDIRKMGGLLKTMPLTSTSLTIGSLALAGMPFLTGFYSKDHIIETANMSYTNAWALSITLIATSLTSAYSTRMILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKMKSPLCTFYFSNMLGFYPSITHRTIPYLGLLTSQNLPLLLLDLTWLEKLLPKTISQHQISTSIITSTQKGMIKLYFLSFFFPLILTLLLIT'
#Sequence('ND5', record.seq[12336:14148] , 12336, 14148)

#ND6 = record.seq[14147:14673]
#true_ND6 = transcribe_2(reverse_complement(ND6))
ND6 = 'MMYALFLLSVGLVMGFVGFSSKPSPIYGGLVLIVSGVVGCVIILNFGGGYMGLMVFLIYLGGMMVVFGYTTAMAIEEYPEAWGSGVEVLVSVLVGLAMEVGLVLWVKEYDGVVVVVNFNSVGSWMIYEGEGSGLIREDPIGAGALYDYGRWLVVVTGWTLFVGVYIVIEIARGN'
#Sequence('ND6', true_ND6 , 14673, 14147)
# T -> reverse_c -> U -> transccribe_2 -> T 

CYTB = 'MTPMRKTNPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGLYYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDLVQWIWGGYSVDSPTLTRFFTFHFILPFIIAALATLHLLFLHETGSNNPLGITSHSDKITFHPYYTIKDALGLLLFLLSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLFAYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLILTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA'
#Sequence('CYTB', record.seq[14746:15887] , 14746, 15887)


# In[9]:


# Heavy tRNA punctuated
# H_tRNA_enforced is a list of sequences (H_1,H_2, etc.)
# translate_all_3frames takes a sequence, translates each frame and returns a long list of AA with * denoting stop codons
# create proteins creates ORF objects by first splitting protein sequences on either side of *. The longest protein is the
# ORF but it also contains a list of 'sister' proteins that are always shorter. MS or ribosome profiling data might match 
# the sisters better so I will need to do a local BLAST and not a global one (for MS). And also perhaps have a consequent
# graphical method for ribosome profiling data (let's see what kind of data I find)

#proteins that do not have a stop codon but that are cleaved because of tRNA will show up and that's ok
#as a measure: of 224 prots in H_tRNA_enforced: 201 have a stop and 23 dont
prots_H_tRNA_mit = []
for H in H_tRNA_enforced:
    frames = translate_all_3frames(H, mit_alt)
    prots_H_tRNA_mit.append(create_proteins(frames, H))
    
# Heavy full strand, no punctuation, 248 proteins
prots_H_full_mit = create_proteins(translate_all_3frames(H_full,mit_alt), H_full)

# Light tRNA punctuated.  285 proteins, 10 dont have a stop codon, 275 do
prots_L_tRNA_mit = []
for L in L_tRNA_enforced:
    frames = translate_all_3frames(L, mit_alt)
    prots_L_tRNA_mit.append(create_proteins(frames, L))

# Light full strand, no punctuation. 361 proteins
prots_L_full_mit = create_proteins(translate_all_3frames(L_full,mit_alt), L_full)


#609 proteins in total. Of course some proteins in tRNA punctuation are not there but most are. Maybe do a final repo
# that has ALL the proteins with a tag to know where it comes from or something like that

Total_prots_mit = prots_L_full_mit + prots_H_full_mit

#--------------------------------------------------------------------------------------------------------------------

#Same but with a cytoplasmic code

# Heavy tRNA punctuated
prots_H_tRNA_cyt = []
for H in H_tRNA_enforced:
    frames = translate_all_3frames(H, cyt)
    prots_H_tRNA_cyt.append(create_proteins(frames, H))

# Heavy full
prots_H_full_cyt = create_proteins(translate_all_3frames(H_full,cyt), H_full)

# Light tRNA punctuated
prots_L_tRNA_cyt = []
for L in L_tRNA_enforced:
    frames = translate_all_3frames(L, cyt)
    prots_L_tRNA_cyt.append(create_proteins(frames, L))
    
prots_L_full_cyt = create_proteins(translate_all_3frames(L_full,cyt), L_full)

# ---------------With AGA and AGG as stop codon and regular mitochondrial code

prots_H_tRNA_mit2 = []
for H in H_tRNA_enforced:
    frames = translate_all_3frames(H, mit)
    prots_H_tRNA_mit2.append(create_proteins(frames, H))
    
# Heavy full strand, no punctuation, 248 proteins
prots_H_full_mit2 = create_proteins(translate_all_3frames(H_full,mit), H_full)

# Light tRNA punctuated.  285 proteins, 10 dont have a stop codon, 275 do
prots_L_tRNA_mit2 = []
for L in L_tRNA_enforced:
    frames = translate_all_3frames(L, mit)
    prots_L_tRNA_mit2.append(create_proteins(frames, L))

# Light full strand, no punctuation. 361 proteins
prots_L_full_mit2 = create_proteins(translate_all_3frames(L_full,mit), L_full)



# # Test
# 
# Ceci est du texte entre le code : 

# In[6]:


#Retrouver les séquences de référence et les altorf connus

ref_H = (HUMANIN, ALTCO1, ALTND4, SHMOOSE, MOTSC, SHLP6,ND1,ND2,ND3,ND4,ND4L,ND5, CO1,CO2,CO3,ATP8,ATP6,CYTB)
ref_L = (ND6, SHLP3, SHLP5, SHLP1, SHLP2, SHLP4)

#modified so that ref_H is a list of protein sequences and not a sequence object
def find_reforf(protein, ref_H):
    # ref_orf is a Sequence object
    # prot is a ORF object
    for prot in protein:
        
        is_ref = 0
    
    #sometimes the actual refPROT is shorter because it's actually in the sister list of our ORF, 
    #hence the 'in' and not '=='
        for ref_prot in ref_H:
            
            if 'L' in prot.mother:
                # [1:] to account for 2 translation tables (first codon sometimes fucks up (I instead of M in ref_prot))
                if ref_prot[1:] in prot.sequence:
                    
                    
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
            #if in H
            else:
                if ref_prot[1:] in prot.sequence:
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
                    
        if is_ref == 0:
            print(prot.mother, ', start: ', prot.start, ', end: ', prot.end, ', Reading frame: ', prot.reading_frame)
            print('sequence: ', prot.sequence)
            if len(prot.list_sisters) == 0:
                print('No sister proteins')
            else:
                print('sisters: ', prot.list_sisters)
            print('')
            


print('____________________________________________ H strand tRNA punctuated ____________________________________________')

for proteins in prots_H_tRNA_mit:
    find_reforf(proteins, ref_H)

print('____________________________________________ L strand tRNA punctuated ____________________________________________')

#because prots_L_tRNA is a list of individual mothers (L1,L2 etc.)
for proteins in prots_L_tRNA_mit:
    find_reforf(proteins, ref_L)
            
print('____________________________________________ H strand full ____________________________________________')

find_reforf(prots_H_full_mit, ref_H)

print('____________________________________________ L strand full ____________________________________________')

find_reforf(prots_L_full_mit, ref_L)

print(' MIT 2 **************************** MIT 2 **************************** MIT 2 **************************** MIT 2')


print('____________________________________________ H strand tRNA punctuated MIT2____________________________________________')

for proteins in prots_H_tRNA_mit:
    find_reforf(proteins, ref_H)

print('____________________________________________ L strand tRNA punctuated MIT2____________________________________________')

#because prots_L_tRNA is a list of individual mothers (L1,L2 etc.)
for proteins in prots_L_tRNA_mit:
    find_reforf(proteins, ref_L)
            
print('____________________________________________ H strand full MIT2____________________________________________')

find_reforf(prots_H_full_mit, ref_H)

print('____________________________________________ L strand full MIT2____________________________________________')

find_reforf(prots_L_full_mit, ref_L)
    


# In[7]:


# Same, find refprots but with cytoplasmic code: 

def find_reforf(protein, ref_H):
    # ref_orf is a Sequence object
    # prot is a ORF object
    
    for prot in protein:
        
        is_ref = 0
    
    #sometimes the actual refPROT is shorter because it's actually in the sister list of our ORF, 
    #hence the 'in' and not '=='
        for ref_prot in ref_H:
            
            if 'L' in prot.mother:
                # [1:] to account for 2 translation tables (first codon sometimes fucks up (I instead of M in ref_prot))
                if ref_prot[1:] in prot.sequence:
                    
                    
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
            #if in H
            else:
                if ref_prot[1:] in prot.sequence:
                    print('***************************************** MATCH: REF PROT ********************************************')
                    print('generated sequence: ', prot.sequence)
                    print('genbank sequence:  ', ref_prot)
                    print('******************************************************************************************************')
                    is_ref += 1
                    
        if is_ref == 0:
            print(prot.mother, ', start: ', prot.start, ', end: ', prot.end, ', Reading frame: ', prot.reading_frame)
            print('sequence: ', prot.sequence)
            if len(prot.list_sisters) == 0:
                print('No sister proteins')
            else:
                print('sisters: ', prot.list_sisters)
            print('')
            


print('______________________________________ H strand tRNA punctuated CYT ______________________________________')

for proteins in prots_H_tRNA_cyt:
    find_reforf(proteins, ref_H)

print('______________________________________ L strand tRNA punctuated CYT_______________________________________')

#because prots_L_tRNA is a list of individual mothers (L1,L2 etc.)
for proteins in prots_L_tRNA_cyt:
    find_reforf(proteins, ref_L)
            
print('____________________________________________ H strand full CYT____________________________________________')

find_reforf(prots_H_full_cyt, ref_H)

print('____________________________________________ L strand full CYT____________________________________________')

find_reforf(prots_L_full_cyt, ref_L)
    


# In[15]:


#turn the whole repo thing into a fasta file for GUI
# Size of files for min_len of 10 AA
#there will be redundance of course but I might make fasta file for all of the tree branchings

#H_full_mit = 140 proteins
fasta_H_full_mit = []
for protein in prots_H_full_mit:
    fasta_string = f">{'mother: '+ protein.mother + ', code: mit, start: '+ str(protein.start) + ', end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
    fasta_H_full_mit.append(fasta_string)

#H_tRNA_mit = 132 proteins
fasta_H_tRNA_mit = []
for prot in prots_H_tRNA_mit:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: mit, start: '+ str(protein.start) + '  end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
        fasta_H_tRNA_mit.append(fasta_string)
        
#H_full_cyt = 80 proteins
fasta_H_full_cyt = []
for protein in prots_H_full_cyt:
    fasta_string = f">{'mother: '+ protein.mother + ', code: cyt, start: '+ str(protein.start) + '  end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
    fasta_H_full_cyt.append(fasta_string)

#H_tRNA_cyt = 76 proteins
fasta_H_tRNA_cyt = []
for prot in prots_H_tRNA_cyt:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: cyt, start: '+ str(protein.start) + '  end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
        fasta_H_tRNA_cyt.append(fasta_string)

#L_full_mit = 233 proteins
fasta_L_full_mit = []
for protein in prots_L_full_mit:
    fasta_string = f">{'mother: '+ protein.mother + ', code: mit, start: '+ str(L_to_H(protein.start)) + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
    fasta_L_full_mit.append(fasta_string)

#L_tRNA_mit = 222 proteins
fasta_L_tRNA_mit = []
for prot in prots_L_tRNA_mit:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: mit, start: '+ str(L_to_H(protein.start))  + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
        fasta_L_tRNA_mit.append(fasta_string)

#L_full_cyt = 133 proteins
fasta_L_full_cyt = []
for protein in prots_L_full_cyt:
    fasta_string = f">{'mother: '+ protein.mother + ', code: cyt, start: '+ str(L_to_H(protein.start))  + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
    fasta_L_full_cyt.append(fasta_string)

#L_tRNA_cyt = 129 proteins
fasta_L_tRNA_cyt = []
for prot in prots_L_tRNA_cyt:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: cyt, start: '+ str(L_to_H(protein.start))  + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
        fasta_L_tRNA_cyt.append(fasta_string)

# MIT 2____________________________________MIT 2_______________________________________MIT 2__________________________________

# 742 proteins in total

fasta_H_full_mit2 = []
for protein in prots_H_full_mit2:
    fasta_string = f">{'mother: '+ protein.mother + ', code: mit2, start: '+ str(protein.start) + ', end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
    fasta_H_full_mit2.append(fasta_string)

#H_tRNA_mit = 132 proteins
fasta_H_tRNA_mit2 = []
for prot in prots_H_tRNA_mit2:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: mit2, start: '+ str(protein.start) + '  end: ' + str(protein.end)+ ' '}\n{protein.sequence}\n"
        fasta_H_tRNA_mit2.append(fasta_string)
        
#L_full_mit = 233 proteins
fasta_L_full_mit2 = []
for protein in prots_L_full_mit2:
    fasta_string = f">{'mother: '+ protein.mother + ', code: mit2, start: '+ str(L_to_H(protein.start)) + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
    fasta_L_full_mit2.append(fasta_string)

#L_tRNA_mit = 222 proteins
fasta_L_tRNA_mit2 = []
for prot in prots_L_tRNA_mit2:
    for protein in prot:
        fasta_string = f">{'mother: '+ protein.mother + ', code: mit2, start: '+ str(L_to_H(protein.start))  + '  end: ' + str(L_to_H(protein.end))+ ' '}\n{protein.sequence}\n"
        fasta_L_tRNA_mit2.append(fasta_string)

total_fasta = fasta_L_tRNA_mit2 + fasta_H_tRNA_mit2 + fasta_L_full_mit2 + fasta_H_full_mit2 + fasta_H_full_mit + fasta_H_full_cyt + fasta_L_full_mit + fasta_L_full_cyt + fasta_H_tRNA_mit + fasta_H_tRNA_cyt + fasta_L_tRNA_mit + fasta_H_tRNA_cyt
# 1834: print(len(total_fasta))

print(total_fasta)


# In[15]:


#CODON USAGE


# In[14]:


#have some sort of reference point of codon usage
#http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606.mitochondrion

# Data
data = {
    'ND1': {'len': 956, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0,  'UAA': 0,  'AUG': 3, 'AUU': 10, 'AUA': 13, 'GUG': 0, 'CUG': 5, 'UGA': 9},
                                     '+1': {'AGA': 2, 'AGG': 0, 'UAG': 13, 'UAA': 10, 'AUG': 3, 'AUU': 1,  'AUA': 1 , 'GUG': 1, 'CUG': 2, 'UGA': 1},
                                     '+2': {'AGA': 3, 'AGG': 9, 'UAG': 1,  'UAA': 1,  'AUG': 6, 'AUU': 8,  'AUA': 7 , 'GUG': 1, 'CUG': 2, 'UGA': 4 }}},
    
    'ND2': {'len': 1042, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0,  'UAA': 0,  'AUG': 1, 'AUU': 8, 'AUA': 23, 'GUG': 0, 'CUG': 4, 'UGA': 10},
                                      '+1': {'AGA': 0, 'AGG': 1, 'UAG': 15, 'UAA': 28, 'AUG': 0, 'AUU': 2, 'AUA': 5 , 'GUG': 1, 'CUG': 1, 'UGA': 1},
                                      '+2': {'AGA': 4, 'AGG': 9, 'UAG': 1,  'UAA': 5,  'AUG': 6, 'AUU': 6, 'AUA': 3 , 'GUG': 0, 'CUG': 4, 'UGA': 1}}},
    
    'ND3': {'len': 346, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0, 'UAA': 0,  'AUG': 1, 'AUU': 5, 'AUA': 7, 'GUG': 0, 'CUG': 2, 'UGA': 4},
                                      '+1': {'AGA': 0, 'AGG': 0, 'UAG': 8, 'UAA': 10, 'AUG': 1, 'AUU': 0, 'AUA': 0, 'GUG': 0, 'CUG': 0, 'UGA': 0}, 
                                      '+2': {'AGA': 3, 'AGG': 1, 'UAG': 0, 'UAA': 1,  'AUG': 1, 'AUU': 6, 'AUA': 0, 'GUG': 2, 'CUG': 1, 'UGA': 2}}},
    
    'CO2': {'len': 684, 'codons': {'RF': {'AGA': 0,  'AGG': 0, 'UAG': 2,  'UAA': 0,  'AUG': 2, 'AUU': 7, 'AUA': 8, 'GUG': 0, 'CUG': 3, 'UGA': 3},
                                      '+1': {'AGA': 1,  'AGG': 2, 'UAG': 10, 'UAA': 11, 'AUG': 4, 'AUU': 0, 'AUA': 1, 'GUG': 1, 'CUG': 2, 'UGA': 0},
                                      '+2': {'AGA': 10, 'AGG': 4, 'UAG': 0,  'UAA': 2,  'AUG': 3, 'AUU': 5, 'AUA': 2, 'GUG': 0, 'CUG': 4, 'UGA': 6}}},

    'CO3': {'len': 784, 'codons': {'RF': {'AGA': 0, 'AGG': 0,  'UAG': 0, 'UAA': 0,  'AUG': 3, 'AUU': 7, 'AUA': 8, 'GUG': 1, 'CUG': 3, 'UGA': 9},
                                  '+1': {'AGA': 0, 'AGG': 0,  'UAG': 9, 'UAA': 13, 'AUG': 4, 'AUU': 3, 'AUA': 2, 'GUG': 0, 'CUG': 1, 'UGA': 2},
                                  '+2': {'AGA': 6, 'AGG': 11, 'UAG': 1, 'UAA': 0,  'AUG': 4, 'AUU': 8, 'AUA': 4, 'GUG': 1, 'CUG': 6, 'UGA': 0}}},
    
    'ND4': {'len': 1378, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0,  'UAA': 0,  'AUG': 3, 'AUU': 16, 'AUA': 24, 'GUG': 1, 'CUG': 4, 'UGA': 12},
                                      '+1': {'AGA': 1, 'AGG': 0, 'UAG': 16, 'UAA': 26, 'AUG': 2, 'AUU': 0,  'AUA': 2, 'GUG': 1, 'CUG': 3, 'UGA': 3},
                                      '+2': {'AGA': 6, 'AGG': 8, 'UAG': 2,  'UAA': 5,  'AUG': 7, 'AUU': 10, 'AUA': 3, 'GUG': 0, 'CUG': 7, 'UGA': 1}}},
    
    'ND4L': {'len': 297, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0,  'UAA': 1, 'AUG': 1, 'AUU': 5, 'AUA': 9, 'GUG': 2, 'CUG': 1, 'UGA': 0},
                                       '+1': {'AGA': 0, 'AGG': 0, 'UAG': 10, 'UAA': 4, 'AUG': 1, 'AUU': 0, 'AUA': 3, 'GUG': 0, 'CUG': 0, 'UGA': 0},
                                       '+2': {'AGA': 2, 'AGG': 2, 'UAG': 0,  'UAA': 1, 'AUG': 2, 'AUU': 1, 'AUA': 2, 'GUG': 0, 'CUG': 1, 'UGA': 0}}},
    
    'CO1': {'len': 1545, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 0,  'UAA': 0,  'AUG': 3, 'AUU': 21, 'AUA': 28, 'GUG': 3, 'CUG': 4, 'UGA': 16},
                                      '+1': {'AGA': 2, 'AGG': 0, 'UAG': 19, 'UAA': 26, 'AUG': 6, 'AUU': 2,  'AUA': 2, 'GUG': 3, 'CUG': 12, 'UGA': 4},
                                      '+2': {'AGA': 9, 'AGG': 10, 'UAG': 3,  'UAA': 2,  'AUG': 6, 'AUU': 11, 'AUA': 4, 'GUG': 2, 'CUG': 8, 'UGA': 1}}},
    
    'ND5': {'len': 1821, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 1,  'UAA': 0,  'AUG': 4, 'AUU': 27, 'AUA': 36, 'GUG': 0, 'CUG': 8, 'UGA': 11},
                                      '+1': {'AGA': 3, 'AGG': 0, 'UAG': 27, 'UAA': 35, 'AUG': 6, 'AUU': 1,  'AUA': 2, 'GUG': 0, 'CUG': 2, 'UGA': 2},
                                      '+2': {'AGA': 9, 'AGG': 12, 'UAG': 3,  'UAA': 7,  'AUG': 7, 'AUU': 15, 'AUA': 7, 'GUG': 3, 'CUG': 8, 'UGA': 2}}},

    'ND6': {'len': 525, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 1,  'UAA': 0,  'AUG': 1, 'AUU': 9, 'AUA': 11, 'GUG': 10, 'CUG': 3, 'UGA': 3},
                                  '+1': {'AGA': 0, 'AGG': 0, 'UAG': 7,  'UAA': 10, 'AUG': 1, 'AUU': 1,  'AUA': 0, 'GUG': 7, 'CUG': 2, 'UGA': 8},
                                  '+2': {'AGA': 5, 'AGG': 2, 'UAG': 1,  'UAA': 0,  'AUG': 3, 'AUU': 5, 'AUA': 2, 'GUG': 2, 'CUG': 1, 'UGA': 7}}},
    
    'CYTB': {'len': 1140, 'codons': {'RF': {'AGA': 0, 'AGG': 0, 'UAG': 1,  'UAA': 0,  'AUG': 2, 'AUU': 14, 'AUA': 20, 'GUG': 0, 'CUG': 2, 'UGA': 11},
                                       '+1': {'AGA': 1, 'AGG': 0, 'UAG': 10, 'UAA': 18, 'AUG': 3, 'AUU': 1,  'AUA': 0, 'GUG': 0, 'CUG': 0, 'UGA': 3},
                                       '+2': {'AGA': 7, 'AGG': 6, 'UAG': 1,  'UAA': 3,  'AUG': 4, 'AUU': 9, 'AUA': 4, 'GUG': 1, 'CUG': 7, 'UGA': 1}}}
    

}


codons_list = ['AGA', 'AGG', 'UAG', 'UAA', 'AUG', 'AUU', 'AUA', 'GUG', 'CUG', 'UGA']

frames_list = ['RF', '+1', '+2']

for codon in codons_list:
    adjusted_frequencies = []

    for frame in frames_list:
        frame_values = np.array([data[transcript]['codons'][frame][codon] / data[transcript]['len'] for transcript in data])
        adjusted_frequencies.append(frame_values)

    mean_values = [np.mean(frame_values) for frame_values in adjusted_frequencies]
    std_values = [np.std(frame_values) for frame_values in adjusted_frequencies]

    # Tukey HSD test
    all_values = np.concatenate(adjusted_frequencies)
    labels = np.repeat(frames_list, len(data))
    tukey_result = pairwise_tukeyhsd(all_values, labels)

    plt.figure()
    plt.title(f'Histogram for codon {codon}')
    plt.bar(frames_list, mean_values, yerr=std_values, capsize=5, alpha=0.5)
    plt.xlabel('Reading Frame')
    plt.ylabel('Mean Adjusted Frequency')

    # Add adjusted p-values to the histogram
    #for i, (low, high) in enumerate(tukey_result.confint):
        #comparison = tukey_result._results_table.data[i + 1][0]
       # p_value = tukey_result._results_table.data[i + 1][4]
       # y = max(low, high)
       # plt.text(i, y, f'p = {p_value:.3f}', ha='center', va='bottom')

plt.show()

   


# In[12]:


#check reading frame according to refORF

def codon_count_3frames(sequences, ref_prot_RF):
    data = {}
    for seq_obj in sequences:
        seq_name = seq_obj.name
        seq_len = len(seq_obj.rna_seq)
        seq = seq_obj.rna_seq

        total_count = {}
        codons = {'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0, 'CUU': 0, 'CUC': 0, 'CUA': 0, 'CUG': 0, 'AUU': 0,
                  'AUC': 0, 'AUA': 0, 'AUG': 0, 'GUU': 0, 'GUC': 0, 'GUA': 0, 'GUG': 0, 'UCU': 0, 'UCC': 0,
                  'UCA': 0, 'UCG': 0, 'CCU': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ACU': 0, 'ACC': 0, 'ACA': 0,
                  'ACG': 0, 'GCU': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0}

        for frame in range(3):
            for i in range(frame, len(seq), 3):
                codon = seq[i:i + 3]

                if codon in codons:
                    codons[codon] += 1

            if ref_prot_RF == frame:
                key = 'RF'
            elif (ref_prot_RF + 1) % 3 == frame:
                key = '+1'
            else:
                key = '+2'

            total_count[key] = codons.copy()
            codons = {k: 0 for k in codons}

        data[seq_name] = {'len': seq_len, 'codons': total_count}

    return data

sequence_objects = [ND1,ND2,ND3,ND4,ND4L,ND5, CO1,CO2,CO3,ATP8,ATP6,CYTB, ND6]  
data = codon_count_3frames(sequence_objects, 0)

codons_list = ['UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC',
               'GUA', 'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG', 'ACU', 'ACC', 'ACA', 'ACG',
               'GCU', 'GCC', 'GCA', 'GCG']

frames_list = ['RF', '+1', '+2']

for codon in codons_list:
    adjusted_frequencies = []

    for frame in frames_list:
        frame_values = np.array([data[transcript]['codons'][frame][codon] / data[transcript]['len'] for transcript in data])
        adjusted_frequencies.append(frame_values)

    mean_values = [np.mean(frame_values) for frame_values in adjusted_frequencies]
    std_values = [np.std(frame_values) for frame_values in adjusted_frequencies]

    # Tukey HSD test
    all_values = np.concatenate(adjusted_frequencies)
    labels = np.repeat(frames_list, len(data))
    tukey_result = pairwise_tukeyhsd(all_values, labels)

    plt.figure()
    plt.title(f'Histogram for codon {codon}')
    plt.bar(frames_list, mean_values, yerr=std_values, capsize=5, alpha=0.5)
    plt.xlabel('Reading Frame')
    plt.ylabel('Mean Adjusted Frequency')

    # Add adjusted p-values to the histogram
    #for i, (low, high) in enumerate(tukey_result.confint):
        #comparison = tukey_result._results_table.data[i + 1][0]
       # p_value = tukey_result._results_table.data[i + 1][4]
       # y = max(low, high)
       # plt.text(i, y, f'p = {p_value:.3f}', ha='center', va='bottom')

plt.show()


# In[13]:


def codon_count_3frames(sequences, ref_prot_RF):
    data = {}
    for seq_obj in sequences:
        seq_name = seq_obj.name
        seq_len = len(seq_obj.rna_seq)
        seq = seq_obj.rna_seq

        total_count = {}
        codons = {'UAU': 0, 'UAC': 0, 'UAA': 0, 'UAG': 0,
                  'CAU': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0,
                  'GAC': 0, 'GAA': 0, 'GAG': 0, 'UGU': 0, 'UGC': 0, 'UGA': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0,
                  'CGA': 0, 'CGG': 0, 'AGU': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 'GGU': 0, 'GGC': 0, 'GGA': 0,
                  'GGG': 0}

        for frame in range(3):
            for i in range(frame, len(seq), 3):
                codon = seq[i:i + 3]

                if codon in codons:
                    codons[codon] += 1

            if ref_prot_RF == frame:
                key = 'RF'
            elif (ref_prot_RF + 1) % 3 == frame:
                key = '+1'
            else:
                key = '+2'

            total_count[key] = codons.copy()
            codons = {k: 0 for k in codons}

        data[seq_name] = {'len': seq_len, 'codons': total_count}

    return data

sequence_objects = [ND1,ND2,ND3,ND4,ND4L,ND5, CO1,CO2,CO3,ATP8,ATP6,CYTB, ND6]  
data = codon_count_3frames(sequence_objects, 0)

codons_list = ['UAU', 'UAC', 'UAA', 'UAG', 'CAU', 'CAC', 'CAA', 'CAG', 'AAU', 'AAC',
               'AAA', 'AAG', 'GAU', 'GAC', 'GAA', 'GAG', 'UGU', 'UGC', 'UGA', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG',
               'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA', 'GGG']

frames_list = ['RF', '+1', '+2']

#make function out of this
for codon in codons_list:
    adjusted_frequencies = []

    for frame in frames_list:
        frame_values = np.array([data[transcript]['codons'][frame][codon] / data[transcript]['len'] for transcript in data])
        adjusted_frequencies.append(frame_values)

    mean_values = [np.mean(frame_values) for frame_values in adjusted_frequencies]
    std_values = [np.std(frame_values) for frame_values in adjusted_frequencies]

    # Tukey HSD test
    all_values = np.concatenate(adjusted_frequencies)
    labels = np.repeat(frames_list, len(data))
    tukey_result = pairwise_tukeyhsd(all_values, labels)

    plt.figure()
    plt.title(f'Histogram for codon {codon}')
    plt.bar(frames_list, mean_values, yerr=std_values, capsize=5, alpha=0.5)
    plt.xlabel('Reading Frame')
    plt.ylabel('Mean Adjusted Frequency')

    # Add adjusted p-values to the histogram
    #for i, (low, high) in enumerate(tukey_result.confint):
        #comparison = tukey_result._results_table.data[i + 1][0]
       # p_value = tukey_result._results_table.data[i + 1][4]
       # y = max(low, high)
       # plt.text(i, y, f'p = {p_value:.3f}', ha='center', va='bottom')

plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





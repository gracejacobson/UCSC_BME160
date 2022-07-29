#!/usr/bin/env python3 
# Name: Grace Jacobson
# Group Members: None

'''
converter  input 

This program can:
Take short amino acid abbreviation input and convert to long amino acid abbreviation.
OR
Take long amino acid abbreviation input and convert to short amino acid abbreviation.
OR
Take RNA or DNA three letter codon input and convert to long amino acid abbreviation.
'''

short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
                    # Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dna_Codon_Table = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

def main():
    """The main function takes the input data and prints the value for the dictionary in which it is a key.
    If it does not match a key in any dictionary, it will print unknown."""

    data = input('Enter RNA codon, DNA codon, or amino acid(1 letter or 3 letter code):')
    key = data.upper() #makes the input uppercase to match keys in all dictionaries, in case input is lowercase
    
    if key in short_AA.keys(): #identifies the key in short amino acid abbreviation dictionary
        answerShortAA = short_AA.get(key) #the answer is the value for the key 
        print(key + " = " + answerShortAA) #prints the input with equal to the converted amino acid abbreviation
        
    if key in long_AA.keys(): #identifies key in the long amino acid abbreviation dictionary
        answerLongAA = long_AA.get(key) #the answer is the value for the key in the dictionary
        print(key + " = " + answerLongAA)#prints the input equal to the converted amino acid abbreviation
        
    if key in dna_Codon_Table.keys(): #identifies the key in the DNA codon table
        answerDNA = dna_Codon_Table.get(key) #answer is the value for the key in the dictionary
        print(key + " = " + answerDNA) #prints the input equal to the converted amino acid long abbreviation 
        
    if key in RNA_codon_table.keys(): #identifies the input as a key in the RNA codon table
        answerRNA = RNA_codon_table.get(key) #answer is the value for the key in the dictionary
        print(key + " = " + answerRNA) #prints input equal to the converted amino acid long abbreviation
    
    #if the input is not a key in any of the dictionaries, it will print unknown
    if (key not in short_AA.keys() and key not in long_AA.keys() and key not in dna_Codon_Table.keys() and key not in RNA_codon_table.keys()):
        print("Unknown")
       
main()
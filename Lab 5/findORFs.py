#!/usr/bin/env python3
# Name: Grace Jacobson (gejacobs)
# Group Members: Ken Nicolay

"""
findORFS.py

This program will import FastAreader and ORFfinder from my sequenceAnalysis
module. 

Using these classes and command line, findORFs will be able to return
possible ORFS for a set of parameters, including several start codon(s),
stop codon(s), if longest gene is desired, and minimum gene length.

findORFs takes input from STDIN, and outputs the frame, beginning,
end, and length of ORFs found to STDOUT. 
"""


########################################################################
# Pseudocode/Design
########################################################################
"""
init:
 this method will capture the options from the command line, and also
 generate the dictionary that will be used in addCodons. i will use a
 dictionary of lists, with lists corresponding to the positions of 
 start codon(s) and stop codon(s). this will help me keep track of frames
 and which start codons i am considering. 

 to get the codons for the reverse strand, i have to find the complement,
 which i will search for in the reverse sequence.

addCodons:
 this method will add the positions of the start codon(s) and stop codon(s)
 i will go by steps 3, adding +1 to iterate through each frame, looking for
 codons stored in the init method. i will add positions to the list corresponding
 to each codon in the dictionary. 

orfGen:
 this method will use the positions gathered in addCodons to find ORFs. 

 i have to consider 5 cases:
 -no start or stop in frame
    in this case, i will return the whole length of the sequence, beginning
    to end.
 -no stop, but a start in frame
    in this case, i count from the start to the end of the sequence(top strand),
    or end to the start (bottom strand)
 -no start, but a stop in frame
    in this case, i will count from beginning of the sequence to the stop (top
    strand), or from the stop to the end of the sequence (bottom strand)
 -starts and stops exist in the frame
    in this case, i will have to get the length from the starts and stops
    counted.
 -dangling stop case
    in this case, there is one start/stop at the end of the sequence, with
    multiple starts. I have to consider each start with the same stop.
 
 Then, i will sort each sequence and print it accordingly

"""

from sequenceAnalysis import FastAreader
from sequenceAnalysis import ORFfinder
import sys

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = [],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = [], nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
        #if starts or stops is not set to a parameter
        if self.args.start == []:
            self.args.start = ['ATG'] #default start
        if self.args.stop == []:
            self.args.stop = ['TAG','TGA','TAA'] #default stops


########################################################################
# Main
########################################################################
def main(inFile = None, options = None):
    '''
    Finds genes using ORFfinder and FastAReader from sequenceAnalysis.   
    '''
    thisCommandLine = CommandLine(options)
    myReader = FastAreader('') 
    myOrf = ORFfinder(thisCommandLine)
    for head, seq in myReader.readFasta():
        myOrf.__init__(thisCommandLine)
        print(head)
        myOrf.addCodons(seq)
        myOrf.orfGen(seq)
    
if __name__ == "__main__":
    main(inFile = '', options = sys.argv[1:]) 
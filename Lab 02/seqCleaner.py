#!/usr/bin/env python3 
# Name: Grace Jacobson  
# Group Members: None

""""
seqCleaner

input: a DNA sequence in upper or lower case with AT MOST one block of ambiguous N characters.
output: a DNA sequence in upper case with a block of N characters replaces with the count of ambiguous bases in brackets.

A regular expression is used to replace one or more -indicated by + sign-  Ns in the sequence

Sources
https://docs.python.org/3/library/re.html 
https://cheatography.com/davechild/cheat-sheets/regular-expressions/ 

""""
import re #regular expression

class DNAstring (str):
    """This class represents the DNA string input that needs to be cleaned. 
    Contains a purify method to replace a block of Ns.
    """
    def length (self):
        return (length(self))
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        numberN = self.count("N") #counts the number of N in a single block
        return  re.sub('N+','{' + str(numberN) + '}', self) #replaces a block of one or more Ns with the number of N counted, within brackets
        
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?').upper() #make all bases uppercase in input
    thisDNA = DNAstring (data) #creates object of DNAstring class with input
    pureData = thisDNA.purify() #applies purify method to data
    print (pureData) #output
    
main()
#!/usr/bin/env python3
# Name: Grace Jacobson (gejacobs)
# Group Members: Kenneth Nicolay

'''
findUnique.py

Takes FASTA sequences of mitochondrial tRNA and returns
unique (does not occur in any other sequence) and essential 
(minimized) tRNA sets.

'''
import sys
import codecs
import io
sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding = 'utf-8') #encoding for stdout
sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding = 'utf-8')

########################################################################
# tRNA
########################################################################

class TRNA:
    '''This class contains methods for building powersets, finding unique psets,
    and finding essential psets.'''
    
    def __init__(self, header, sequence):
        '''Creates tRNA object'''
        pass
    
    def _buildPset(self, name, inSeq):
        '''This method finds all possible substrings and stores them in a set'''
        pSet = set()
        x = len(inSeq)
        for i in range(0, x): 
            for j in range(i+1, x+1):
                pSet.add(inSeq[i:j]) #adds to list from i to j
        return pSet
    
    def buildUnique(self, thisAllpSets):
        '''This method finds all the unique tRNA sets in each powerset'''
        unionpSet = {}
        for header1,pSet1 in thisAllpSets.items(): 
            tempUnion = set() #tempunion to be removed
            for header2,pSet2 in thisAllpSets.items():
                if header1!=header2:
                    tempUnion = tempUnion.union(pSet2) #keeps same elements from second pset
            unionpSet[header1] = pSet1.difference(tempUnion) #removes common tRNA sets from first pset
        return unionpSet 
                    
    def buildEssential(self, thisAllUnique):
        '''This method keeps essential tRNA sets by discarding extensions of a string'''
        for header in thisAllUnique.keys():
            tempEssential = thisAllUnique[header].copy()#by discarding while iterating, the set size changes. i copy each time to prevent an error
            for string1 in thisAllUnique[header]:
                for string2 in thisAllUnique[header]:
                    if string1 in string2 and string1!=string2: #if string1 is in string 2 and they are not the same
                        tempEssential.discard(string2) #discards unessential string
            thisAllUnique[header] = tempEssential #temp set with discarded string is set to the input set
        return thisAllUnique

########################################################################
# FastAreader
########################################################################

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return codecs.open(self.fname, encoding='utf-8') #codec encoding
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

########################################################################
# Main
########################################################################

def main(inCL=None):
    '''Finds and prints unique and essential pSets by calling methods of TRNA class'''
    allpSets = {}
    allUnique = {}
    allEssential = {}
    allSequences = {}
    myReader =  FastAreader('')
   
    for header, sequence in myReader.readFasta():
        sequence = sequence.replace('.','') #remove unwanted characters in sequence and headers
        sequence = sequence.replace('_','') 
        sequence = sequence.replace('-','')
        header = header.replace(" ","")
        allSequences[header] = sequence #store header and sequences
        thistRNA = TRNA(header, sequence) 
        allpSets[header] = thistRNA._buildPset(header, sequence) #creates powersets
        
    allUnique = thistRNA.buildUnique(allpSets) #keeps unique psets
    allEssential = thistRNA.buildEssential(allUnique) #keeps unique and essential psets
    
    for header in sorted(allSequences.keys()): #sort by header
        print(header)
        print(allSequences[header])
        for sequence in allSequences.values(): 
            outputDict = {} #create dictionary for output
            for i in range (0, len(sequence)):
                for string1 in allEssential[header]:
                    x = sequence.find(string1, i) #searches for string at each index position
                    if (x == i): #if the index positions are the same
                        outputDict[x] = string1 #adds to output dictionary with index as key
            sortedKeys = sorted(outputDict.keys()) #create sorted list of index positions
            for a in sortedKeys:
                print("."*a + outputDict[a]) #prints periods by index position with string
    pass

if __name__ == "__main__":
    main()  

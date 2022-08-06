import sys

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
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
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
        
class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    rnaCount = {
        '-': {'UAA' : 0, 'UAG' : 0, 'UGA' : 0},
        'A': {'GCU' : 0, 'GCC' : 0, 'GCA' : 0, 'GCG' : 0},
        'G': {'GGU' : 0, 'GGC' : 0, 'GGA' : 0, 'GGG' : 0},
        'M': {'AUG' : 0},
        'S': {'AGU' : 0, 'AGC' : 0, 'UCU' : 0, 'UCC' : 0, 'UCA' : 0, 'UCG' : 0},
        'C': {'UGU' : 0, 'UGC' : 0},
        'H': {'CAU' : 0, 'CAC' : 0},  
        'N': {'AAU' : 0, 'AAC' : 0},  
        'T': {'ACU' : 0, 'ACC' : 0, 'ACA' : 0, 'ACG' : 0},
        'D': {'GAU' : 0, 'GAC' : 0},
        'I': {'AUU' : 0, 'AUC' : 0, 'AUA' : 0},
        'P': {'CCU' : 0, 'CCC' : 0, 'CCA' : 0, 'CCG' : 0},
        'V': {'GUU' : 0, 'GUC' : 0, 'GUA' : 0, 'GUG' : 0},
        'E': {'GAA' : 0, 'GAG' : 0},
        'K': {'AAA' : 0, 'AAG' : 0}, 
        'Q': {'CAA' : 0, 'CAG' : 0},
        'W': {'UGG' : 0},
        'F': {'UUU' : 0, 'UUC' : 0}, 
        'L': {'UUA' : 0, 'UUG' : 0, 'CUU' : 0, 'CUC' : 0, 'CUA': 0, 'CUG' : 0},
        'R': {'CGU' : 0, 'CGC' : 0, 'CGA' : 0, 'CGG' : 0, 'AGA': 0, 'AGG' : 0},
        'Y': {'UAU' : 0, 'UAC' : 0},
    }
    
    
    nucComp = { 'A': 0, 'T' : 0, 'C' : 0, 'G': 0, 'N' : 0 , 'U' : 0}

    def __init__ (self, inString=''):
        self.addSequence(inString)
        return
     
    def addSequence (self, inSeq):
        for i in range (0, len(inSeq), 3):
            codon = inSeq[i:i+3]
            rnaCodon = codon.replace("T", "U")
            if rnaCodon in (self.rnaCodonTable.keys()):
                if rnaCodon in (self.rnaCount[self.dnaCodonTable[codon]].keys()):
                    self.rnaCount[self.rnaCodonTable[rnaCodon]][rnaCodon] += 1         
        for i in inSeq:
            self.nucComp[i] += 1
        return
    
    def aaComposition(self):
        aaComp = {}
        for aa in self.rnaCount.keys():
            aaComp[aa] = sum(self.rnaCount[aa].values())
        return aaComp
        
    def nucComposition(self):
        return self.nucComp
    
    def codonComposition(self):
        codonComp = {}
        for codondict in self.rnaCount.values():
            for codon, value in codondict.items():
                codonComp[codon] = value    
        return codonComp
    
    def nucCount(self):
        nucCount = sum(self.nucComp.values())
        return nucCount
    
def main ():
    myReader = FastAreader('testGenome2.fa')# make sure to change this to use stdin
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    
    #print sequence length
    nucCount = (myNuc.nucCount()) / 1000000 #convert to Mb
    print("sequence length = {0:.2f} Mb".format(nucCount)) #print with 2 decimals
    
    #print GC content
    GCcontent = ((myNuc.nucComp['G']+ myNuc.nucComp['C']) / (myNuc.nucCount())) * 100
    print('\nGC content = {0:.1f}%\n '.format(GCcontent))
     
    # sort codons in alpha order, by Amino Acid
    for aa in sorted(myNuc.rnaCount.keys()):
        denom = sum(myNuc.rnaCount[aa].values())
        for codon in sorted(myNuc.rnaCount[aa].keys()):
            if denom == 0:
                val = 0
            else:
                val = myNuc.rnaCount[aa][codon] / denom
            print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val*100, myNuc.rnaCount[aa][codon]))

if __name__ == "__main__":
    main()
    
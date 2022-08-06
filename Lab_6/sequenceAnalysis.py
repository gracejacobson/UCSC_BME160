import sys 

########################################################################
# ORFfinder
######################################################################## 
class ORFfinder:
    """
    Takes FASTA DNA sequence and finds possible open reading frames (ORFS).
    Prints frame (+ for top strand, - for bottom strand), beginning position
    
    (first base of start codon for top strand, last base of stop codon for
    bottom strand), end position (last base of stop codon for top strand,
    first base of start codon for bottom strand), and length of the ORF.

    """
    startCodon =  []
    stopCodon = []
    startCodonNeg =[] #i use neg to refer to bottom strand (- frames)
    stopCodonNeg = []
    codonDict = {}        
    minGene = 0
    longestGene = False

    def __init__ (self, commandLine='', inString =''):
        """This method initializes command line parameters and adds desired
        codons for each frame in the dictionary"""
        self.codonDict = {
            '+1':{},  
            '-1':{},    
            '+2':{},  
            '-2':{},
            '+3':{},  
            '-3':{}}
        self.startCodon =  []
        self.stopCodon = []
        self.startCodonNeg =[]
        self.stopCodonNeg = []
        self.startCodon = commandLine.args.start #adds start codons from command line to list
        print(self.startCodon)
        self.stopCodon = commandLine.args.stop #adds stop codons from command line to list
        self.minGene = commandLine.args.minGene #creates minGene parameter
        self.longestGene = commandLine.args.longestGene #creates longestgene parameter

        #positive frames
        for codon in self.startCodon:
            for frame in ['+1', '+2', '+3']:
                self.codonDict[frame][codon] =[] #adds start codon(s) to dictionary
        for codon in self.stopCodon:
            for frame in ['+1', '+2', '+3']:
                self.codonDict[frame][codon] =[] #adds stop codon(s) to dictionary
        #negative frames
        for rawCodon in self.startCodon: #creates complement to start codon(s)
            codon = rawCodon.replace('A', 'X')
            codon = codon.replace('C', 'Y')
            codon = codon.replace('T', 'A')
            codon = codon.replace('G', 'C')
            codon = codon.replace('X', 'T')
            codon = codon.replace('Y', 'G')
            self.startCodonNeg.append(codon) #adds to list of bottom strand start codons       
            for frame in ['-1', '-2', '-3']:
                self.codonDict[frame][codon] =[]
        for rawCodon in self.stopCodon: #creates complement to stop codon(s)
            codon = rawCodon.replace('A', 'X')
            codon = codon.replace('C', 'Y')
            codon = codon.replace('T', 'A')
            codon = codon.replace('G', 'C')
            codon = codon.replace('X', 'T')
            codon = codon.replace('Y', 'G')
            self.stopCodonNeg.append(codon) #adds to list of bottom strand stop codons
            for frame in ['-1', '-2', '-3']:
                self.codonDict[frame][codon] =[] #adds to dictionary
        return
    
    def addCodons(self, inSeq):
     """Adds positions of desired start(s) and stop codon(s) to codonDict"""

     inSeqReverse = inSeq[::-1] #taking reverse of sequence to look for codons in bottom strand
        
     for i in range (0, len(inSeq), 3):

        #top strand codons
        frame1Codon = inSeq[i:i+3]
        frame2Codon = inSeq[i+1:i+1+3]
        frame3Codon = inSeq[i+2:i+2+3]
           
        #Start Codon(s) for top strand
        for codon in self.startCodon: #searches for codons matching start codons in list
            if (frame1Codon == codon):
                self.codonDict['+1'][codon].append(i+1) #adds to dictionary
            elif (frame2Codon == codon):
                self.codonDict['+2'][codon].append(i+2)
            elif (frame3Codon == codon):
                self.codonDict['+3'][codon].append(i+3)
        for codon in self.stopCodon: #searches for stop codons matching stop codon in list
            if (frame1Codon == codon):
                self.codonDict['+1'][codon].append(i+3) #adds to dictionary
            elif (frame2Codon == codon):
                self.codonDict['+2'][codon].append(i+4)
            elif (frame3Codon == codon):
                self.codonDict['+3'][codon].append(i+5)            

        #bottom strand codons
        frame1CodonNeg = inSeqReverse[i:i+3]
        frame2CodonNeg = inSeqReverse[i+1:i+1+3]
        frame3CodonNeg = inSeqReverse[i+2:i+2+3]                 
                
        #Start Codon(s) for bottom strand
        for codon in self.startCodonNeg: #searches for codons matching complement start codon in list
            if (frame1CodonNeg == codon):
                self.codonDict['-1'][codon].append(len(inSeq)-i) #subtract from i because its reverse
            elif (frame2CodonNeg == codon):
                self.codonDict['-2'][codon].append(len(inSeq)-i-1)
            elif (frame3CodonNeg == codon):
                self.codonDict['-3'][codon].append(len(inSeq)-i-2)
        
        #Stop Codon(s)  for bottom strand
        for codon in self.stopCodonNeg:
            if (frame1CodonNeg == codon):
                self.codonDict['-1'][codon].append(len(inSeq)-i-2) 
            elif (frame2CodonNeg == codon):
                self.codonDict['-2'][codon].append(len(inSeq)-i-3)
            elif (frame3CodonNeg == codon):
                self.codonDict['-3'][codon].append(len(inSeq)-i-4)   
     return
 
    def orfGen (self, inSeq):
     """This method uses positions from addCodons to find ORFs in sequence."""
     
     orfData = [] #this is the list i will be using to compile the data to print for the ORF
      
     for frame in self.codonDict.keys():
         startsPos = [] #list of start positions for top strand
         stopsPos = [] #list of stop positions for top strand
         startsNeg = [] #list of start positions for bottom strand
         stopsNeg = [] #list of stop positions for bottom strand

         #finding ORFS in top strand
         if '+' in frame:
            
             #adds positions for start codons in frame
             for codon in self.startCodon:
                 for element in self.codonDict[frame][codon]:
                     startsPos.append(element) #adds to list
             #adds stop positions for stop codons in frame
             for codon in self.stopCodon:
                 for element in self.codonDict[frame][codon]:
                     stopsPos.append(element)#adds to list

             startsPos.sort() #sort by size
             stopsPos.sort()

            #dangling stop with starts
            #i don't need to consider the +1 frame because i consider all starts and stops later
             if (frame=='+2') or (frame=='+3'): 
                 if (len(stopsPos)==1) and (len(startsPos)>0): #check that there is one stop, with multiple starts
                     for x in range(len(stopsPos)):
                         dangleStop = len(inSeq) - stopsPos[x] #finds leftover sequence after dangle stop
                         orfData.append((int(frame), (int(frame)), 1, stopsPos[x], len(inSeq)-dangleStop))
                         #i record the frame twice, so i can sort by + and - frames later
            
            #no starts and stops
             if not (startsPos) and not (stopsPos):
                 #i simply take the length of the entire sequence since there is no start and stop
                 orfData.append((int(frame), int(frame), 1, len(inSeq), len(inSeq)))

            #no starts
             if not (startsPos) and (len(stopsPos)>0): #check there are no starts and at least 1 stop
                 for x in range(len(stopsPos)):
                     for y in stopsPos:
                         noStartLenPos = y #length of sequence is from start to stop
                         orfData.append((int(frame),int(frame), 1, stopsPos[x] , noStartLenPos))
            
            #no stops
             if not (stopsPos) and (len(startsPos)>0): #check there are starts and no stops
                 for x in range(len(startsPos)):
                     for y in startsPos:
                        noStopLenPos = len(inSeq) - y #length is from end to start
                        #i return this length and set the start position to the start to end of sequuence
                        orfData.append((int(frame),int(frame), startsPos[x], len(inSeq), noStopLenPos))
            
            #starts and stops
             if (len(startsPos)>0) and (len(stopsPos)>0):
                 #if longest gene is not desired
                 if self.longestGene == False:
                    for x in range(len(stopsPos)):
                         for y in range(len(startsPos)):
                             #the length is all possible stops positions minus start positions
                             lengthsPos = stopsPos[x] - startsPos[y]
                             if lengthsPos>0: #only return positive lengths
                                 # i add to list with 
                                 orfData.append((int(frame),int(frame), startsPos[y], stopsPos[x],lengthsPos+1))

                 #if longest gene IS desired
                 else:

                     smallestStart=[]
                     prevStop = 0
                     # iterate through each stop, looking for the earliest start after the last stop
                     for stop in stopsPos:
                        for start in startsPos:
                            if (prevStop < start < stop): #if start is between stops
                                smallestStart.append(start)
                        if(len(smallestStart)>0):
                            minStart = min(smallestStart) #take smallest start for maximum length
                            lengthsPos = stop - minStart + 1
                            orfData.append((int(frame)*-4,int(frame), minStart, stop, lengthsPos))
                            smallestStart=[] #reset to move onto next sotp
                            prevStop = stop #set current stop to previous to move onto next stop
                        else:
                            # no start or starts found to get the min of, so move on to next stop
                            next

         #finding ORFs in bottom strand
         if '-' in frame:      
             #add positions to list for start codon(s)
             for codon in self.startCodonNeg:
                 for element in self.codonDict[frame][codon]:
                     startsNeg.append(element)
             #adds positions to list for stop codons
             for codon in self.stopCodonNeg:
                 for element in self.codonDict[frame][codon]:
                     stopsNeg.append(element)

             stopsNeg.sort() #sort by size
             startsNeg.sort()

            #no starts and stops
             if not (startsNeg) and not (stopsNeg): #check that there are no start and stops in list
                 #i return the sequence from start to end with the whole length
                 #here i multiply the frame by -4, this is to make the most negative frames largest
                 #i do this so i can sort by reverse frames, having +1 first and -3 (now 12) last.
                 orfData.append(((int(frame)*-4),int(frame), 1, len(inSeq), len(inSeq)))
                
            #no starts
             if not (startsNeg) and (len(stopsNeg)>0): #check that there are no starts, and at least 1 stop
                 for x in range(len(stopsNeg)):
                     for y in stopsNeg:
                         noStartLenNeg = len(inSeq) - y  #length is the sequence minus the start position
                         orfData.append(((int(frame)*-4),int(frame), stopsNeg[x], len(inSeq), noStartLenNeg+1))
                 
            #dangling stop with starts
             if (frame=='-2') or (frame=='-3'): 
             #again, i don't consider the -1 frame because i consider all starts and stops later
                 if (len(stopsNeg)==1) and (len(startsNeg)>0): #check there is one stop and at least 1 start
                     for x in range(len(stopsNeg)):
                         dangleStop = len(inSeq) - stopsNeg[x] +1 #length is sequence minus stop position
                         orfData.append(((int(frame)*-4), (int(frame)), stopsNeg[x], len(inSeq), dangleStop))
            
        
            #no stops
             if not (stopsNeg) and (len(startsNeg)>0): #check there are no stops, and at least 1 start
                 for x in range(len(startsNeg)):
                     for y in startsNeg:
                         noStopLenNeg = len(inSeq) - y  #length is sequence length minus the start position
                         orfData.append(((int(frame)*-4),int(frame), len(inSeq), startsNeg[x], noStopLenNeg+1))
                
                
            #starts and stops
             if (len(startsNeg)>0):
                 #if longest gene is not desired
                 if self.longestGene == False:
                     for x in range(len(stopsNeg)):
                         for y in range(len(startsNeg)):
                             #i find all possible lengths by subtracting possible stops from possible starts
                             lengthsNeg = startsNeg[y] - stopsNeg[x]  
                             if (len(stopsNeg)>0): #only return postive lengths that exist
                                 #i append to the orfData list
                                 orfData.append(((int(frame)*-4),int(frame), stopsNeg[x], startsNeg[y], lengthsNeg+1))
                 #if longest gene is desired
                 else:
                     biggestStart=[]
                     prevStop = 0
                     for stop in stopsNeg:
                        for start in startsNeg:
                            if (prevStop < start < stop): #look for start between stops
                                biggestStart.append(start)
                        if(len(biggestStart)>0):
                            maxStart = max(biggestStart) #take maximum start for longest length
                            lengthsNeg = maxStart -prevStop + 1
                            orfData.append((int(frame)*-4,int(frame), prevStop, maxStart, lengthsNeg)) #add to orfData
                        else:
                            # no start or starts found to get the min of, so move on to next stop
                            next
                         #reset because we are going to the next stop
                        biggestStart=[]
                        prevStop = stop

         #here i sort my list ORFS, each containing a list of frames, starts, ends, and lengths
         #i sort by l[4] which is sequence lengths first, then by start positions which is l[2]
         #then finally i sort by frames in reverse, where -3 frame is 12, -2 frame is 8, etc.
         orfData.sort(reverse=True, key = lambda l: (l[4], -l[2], -l[0]))

    #print statements

     #for longest gene desired
     if self.longestGene == True:
         for x in orfData:
             #print only if length is larger or equal to minimum gene 
             if x[4] > self.minGene:
                 #print with positive sign for frames
                 if x[1]>=0:
                    print(('{0:+}'.format(x[1]))+ "\t" + str(x[2]) + "..\t" + str(x[3])+ "\t" +str(x[4]))
                 #print statement for reverse frames
                 if x[1]<0:
                    print(str(x[1])+ "\t" +str(x[2])+ "..\t" +  str(x[3])+ "\t" + str(x[4]))

     #if longest gene not desired
     if self.longestGene == False:  
         for x in orfData:
             #print only if length is larger or equal to minimum gene 
             if x[4] > self.minGene:
                 #print with positive sign for frames
                 if x[1]>=0:
                    print(('{0:+}'.format(x[1]))+ "\t" + str(x[2]) + "..\t" + str(x[3])+ "\t" +str(x[4]))
                 #print statement for reverse frames
                 if x[1]<0:
                    print(str(x[1])+ "\t" +str(x[2])+ "..\t" +  str(x[3])+ "\t" + str(x[4]))


     return 

########################################################################
# NucParams
########################################################################                 
class NucParams:
    """
    Takes DNA or RNA sequence string and keeps dictionary of amino acid count, 
    codon count, and nucleotide (A,T,C,G,N, and U) counts. 
    
    Prints sequence length in Mb, GC content %, and codon with corresponding 
    amino acid with relative frequency and count.
    """
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
    } #i used a nested dictionary with the corresponding codons for each amino acid
    #this makes it simpler to return codon composition and codon composition
    
    #nucleotide dictionary
    nucComp = { 'A': 0, 'T' : 0, 'C' : 0, 'G': 0, 'N' : 0 , 'U' : 0}

    def __init__ (self, inString=''):
        """Calls addSequence"""
        self.addSequence(inString)
        return
     
    def addSequence (self, inSeq):
        """Accepts sequence and counts codons and nucleotides"""
        for i in range (0, len(inSeq), 3):
            codon = inSeq[i:i+3]
            rnaCodon = codon.replace("T", "U") #converts to RNA if input was in DNA
            if rnaCodon in (self.rnaCodonTable.keys()):
                if rnaCodon in (self.rnaCount[self.rnaCodonTable[rnaCodon]].keys()): #verify codon exists in dictionary
                    self.rnaCount[self.rnaCodonTable[rnaCodon]][rnaCodon] += 1 #counts codons in nested dictionary corresponding to amino acid      
        for i in inSeq:
            self.nucComp[i] += 1 #counts nucleotides
        return
    
    def aaComposition(self):
        """Returns dictionary of counts for amino acid and stop codons"""
        aaComp = {}
        for aa in self.rnaCount.keys():
            aaComp[aa] = sum(self.rnaCount[aa].values()) #creates dictionary with key as amino acid and value as sum of codon values from nested rnaCount dict
        return aaComp
        
    def nucComposition(self):
        """ Returns dictionary of nucleotide counts"""
        return self.nucComp
    
    def codonComposition(self):
        """Returns dictionary of codon counts"""
        codonComp = {}
        for codondict in self.rnaCount.values(): 
            for codon, value in codondict.items():
                codonComp[codon] = value #creates dictionary of each codon with corresponding value from rnaCount dict
        return codonComp
    
    def nucCount(self):
        """ Returns count of every valid nucleotide found"""
        nucCount = sum(self.nucComp.values()) #sums values for nucComp dict
        return nucCount

########################################################################
# ProteinParam
######################################################################## 
class ProteinParam :
    """
    This program calculcates the number of amino acids, molecular weight, molar extinction coefficient,
    mass extinction coefficient, theoretical pI, and amino acid composition of an amino acid sequence.

    If no amino acid is entered, the program will not run.
    If a sequence is entered and it contains a character that is not an amino acid, a warning message
    appears and the program continues.


    """
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
    
    #aaComp dictionary keeps a count of the number of each amino acid in the sequence with an intitial value of 0
    aaComp = {
        'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
        'H': 0,  'N': 0,  'T': 0, 'D': 0, 'I': 0,
        'P': 0,  'V': 0,  'E': 0, 'K': 0, 'Q': 0,
        'W': 0,  'F': 0,  'L': 0, 'R':0,  'Y': 0  
        }
    
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    e = 0.01
    
    def __init__ (self, protein):
        """This method computes and saves the aaComposition dictionary from the input amino acid string, all amino acids have an intiial value of 0"""

        for aa in protein: #iterates over protein sequence
            if aa not in self.aaComp.keys(): 
                print ("Warning: Sequence contains invalid amino acid.") #if something is entered that is not an AA, a warning message is printed
            else: self.aaComp [aa] += 1 #adds one value to aaComp dictionary for each amino acid counted
        
        
        
    def aaCount (self):
        """This method computes the number of amino acids in the input string"""
        aaNum = sum(self.aaComp.values()) #creates aaNum object which is the sum of aaComp dictionary values for each AA
        return aaNum #returns number of AA

    def pI (self):
        """This method prints the pH which gives a neutral net charge."""
        
        for pH in range(0,1400+1): #includes pH of 1400
            charge = self._charge_(pH/100) #calls charge method

            if (abs(charge) < self.e): #e for 2 decimal accuracy of pH, close to 0
                return pH/100 #returns pH
            
        return -1 #returns -1 if no pH is within precision of e
            
    def aaComposition (self) :
        """This returns the dictionary of each amino acid and its count in the input sequence"""
        return self.aaComp #returns aaComp dictionary with each AA counted

    def _charge_ (self, pH): #added pH
        """This method calculates the net charge on a protein at a specific pH"""
        
        #positive charges of the netcharge equation
        kCharge = self.aaComp['K'] * ((10 ** self.aa2chargePos['K']) / ((10 ** self.aa2chargePos['K']) + 10 ** pH))
        rCharge = self.aaComp['R'] * ((10 ** self.aa2chargePos['R']) / ((10 ** self.aa2chargePos['R']) + 10 ** pH))
        hCharge = self.aaComp['H'] * ((10 ** self.aa2chargePos['H']) / ((10 ** self.aa2chargePos['H']) + 10 ** pH))
        nTermCharge = (10 ** self.aaNterm) / ((10 ** self.aaNterm ) + (10 ** pH))
        
        #adding positive charges to get the total positive charges
        posCharge = kCharge + rCharge + hCharge + nTermCharge
        
        #negative charges of the netcharge equation
        dCharge = self.aaComp['D'] * ((10 ** pH)/ ((10 ** self.aa2chargeNeg['D']) + 10 ** pH))
        eCharge = self.aaComp['E'] * ((10 ** pH)/ ((10 ** self.aa2chargeNeg['E']) + 10 ** pH))
        cCharge = self.aaComp['C'] * ((10 ** pH)/ ((10 ** self.aa2chargeNeg['C']) + 10 ** pH))
        yCharge = self.aaComp['Y'] * ((10 ** pH)/ ((10 ** self.aa2chargeNeg['Y']) + 10 ** pH))
        cTermCharge = (10 ** pH) / ((10 ** self.aaCterm) + (10 ** pH))
        
        #adding negative charges up to get the total negative charges
        negCharge = dCharge+ eCharge + cCharge + yCharge + cTermCharge

        #create netCharge object which returns the answer
        netCharge = posCharge - negCharge
        
        return netCharge #returns the answer
        
    def molarExtinction (self):
        """This method calculates the molar extinction coefficient at 280nm, this indicates how much light a protein absorbs"""
        eY = self.aaComp['Y'] * self.aa2abs280['Y'] #extinction coefficient for tyrosine multipled by the number in the sequence
        eW = self.aaComp['W'] * self.aa2abs280['W'] #extinction coefficient for tryptophan multipled by the number of tryptophan in the sequence
        eC = self.aaComp['C'] * self.aa2abs280['C'] #extinction coefficient for cysteine multiplied by the number of cysteines in the sequence
        return eY + eW + eC #returns the molar extinction coefficient

    def massExtinction (self):
        """This method calculates the mass extinction coefficient from the molar extinction coefficient 
        by dividing by the molecular weight of the protein sequence"""
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        """This method calculates the molecular weight of the amino acid sequence"""
        aaNum = sum(self.aaComp.values()) #create aaNum object for molecular weight calculation 
        mwAA = sum(self.aaComp[aa] * self.aa2mw[aa] for aa in self.aaComp) #creates mw object which is the sum of the values in aaComp times their molecular weight for amino acids in aaComp
        return mwAA - ((aaNum-1) * self.mwH2O) #returns mw for the protein sequence 

########################################################################
# FastAReader
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

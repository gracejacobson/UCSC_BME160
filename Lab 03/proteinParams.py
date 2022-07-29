#!/usr/bin/env python3
# Name: Grace Jacobson (gejacobs)
# Group Members: None
"""
ProteinParam

This program calculcates the number of amino acids, molecular weight, molar extinction coefficient,
mass extinction coefficient, theoretical pI, and amino acid composition of an amino acid sequence.

If no amino acid is entered, the program will not run.
If a sequence is entered and it contains a character that is not an amino acid, a warning message
appears and the program continues.


"""
class ProteinParam :
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

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()
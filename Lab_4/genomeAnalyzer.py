from sequenceAnalysis import NucParams
from sequenceAnalysis import FastAreader

def main ():
    """Calls FastAreader and NucParams, prints output"""
    myReader = FastAreader()
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    
    #print sequence length
    nucCount = (myNuc.nucCount()) / 1000000 #convert to Mb
    print("sequence length = {0:.2f} Mb".format(nucCount)) #prints with 2 decimals
    
    #print GC content
    GCcontent = ((myNuc.nucComp['G']+ myNuc.nucComp['C']) / (myNuc.nucCount())) * 100 #calculates percentage of G and C by getting value from nucComp dict
    print('\nGC content = {0:.1f}%\n '.format(GCcontent)) #prints GC content
     
    #sort and print codons in alphabetical order by amino acid
    for aa in sorted(myNuc.rnaCount.keys()):
        denom = sum(myNuc.rnaCount[aa].values()) #creates denominator for division for relative frequency of codon for amino acid
        for codon in sorted(myNuc.rnaCount[aa].keys()):
            if denom == 0:
                val = 0 #avoids division by zero and creating error
            else:
                val = myNuc.rnaCount[aa][codon] / denom #finds relative frequency of codon for an amino acid
            print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val*100, myNuc.rnaCount[aa][codon]))

if __name__ == "__main__":
    main()
    
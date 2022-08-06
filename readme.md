# BME 160: Research Programming in the Life Sciences
This repo contains assignments from the BME160 course at UCSC which focused on programming in Python to manipula te biological data. I took this course in Spring 2021.

## Lab 2:
1. **converter.py**: Convert from amino acid abbreviations to DNA/RNA codons and vice versa.
2. **coordinateMathSoln.py**: Calculate bond length and angles between C-Ca-N.
3. **fastqParse.py**: Filter out run data from FASTQ files.
4. **seqCleaner.py**: Returns sequence with uppercase and replaces block of ambiguous bases with their count.

## Lab 3:
1. **proteinparams.py**: Calculates protein properties such as molecular weight, amino acid composition, theoretical PI, etc.

## Lab 4:
1. **FastAreader.py**: Reads FastA files for analysis.
2. **genomeAnalyzer.py**: Calls FastAreader and nucparams and formats their outputs.
3. **nucparams.py**: Calculates sequence properties such as amino acid composition and codon usage.
4. **sequenceAnalysis.py**: Contains NucParams,FastAreader, and ProteinParams classes.

## Lab 5:
1. **findORFs.py**: findORFs will be able to return possible ORFS for a set of parameters, including several start codon(s), stop codon(s), if longest gene is desired, and minimum gene length. Calls fastAreader and ORFfinder from sequenceAnalysis.py.
2. **sequenceAnalysis.py**: Updated to include ORFfinder class to find all possible open reading frames.

## Lab 6:
1. **findUnique.py**: Takes FASTA sequences of mitochondrial tRNA and returns unique (does not occur in any other sequence) and essential (minimized) tRNA sets.
2. **sequenceAnalysis.py**: Inludes ORFfinder, NucParams, FastAreader, and ProteinParams classes.
   


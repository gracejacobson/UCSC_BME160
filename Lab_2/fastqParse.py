#!/usr/bin/env python3 
# Name: Grace Jacobson
# Group Members: None

'''
fastqParse

input: a FASTQ seqname line starting with @, with seven fields separated by :
output: run information with each field on a new line, with an = sign


'''

class FastqString (str):
    ''' A class to represent a FASTQ string containing run information. Contains parse method to separate each field.'''
    
    def parse(self):
        ''' This method will separate each field and print the run information.'''
        parseList = self.split(":") #forms a list of each field, which were separated by : . the first item is a blank space therefore the program starts printing at item 1
        print("Instrument = "+parseList[1]) #prints the second item in the list, the instrument
        print("Run ID = "+parseList[2]) #prints the third item in the list, the Run ID
        print("Flow Cell ID = "+parseList[3]) #prints the fourth item in the list, the flow cell ID
        print("Flow Cell Lane = "+parseList[4]) #prints the fifth item in the list, the flwo cell lane
        print("Tile Number = "+parseList[5]) #prints the sixth item in the list, the tile number
        print("X-coord = "+parseList[6]) #prints the seventh item in the list, the x-coordinate
        print("Y-coord = "+parseList[7]) #prints the eighth item in the list, the y-coordinate
    
def main():
    ''' The main function will take an FASTQ seqname line and separate it from the @ symbol and direct to the parse method '''
    data = input("Enter FASTQ seqname line:") 
    fastqList = data.split('@') #separates the seven fields from the @ symbol, creating a list with two items of an empty space and the rest of the seqname line
    fastqString = ":".join(fastqList) #creates a string, joining the space to the rest of the seven fields with a :
    FastqString.parse(fastqString) #applies the parse method of the FastqString object to the string
    
main()
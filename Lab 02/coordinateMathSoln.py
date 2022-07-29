#!/usr/bin/env python3 
# Name: Grace Jacobson
# Group Members: None

'''
coordinateMathSoln

Takes coordinates contained in parantheses for C, N, and Ca atoms with an '=' sign in one input line. 
Prints N-C and N-Ca bond lengths and C-N-Ca bond angle in degrees as its output.

'''

import math
import re

class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
        
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

def main():
    ''' 
    This function uses regular expressions to find three coordinates for C, N, and Ca. Using a regular
    expression allows the program to find the coordinates for each point no matter the number of spaces between
    each atom and '=' sign. After forming the triad, the function prints each bond length and bond angle. 
    The bond angle is converted to degrees.
    
    Resources:
            https://www.regular-expressions.info/floatingpoint.html for regular expressions to find coordinates
            https://regexr.com/ for editing the regular expression
            https://www.geeksforgeeks.org/python-convert-list-of-tuples-into-list/ for converting tuples into a list
            
    '''
    coordinates = input("Enter Coordinates: ")
    
    #creating list of p coordinates using a regular expression
    pTuples = re.findall(".*C+\s*=\s*\([-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)", coordinates)
    p = [] #creates p as a list of coordinates
    for x in pTuples: #converts the tuple into a mutable list
            for t in x:
                p.append(float(t))
                
    #creating a list of q coordinates using a regular expression
    qTuples = re.findall(".*N+\s*=\s*\([-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)", coordinates)
    q = [] #creates q as a list of coordinates
    for x in qTuples: #converts the tuple into a mutable list
            for t in x:
                q.append(float(t))
    
    #creating a list of r coordinates using a regular expression
    rTuples = re.findall(".*Ca+\s*=\s*\([-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)\s*,\s*[-+]?([0-9]*\.[0-9]+|[0-9]+)", coordinates)
    r = [] #creates r as a list of coordinates
    for x in rTuples: #converts the tuple into a mutable list
            for t in x:
                r.append(float(t))
    
    p1 = Triad(p,q,r) #forms triad
    print("N-C bond length = {0:0.2f}".format(p1.dPQ())) #prints the bond length with 2 decimal places
    print("N-Ca bond length = {0:0.2f}".format(p1.dQR())) #prints the bond length with 2 decimal places
    
    angleDegrees = math.degrees(p1.angleQ()) #i used math.degrees to convert angleQ to degrees
    print("C-N-Ca bond angle = {0:0.1f}".format(angleDegrees)) #prints the bond angle in degrees with one decimal place
    
main()
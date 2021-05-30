#!/usr/bin/python3
import sys
from biquad.biquadfunc import BiquadField as BF

#the script has 3 input parameters p, q, max_trace
#if these are not provided, takes the default ones below
args = sys.argv
if len(args)<4:
   p = 13
   q = 5
   max_trace = 65
else:
   p = int(args[1])
   q = int(args[2])
   max_trace = int(args[3])

print("Start of computation.")
K = BF.createField(p,q) #Creates the field class instance
#Computes all sums of squares up to max_trace and their lengths
#number of elements is stored in s, elemenets in squares (for possible further use)
s, squares = K.lengths_up_to_trace(max_trace,return_squares = True)

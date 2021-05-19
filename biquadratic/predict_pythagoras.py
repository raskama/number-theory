import sys
import biquad
from biquad import BiquadField as BF


args = sys.argv
#the script has 3 input parameters p, q, max_trace
#i.e. can be called like python predict_pythagoras.py p q max_trace
#if these are not provided, takes the default ones below
if len(args)<4:
   p = 13
   q = 5
   max_trace = 5
else:
   p = int(args[1])
   q = int(args[2])
   max_trace = int(args[3])



print("Start of computation.")
K = BF.createField(p,q)
K.update_elements(max_trace)
K.compare_together(max_trace)

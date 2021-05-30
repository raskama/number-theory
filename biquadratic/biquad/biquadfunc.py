#!/usr/bin/python3

import numpy as np
import sympy #for prime generator
from sympy.solvers.diophantine.diophantine import square_factor
import sys
from itertools import repeat
import multiprocessing
from multiprocessing import freeze_support
from abc import ABC, abstractmethod
from math import gcd

def isSquareFree(n):
   if square_factor(n)>1:
      return False
   else:
      return True

#sqrt working even for high values, for high values gives floor of sqrt 
def isqrt(x):
    if x < 0:
        raise ValueError('Square root is not defined for negative numbers')
    if x < (1 << 50):
        return np.sqrt(x)  # use math's sqrt() for small parameters
    n = int(x)
    if n <= 1:
        return n  # handle sqrt(0)==0, sqrt(1)==1
    
    r = 1 << ((n.bit_length() + 1) >> 1)
    while True:
        newr = (r + n // r) >> 1  # next estimate by Newton-Raphson
        if newr >= r:
            return r
        r = newr


#Checks if floor((a*sqrt(D)+b)/c) == x, returns -1 if floor is less, 0 if correct, +1 if more
def checkFloor(D,a,b,c,x):
   L = (x*c-b)*(x*c-b)
   R = (x*c-b+c)*(x*c-b+c)
   if a>=0 and c>=0:
      if c*x-b>=0 and (L > D*a*a):
         return 1 
      if c*x+c-b < 0 or (a*a*D >= R):
         return -1
   elif a >=0 and c<0:
      if c*x-b < 0 or (a*a*D > L):
         return 1
      if c*x-b+c>=0 and (R >= D*a*a):
         return -1
   elif a<0 and c>=0:
      if c*x-b>0 or (L < D*a*a):
         return 1
      if c*x+c-b < 0 and (a*a*D <= R):
         return -1
   else:
      if c*x-b<0 and (L < D*a*a):
         return 1
      if c*x+c-b > 0 or (a*a*D >= R):
         return
   return 0   

#Safely counts floor of (a*sqrt(D)+b)/c
def safeFloor(D,a,b,c):
   if a == 0:
      return b // c
   x = int(((a*isqrt(D)+b)//c)) #if there is overflow python will probably give warning anyway
   e = checkFloor(D,a,b,c,x) #checks if floor was computed correctly
   first_x = x
   exp = 0
   last_e = 0
   rise = True
   while e != 0: #corrects approximation error
      if last_e*e < 0:
         rise = False
      if rise:
         exp = max(exp*2,1)
      else:
         exp = max(exp // 2,1)      
      x+= -e*exp
      last_e = e
      e = checkFloor(D,a,b,c,x)
   if first_x != x: 
      print("Approximation error was fixed. The problem occured for x = floor((a*sqrt(D)+b)/c), the original result was D = {0}, a={1}, b={2},c={3}, x ={4}, the correct one x = {5}.".format(D,a,b,c,first_x,x))
   return x

#generate all k-tuples (actually lists) with sum n, meaning (a1,a2,...,ak), a1 >= a2 >= a3... such that a1+a2+...+ak = n, such that a1<=l
def combinations(n,k,l):
   if l < (n // k):
      return
   if k == 1:
      if n>l:
         return
      else:
         yield [n]
         return
   for i in range(0,min(l+1,n+1)):
      for comb in combinations(n-i,k-1,i):
         comb.insert(0,i)
         yield comb

def take_one_from_each(elements, indices,k):
   if k == len(indices) -1:
      for item in elements[indices[k]]:
         yield [item]
      return
   
   for subseq in take_one_from_each(elements, indices,k+1):
      for item in elements[indices[k]]:
         subseq.append(item)
         yield subseq
         subseq.pop()

#creates all sums of two elements -- one from elements1, one from elements2; elements are quadruples
def create_sums_of_2(elements1,elements2,max_ind):      
   res = [set() for x in range(max_ind+1)]
   for i in range(min(max_ind+1,len(elements1))):
      for j in range(min(max_ind+1-i,len(elements2))):
         for item1 in elements1[i]:
            for item2 in elements2[j]:
               res[i+j].add((item1[0]+item2[0],item1[1]+item2[1],item1[2]+item2[2],item1[3]+item2[3]))
   return res



## Basic abstract class for biquadratic fields, containing most of the function
## Specific classes for fields of diferent types and type specific functions (trace, square,..) are defined later
class BiquadField(ABC):

   def __init__(self, p,q):
      self.p = p
      self.q = q
      self.r = p*q // (gcd(p,q)*gcd(p,q))
      self.r0 = gcd(p,q)
      self.p0 = gcd(self.r,q)
      self.q0 = gcd(p,self.r)
      self.elem_cap = -1
      self.elements = []

   #Creates the field given squarefree p, q. 
   #Chooses correct type of field and returns instance of that type.
   @staticmethod
   def createField(p,q):
      if not isSquareFree(p) or not isSquareFree(p):
         return None
      r = p*q // (gcd(p,q)*gcd(p,q))
      t = [p,q,r]
      s = [p % 4, q%4, r %4]
      if 2 in s:
         if 3 in s:
            return BiquadFieldI(t[s.index(2)],t[s.index(3)])
         else:
            return BiquadFieldII_III(t[s.index(2)],t[s.index(1)])
      elif 3 in s:
         return BiquadFieldII_III(t[s.index(3)],t[s.index(1)])
      elif gcd(p,q) % 4 == 1:
         return BiquadFieldIVa(p,q)
      else:
         return BiquadFieldIVb(p,q)
      
   #computes trace
   @abstractmethod
   def trace(self,tple):
      pass   

   #computes trace of square
   @abstractmethod
   def square_trace(self,tple):
      pass
   
   #computes the tuple of square
   @abstractmethod
   def tuple_of_square(self,tple):
      pass

   #computes the necessary bounds for a,b,c,d needed for square of (a,b,c,d) to have small trace
   @abstractmethod
   def coeff_bounds(self, trace):
      pass

   #computes all elements with trace smaller than max_trace
   def elements_with_small_trace(self,max_trace): 
      elements = [set() for i in range(max_trace+1)]
      max_a, max_b, max_c, max_d = self.coeff_bounds(max_trace)
      for a in range(max_a+1):
         for b in range(-max_b,max_b+1):
            for c in range(-max_c,max_c+1):
               for d in range(-max_d,max_d+1):
                  trace = self.square_trace((a,b,c,d))
                  if trace <= max_trace:
                     elements[trace].add(self.tuple_of_square((a,b,c,d)))
      return elements;

   #finds (only one) square root of tuple
   def square_root(self,tple):
      trace_cap = self.trace(tple)
      max_a, max_b, max_c, max_d = self.coeff_bounds(trace_cap)
      for a in range(max_a+1):
         for b in range(-max_b,max_b+1):
            for c in range(-max_c,max_c+1):
               for d in range(-max_d,max_d+1):
                  if self.tuple_of_square((a,b,c,d)) == tple:
                     return (a,b,c,d)
      return None
   #checks if all elements with small trace are computed 
   def update_elements(self,trace):
      if self.elem_cap < trace:
         self.elements = self.elements_with_small_trace(trace)
         self.elem_cap = trace   


     
   #Computes number of elements of lengths min_k - max_k up to trace max_trace
   #If return_squares, returns the computed squares as well
   def lengths_up_to_trace(self,max_trace, max_k = 7, return_squares = False,print_result = True):
      self.update_elements(max_trace)
      squares = [[set() for j in range(max_trace+1)] for i in range(max_k+1)]
      squares[1] = self.elements
      i = 2
      for i in range(2,max_k+1):
         #print("Pocitam pro i = {0}".format(i));
         sums =  create_sums_of_2(squares[i-1],squares[1],max_trace)
         for j in range(max_trace+1):
            for item in sums[j]:
               ok = True
               for k in range(1,i):
                     if item in squares[k][j]:
                        ok = False
                        break         
               if ok:
                  squares[i][j].add(item)
      s = [sum([len(squares[i][j]) for j in range(max_trace+1)]) for i in range(0,max_k+1)]
      for i in range(max_k,7):
         s.append("?")
      if print_result:
         print("Up to trace {0}, (p,q,r) = ({8},{9},{10}), number of elements with length 1-7: {1}, {2}, {3}, {4}, {5}, {6}, {7}. Total {11}.".format(max_trace,s[1],s[2],s[3],s[4],s[5],s[6],s[7],self.p,self.q, self.r, sum(s[0:max_k+1]))) 
         sys.stdout.flush()
      if not return_squares:
         return s
      else:
         return s, squares
   
   #Computes if element is sum of k squares
   def is_sum_of_k_squares(self, k, tple):
      trace = self.trace(tple)
      self.update_elements(trace)
      #Corner cases
      if k < 1:
         return False
      if k == 1:
         if tple in self.elements[trace]:
            return True
         else:
            return False

      #k = k1+k2
      k1 =  k // 2
      k2 = k-k1
      s, squares = self.lengths_up_to_trace(trace,max_k = max(k1,k2),return_squares = True,print_result = False)
      
      if k1 > k2:
         k1,k2 = k2,k1
      #creates sets of elements of length <=k1, resp. k2
      allk1 = [set() for j in range(trace+1)]
      allk2 = [set() for j in range(trace+1)]
      for i in range(1,k1+1):
         for j in range(trace+1):
            allk1[j] = allk1[j].union(squares[i][j])
      if k1 == k2:
         allk2 = allk1
      else:
         for j in range(trace+1):
            allk2[j] = allk1[j].union(squares[k2][j])
      #for elements with length<=k1 find if there exists complementary element of length <=k2
      for i in range(trace+1):
         for item in allk1[i]:
            x = (tple[0]-item[0],tple[1]-item[1],tple[2]-item[2],tple[3]-item[3])
            if x in allk2[trace-i]:
               return True
      return False

####### The following functions deal with elements of concrete trace by computing different #######
#### combinations how to get to it. They are therefore quite slow compared to the ones above. #####

   #for tuple finds all possible sums of squares forming it, k>0
   def find_squares_forming_sum_by_trace(self,k,tple):
      trace = self.trace(tple)      
      self.update_elements(trace)
      sums = []
      for comb in combinations(trace,k,trace):
         for squares in take_one_from_each(self.elements,comb,0):
            final_sum = tuple(map(sum,zip(*squares)))
            if final_sum == tple:
               if squares not in sums:
                  sums.append(list(squares)) 
      return sums

   #computes all elements with given trace which are sums of k squares
   def sums_of_k_squares_by_trace(self,trace,k,length=True):
      self.update_elements(trace)
      sums = {};
      for comb in combinations(trace,k,trace):
         for squares in take_one_from_each(self.elements,comb,0):
            final_sum = tuple(map(sum,zip(*squares)))
            if final_sum not in sums:
               sums[final_sum] = [list(squares)]
            else:
               if list(squares) not in sums[final_sum]:
                  sums[final_sum].append(list(squares))
            if final_sum == (-1,-1,-1,-1):
               sums.add(squares)
      if length:
         return len(sums)
      else:
         return len(sums),sums

   #for a given trace compares number of elements which are sums of min_k-max_k squares
   def compare_by_trace(self,trace, min_k = 1, max_k = 7):
      self.update_elements(trace)
      s = []
      for i in range(1,8):
         if min_k<= i and i <=max_k:
            s.append(self.sums_of_k_squares_by_trace(trace,i))
         else:
            s.append("?")
      print("For trace {0}, p = {8}, q = {9}, r = {10} there are this many elements: {1}, {2}, {3}, {4}, {5}, {6}, {7}.".format(trace,s[0],s[1],s[2],s[3],s[4],s[5],s[6],self.p,self.q, self.r))
      sys.stdout.flush()
      return s

   #For a given trace and k finds elements which are sums of k squares but not sums of k-1 squares, returns tuples of square roots
   def find_culprits_by_trace(self,trace,k, known_roots = {}):
      s_less = self.sums_of_k_squares_by_trace(trace,k-1,False)
      s_k = self.sums_of_k_squares_by_trace(trace,k,False)
      result = []
      if s_less[0] == s_k[0]:
         return result
      for item in s_k[1]:
         if item not in s_less[1]:
            int_sums = s_k[1][item]
            desquared_sums = [item]
            for sum in int_sums:
               new_sum = []
               for square in sum:
                  if square in known_roots:
                     root = known_roots[square] 
                  else:
                     root = self.square_root(square)
                     known_roots[square] = root
                  new_sum.append(root)
               desquared_sums.append(new_sum)
            result.append(desquared_sums)
      return result

   
### IMPLEMENTATION OF TYPE OF FIELD SPECIFIC FUNCTIONS (as well as classes for these fields ###

class BiquadFieldI(BiquadField):

   def trace(self,tple):
      return 4*tple[0]
   
   def square_trace(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      return 4*a*a+self.p*(2*b+d)*(2*b+d)+4*self.q*c*c+self.r*d*d


   def coeff_bounds(self,trace):
      max_d = safeFloor((trace//self.r)+1,1,0,1)
      max_c = safeFloor((trace//self.q) +1,1,0,2)
      max_b = safeFloor((trace//self.p) +1,1,max_d,2)
      max_a = safeFloor(trace,1,0,2)
      return max_a,max_b,max_c,max_d

   def tuple_of_square(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      D = ((2*b+d)*(2*c)*self.r0+2*a*d)
      B = ((2*c)*d*self.p0+(2*a)*(2*b+d)-D)//2
      C = ((2*b+d)*d*self.q0+(2*a)*(2*c))//2
      A = self.square_trace(tple)//4
      return (A,B,C,D)
   

class BiquadFieldII_III(BiquadField):

   def trace(self,tple):
      return 4*tple[0]+2*tple[2]
   
   def square_trace(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      return (2*a+c)*(2*a+c)+self.p*(2*b+d)*(2*b+d)+self.q*c*c+self.r*d*d


   def coeff_bounds(self,trace):
      max_d = safeFloor(((trace+1)//self.r)+1,1,0,1)
      max_c = safeFloor(((trace+1)//self.q) +1,1,0,1)
      max_b = safeFloor(((trace+1)//self.p) +1,1,max_d,2)
      max_a = safeFloor((trace+1),1,max_c,2)
      return max_a,max_b,max_c,max_d

   def tuple_of_square(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      D = ((2*b+d)*(c)*self.r0+(2*a+c)*d)
      B = (c*d*self.p0+(2*a+c)*(2*b+d)-D)//2
      C = ((2*b+d)*d*self.q0+(2*a+c)*(c))
      A = (self.square_trace(tple)-2*C)//4
      return (A,B,C,D)
   
class BiquadFieldIVa(BiquadField):

   def trace(self,tple):
      return 4*tple[0]+2*tple[1]+2*tple[2]+tple[3]
   
   def square_trace(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      return (2*a+b+c)*(2*a+b+c)+self.p*b*b+self.q*c*c+d*(2*a+(1+self.p)*b+(1+self.q)*c+(d+self.p*d+self.q*d+self.r*d)//4)

   #Could be perhaps best seen, when looking at numbers in a way x = a' + b'sqrt(p)+c'sqrt(q)+d'sqrt(r)
   def coeff_bounds(self,trace):
      max_d = safeFloor((4*(trace+1)//self.r)+1,1,0,1)
      max_b = safeFloor(((trace+1)//self.p) +1,2,max_d,2)
      max_c = safeFloor(((trace+1)//self.q) +1,2,max_d,2)
      max_a = safeFloor((trace+1),2,max_d+2*max_c+2*max_b,4)
      return max_a,max_b,max_c,max_d

   def tuple_of_square(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      D = ((2*b+d)*(2*c+d)*self.r0+(4*a+2*b+2*c+d)*d)//2
      B = ((2*c+d)*d*self.p0+(4*a+2*b+2*c+d)*(2*b+d)-2*D)//4
      C = ((2*b+d)*d*self.q0+(4*a+2*b+2*c+d)*(2*c+d)-2*D)//4
      A = (self.square_trace(tple)-2*B-2*C-D)//4
      return (A,B,C,D)

class BiquadFieldIVb(BiquadField):

   def trace(self,tple):
      return 4*tple[0]+2*tple[1]+2*tple[2]+tple[3]
   
   def square_trace(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      return (2*a+b+c)*(2*a+b+c)+self.p*b*b+self.q*c*c+d*(2*a+(1-self.p)*b+(1+self.q)*c+(d+self.p*d+self.q*d+self.r*d)//4)

   def coeff_bounds(self,trace):
      max_d = safeFloor((4*(trace+1)//self.r)+1,1,0,1)
      max_b = safeFloor(((trace+1)//self.p) +1,2,max_d,2)
      max_c = safeFloor(((trace+1)//self.q) +1,2,max_d,2)
      max_a = safeFloor((trace+1),2,max_d+2*max_c+2*max_b,4)
      return max_a,max_b,max_c,max_d

   def tuple_of_square(self,tple):
      a = tple[0];b=tple[1];c=tple[2];d=tple[3];
      D = ((2*b-d)*(2*c+d)*self.r0+(4*a+2*b+2*c+d)*d)//2
      B = ((2*c+d)*d*self.p0+(4*a+2*b+2*c+d)*(2*b-d)+2*D)//4
      C = ((2*b-d)*d*self.q0+(4*a+2*b+2*c+d)*(2*c+d)-2*D)//4
      A = (self.square_trace(tple)-2*B-2*C-D)//4
      return (A,B,C,D)


##############################################################################################

#Example piece of code and usage:
def main(args):
   # The code below works with field Q(sqrt(13),sqrt(17)) (or any other given by input arguments)
   # For each trace up to max_trace computes the number of elements with that trace which can be represented as a sum of 1,..,7 squares.
   # Unused version in the end uses multiprocessing to parallelize the computations
   # This is certainly NOT the fastest ways for predictiong Pythagoras number (for that see functions in predict_pythagoras.py),
   # but it can certainly give some good insights and offers time-memory tradeoff
   
   if len(args)<2:
      p = 13
      q = 17
   else:
      p = int(args[1])
      q = int(args[2])
      
   K = BiquadField.createField(p,q)
   max_trace = 60
   K.update_elements(max_trace)

   if True:
      for i in range(0,max_trace):
         K.compare_by_trace(i);
   else:
      #multiprocessing implementation of the code above
      with multiprocessing.Pool(processes=8) as pool:
         pool.starmap(K.compare, zip(range(0,max_trace,1)))
   
   

if __name__ == "__main__":
   main(sys.argv)
   

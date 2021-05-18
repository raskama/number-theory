import numpy as np
import json
import sympy #for primes generator
from sympy.ntheory.residue_ntheory import sqrt_mod_iter
from sympy.solvers.diophantine.diophantine import square_factor
import sys
from itertools import repeat
import multiprocessing
from multiprocessing import freeze_support
from abc import ABC, abstractmethod

def isSquareFree(n):
	if square_factor(n)>1:
		return False
	else:
		return True

#sqrt working even for high values, however for high values gives floor of sqrt (but whatever, thats probably what i want and it gets checked later either way)
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


#Checks if floor of (a*sqrt(D)+b)/c is x, returns -1 if it is less, 0 if correct, +1 if more
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
	x = int(((a*isqrt(D)+b)//c)) #ono to pro moc velká D asi zařve samo
	e = checkFloor(D,a,b,c,x)
	first_x = x
	exp = 0
	last_e = 0
	rise = True
	while e != 0:
		#print("Nove x = {0}, err = {1}, exp = {2}".format(x,e,exp))	
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
		print("Opravena zaukrouhlovaci chyba. Skapalo to pro x = floor((a*sqrt(D)+b)/c), kde vyslo D = {0}, a={1}, b={2},c={3}, x ={4}, spravne x = {5}.".format(D,a,b,c,first_x,x))
	return x



def my_gcd(a,b):
	A = a
	B = b
	if a < 0:
		a = -a
	if b < 0:
		b = -b
	if a < b:
		c = b
		b = a
		a = c
	while b>0:
		c = a % b
		a = b
		b = c

	return a


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

def create_sums_of_2(elements1,elements2,max_ind):		
	res = [set() for x in range(max_ind+1)]
	for i in range(max_ind+1):
		for j in range(max_ind+1-i):
			for item1 in elements1[i]:
				for item2 in elements2[j]:
					res[i+j].add((item1[0]+item2[0],item1[1]+item2[1],item1[2]+item2[2],item1[3]+item2[3]))
	return res



def do_analysis():
	primes = list(sympy.sieve.primerange(0, 2000))
	good  = {}
	known = {}
	### For q = 5
	known[11] = [140,148,152,160,164,168,172,176,180,184]
	known[19] = [224,228,232,236,244,248,252,256,260,264]
	known[23] = [260,264,268,276,280,284,288,292,296,300]
	known[31] = [348,352,356,360,364,368,372,376,380,384]
	### For p = 2
	#known[13] = [112,144]
	#known[17] = [120,128,136,160,168]
	for i in range(1):
		p = 19
		if p % 5 == 0:
			continue
		if not isSquareFree(p,primes):
			continue
		q = 437
		print("Pocitam prvky pro p = {0}, q = {1}".format(p,q))
		r = p*q // (my_gcd(p,q)*my_gcd(p,q))
		trace_cap = 400
		elements = elements_with_small_trace_II_III(p,q,r,trace_cap)	
		sstdout = sys.stdout
		known_roots = {}
		if True:
			t1 = (1,0,0,0); t2 = (1,0,0,0)
			t3 = (2,0,0,0); t4 = (3,0,0,0)
			t5 = (1,1,0,0)
			tup1 = tuple_of_square_II_III(p,q,r,t1)
			tup2 = tuple_of_square_II_III(p,q,r,t2)
			tup3 = tuple_of_square_II_III(p,q,r,t3)
			tup4 = tuple_of_square_II_III(p,q,r,t4)
			tup5 = tuple_of_square_II_III(p,q,r,t5)
			tups = [tup1,tup2,tup3,tup4,tup5]
			res = tuple(map(sum,zip(*tups)))
			#res = (4*p+2,2,2*p,0)
			print(res)
			if (4*res[0]+2*res[2])>trace_cap:
				elements = elements_with_small_trace_II_III(p,q,r,4*res[0]+2*res[2]+5)
			print("Predpocitano")
			#tryhard = is_sum_of_k_squares(p,q,r,elements,4*res[0]+2*res[2],4,res)
			tryhard = is_sum_of_4_squares(elements,4*res[0]+2*res[2],res)
			if tryhard:
				print("Je souctem :{ UAAAAAAUUU POZOOOOOR")
				print(find_squares_forming_sum(p,q,r,elements,4*res[0]+2*res[2],3,res))
				for myset in find_squares_forming_sum(p,q,r,elements,4*res[0]+2*res[2],3,res):
					for square in myset:
						print(square_root_II_III(p,q,r,4*square[0]+2*square[2],square))					
			else:
				print("Juch, neni souctem ctvercu.")
			continue

		found = 0;
		if p in known:
			my_it = iter(known[p])
		else:
			my_it = range(0,trace_cap,2)
		#my_it = iter([16*p+20])
		for tr in my_it:
			if found == 5:
				break
			res = find_culprits(p,q,r,elements,tr,5,known_roots)
			if res!=[]:
				found+=1
				for item in res:
					for tuples in item:
						t_item = tuple(tuples)
						if t_item in good:
							good[t_item].append(p)
						else:
							good[t_item] = [p]
				print("Nasel jsem prvky se stopou {0}".format(tr))
	for item in good:
		print("{0} cisel ( {1} ) pro {2}".format(len(good[item]),good[item],item))

class BiquadField(ABC):

	def __init__(self, p,q):
		self.p = p
		self.q = q
		self.r = p*q // (my_gcd(p,q)*my_gcd(p,q))
		self.r0 = my_gcd(p,q)
		self.p0 = my_gcd(self.r,q)
		self.q0 = my_gcd(p,self.r)
		self.elem_cap = -1
		self.elements = []


	@staticmethod
	def createField(p,q):
		if not isSquareFree(p) or not isSquareFree(p):
			print("DOPLNIT ERROR, nectvrcova")
			return None
		r = p*q // (my_gcd(p,q)*my_gcd(p,q))
		t = [p,q,r]
		s = [p % 4, q%4, r %4]
		if 2 in s:
			if 3 in s:
				return BiquadFieldI(t[s.index(2)],t[s.index(3)])
			else:
				return BiquadFieldII_III(t[s.index(2)],t[s.index(1)])
		elif 3 in s:
			return BiquadFieldII_III(t[s.index(3)],t[s.index(1)])
		elif my_gcd(p,q) % 4 == 1:
			return BiquadFieldIVa(p,q)
		else:
			print("FUUUUUUCK")
			return None
		

	@abstractmethod
	def trace(self,tple):
		pass	

	@abstractmethod
	def square_trace(self,tple):
		pass
	
	@abstractmethod
	def tuple_of_square(self,tple):
		pass

	@abstractmethod
	def coeff_bounds(self, trace):
		pass

	def elements_with_small_trace(self,trace_cap): 
		elements = [set() for i in range(trace_cap+1)]
		max_a, max_b, max_c, max_d = self.coeff_bounds(trace_cap)
		for a in range(max_a+1):
			for b in range(-max_b,max_b+1):
				for c in range(-max_c,max_c+1):
					for d in range(-max_d,max_d+1):
						trace = self.square_trace((a,b,c,d))
						if trace <= trace_cap:
							elements[trace].add(self.tuple_of_square((a,b,c,d)))
		return elements;

	def square_root(self,tple):
		trace_cap = self.trace(tple)
		max_a, max_b, max_c, max_d = self.coeff_bounds(trace_cap)
		for a in range(max_a+1):
			for b in range(-max_b,max_b+1):
				for c in range(-max_c,max_c+1):
					for d in range(-max_d,max_d+1):
						if self.tuple_of_square((a,b,c,d)) == tple:
							return (a,b,c,d)

	def update_elements(self,trace):
		if self.elem_cap < trace:
			self.elements = self.elements_with_small_trace(trace)
			self.elem_cap = trace	


	def compare_together(self,max_trace, min_k = 1, max_k = 7):
		self.update_elements(max_trace)
		squares = [[] for i in range(max_k+1)]
		squares[1] = self.elements
		i = 2
		for i in range(2,max_k+1):
			squares[i] =  create_sums_of_2(squares[i-1],squares[1],max_trace)
		
		s = [sum([len(squares[i][j]) for j in range(max_trace+1)]) for i in range(1,max_k+1)]
		for i in range(max_k,7):
			s.append("?")
		print("Up to trace {0}, p = {8}, q = {9}, r = {10} there are this many elements: {1}, {2}, {3}, {4}, {5}, {6}, {7}.".format(max_trace,s[0],s[1],s[2],s[3],s[4],s[5],s[6],self.p,self.q, self.r))
		sys.stdout.flush()
		return s
		
		

	def is_sum_of_4_squares(self,tple):
		trace = self.trace(tple)
		self.update_elements(trace)
		squares = create_sums_of_2(self.elements,self.elements,trace)
		for i in range(trace+1):
			for item in squares[i]:
				x = (tple[0]-item[0],tple[1]-item[1],tple[2]-item[2],tple[3]-item[3])
				if x in squares[trace-i]:
					return True
		return False

	def is_sum_of_k_squares(self,k,tple):
		trace = self.trace(tple)
		self.update_elements(trace)
		if k == 4:
			return self.is_sum_of_4_squares(tple)
		for comb in combinations(trace,k,trace):
			for squares in take_one_from_each(self.elements,comb,0):
				final_sum = tuple(map(sum,zip(*squares)))
				if final_sum == tple:
					return True
		return False

	def find_squares_forming_sum(self,k,tple):
		trace = self.trace(tple)		
		self.update_elements(trace)
		sums = []
		for comb in combinations(trace,k,trace):
			for squares in take_one_from_each(elements,comb,0):
				final_sum = tuple(map(sum,zip(*squares)))
				if final_sum == tple:
					if squares not in sums:
						sums.append(list(squares)) #wtf, proc tu ten list musi byt, asi referencovane promenne?
		return sums

	def sums_of_k_squares(self,trace,k,length=True):
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

	def compare(self,trace, min_k = 1, max_k = 7):
		self.update_elements(trace)
		s = []
		for i in range(1,8):
			if min_k<= i and i <=max_k:
				s.append(self.sums_of_k_squares(trace,i))
			else:
				s.append("?")
		print("For trace {0}, p = {8}, q = {9}, r = {10} there are this many elements: {1}, {2}, {3}, {4}, {5}, {6}, {7}.".format(trace,s[0],s[1],s[2],s[3],s[4],s[5],s[6],self.p,self.q, self.r))
		sys.stdout.flush()
		return s

	def find_culprits(self,trace,k, known_roots = {}):
		s_less = self.sums_of_k_squares(trace,k-1,False)
		s_k = self.sums_of_k_squares(trace,k,False)
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

	def coeff_bounds(self,trace):
		max_d = safeFloor((4*(trace+1)//self.r)+1,1,0,1)
		max_b = safeFloor(((trace+1)//self.p) +1,1,max_d,2)
		max_c = safeFloor(((trace+1)//self.q) +1,1,max_d,2)
		max_a = safeFloor((trace+1),2,max_d+2*max_c+2*max_b,4)
		return max_a,max_b,max_c,max_d

	def tuple_of_square(self,tple):
		a = tple[0];b=tple[1];c=tple[2];d=tple[3];
		D = ((2*b+d)*(2*c+d)*self.r0+(4*a+2*b+2*c+d)*d)//2
		B = ((2*c+d)*d*self.p0+(4*a+2*b+2*c+d)*(2*b+d)-2*D)//4
		C = ((2*b+d)*d*self.q0+(4*a+2*b+2*c+d)*(2*c+d)-2*D)//4
		A = (self.square_trace(tple)-2*B-2*C-D)//4
		return (A,B,C,D)

def main(args):
	if len(args)<2:
		p = 13
		q = 5
	else:
		p = int(args[1])
		q = int(args[2])
		
	K = BiquadField.createField(p,q)
	trace_cap = 400
	K.update_elements(trace_cap)
	with multiprocessing.Pool(processes=8) as pool:
		pool.starmap(K.compare, zip(range(0,trace_cap,1)))

if __name__ == "__main__":
	main(sys.argv)
	

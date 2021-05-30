#!/usr/bin/python3

import unittest
import biquad.biquadfunc
from abc import ABC

class RationalIntegers(biquad.biquadfunc.BiquadField):

   def __init__(self):
      self.elem_cap = -1
      self.elements = []
   
   def trace(self,tple):
      return tple[0]

   def square_trace(self, tple):
      return tple[0]*tple[0]

   def coeff_bounds(self,trace):
      max_d = 0; max_c = 0; max_b = 0;
      max_a = biquad.biquadfunc.safeFloor(trace,1,0,1)
      return max_a,max_b,max_c,max_d

   def tuple_of_square(self, tple):
      return (tple[0]*tple[0],0,0,0)


class Test_Rational_Integers(unittest.TestCase):

   def setUp(self):
      self.F = RationalIntegers()

   def test_elements_with_small_trace(self):
      res = [[] for _ in range(7)]; true_res = [[] for _ in range(7)];
      
      res[0] = self.F.elements_with_small_trace(1)
      res[1] = self.F.elements_with_small_trace(10)
      res[2] = self.F.elements_with_small_trace(15)
      res[3] = self.F.elements_with_small_trace(49)
      res[4] = self.F.elements_with_small_trace(0)
      res[5] = self.F.elements_with_small_trace(55)
      res[6] = self.F.elements_with_small_trace(83)

      true_res[0] = [{(0,0,0,0)},{(1,0,0,0)}]
      true_res[1] = [set() for _ in range(11)]
      true_res[2] = [set() for _ in range(16)]
      true_res[3] = [set() for _ in range(50)]
      true_res[4] = [set() for _ in range(1)]
      true_res[5] = [set() for _ in range(56)]
      true_res[6] = [set() for _ in range(84)]

      true_res[1][0].add((0,0,0,0))
      true_res[2][0].add((0,0,0,0))
      true_res[3][0].add((0,0,0,0))
      true_res[4][0].add((0,0,0,0))
      true_res[5][0].add((0,0,0,0))
      true_res[6][0].add((0,0,0,0)) 
      true_res[1][1].add((1,0,0,0))
      true_res[2][1].add((1,0,0,0))
      true_res[3][1].add((1,0,0,0))
      true_res[5][1].add((1,0,0,0))
      true_res[6][1].add((1,0,0,0))
      true_res[1][4].add((4,0,0,0))
      true_res[2][4].add((4,0,0,0))
      true_res[3][4].add((4,0,0,0))
      true_res[5][4].add((4,0,0,0))
      true_res[6][4].add((4,0,0,0))
      true_res[1][9].add((9,0,0,0))
      true_res[2][9].add((9,0,0,0))
      true_res[3][9].add((9,0,0,0))
      true_res[5][9].add((9,0,0,0))
      true_res[6][9].add((9,0,0,0))
      true_res[3][16].add((16,0,0,0))
      true_res[5][16].add((16,0,0,0))
      true_res[6][16].add((16,0,0,0))
      true_res[3][25].add((25,0,0,0))
      true_res[5][25].add((25,0,0,0))
      true_res[6][25].add((25,0,0,0))
      true_res[3][36].add((36,0,0,0))
      true_res[5][36].add((36,0,0,0))
      true_res[6][36].add((36,0,0,0))
      true_res[3][49].add((49,0,0,0))
      true_res[5][49].add((49,0,0,0))
      true_res[6][49].add((49,0,0,0))
      true_res[6][64].add((64,0,0,0))
      true_res[6][81].add((81,0,0,0))

      assert(res[0] == true_res[0])
      assert(res[1] == true_res[1])
      assert(res[2] == true_res[2])
      assert(res[3] == true_res[3])
      assert(res[4] == true_res[4])
      assert(res[5] == true_res[5])
      assert(res[6] == true_res[6])


   def test_square_root(self):
      assert(self.F.square_root((0,0,0,0)) == (0,0,0,0))
      assert(self.F.square_root((1,0,0,0)) == (1,0,0,0))
      assert(self.F.square_root((49,0,0,0)) == (7,0,0,0))
      assert(self.F.square_root((123*123,0,0,0)) == (123,0,0,0))

      self.assertRaises(ValueError,self.F.square_root,(-1,0,0,0))
      assert(self.F.square_root((8,0,0,0)) == None)
      assert(self.F.square_root((1,1,0,0)) == None)

   def test_update_elements(self):
      self.F.update_elements(1)
      assert(self.F.elements == self.F.elements_with_small_trace(1))

      self.F.update_elements(16)
      assert(self.F.elements == self.F.elements_with_small_trace(16))

      self.F.update_elements(7)
      assert(self.F.elements == self.F.elements_with_small_trace(16))

   def test_lengths_up_to_trace(self):
      x1 = self.F.lengths_up_to_trace(1,print_result = False)
      x5 = self.F.lengths_up_to_trace(5,print_result = False)
      x20 = self.F.lengths_up_to_trace(20,print_result = False)
      x25 = self.F.lengths_up_to_trace(25,print_result = False)
      x25_3 = self.F.lengths_up_to_trace(25,print_result = False, max_k = 3)
      x25_2 = self.F.lengths_up_to_trace(25,print_result = False, max_k = 2)
      x25_1 = self.F.lengths_up_to_trace(25,print_result = False, max_k = 1)

      assert(x1 == [0,2,0,0,0,0,0,0])
      assert(x5 == [0,3,2,1,0,0,0,0])
      assert(x20 == [0,5,8,6,2,0,0,0])
      assert(x25 == [0,6,8,9,3,0,0,0])
      assert(x25_3 == [0,6,8,9,"?","?","?","?"])
      assert(x25_2 == [0,6,8,"?","?","?","?","?"])
      assert(x25_1 == [0,6,"?","?","?","?","?","?"])      

   def test_is_sum_of_k_squares(self):
      assert(self.F.is_sum_of_k_squares(6,(1,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(3,(1,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(2,(1,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(1,(1,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(0,(1,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(100,(1,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(100,(0,1,0,0)) == False)

      assert(self.F.is_sum_of_k_squares(2,(5,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(3,(5,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(1,(5,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(3,(7,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(4,(7,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(3,(15,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(4,(15,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(4,(4,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(2,(22,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(2,(23,0,0,0)) == False)
      assert(self.F.is_sum_of_k_squares(2,(17,0,0,0)) == True)
      assert(self.F.is_sum_of_k_squares(3,(22,0,0,0)) == True)

   def test_find_squares_forming_sum_by_trace(self):

      res23_4 = self.F.find_squares_forming_sum_by_trace(4,(23,0,0,0))
      true23_4 = [[(1,0,0,0),(4,0,0,0),(9,0,0,0),(9,0,0,0)]]
      assert(res23_4 == true23_4)

      res23_3 = self.F.find_squares_forming_sum_by_trace(3,(23,0,0,0))
      true23_3 = []
      assert(res23_3 == true23_3)

      res23_5 = self.F.find_squares_forming_sum_by_trace(5,(23,0,0,0))
      assert([(0,0,0,0),(1,0,0,0),(4,0,0,0),(9,0,0,0),(9,0,0,0)] in res23_5)
      assert([(1,0,0,0),(1,0,0,0),(1,0,0,0),(4,0,0,0),(16,0,0,0)] in res23_5)
      assert(len(res23_5) == 2)

      res20_5 = self.F.find_squares_forming_sum_by_trace(5,(20,0,0,0))
      assert([(1,0,0,0),(1,0,0,0),(1,0,0,0),(1,0,0,0),(16,0,0,0)] in res20_5)
      assert([(0,0,0,0),(0,0,0,0),(0,0,0,0),(4,0,0,0),(16,0,0,0)] in res20_5)
      assert([(0,0,0,0),(1,0,0,0),(1,0,0,0),(9,0,0,0),(9,0,0,0)] in res20_5)
      assert([(4,0,0,0),(4,0,0,0),(4,0,0,0),(4,0,0,0),(4,0,0,0)] in res20_5)
      assert(len(res20_5) == 4)

      res25_2 = self.F.find_squares_forming_sum_by_trace(2,(25,0,0,0))
      assert([(0,0,0,0),(25,0,0,0)] in res25_2)
      assert([(9,0,0,0),(16,0,0,0)] in res25_2)
      assert(len(res25_2) == 2)

      res9_3 = self.F.find_squares_forming_sum_by_trace(3,(9,0,0,0))
      assert([(0,0,0,0),(0,0,0,0),(9,0,0,0)] in res9_3)
      assert([(1,0,0,0),(4,0,0,0),(4,0,0,0)] in res9_3)
      assert(len(res9_3) == 2)

      res9_2 = self.F.find_squares_forming_sum_by_trace(2,(9,0,0,0))
      true9_2 = [[(0,0,0,0),(9,0,0,0)]]
      assert(res9_2 == true9_2)


      res8_2 = self.F.find_squares_forming_sum_by_trace(2,(8,0,0,0))
      true8_2 = [[(4,0,0,0),(4,0,0,0)]]
      assert(res8_2 == true8_2)

      res1_1 = self.F.find_squares_forming_sum_by_trace(1,(1,0,0,0))
      true1_1 = [[(1,0,0,0)]]
      assert(res1_1 == true1_1)

      res1_2 = self.F.find_squares_forming_sum_by_trace(2,(1,0,0,0))
      true1_2 = [[(0,0,0,0),(1,0,0,0)]]
      assert(res1_2 == true1_2)


   def test_sums_of_k_squares_by_trace(self):
      assert(self.F.sums_of_k_squares_by_trace(15,1) == 0)
      assert(self.F.sums_of_k_squares_by_trace(15,2) == 0)
      assert(self.F.sums_of_k_squares_by_trace(15,3) == 0)
      assert(self.F.sums_of_k_squares_by_trace(15,4) == 1)
      assert(self.F.sums_of_k_squares_by_trace(15,5) == 1)
 
      assert(self.F.sums_of_k_squares_by_trace(1,1) == 1)
      assert(self.F.sums_of_k_squares_by_trace(1,2) == 1)
      assert(self.F.sums_of_k_squares_by_trace(1,3) == 1)
      assert(self.F.sums_of_k_squares_by_trace(1,4) == 1)
      assert(self.F.sums_of_k_squares_by_trace(1,5) == 1)

      assert(self.F.sums_of_k_squares_by_trace(2,1) == 0)
      assert(self.F.sums_of_k_squares_by_trace(2,2) == 1)
      assert(self.F.sums_of_k_squares_by_trace(2,3) == 1)
      assert(self.F.sums_of_k_squares_by_trace(2,4) == 1)
      assert(self.F.sums_of_k_squares_by_trace(2,5) == 1)
 
      assert(self.F.sums_of_k_squares_by_trace(3,1) == 0)
      assert(self.F.sums_of_k_squares_by_trace(3,2) == 0)
      assert(self.F.sums_of_k_squares_by_trace(3,3) == 1)
      assert(self.F.sums_of_k_squares_by_trace(3,4) == 1)
      assert(self.F.sums_of_k_squares_by_trace(3,5) == 1)
 
      assert(self.F.sums_of_k_squares_by_trace(0,1) == 1)
      assert(self.F.sums_of_k_squares_by_trace(0,2) == 1)
      assert(self.F.sums_of_k_squares_by_trace(0,3) == 1)
      assert(self.F.sums_of_k_squares_by_trace(0,4) == 1)
      assert(self.F.sums_of_k_squares_by_trace(0,5) == 1)

def test_find_culprits_by_trace(self):

      assert(self.F.find_culprits_by_trace(15,5) == [])
      assert(self.F.find_culprits_by_trace(15,4) == [(15,0,0,0),[(1,0,0,0),(1,0,0,0),(4,0,0,0),(9,0,0,0)]])
      assert(self.F.find_culprits_by_trace(15,3) == [])
 

if __name__ == '__main__':
    unittest.main()

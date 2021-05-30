#!/usr/bin/python3

import unittest
from biquad.biquadfunc import BiquadField as BF
from abc import ABC

class BiquadFieldType(ABC):

   def test_coeff_bounds(self):
      l1 = 20; l2 = 40; l3 = 100; l4 = 200; l5 = 400
      l6 = 1; l7 = 10; l8 = 150; l9 = 1000; l10 = 2000
      
      big = self.F.elements_with_small_trace(3000)
      s1 = self.F.elements_with_small_trace(l1)
      s2 = self.F.elements_with_small_trace(l2)
      s3 = self.F.elements_with_small_trace(l3)
      s4 = self.F.elements_with_small_trace(l4)
      s5 = self.F.elements_with_small_trace(l5)
      s6 = self.F.elements_with_small_trace(l6)
      s7 = self.F.elements_with_small_trace(l7)
      s8 = self.F.elements_with_small_trace(l8)
      s9 = self.F.elements_with_small_trace(l9)
      s10 = self.F.elements_with_small_trace(l10)

      assert(sum([len(big[i]) for i in range(0,l1+1)]) == sum([len(s1[i]) for i in range(0,l1+1)]))
      assert(sum([len(big[i]) for i in range(0,l2+1)]) == sum([len(s2[i]) for i in range(0,l2+1)]))
      assert(sum([len(big[i]) for i in range(0,l3+1)]) == sum([len(s3[i]) for i in range(0,l3+1)]))
      assert(sum([len(big[i]) for i in range(0,l4+1)]) == sum([len(s4[i]) for i in range(0,l4+1)]))
      assert(sum([len(big[i]) for i in range(0,l5+1)]) == sum([len(s5[i]) for i in range(0,l5+1)]))
      assert(sum([len(big[i]) for i in range(0,l6+1)]) == sum([len(s6[i]) for i in range(0,l6+1)]))
      assert(sum([len(big[i]) for i in range(0,l7+1)]) == sum([len(s7[i]) for i in range(0,l7+1)]))
      assert(sum([len(big[i]) for i in range(0,l8+1)]) == sum([len(s8[i]) for i in range(0,l8+1)]))
      assert(sum([len(big[i]) for i in range(0,l9+1)]) == sum([len(s9[i]) for i in range(0,l9+1)]))
      assert(sum([len(big[i]) for i in range(0,l10+1)]) == sum([len(s10[i]) for i in range(0,l10+1)]))




class BiquadFieldTypeI_2_3(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(2,3)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 0)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == 0)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 4)
      assert(self.F.trace(self.t10) == 16)
      assert(self.F.trace(self.t11) == 172)
      assert(self.F.trace(self.t12) == 96)
      assert(self.F.trace(self.t13) == 388)
      assert(self.F.trace(self.t14) == 4)
      assert(self.F.trace(self.t15) == -16)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 8)
      assert(self.F.square_trace(self.t3) == 12)
      assert(self.F.square_trace(self.t4) == 8)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 40)
      assert(self.F.square_trace(self.t10) == 54*4)
      assert(self.F.square_trace(self.t11) == 3810*4)
      assert(self.F.square_trace(self.t12) == 1124*4)
      assert(self.F.square_trace(self.t13) == 10650*4)
      assert(self.F.square_trace(self.t14) == 10*4)
      assert(self.F.square_trace(self.t15) == 297*4)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (2,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (3,0,0,0))
      assert(self.F.tuple_of_square(self.t4) == (2,0,1,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (10,2,5,8))
      assert(self.F.tuple_of_square(self.t10) == (54,16,23,36))
      assert(self.F.tuple_of_square(self.t11) == (3810,1170,1995,1256*2))
      assert(self.F.tuple_of_square(self.t12) == (1124,632,224,320))
      assert(self.F.tuple_of_square(self.t13) == (10650,5140,-510,-4204))
      assert(self.F.tuple_of_square(self.t14) == (10,-2,5,-8))
      assert(self.F.tuple_of_square(self.t15) == (297,-170,-11,2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeI_6_103(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(6,103)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 0)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == 0)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 4)
      assert(self.F.trace(self.t10) == 16)
      assert(self.F.trace(self.t11) == 172)
      assert(self.F.trace(self.t12) == 96)
      assert(self.F.trace(self.t13) == 388)
      assert(self.F.trace(self.t14) == 4)
      assert(self.F.trace(self.t15) == -16)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 4*6)
      assert(self.F.square_trace(self.t3) == 4*103)
      assert(self.F.square_trace(self.t4) == 156*4)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 272*4)
      assert(self.F.square_trace(self.t10) == 656*4)
      assert(self.F.square_trace(self.t11) == 60248*4)
      assert(self.F.square_trace(self.t12) == 4996*4)
      assert(self.F.square_trace(self.t13) == 85826*4)
      assert(self.F.square_trace(self.t14) == 272*4)
      assert(self.F.square_trace(self.t15) == 17615*4)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (6,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (103,0,0,0))
      assert(self.F.tuple_of_square(self.t4) == (156,0,3,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (272,102,11,4*2))
      assert(self.F.tuple_of_square(self.t10) == (656,216,37,18*2))
      assert(self.F.tuple_of_square(self.t11) == (60248,23270,3061,1256*2))
      assert(self.F.tuple_of_square(self.t12) == (4996,1432,480,320))
      assert(self.F.tuple_of_square(self.t13) == (85826,2940,-1918,-4204))
      assert(self.F.tuple_of_square(self.t14) == (272,-102,11,-8))
      assert(self.F.tuple_of_square(self.t15) == (17615,-6470,79,2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeI_62_67(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(62,67)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 0)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == 0)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 4)
      assert(self.F.trace(self.t10) == 16)
      assert(self.F.trace(self.t11) == 172)
      assert(self.F.trace(self.t12) == 96)
      assert(self.F.trace(self.t13) == 388)
      assert(self.F.trace(self.t14) == 4)
      assert(self.F.trace(self.t15) == -16)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 62*4)
      assert(self.F.square_trace(self.t3) == 67*4)
      assert(self.F.square_trace(self.t4) == 1054*4)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 1246*4)
      assert(self.F.square_trace(self.t10) == 2082*4)
      assert(self.F.square_trace(self.t11) == 222774*4)
      assert(self.F.square_trace(self.t12) == 33332*4)
      assert(self.F.square_trace(self.t13) == 527982*4)
      assert(self.F.square_trace(self.t14) == 1246*4)
      assert(self.F.square_trace(self.t15) == 87805*4)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (62,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (67,0,0,0))
      assert(self.F.tuple_of_square(self.t4) == (1054,0,31,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (1246,66,95,8))
      assert(self.F.tuple_of_square(self.t10) == (2082,144,233,36))
      assert(self.F.tuple_of_square(self.t11) == (222774,15314,17985,1256*2))
      assert(self.F.tuple_of_square(self.t12) == (33332,1144,4064,320))
      assert(self.F.tuple_of_square(self.t13) == (527982,3732,-21630,-4204))
      assert(self.F.tuple_of_square(self.t14) == (1246,-66,95,-8))
      assert(self.F.tuple_of_square(self.t15) == (87805,-4202,1339,2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()

class BiquadFieldTypeI_30_39(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(30,39)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 0)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == 0)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 4)
      assert(self.F.trace(self.t10) == 16)
      assert(self.F.trace(self.t11) == 172)
      assert(self.F.trace(self.t12) == 96)
      assert(self.F.trace(self.t13) == 388)
      assert(self.F.trace(self.t14) == 4)
      assert(self.F.trace(self.t15) == -16)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 30*4)
      assert(self.F.square_trace(self.t3) == 39*4)
      assert(self.F.square_trace(self.t4) == 40*4)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 140*4)
      assert(self.F.square_trace(self.t10) == 572*4)
      assert(self.F.square_trace(self.t11) == 31220*4)
      assert(self.F.square_trace(self.t12) == 8932*4)
      assert(self.F.square_trace(self.t13) == 32858*4)
      assert(self.F.square_trace(self.t14) == 140*4)
      assert(self.F.square_trace(self.t15) == 4747*4)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (30,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (39,0,0,0))
      assert(self.F.tuple_of_square(self.t4) == (40,0,5,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (140,6,17,10*2))
      assert(self.F.tuple_of_square(self.t10) == (572,8,51,46*2))
      assert(self.F.tuple_of_square(self.t11) == (31220,1986,4127,2650*2))
      assert(self.F.tuple_of_square(self.t12) == (8932,584,736,288*2))
      assert(self.F.tuple_of_square(self.t13) == (32858,4856,-3326,-2038*2))
      assert(self.F.tuple_of_square(self.t14) == (140,-6,17,-20))
      assert(self.F.tuple_of_square(self.t15) == (4747,-730,169,-69*2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


############################################## FIELD TYPE B2,3 ################################

class BiquadFieldTypeII_2_5(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(2,5)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 6)
      assert(self.F.trace(self.t10) == 20)
      assert(self.F.trace(self.t11) == 206)
      assert(self.F.trace(self.t12) == 100)
      assert(self.F.trace(self.t13) == 390)
      assert(self.F.trace(self.t14) == 6)
      assert(self.F.trace(self.t15) == -2)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 8)
      assert(self.F.square_trace(self.t3) == 6)
      assert(self.F.square_trace(self.t4) == 12)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 42)
      assert(self.F.square_trace(self.t10) == 57*4)
      assert(self.F.square_trace(self.t11) == 8553*2)
      assert(self.F.square_trace(self.t12) == 1182*4)
      assert(self.F.square_trace(self.t13) == 22459*2)
      assert(self.F.square_trace(self.t14) == 21*2)
      assert(self.F.square_trace(self.t15) == 553*2)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (2,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (1,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (2,0,2,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (6,4,9,6))
      assert(self.F.tuple_of_square(self.t10) == (40,28,34,24))
      assert(self.F.tuple_of_square(self.t11) == (5736//2,1646,2817,1018*2))
      assert(self.F.tuple_of_square(self.t12) == (1004,688,178*2,132*2))
      assert(self.F.tuple_of_square(self.t13) == (23672//2,5194,-1213,-2129*2))
      assert(self.F.tuple_of_square(self.t14) == (6,-4,9,-6))
      assert(self.F.tuple_of_square(self.t15) == (470//2,-142,83,-26))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeIII_7_21(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(7,21)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 6)
      assert(self.F.trace(self.t10) == 20)
      assert(self.F.trace(self.t11) == 206)
      assert(self.F.trace(self.t12) == 100)
      assert(self.F.trace(self.t13) == 390)
      assert(self.F.trace(self.t14) == 6)
      assert(self.F.trace(self.t15) == -2)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 4*7)
      assert(self.F.square_trace(self.t3) == 22)
      assert(self.F.square_trace(self.t4) == 10)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 24*4)
      assert(self.F.square_trace(self.t10) == 265*2)
      assert(self.F.square_trace(self.t11) == 7238*4)
      assert(self.F.square_trace(self.t12) == 2450*4)
      assert(self.F.square_trace(self.t13) == 23333*2)
      assert(self.F.square_trace(self.t14) == 24*4)
      assert(self.F.square_trace(self.t15) == 362*4)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (7,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (5,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (2,0,1,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (21,-6,6,12*2))
      assert(self.F.tuple_of_square(self.t10) == (238//2,-16,27,54*2))
      assert(self.F.tuple_of_square(self.t11) == (6096,-666,1142*2,3109*2))
      assert(self.F.tuple_of_square(self.t12) == (2336,488,114*2,324*2))
      assert(self.F.tuple_of_square(self.t13) == (23842//2,5120,-509,-2033*2))
      assert(self.F.tuple_of_square(self.t14) == (21,6,6,-24))
      assert(self.F.tuple_of_square(self.t15) == (343,26,19*2,-118*2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeII_22_77(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(22,77)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 0)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 0)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == 0)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == 0)
      assert(self.F.trace(self.t9) == 6)
      assert(self.F.trace(self.t10) == 20)
      assert(self.F.trace(self.t11) == 206)
      assert(self.F.trace(self.t12) == 100)
      assert(self.F.trace(self.t13) == 390)
      assert(self.F.trace(self.t14) == 6)
      assert(self.F.trace(self.t15) == -2)


   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 22*4)
      assert(self.F.square_trace(self.t3) == 39*2)
      assert(self.F.square_trace(self.t4) == 9*4)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 149*2)
      assert(self.F.square_trace(self.t10) == 375*4)
      assert(self.F.square_trace(self.t11) == 36105*2)
      assert(self.F.square_trace(self.t12) == 6390*4)
      assert(self.F.square_trace(self.t13) == 33703*2)
      assert(self.F.square_trace(self.t14) == 149*2)
      assert(self.F.square_trace(self.t15) == 2729*2)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (22,0,0,0))
      assert(self.F.tuple_of_square(self.t3) == (19,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (8,0,2,0))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (70,-10,9,18*2))
      assert(self.F.tuple_of_square(self.t10) == (358,-40,34,82*2))
      assert(self.F.tuple_of_square(self.t11) == (33288//2,-1618,2817,4503*2))
      assert(self.F.tuple_of_square(self.t12) == (6212,376,178*2,452*2))
      assert(self.F.tuple_of_square(self.t13) == (34916//2,5012,-1213,-1969*2))
      assert(self.F.tuple_of_square(self.t14) == (70,10,9,-18*2))
      assert(self.F.tuple_of_square(self.t15) == (2646//2,-30,83,-188*2))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()

################################## FIELDS OF TYPE IVa ######################################

class BiquadFieldTypeIVa_5_13(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(5,13)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 6)
      assert(self.F.square_trace(self.t3) == 14)
      assert(self.F.square_trace(self.t4) == 21)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 77)
      assert(self.F.square_trace(self.t10) == 341)
      assert(self.F.square_trace(self.t11) == 27279)
      assert(self.F.square_trace(self.t12) == 1526*4)
      assert(self.F.square_trace(self.t13) == 27483*2)
      assert(self.F.square_trace(self.t14) == 33)
      assert(self.F.square_trace(self.t15) == 1441)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (1,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (3,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (3,3,1,1))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (8,12,6,9))
      assert(self.F.tuple_of_square(self.t10) == (40,96//2,54//2,31))
      assert(self.F.tuple_of_square(self.t11) == ((27279-9035-7137+2569)//4,(9035-2569)//2,(7137-2569)//2,2569))
      assert(self.F.tuple_of_square(self.t12) == (832,482*2,114*2,98*4))
      assert(self.F.tuple_of_square(self.t13) == ((27483-2403*2+1495*2-2481)//2,2403*2+2481,-1495*2+2481,-2481*2))
      assert(self.F.tuple_of_square(self.t14) == (8,-4,6,-3))
      assert(self.F.tuple_of_square(self.t15) == (1640//4,-294//2,76//2,19))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()

class BiquadFieldTypeIVa_65_85(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(65,85)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 33*2)
      assert(self.F.square_trace(self.t3) == 43*2)
      assert(self.F.square_trace(self.t4) == 93)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 413)
      assert(self.F.square_trace(self.t10) == 1565)
      assert(self.F.square_trace(self.t11) == 98847)
      assert(self.F.square_trace(self.t12) == 5810*4)
      assert(self.F.square_trace(self.t13) == 48201*2)
      assert(self.F.square_trace(self.t14) == 225)
      assert(self.F.square_trace(self.t15) == 5425)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (16,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (21,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (20,3,2,3))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (368//4,6,3,27))
      assert(self.F.tuple_of_square(self.t10) == (1416//4,18,6,101))
      assert(self.F.tuple_of_square(self.t11) == ((98847-10257-9269+6423)//4,(10257-6423)//2,(9269-6423)//2,6423))
      assert(self.F.tuple_of_square(self.t12) == (5100,370*2,114*2,226*4))
      assert(self.F.tuple_of_square(self.t13) == ((48201-2623*2+2199*2-3121)//2,2623*2+3121,-2199*2+3121,-3121*2))
      assert(self.F.tuple_of_square(self.t14) == (208//4,-2,15,-9))
      assert(self.F.tuple_of_square(self.t15) == (5484//4,-167,153,-31))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeIVa_101_13(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(101,13)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 102)
      assert(self.F.square_trace(self.t3) == 14)
      assert(self.F.square_trace(self.t4) == 357)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 605)
      assert(self.F.square_trace(self.t10) == 1829)
      assert(self.F.square_trace(self.t11) == 120351)
      assert(self.F.square_trace(self.t12) == 8918*4)
      assert(self.F.square_trace(self.t13) == 115275*2)
      assert(self.F.square_trace(self.t14) == 561)
      assert(self.F.square_trace(self.t15) == 27313)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (25,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (3,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == ((357-51-7+1)//4,3,25,1))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == ((605-165-33+9)//4,12,78,9))
      assert(self.F.tuple_of_square(self.t10) == ((1829-421-127+31)//4,48,195,31))
      assert(self.F.tuple_of_square(self.t11) == ((120351-32721-9035+2569)//4,(9035-2569)//2,(32721-2569)//2,2569))
      assert(self.F.tuple_of_square(self.t12) == (8918-1748-580+98,482*2,1650*2,98*4))
      assert(self.F.tuple_of_square(self.t13) == ((115275+2*9943-2*2403-2481)//2,2*2403+2481,-2*9943+2481,-2481*2))
      assert(self.F.tuple_of_square(self.t14) == (416//4,-4,78,-3))
      assert(self.F.tuple_of_square(self.t15) == ((27313-2255+275+19)//4,-147,2236//2,19))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()

################################## FIELDS OF TYPE IVb ######################################

class BiquadFieldTypeIVb_21_33(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(21,33)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 22)
      assert(self.F.square_trace(self.t3) == 34)
      assert(self.F.square_trace(self.t4) == 33)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 119)
      assert(self.F.square_trace(self.t10) == 539)
      assert(self.F.square_trace(self.t11) == 37911)
      assert(self.F.square_trace(self.t12) == 2054*4)
      assert(self.F.square_trace(self.t13) == 43731*2)
      assert(self.F.square_trace(self.t14) == 35)
      assert(self.F.square_trace(self.t15) == 2665)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (5,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (8,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (8,2,-1,-1))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (18,15,4,9))
      assert(self.F.tuple_of_square(self.t10) == (77,146//2,34//2,51))
      assert(self.F.tuple_of_square(self.t11) == ((37911-5213-6487-2663)//4,(5213+2663)//2,(6487-2663)//2,2663))
      assert(self.F.tuple_of_square(self.t12) == (2054-440-216-138,578*2,78*2,138*4))
      assert(self.F.tuple_of_square(self.t13) == ((43731-4614*2+2518*2+3461)//2,4614*2-3461,-2518*2+3461,-3461*2))
      assert(self.F.tuple_of_square(self.t14) == (10,-5,4,-3))
      assert(self.F.tuple_of_square(self.t15) == ((2665+293+427-129)//4,-164//2,-556//2,129))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()

class BiquadFieldTypeIVb_77_105(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(77,105)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 39*2)
      assert(self.F.square_trace(self.t3) == 53*2)
      assert(self.F.square_trace(self.t4) == 87)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 317)
      assert(self.F.square_trace(self.t10) == 1361)
      assert(self.F.square_trace(self.t11) == 84541)
      assert(self.F.square_trace(self.t12) == 4446*4)
      assert(self.F.square_trace(self.t13) == 93087*2)
      assert(self.F.square_trace(self.t14) == 89)
      assert(self.F.square_trace(self.t15) == 7263)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (19,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (26,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (22,2,-1,-3))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == (256//4,21,2,15))
      assert(self.F.tuple_of_square(self.t10) == ((1361-105-95-101)//4,103,-3,101))
      assert(self.F.tuple_of_square(self.t11) == ((84541-6435-6877-4073)//4,(6435+4073)//2,(6877-4073)//2,4073))
      assert(self.F.tuple_of_square(self.t12) == (4446-456-264-234,690*2,30*2,234*4))
      assert(self.F.tuple_of_square(self.t13) == ((93087-4834*2+3354*2+4981)//2,4834*2-4981,-3354*2+4981,-4981*2))
      assert(self.F.tuple_of_square(self.t14) == (24,-7,6,-5))
      assert(self.F.tuple_of_square(self.t15) == ((7263+383+661-259)//4,-(124)//2,(-661-259)//2,259))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()


class BiquadFieldTypeIVb_141_21(BiquadFieldType,unittest.TestCase):

   def setUp(self):
      self.F = BF.createField(141,21)
      self.t1 = (1,0,0,0)
      self.t2 = (0,1,0,0)
      self.t3 = (0,0,1,0)
      self.t4 = (0,0,0,1)
      self.t5 = (-1,0,0,0)
      self.t6 = (0,-1,0,0)
      self.t7 = (0,0,-1,0)
      self.t8 = (0,0,0,-1)
      self.t9 = (1,1,1,1)
      self.t10 = (4,3,2,1)
      self.t11 = (43,14,17,13)
      self.t12 = (24,14,2,4)
      self.t13 = (97,27,1,-22)
      self.t14 = (1,-1,1,-1)
      self.t15 = (-4,2,7,-9)

   def test_trace(self):
      assert(self.F.trace(self.t1) == 4)
      assert(self.F.trace(self.t2) == 2)
      assert(self.F.trace(self.t3) == 2)
      assert(self.F.trace(self.t4) == 1)
      assert(self.F.trace(self.t5) == -4)
      assert(self.F.trace(self.t6) == -2)
      assert(self.F.trace(self.t7) == -2)
      assert(self.F.trace(self.t8) == -1)
      assert(self.F.trace(self.t9) == 9)
      assert(self.F.trace(self.t10) == 27)
      assert(self.F.trace(self.t11) == 247)
      assert(self.F.trace(self.t12) == 132)
      assert(self.F.trace(self.t13) == 422)
      assert(self.F.trace(self.t14) == 3)
      assert(self.F.trace(self.t15) == -7)

   def test_square_trace(self):
      assert(self.F.square_trace(self.t1) == 4)
      assert(self.F.square_trace(self.t2) == 71*2)
      assert(self.F.square_trace(self.t3) == 11*2)
      assert(self.F.square_trace(self.t4) == 123)
      assert(self.F.square_trace(self.t5) == self.F.square_trace(self.t1))
      assert(self.F.square_trace(self.t6) == self.F.square_trace(self.t2))
      assert(self.F.square_trace(self.t7) == self.F.square_trace(self.t3))
      assert(self.F.square_trace(self.t8) == self.F.square_trace(self.t4))
      assert(self.F.square_trace(self.t9) == 185)
      assert(self.F.square_trace(self.t10) == 1277)
      assert(self.F.square_trace(self.t11) == 48681)
      assert(self.F.square_trace(self.t12) == 6578*4)
      assert(self.F.square_trace(self.t13) == 145017*2)
      assert(self.F.square_trace(self.t14) == 125)
      assert(self.F.square_trace(self.t15) == 12763)

   def test_tuple_of_square(self):
      assert(self.F.tuple_of_square(self.t1) == (1,0,0,0))
      assert(self.F.tuple_of_square(self.t2) == (35,1,0,0))
      assert(self.F.tuple_of_square(self.t3) == (5,0,1,0))
      assert(self.F.tuple_of_square(self.t4) == (36,1,-11,-1))
      assert(self.F.tuple_of_square(self.t5) == self.F.tuple_of_square(self.t1))
      assert(self.F.tuple_of_square(self.t6) == self.F.tuple_of_square(self.t2))
      assert(self.F.tuple_of_square(self.t7) == self.F.tuple_of_square(self.t3))
      assert(self.F.tuple_of_square(self.t8) == self.F.tuple_of_square(self.t4))
      assert(self.F.tuple_of_square(self.t9) == ((185-37-15-9)//4,12,14,9))
      assert(self.F.tuple_of_square(self.t10) == ((1277-185-85-51)//4,136//2,134//2,51))
      assert(self.F.tuple_of_square(self.t11) == ((48681-10387-3991-2663)//4,(3991+2663)//2,(10387-2663)//2,2663))
      assert(self.F.tuple_of_square(self.t12) == (6578-696-424-138,562*2,558*2,138*4))
      assert(self.F.tuple_of_square(self.t13) == ((145017+2*10878-2*4394+3461)//2,2*4394-3461,-2*10878+3461,-3461*2))
      assert(self.F.tuple_of_square(self.t14) == (27,-4,14,-3))
      assert(self.F.tuple_of_square(self.t15) == ((12763+2767+203-129)//4,(-203+129)//2,-2896//2,129))


   def test_coeff_bounds(self):
      super().test_coeff_bounds()




if __name__ == '__main__':
    unittest.main()



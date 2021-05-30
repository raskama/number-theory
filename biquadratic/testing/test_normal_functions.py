#!/usr/bin/python3

import unittest
import biquad.biquadfunc
import copy
from abc import ABC


class Test_Normal_Functions(unittest.TestCase):

   def test_safe_floor(self):
      assert(biquad.biquadfunc.safeFloor(5,1,0,1) == 2)
      assert(biquad.biquadfunc.safeFloor(13,1,13,12) == 1)
      assert(biquad.biquadfunc.safeFloor(2<<999,1,-1,1) == ((2<<499) -1))

      t = 99876501662822531623755033154866851165779327101164713298781768659855791732817126707721912969349429671575017847920134313493405850907283733235066629502180206991759018980830247980281373951132916451211308874728315332882047072951690316853513613124748415498408350378170331574728985154681394142345592878415462173103829493531232129810064360798164645509918721215517018150880048267736276179413883490772837808938946935412173930516752935893864948545545868005791779846819027014636810002447321066185432308104873353344400014681790302494260673000134249873572730465752777542605018416497989520044102530904932349141264554332560206044856679990179989515957432998831958253330235056494361950101503080234961024504319076448334775590000822443543393550016796178565112935343770946263570248213514219009555572786880393927363300521003058872074282857982178436005611261797830841660839540377800511575920273942048270905730515416217024036617670233696896801870130605668939074043329717175218364885931690997714331070902502079812227023485172989353422846078647119013116198848075353379456395956908903759472755994029942865113834520776069191360689106781672033880910685382006310073332806811166661351614184882343334624628850242564492107044560855227463677838045756594927696444459828570574997338909356825125636857409715062688641715531258155526220922355629058568936969881106558806808858415409488454086324834848149325138127852026987637509099160790456519071204670795663936707133648425766945282541283150039810785370014446991688350177707621802516694727055266
      assert(biquad.biquadfunc.safeFloor(2<<10000,1,0,2) == t)
      assert(biquad.biquadfunc.safeFloor(0,12,3,2) == 1)
      assert(biquad.biquadfunc.safeFloor(101,13,89,7) == 31)
      assert(biquad.biquadfunc.safeFloor(37,1,0,1) == 6)
      assert(biquad.biquadfunc.safeFloor(36,1,0,1) == 6)
      assert(biquad.biquadfunc.safeFloor(35,1,0,1) == 5)

   def test_combinations(self):
      assert(list(biquad.biquadfunc.combinations(5,5,1)) == [[1,1,1,1,1]])
      assert(list(biquad.biquadfunc.combinations(5,5,0)) == [])
      assert(list(biquad.biquadfunc.combinations(5,5,2)) == [[1,1,1,1,1],[2,1,1,1,0],[2,2,1,0,0]])
      assert(list(biquad.biquadfunc.combinations(5,3,6)) == [[2,2,1],[3,1,1],[3,2,0],[4,1,0],[5,0,0]])
      assert(list(biquad.biquadfunc.combinations(3,3,3)) == [[1,1,1],[2,1,0],[3,0,0]])
      assert(list(biquad.biquadfunc.combinations(1,7,7)) == [[1,0,0,0,0,0,0]])
      assert(list(biquad.biquadfunc.combinations(3,1,3)) == [[3]])
      assert(list(biquad.biquadfunc.combinations(3,1,2)) == [])

   def test_take_one_from_each(self):
      a = [[1,11,111],[2,22],[3],[],[5,55,555],[6],[7]]

      l025_0 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[0,2,5],0):
        l025_0.append(copy.deepcopy(item))
      assert(l025_0 == [[6,3,1],[6,3,11],[6,3,111]])

      l421_0 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[4,2,1],0):
        l421_0.append(copy.deepcopy(item))
      assert(l421_0 == [[2,3,5],[2,3,55],[2,3,555],[22,3,5],[22,3,55],[22,3,555]])

      l04_1 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[0,4],1):
        l04_1.append(copy.deepcopy(item))
      assert(l04_1 == [[5],[55],[555]])

      l0_0 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[0],0):
        l0_0.append(copy.deepcopy(item))
      assert(l0_0 == [[1],[11],[111]])

      l035_0 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[0,3,5],0):
        l035_0.append(copy.deepcopy(item))
      assert(l035_0 == [])

      l65_0 = []
      for item in biquad.biquadfunc.take_one_from_each(a,[6,5],0):
         l65_0.append(copy.deepcopy(item))
      assert(l65_0 == [[6,7]])


   def test_create_sums_of_2(self):
      a = [[(1,2,3,0),(0,0,4,0)],[],[(1111,2,0,2)]]
      b = [[],[(1,1,1,1)],[],[(1,2,3,3),(5,5,5,3)]]
   
      assert(biquad.biquadfunc.create_sums_of_2(a,b,7) == [set(),{(2,3,4,1),(1,1,5,1)},set(),{(1112,3,1,3),(2,4,6,3),(6,7,8,3),(1,2,7,3),(5,5,9,3)},set(),{(1112,4,3,5),(1116,7,5,5)},set(),set()])

      assert(biquad.biquadfunc.create_sums_of_2(a,b,3) == [set(),{(2,3,4,1),(1,1,5,1)},set(),{(1112,3,1,3),(2,4,6,3),(6,7,8,3),(1,2,7,3),(5,5,9,3)}])



if __name__ == '__main__':
    unittest.main()

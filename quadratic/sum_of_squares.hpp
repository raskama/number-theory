#ifndef SUM_OF_SQUARES_H
#define SUM_OF_SQUARES_H

#include <NTL/ZZ.h>
#include <vector>

namespace SumOfSquares
{

   /*Checks if n is squarefree
    *WARNING: requires n to be relatively small, so all primes^2 can fit in LONG_MAX
    */
   long isSquarefree(const NTL::ZZ& n);

   /* Computes the next coefficient in continued fraction. 
    * Format: residue = (a*sqrt(d)+b)/c
    * coeff = floor(residue);
    * updates a,b,c as the new coefficients of (residue - coeff), according to the residue format
    */
   void get_next_CF(NTL::ZZ& coeff, const NTL::ZZ& D, NTL::ZZ& a, NTL::ZZ& b, NTL::ZZ& c);


   //Finds all indecomposables in O_K up to conjugation and multiplication by epsilon^2
   // In form x+y*omega_D 
   //x in vector indecomp[0], y in vector indecomp[1]
   void find_indecomp(const NTL::ZZ& D, std::vector<NTL::ZZ> * indecomp);


   /*Checks if element m*(x+y(1+\sqrt{D})/2) for D \equiv 1 (mod 4)
    * satisfies condition in Peters theorem, e.g. if it is a sum of squares
    * Return: 1 if True, 0 if False
    */
   long test_peters_mod1(const NTL::ZZ& D, const NTL::ZZ& m, const NTL::ZZ& x, const NTL::ZZ& y);


   /*Checks if element m*(x+y*\sqrt{D}) for D \equiv 2,3 (mod 4)
    * satisfies condition in Peters theorem, e.g. if it is a sum of squares
    * Return: 1 if True, 0 if False
    */
   long test_peters_mod23(const NTL::ZZ& D, const NTL::ZZ& m, const NTL::ZZ& x, const NTL::ZZ& y);


   /*Solves the problem for D \equiv 1 (mod 4) for min_m <= m <= max_m
    * -for each D up to given bound finds all unique indecomposables
    * -checks if each indecomposable is sum of squares using Peters
    * Results for m is stored as vector in res[m]
    */
   void solve_mod1(int min_m, int max_m, std::vector<std::vector<NTL::ZZ>> &res);


   /*Solves the problem for D \equiv 23 (mod 4) for even m in  min_m <= m <= max_m
    * -for each D up to given bound finds all unique indecomposables
    * -checks if each indecomposable is sum of squares using Peters
    * Results for 2*m is stored as vector in res[m]
    */
   void solve_mod23(int min_m, int max_m, std::vector<std::vector<NTL::ZZ>> &res);

   /* Checks if all elements of mO+ are sum of squares in Q(sqrt(D))
   * returns 1 if yes, 0 otherwise
   */
   long all_elements_sum_of_squares(const NTL::ZZ& D, const NTL::ZZ& m);

   /* For fixed m finds all D such that all elements of mO+ are sums of squares
   * Stores result in res
   */
   void all_D_for_fixed_m(const NTL::ZZ& m, std::vector<NTL::ZZ>& res);

 
   /* Dumps data into JSON data format for further easy use in Python scripts
    * Format: [array of arrays of results, max_m+1 = len(res)]
    * E. g. for max_m = 3, m equiv 1 (mod 4)
    * [[[],[5], [5], [5,13,17,21]],4]
    */
   void dump_res_into_json(std::vector<std::vector<NTL::ZZ>>& res, std::string filename);


}

#endif

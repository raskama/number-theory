#include "sum_of_squares.hpp"
#include <fstream>
#include <climits>
#include <algorithm>

using namespace std;
using namespace NTL;

namespace SumOfSquares
{
   /*Checks if n is squarefree
    *WARNING: requires n to be relatively small, so all primes^2 can fit in LONG_MAX
    */
   long isSquarefree(const ZZ& n)
   {
      PrimeSeq s;
      long p,p2;
      
      if(n==0) return 0;
      p = s.next();
      p2 = p*p;
      while (p2 <= n){
	 if (divide(n,p2) == 1){
	    return 0;
	 }
	 p = s.next();
	 p2 = p*p;
      }
      return 1;
   }


   /* Computes the next coefficient in continued fraction. 
    * Format: residue = (a*sqrt(d)+b)/c
    * coeff = floor(residue);
    * updates a,b,c as the new coefficients of (residue - coeff), according to the residue format
    */
   void get_next_CF(ZZ& coeff, const ZZ& D, ZZ& a, ZZ& b, ZZ& c)
   {
      ZZ aa, bb, cc, gcd,temp1,temp2,temp3;
      mul(temp1,a,a); mul(temp1,temp1,D); //temp1 = a*a*D
      SqrRoot(coeff,temp1); add(coeff,coeff,b);
      div(coeff,coeff,c);

      mul(aa,a,c); //aa = a*c
      mul(temp2,c,c); mul(temp2,temp2,coeff); mul(temp3,b,c);
      sub(bb,temp2,temp3); //bb = c*c*coeff-b*c;
      mul(temp3,temp3,2); sub(temp3,temp2,temp3); mul(temp3,temp3,coeff);
      mul(temp2,b,b); add(temp3,temp2,temp3);
      sub(cc,temp1,temp3); //cc = a*a*D-b*b-c*c*coeff*coeff+2*b*c*coeff;
      
      GCD(gcd,bb,cc); GCD(gcd,gcd,aa);
      div(a,aa,gcd); div(b,bb,gcd); div(c,cc,gcd);
   }


   //Finds all indecomposables in O_K up to conjugation and multiplication by epsilon^2
   // In form x+y*omega_D 
   //x in vector indecomp[0], y in vector indecomp[1]
   void find_indecomp(const ZZ& D, vector<ZZ> * indecomp)
   {
      vector<ZZ> pi, qi;
      ZZ a,b,c,coeff,reference,temp1, temp2;
      long i,n,t1;
      i = 0; n = -100; //just some random number <-3
      indecomp[0] = {}; indecomp[1] = {};

      pi = {(ZZ)0,(ZZ)1}; qi = {(ZZ)1,(ZZ)0};
      if (rem(D,4) == 1) {
	 t1 = 1;
	 a = 1; b = 1; c = 2;
	 SqrRoot(reference,D);
	 reference = (reference+1+2*((reference-1)/2))/2;
      } else {
	 t1 = 0;
	 a = 1; b = 0; c = 1;
	 SqrRoot(reference,D);
	 mul(reference,reference,2);
      }

      while (i!= n+1) {
	 get_next_CF(coeff,D,a,b,c);
	 mul(temp1,pi.end()[-1],coeff); add(temp1, temp1, pi.end()[-2]);
	 mul(temp2,qi.end()[-1],coeff); add(temp2, temp2, qi.end()[-2]);
	 pi.push_back(temp1);
	 qi.push_back(temp2);
	 
	 if((i % 2) == 1) {
	    temp1 = qi.end()[-3]; temp2 = pi.end()[-3];
	    for(long k = 0; k<coeff;k++){
	       indecomp[1].push_back(temp1);
	       if(t1 == 0){
		  indecomp[0].push_back(temp2);
	       }
	       else{
		  indecomp[0].push_back(temp2-temp1);
	       }
	       add(temp1,temp1,qi.end()[-2]);
	       add(temp2,temp2,pi.end()[-2]);
	    }
	 }
	 
	 if(i>0 && coeff == reference && n<0){
	    n = 2*i;
	 }
	 i+=1;
      }
   }


   /*Checks if element m*(x+y(1+\sqrt{D})/2) for D \equiv 1 (mod 4)
    * satisfies condition in Peters theorem, e.g. if it is a sum of squares
    * Return: 1 if True, 0 if False
    */
   long test_peters_mod1(const ZZ& D, const ZZ& m, const ZZ& x, const ZZ& y)
   {
      ZZ norm, low, high, temp1,temp2,temp3;

      mul(temp1,x,2); add(temp1,temp1,y); //temp1 = 2x+y
      mul(temp2,y,y); mul(temp2,temp2,D);
      mul(temp3,temp1,temp1); sub(temp2,temp3,temp2);
      div(norm,temp2,4); //norm = ((2*x+y)*(2*x+y)-D*y*y)/4;

      mul(temp1,temp1,m); 
      temp2 = (ZZ) (4*m); mul(temp2,temp2,m); mul(temp2,temp2,norm);
      SqrRoot(temp2,temp2);
      sub(low,temp1,temp2); // low = m*(2*x+y) - SqrRoot(4*m*m*norm);
      add(high,temp1,temp2); //high = m*(2*x+y) + SqrRoot(4*m*m*norm);
      if(IsOdd(m*y)){
	      add(high,high,D); add(low,low,D);
      }
      mul(temp1,2,D);//temp1 = 2*D
      rem(temp2,low,temp1); rem(temp3,high,temp1); //temp2 = low % 2D, temp3 = high % 2D
      if((temp2 >= temp3) || (temp2 == 0) || (high-low >= temp1)){
	      return 1;
      }
      return 0;
   }


   /*Checks if element m*(x+y*\sqrt{D}) for D \equiv 2,3 (mod 4)
    * satisfies condition in Peters theorem, e.g. if it is a sum of squares
    * Return: 1 if True, 0 if False
    */
   long test_peters_mod23(const ZZ& D, const ZZ& m, const ZZ& x, const ZZ& y)
   {
      ZZ norm, low, high, temp1,temp2,temp3;
      
      mul(temp1,m,y);
      if(IsOdd(temp1)){
	      return 0;
      }
      mul(temp1,x,x); 
      mul(temp2,D,y); mul(temp2,temp2,y);
      sub(norm,temp1,temp2); //norm = x*x-D*y*y

      mul(temp1,norm,m); mul(temp1, temp1,m);
      SqrRoot(temp1,temp1);
      mul(temp2,m,x);
      sub(low,temp2,temp1); // low = m*x - SqrRoot(m*m*norm);
      add(high,temp2,temp1); //high = m*x + SqrRoot(m*m*norm);

      mul(temp1,2,D);//temp1 = 2*D
      rem(temp2,low,temp1); rem(temp3,high,temp1); //temp2 = low % 2D, temp3 = high % 2D
      if((temp2 >= temp3) || (temp2 == 0) || (high-low >= temp1)){
	      return 1;
      }
      return 0;
   }


      
   /*Solves the problem for D \equiv 1 (mod 4) for min_m <= m <= max_m
    * -for each D up to given bound finds all unique indecomposables
    * -checks if each indecomposable is sum of squares using Peters
    * Results for m is stored as vector in res[m]
    */
   void solve_mod1(int min_m, int max_m, vector<vector<ZZ>> &res)
   {
      ZZ bound_d,D,temp1,temp2, temp3, l1,l2,l3,r1,r2;
      ZZ c0,c4,c9,c16,c25,c36;
      ZZ cs1120,cs10240,cs17920,cs51840,cs280,cs12960,cs22680,cs100000;
      long lbound[3],rbound[3],ok;

      vector<ZZ> v;
      vector<ZZ> indecomposables[2];

      for(int i = res.size(); i<=max_m;i++){
	 res.push_back(v);
      }
      bound_d = max_m*max_m+8*max_m+16;
      D = 5;
      //Precomputing some constants for faster future computations
      c0 = (ZZ) 0; c4 = (ZZ) 4; c9 = (ZZ) 9; c16 = (ZZ) 16; c25 = (ZZ) 25; c36 = (ZZ) 36; 
      cs1120 = (ZZ) 1120; SqrRoot(cs1120,cs1120); add(cs1120,cs1120,1);
      cs10240 = (ZZ) 10240; SqrRoot(cs10240,cs10240); 
      cs17920 = (ZZ) 17920; SqrRoot(cs17920,cs17920); add(cs17920,cs17920,1);
      cs51840 = (ZZ) 51840; SqrRoot(cs51840,cs51840);
      cs280 = (ZZ) 280; SqrRoot(cs280,cs280); add(cs280,cs280,1);   
      cs12960 = (ZZ) 12960; SqrRoot(cs12960,cs12960);
      cs22680 = (ZZ) 22680; SqrRoot(cs22680,cs22680); add(cs22680,cs22680,1);
      cs100000 = (ZZ) 100000; SqrRoot(cs100000,cs100000); 

      while(D < bound_d){
	if(isSquarefree(D) == 1){
	   find_indecomp(D,indecomposables);
	   //even m
	   mul(temp1,c4,D); SqrRoot(temp1,temp1);
	   sub(l1,temp1,c16); add(r1,temp1,cs1120);
	   if (l1 <= max_m){
	      mul(temp2,c16,D); SqrRoot(temp2,temp2);
	      mul(temp3,c36,D); SqrRoot(temp3,temp3);
	      sub(l2,temp2,cs10240); add(r2,temp2,cs17920);
	      sub(l3,temp3,cs51840);
	      
	      rbound[0] = -1; rbound[1] = -1; rbound[2] =-1;
	      lbound[0] = 0; lbound[1] = 0; lbound[2] =0;
	   
	      if(l1 <= c0){
		   lbound[0]=2; rbound[0]=max_m;
	      }
	      else if(r1 >= l2 || r1 >= max_m)
	      {
		   conv(lbound[0],l1); rbound[0]=max_m;
	      }
	      else if (r2 >=l3 || r2 >= max_m)
	      {
		   conv(lbound[0],l1); conv(rbound[0],r1);
		   conv(lbound[1],l2); rbound[1]=max_m;
	      }
	      else
	      {
		   conv(lbound[0],l1); conv(rbound[0],r1);
		   conv(lbound[1],l2); conv(rbound[1],r2);
		   conv(lbound[2],l3); rbound[2] = max_m;
	      }

	      for(int i=0;i<3;i++){
               lbound[i] = max(lbound[i],min_m);
		   lbound[i]+=lbound[i] % 2;
		   for(int m = lbound[i]; m <=rbound[i]; m+=2){
		      ok = 1;
		      for(int i = 0; i<indecomposables[0].size(); i++){
		         if (!test_peters_mod1(D,(ZZ) m,indecomposables[0][i],indecomposables[1][i])){
			      ok = 0;
			      break;
		         }
		      }
		      if (ok){
		         res[m].push_back(D);
		      }
		   }
	      }
	      
	   }
	   //odd m
	   SqrRoot(temp1,D);
	   mul(temp2,c9,D); SqrRoot(temp2,temp2);
	   mul(temp3,c25,D); SqrRoot(temp3,temp3);
	   sub(l1,temp1,c4); add(r1,temp1,cs280);
	   sub(l2,temp2,cs12960); add(r2,temp2,cs22680);
	   sub(l3,temp3,cs100000);

	   rbound[0] = -1; rbound[1] = -1; rbound[2] =-1;
	   lbound[0] = 0; lbound[1] = 0; lbound[2] =0;
	   
	   if(l1 <= c0){
	      lbound[0]=1; rbound[0] = max_m;
	   }
	   else if(r1 >= l2 || r1 >=max_m)
	   {
	      conv(lbound[0],l1); rbound[0] = max_m;
	   }
	   else if (r2 >=l3 || r2 >= max_m)
	   {
	      conv(lbound[0],l1); conv(rbound[0],r1);
	      conv(lbound[1],l2); rbound[1] = max_m;
	   }
	   else
	   {
	      conv(lbound[0],l1); conv(rbound[0],r1);
	      conv(lbound[1],l2); conv(rbound[1],r2);
	      conv(lbound[2],l3); rbound[2] = max_m;
	   }
	   
	   for(int i=0;i<3;i++){
            lbound[i] = max(lbound[i],min_m);
	      lbound[i]+= 1 - (lbound[i] %2);
	      for(int m = lbound[i]; m <=rbound[i]; m+=2){
		 ok = 1;
		 for(int i = 0; i<indecomposables[0].size(); i++){
		    if (!test_peters_mod1(D,(ZZ)m,indecomposables[0][i],indecomposables[1][i])){
		       ok = 0;
		       break;
		    }
		 }
		 if (ok){
		    res[m].push_back(D);
		 }
	      }
	   }
	}
	add(D,D,4);
	indecomposables[0].clear(); indecomposables[1].clear();
      }
   }


   /*Solves the problem for D \equiv 23 (mod 4) for even m in  min_m <= m <= max_m
    * -for each D up to given bound finds all unique indecomposables
    * -checks if each indecomposable is sum of squares using Peters
    * Results for 2*m is stored as vector in res[m]
    */
   void solve_mod23(int min_m, int max_m, vector<vector<ZZ>> &res)
   {
      ZZ bound_d,D,mm,temp1,temp2, temp3;
      ZZ l1,l2,l3,r1,r2;
      ZZ c0,c4,c9;
      ZZ cs70,cs640,cs1120,cs3240;
      long lbound[3],rbound[3],ok,Dtype;

      vector<ZZ> v;
      vector<ZZ> indecomposables[2];
      
      max_m = max_m / 2;
      min_m = (min_m % 2) + min_m/2;

      for(int i = res.size(); i<=max_m;i++){
	 res.push_back(v);
      }
      bound_d = max_m*max_m+8*max_m+16;
      D = 2; Dtype = 0;
      //Precomputing some constants for faster future computations
      c0 = (ZZ) 0; c4 = (ZZ) 4; c9 = (ZZ) 9;
      
      cs70 = (ZZ) 70; SqrRoot(cs70,cs70); add(cs70,cs70,1);
      cs640 = (ZZ) 640; SqrRoot(cs640,cs640);
      cs1120 = (ZZ) 1120; SqrRoot(cs1120,cs1120); add(cs1120,cs1120,1);
      cs3240 = (ZZ) 3240; SqrRoot(cs3240,cs3240);
      
      while(D < bound_d){
	if(isSquarefree(D) == 1){
	   find_indecomp(D,indecomposables);

	   SqrRoot(temp1,D);
	   mul(temp2,c4,D); SqrRoot(temp2,temp2); 
	   mul(temp3,c9,D); SqrRoot(temp3,temp3);
       
	   sub(l1,temp1,c4); add(r1,temp1,cs70);
	   sub(l2,temp2,cs640); add(r2,temp2,cs1120);
	   sub(l3,temp3,cs3240);
	 
	   rbound[0] = -1; rbound[1] = -1; rbound[2] =-1;
	   lbound[0] = 0; lbound[1] = 0; lbound[2] =0;

	   if(l1 <= c0){
	      lbound[0]=1; rbound[0] = max_m;
	   }
	   else if(r1 >= l2 || r1 >=max_m)
	   {
	      conv(lbound[0],l1); rbound[0] = max_m;
	   }
	   else if (r2 >=l3 || r2 >= max_m)
	   {
	      conv(lbound[0],l1); conv(rbound[0],r1);
	      conv(lbound[1],l2); conv(rbound[1],max_m);

	   }
	   else
	   {
	      conv(lbound[0],l1); conv(rbound[0],r1);
	      conv(lbound[1],l2); conv(rbound[1],r2);
	      conv(lbound[2],l3); rbound[2] =max_m;
	   }	
	   for(int i=0;i<3;i++){
            lbound[i] = min(lbound[i],min_m);
	      for(int m = lbound[i]; m <=rbound[i]; m+=1){
		   ok = 1;
		   for(int i = 0; i<indecomposables[0].size(); i++){
		      if (!test_peters_mod23(D,(ZZ) (m*2),indecomposables[0][i],indecomposables[1][i])){
		         ok = 0;
		         break;
		      }
		   }
		   if (ok){
		      res[m].push_back(D);
		   }
	      }
	   }
	}
	//Add 1 if D \equiv 2 (mod 4), else add 3.
	if(Dtype == 0){
	   add(D,D,1);
	   Dtype = 1;
	}
	else
	{
	   add(D,D,3);
	   Dtype = 0;
	}
	indecomposables[0].clear(); indecomposables[1].clear();
      }
   }


   /* Checks if all elements of mO+ are sum of squares in Q(sqrt(D))
   * returns 1 if yes, 0 otherwise
   */
   long all_elements_sum_of_squares(const ZZ& D, const ZZ& m)
   {
      long ok;
      vector<ZZ> indecomposables[2];
      
      if(!isSquarefree(D))
      {
         return 0;
      }
      ok = 1;
      find_indecomp(D,indecomposables);
      if(rem(D,4) == 1)
      {
         for(int i = 0; i<indecomposables[0].size(); i++){
            if (!test_peters_mod1(D,m,indecomposables[0][i],indecomposables[1][i])){
               ok = 0;
               break;
            }
          }
      }else{
         for(int i = 0; i<indecomposables[0].size(); i++){
            if (!test_peters_mod23(D, m,indecomposables[0][i],indecomposables[1][i])){
               ok = 0;
               break;
            }
          }
      }
      if (ok){
         return 1;
      }
      return 0;
   }

   /* For fixed m finds all D such that all elements of mO+ are sums of squares
   * Stores result in res
   */
   void all_D_for_fixed_m(const ZZ& m, vector<ZZ>& res)
   {
      ZZ l,r,D,prev,mm,m2,mm2,temp1, temp2, temp3, temp4;
      ZZ c40,c70, c160,c280,c640,c1120;
      c40 = (ZZ) 40; c70 = (ZZ) 70; c160 = (ZZ) 160; c280 = (ZZ) 280;
      c640 = (ZZ) 640; c1120 = (ZZ) 1120;
      long i,isqr,j;
      vector<ZZ> result = {};
      mul(mm,m,m);
      if(IsOdd(m)){
         //odd m, D \equiv 1
         mul(temp1,m,8); add(temp1,temp1,16); add(r, temp1, mm);
         mul(temp2,mm,c1120); SqrRoot(temp1,temp2); 
         add(l, mm,c280); sub(l, l, temp1);
         if(m < SqrRoot(280)){
          l = 0;
         }
         add(prev,r,1);
         i = 1;
         while(i>0){
            //check corner cases
            if( r>= prev){
               r = prev-1;
               l = 5;
               i = -1;
            }
            if (l <= 5){
               l = 5;
               i = -1;
            }

            rem(D,r,(ZZ)4); if (D == 0) {add(D,D,4);}
            sub(D,r,D); add(D,D,1);
            while (D >= l){
               if (all_elements_sum_of_squares(D,m)){
                  result.push_back(D);
               }
               sub(D,D,4);
            }
            
            //Compute next interval bounds
            if(i >0){
               prev = l;
               isqr = (2*i+1)*(2*i+1);
               div(temp1,mm,isqr);mul(temp2,isqr,c160); 
               mul(temp3,mm,c640); SqrRoot(temp3,temp3);
               add(r,temp1,temp2); add(r,r,temp3); add(r,r,1);
               mul(temp2,isqr,c280); SqrRoot(temp4,temp2);
               if ((m / (2*i-1)) <= temp4){
                  l =0;
               }
               else{
                  mul(temp3,mm,c1120); SqrRoot(temp3,temp3);
                  add(l,temp1,temp2); sub(l,l,temp3);
               }
               i+=1;
            }
         }
         
      }else{
         //even m, D \equiv 1
         div(m2,m,2); mul(mm2,m2,m2);
         mul(temp1,m,8); add(temp1,temp1,64); add(r, temp1, mm2);
         mul(temp2,mm,c280); SqrRoot(temp1,temp2); 
         add(l, mm2,c280); sub(l, l, temp1);
         if(m2 < SqrRoot(280)){
          l = 0;
         }
         add(prev,r,1);
         i = 2;
         while(i>0){
            //check corner cases
            if( r>= prev){
               r = prev-1;
               l = 5;
               i = -1;
            }
            if (l <= 5){
               l = 5;
               i = -1;
            }

            rem(D,r,(ZZ)4); if (D == 0) {add(D,D,4);}
            sub(D,r,D); add(D,D,1);
            while (D >= l){
               if (all_elements_sum_of_squares(D,m)){
                  result.push_back(D);
               }
               sub(D,D,4);
            }
            
            //Compute next interval bounds
            if(i >0){
               prev = l;
               isqr = i*i;
               div(temp1,mm2,isqr);mul(temp2,isqr,c160); 
               mul(temp3,mm,c160); SqrRoot(temp3,temp3);
               add(r,temp1,temp2); add(r,r,temp3); add(r,r,1);
               mul(temp2,isqr,c280); SqrRoot(temp4,temp2);
               if ((m2 / i) <= temp4){
                  l =0;
               }
               else{
                  mul(temp3,mm,c280); SqrRoot(temp3,temp3);
                  add(l,temp1,temp2); sub(l,l,temp3);
               }
               i+=1;
            }
         }
         //even m, D \equiv 2,3
         mul(temp1,m,4); add(temp1,temp1,16); add(r, temp1, mm2);
         mul(temp2,mm,c70); SqrRoot(temp1,temp2); 
         add(l, mm2,c70); sub(l, l, temp1);
         if(m2 < SqrRoot(70)){
          l = 0;
         }
         add(prev,r,1);
         i = 2;
         while(i>0){
            //check corner cases
            if( r>= prev){
               r = prev-1;
               l = 2;
               i = -1;
            }
            if (l <= 2){
               l = 2;
               i = -1;
            }
            
            j =0;
            rem(D,r,(ZZ)4); 
            if (D == 2) {j = 1;}
            else if ( D < 2){
            sub(r,r,D); sub(r,r,1);}
            D = r;

            while (D >= l){
               if (all_elements_sum_of_squares(D,m)){
                  result.push_back(D);
               }
               if (j == 0){
                  sub(D,D,1);
                  j=1;
               }
               else{
                  sub(D,D,3);
                  j=0;
               }
            }
            
            //Compute next interval bounds
            if(i >0){
               prev = l;
               isqr = i*i;
               div(temp1,mm2,isqr);mul(temp2,isqr,c40); 
               mul(temp3,mm,c40); SqrRoot(temp3,temp3);
               add(r,temp1,temp2); add(r,r,temp3); add(r,r,1);
               mul(temp2,isqr,c70); SqrRoot(temp4,temp2);
               if ((m2 / i) <= temp4){
                  l =0;
               }
               else{
                  mul(temp3,mm,c70); SqrRoot(temp3,temp3);
                  add(l,temp1,temp2); sub(l,l,temp3);
               }
               i+=1;
            }
         }
      }
      reverse(result.begin(),result.end());
      res = result;

   }


   /* Dumps data into JSON data format for further easy use in Python scripts
    * Format: [array of results, max_m+1]
    * E. g. for max_m = 3, m equiv 1 (mod 4)
    * [[[],[5], [5], [5,13,17,21]],4]
    */
   void dump_res_into_json(vector<vector<ZZ>>& res, string filename)
   {
      ofstream mfile;
      mfile.open(filename);
      mfile << "[[";

      for(int i=0;i<res.size();i++){
	 mfile << "[";
	 for (int j = 0; j < res[i].size(); j++){
	    mfile << res[i][j];
	    if (j < res[i].size() -1){
	       mfile<< ", ";
	    }
	 }
	 mfile << "]";
	 if (i < res.size()-1){
	    mfile << ", ";
	 }
      }

      mfile << "], ";
      mfile << res.size();
      mfile << "]\n";
      mfile.close();
   }
}

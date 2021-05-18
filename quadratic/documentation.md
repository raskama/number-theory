# Sum of squares -- Documentation
<a name="line-4"></a>
## Description
For general information about the project see [README](README.md). Used algorithms are described in referenced article.
The goal is focused on determining whether for given _m_ and _D_ all _m_-multiples of all totally positive integers in _Q(sqrt(D)_ can be represented as the sum of squares.

- For determining this for single pair _(m,D)_, the function _all_elements_sum_of_squares_ is suitable.
 - For finding all such _D_ for a single fixed _m_ use _all_D_for_fixed_m_. 
 - If one wants to compute this for multiple _m_ in a given interval, the fastest is _solve_mod1_ and _solve_mod23_ (the cases D \equiv 1 (mod 4) and D \equiv 23 (mod 4) are handled separately). Unline reruning the previous function for multiple _m_, this does not recompute indecomposable elements for each D multiple times.


## Example
The following code computes the task for all m <= 500, D \equiv 1 (mod 4) and saves the results in `results_500_mod1.json` in JSON format.

```cpp
include "sum_of_squares.hpp"

using namespace NTL;
using namespace SumOfSquares;
using namespace std;

int main(){
   string filename;
   vector<vector<ZZ>> results = {};
   filename = "results_500_mod1.json";
   
   solve_mod1(1,500,results);
   dump_res_into_json(results,filename);
}
```


## Functions
<a name="line-76"></a>
### 🔷 isSquarefree

```cpp
long isSquarefree(const NTL::ZZ& n)
```

Checks if non-negative number _n_ is squarefree. 
Returns 1 if yes, 0 otherwise.
WARNING: requires _n_ to be relatively small, so all primes^2 can fit in LONG_MAX

##### Parameters:

- `n` -- the number one wants to check

### 🔷 get_next_CF

```cpp
void get_next_CF(NTL::ZZ& coeff, const NTL::ZZ& D, NTL::ZZ& a, NTL::ZZ& b, NTL::ZZ& c)
```
For x  in format (a*sqrt(D)+b)/c computes coeff = floor(x). S
Stores the new parameters of the residue = x - coeff in a,b,c (residue in the same format as x above).

- can be used repeatedly to compute coefficients in continued fraction of (a*sqrt(D)+b)/c


##### Parameters:

- `a,b,c,D` -- parameters of the input number x, after running updated to parameters of x-floor(x)
- `coeff` -- stores floor(x) after computations

##### Example use:

```cpp
#include "sum_of_squares.hpp"

using namespace NTL;
using namespace SumOfSquares;
using namespace std;

int main(){
   ZZ a,b,c,D,coeff;

   a = 1; b = 17; c = 3; D = 23; //represents number x = (1*sqrt(23)+17)/3

   cout << "First 15 coefficients in continued fraction of x: ";
   for(int i = 0; i<15;i++){
      get_next_CF(coeff,D,a,b,c);
      cout << coeff << " ";
   }
   cout << "\n";
}
```
Output: ` First 15 coefficients in continued fraction of x: 7 3 1 3 2 1 13 1 2 3 1 3 2 1 13 `


###  🔷 find_indecomp

```cpp
void find_indecomp(const NTL::ZZ& D, std::vector<NTL::ZZ> * indecomp)
```

Finds all unique indecomposable elements of O_K for K = Q(sqrt(D)) (unique up to conjugation and multiplication by epsilon^2).
The indecomposables in form x+y*omega_D are stored in array of 2 vectors `indecomp = [xx,yy]`, where vector xx contains x's and yy contains y's.

##### Parameters:

- `D` -- parameter of the quadratic field
- `indecomp` -- array of vectors for storing the indecomposables


##### Example use:

```
#include "sum_of_squares.hpp"

using namespace NTL;
using namespace SumOfSquares;
using namespace std;

int main(){
   ZZ D;
   vector<ZZ> indecomp[2];

   D = 6;
   find_indecomp(D,indecomp);

   cout << "Unique indecomposable elements for D = 6:\n";
   for(int i = 0; i<indecomp[0].size();i++){
      cout << indecomp[0][i] << " + " << indecomp[1][i] << " * sqrt(6) \n";
   }
}
```
Output:
```
Unique indecomposable elements for D = 6:
1 + 0 * sqrt(6) 
3 + 1 * sqrt(6) 
5 + 2 * sqrt(6) 
27 + 11 * sqrt(6) 
```

###  🔷 test_peters_mod1

```cpp
long test_peters_mod1(const NTL::ZZ& D, const NTL::ZZ& m, const NTL::ZZ& x, const NTL::ZZ& y)
```

For D \equiv 1 (mod 4) tests if element m*(x+y(1+\sqrt{D})/2) is sum of squares using Peters theorem
Returns 1 if true, 0 otherwise.

##### Parameters:

- `D` -- parameter of the quadratic field
- `m,x,y` -- parameters of the tested element

###  🔷 test_peters_mod23

```cpp
long test_peters_mod23(const NTL::ZZ& D, const NTL::ZZ& m, const NTL::ZZ& x, const NTL::ZZ& y)
```

Analogy of the previous function for D \equiv 2,3 (mod 4).

###  🔷 solve_mod1

```cpp
void solve_mod1(int min_m, int max_m, std::vector<std::vector<NTL::ZZ>> &res)
```

For D \equiv 1 (mod 4) and  min_m <= m <= max_m solves the problem described above. For each _m_, corresponding _D_ are stored as a vector in _res_.

Iterates through all possible D (maximal D is computed using known bounds). For each D computes indecomposable elements and for all m in interval (further restricted by the first 3 known bounds in Theorem 2) checks all of them using _test_peters_mod1_.

##### Parameters:

- `min_m, max_m` -- bounds for m
- `res` -- vector of vectors for storing the results

###  🔷 solve_mod23

```cpp
void solve_mod23(int min_m, int max_m, std::vector<std::vector<NTL::ZZ>> &res)
```
Analogy of the previous function for D \equiv 2,3 (mod 4).
WARNING: Since in this case, we need to consider only even _m_, results ignore odd _m_. Therefore results for _m = 2t_ can be found in vector _res[t]_.

###  🔷 all_elements_sum_of_squares

```cpp
long all_elements_sum_of_squares(const NTL::ZZ& D, const NTL::ZZ& m);
```
Checks if all elements of mO+ are sum of squares in Q(sqrt(D))
Returns 1 if yes, 0 otherwise.
##### Parameters:

- `D` -- parameter of the quadratic field
- `m` -- positive integer

###  🔷 all_D_for_fixed_m

```cpp
void all_D_for_fixed_m(const NTL::ZZ& m, std::vector<NTL::ZZ>& res)
```
For fixed m finds all D such that all elements of mO+ are sums of squares
##### Parameters:

- `m` -- positive integer
- `res` -- vector for storing the results


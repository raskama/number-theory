# Representing multiples of _M_ in real quadratic fields as sums of squares

## General Information
- This project contains the code in C++ for computational parts of my number theory article.
- The article can be found at [_arXiv_](https://www.example.com).
- In a nutshell, for positive integer _m_, quadratic field _Q(sqrt(D))_ and its ring of integers, we are interested if all _m_-multiples of totally positive integers can be represented as sums of squares.
- Outputs of these algorithms, created graphs, as well as interactive web browser graphs can be found at <http://www.example.com>.

## Setup
The project uses [NTL library](https://libntl.org/) for fast computations with unlimited integers (for installation and usage see the linked documentation).

The functions can be included from header file _sum\_of\_squares.hpp_, the implementation is in _sum\_of\_squares.cpp_. Functions are wrapped in namespace _SumOfSquares_.

The code was written in C++14 but does not intentionally use any special features of this version,
 so it should work in any other version suitable for NTL.


## Usage
One very basic example of use is written below (it is informational example on how to include and run written code). 


For more complex examples as well as for the description of intended use, see [documentation](documentation.md).

```cpp
#include "sum_of_squares.hpp" //imports project functions

using namespace NTL;
using namespace SumOfSquares; 
//otherwise functions would have to be called in format SumOfSquares::ImportedFunction

int main(){
   ZZ a;
   a = (ZZ) 3;
   
   //prints 1 if a is squarefree, 0 otherwise
   std::cout << isSquarefree(a);
}
````

To run the code in UNIX terminal one could use the following commands (given the code is in file _example.cpp_ in the same folder as _sum_of_squares.hpp/cpp_):

```bash
% g++ -g -O2 -pthread -march=native sum_of_squares.cpp example.cpp -o example -lntl -lgmp -lm
% ./example
```
Most of the parameters are needed for fast usage of NTL library (for further configuration on specific machine see [here](https://libntl.org/doc/tour-unix.html).

The necessary inclusion of _sum\_of\_squares.cpp_ among compiled files should be apparent.

## Room for Improvement
Multithreading can be definitely used to boost speed.


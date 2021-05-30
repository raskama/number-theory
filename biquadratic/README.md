# Pythagoras numbers of orders in biquadratic fields

## General Information
- This project contains the code in Python for computational parts of number theory article with the same name as the title by Jakub Krásenský, Martin Raška and Ester Sgallová.
- The article is available at [_arXiv_](https://arxiv.org/abs/2105.08860).
- In a nutshell, for biquadratic field, the algorithms take all elements up to a given trace and find out their lengths (i.e. how many squares are needed for a representation by the sum of squares). This gives a good lower bound and estimation for Pythagoras number of the field.
- The project includes class for biquadratic fields (with functionality around sums of squares) and some associated functions


## Setup
The code is written in Python 3.8.5 but it should work with any other newer version (and possibly also older, as long as it is Python 3).

There are mostly two possibilities for using the code. Either simply including it as a file or installing it as a package.

For the first option, all necessary code is in file `biquad/biquadfunc.py`. Then you can simply place your Python file in the same folder and include the classes and functions using `import biquadfunc` (in the example below, the file is located in one directory above).

For the rest of this section, we will talk about installing this as a package on UNIX systems.
The written module uses some packages listed in `requirements.txt`. All of them can be installed for example using command `python3 -m pip install -r requirements.txt`. 
For installation of the package then go the the main folder and run `python3 setup.py install` (could possibly require the `sudo` command). By doing this, the package is installed on your device and could be imported from everywhere using `import biquad.biquadfunc`.



## Usage
One basic example of use it the following (can be found in _example_predict_pythagoras.py_). It finds all elements in field _Q(sqrt(_p_), sqrt(_q_))_ up to a given trace _max_trace_, which are sums of squares and computes their lengths.

It can be called from command line using `python3 predict_pythagoras.py p q max_trace`.
For more complex examples as well as for the description of intended use, see [documentation](documentation.md).

```python
#!/usr/bin/python3
import sys
from biquad.biquadfunc import BiquadField as BF

#the script has 3 input parameters p, q, max_trace
#if these are not provided, takes the default ones below
args = sys.argv
if len(args)<4:
   p = 13
   q = 5
   max_trace = 65
else:
   p = int(args[1])
   q = int(args[2])
   max_trace = int(args[3])

print("Start of computation.")
K = BF.createField(p,q) #Creates the field class instance
#Computes all sums of squares up to max_trace and their lengths
#number of elements is stored in s, elemenets in squares (for possible further use)
s, squares = K.lengths_up_to_trace(max_trace,return_squares = True)
```




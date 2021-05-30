# Pythagoras numbers of orders in biquadratic fields -- documentation
<a name="line-4"></a>
## Description
For general information about the project see [README](README.md). Used algorithms are described in referenced article.
The goal is focused on finding Pythagoras number of the ring of integers in biquadratic fields.
Proving that all that all elements can be represented as a sum of at most _k_ squares is obviously a difficult task, since there are infinitely many elements. However one can get a lower bound by finding an element with big length (i.e. requires a lot of squares) and decent estimates of the Pythagoras number by looking at all elements up to a given trace. So these are the two main goals fo the following code.

##### Class structure

For working with biquadratic fields, we use classes `BiquadField{I/II_III/IVa/IVb}`, each representing one of the four different types of biquadratic fields (see article preliminaries). They all inherit from abstract class `BiquadField`, which implements most of the algorithms regarding sum of squares, while these non-abstract classes implement the type-of-field-specific functions such as trace, square of the integer, etc.

Integral elements are represented as quadruples `(a,b,c,d)` with _a,b,c,d_ rational integers and this quadruple correspond to the integral basis of the specific field (each of the fields of type BI, BII/III, BIVa, BIVb have differeny type of integral basis). We will work with the maximal order (for which we have implemented the 4 different derived classes), however, one can implement this kind of derived class for non-maximal order as well (or for quadratic fields), as long as integral elements are represented as quadruples and functions for trace, tuple of square, etc. are implemented.

##### Type of methods
For used algorithms see the mentioned article. The fastest methods generate step by step elements with length 1, 2, 3, etc. up to a given trace.

A slower alternative (which on the other hand requires less memory) is to look at elements of one fixed trace and generate them combining elements of smaller traces. This is less time efficient, but could perhaps provide additional insight and is usefull when high memory use is a problem. These methods are indicated by suffix _by_trace_.

##### Warning
Most of the functions are not written to be foolproof, so wherever user input is used, we require some restrictions.
E.g. when we talk about sum of _k_ squares, we mean _k > 0_. Therefore, e.g. if some function checks if an element is a sum of _k_ squares, we assume the argument _k_ to be positive. The same goes for trace in the most cases, since most of the functions regarding sums of squares are designed to work only with totally positive integers. All of these and other assumptions should be pretty intuitive and are written later in the documentation. You have been warned.


## Abstract class BiquadField

 As was mentioned, from this class are derived all the used classes _BiquadFieldI_, _BiquadFieldII_III_, _BiquadFieldIVa_, _BiquadFieldIVb_, which just implements some of the methods marked here as abstract.

### Attributes

 - `p,q,r` -- squarefree positive integers, corresponds to the integral basis
 -  `p0,q0,r0` -- positive integers, gcd of the other two variables, e.g. _p0 = gcd(q,r)_
 -  `elem_cap` -- integer
 -  `elements` -- list of sets, set _elements[i]_ contains all squares with trace _i_, is counted up to _elem_cap_ and is automatically updated to a needed cap everytime some of the methods uses it

### Constructor

##### `__init__(p,q)`
- automatically computes attributes _r_ and _p0,q0,r0_

### Methods
<a name="line-76"></a>
#### 🔷 createField

```python
createField(p,q)
```
Chooses the correct derived class of biquadratic field (depending on the congruence classes of _p,q_) and returns created instance of this field.
Tries to assign values _p,q_ to the corresponding attributes _p,q_ in the biquadratic field if possible but is sometimes forced to change the order due to the integral basis, e.g. _createField(3,7)_ will return _BiquadFieldI(3,21)_.

If _p,q_ are not squarefree returns `None`.
##### Parameters
- `p,q` -- squarefree positive integers

#### 🔷 trace (abstract method)

```python
trace(tple)
```
Return trace of the integral element represented by _tple_.
##### Parameters
- `tple` -- tuple (a,b,c,d) representing the element

#### 🔷 square_trace (abstract method)

```python
square_trace(tple)
```
Return trace of the square of the integral element represented by _tple_.
##### Parameters
- `tple` -- tuple (a,b,c,d) representing the element

#### 🔷 tuple_of_square (abstract method)

```python
tuple_of_square(tple)
```
Return tuple of the square of the integral element represented by _tple_.
##### Parameters
- `tple` -- tuple (a,b,c,d) representing the element

#### 🔷 coeff_bounds (abstract method)

```python
coeff_bounds(trace)
```
Returns bounds `max_a, max_b, max_c, max_d` such that for all elements _(a,b,c,d)_, such that the trace of their square is at most _trace_, holds _a<=max_a_, _b <= max_b_, etc.
##### Parameters
- `trace` -- positive int

#### 🔷 elements_with_small_trace

```python
elements_with_small_trace(max_trace)
```
Returns all squares with trace up to max_trace in the same format as the attribute elements.
##### Parameters
- `max_trace` -- positive int

#### 🔷 square_root

```python
square_root(tple)
```
Return (only) one element _(a,b,c,d)_ such that _(a,b,c,d)^2 = tple_ (as integral elements), if it exists. Otherwise returns None.
##### Parameters
- `tple` -- tuple representing the squared element

#### 🔷 update_elements

```python
update_elements(trace)
```
Updates the attribute `elements` to contain all the elements with trace up to at least _trace_. 
##### Parameters
- `trace` --  positive int

#### 🔷 lengths_up_to_trace

```python
lengths_up_to_trace(max_trace,max_k = 7, return_squares = False, print_results = True)
```
Finds all elements which are sums of squares up to a given `max_trace` and computes their lengths (computes only elements with length <= max_k).
In the basic settings returns array `s` containing number of elements with the given length, i.e. `s[2]` contains number of elements with length 2.

##### Parameters
- `max_trace` -- positive int
- `max_k` -- positive int, optional, we only compute elements with lengths <= max_k
- `return_squares` -- boolean, optional, if true returns the computed elements as well (the function returns _(s, squares)_), where _squares_ contains the computead elements in the same format as the attribute _elements_
- `print_results` -- boolean, optional, if true, prints the results into the standart output

#### 🔷 is_sum_of_k_squares

```python
is_sum_of_k_squares(k,tple)
```
Computes if elements `tple` can be represented as the sum of `k` squares. Returns boolean value.

##### Parameters
- `k` -- positive int
- `tple` -- tuple of the corresponding integral elements

#### 🔷 find_squares_forming_sum_by_trace

```python
find_squares_forming_sum_by_trace(k,tple)
```
Finds all the decompositions of `tple` into sum of `k` squares and returns list of lists of possible decomposition. Can contain sum decomposition multiple times, just in different ordering.   

##### Parameters
- `k` -- positive int
- `tple` -- tuple of the corresponding integral elements

#### 🔷 sums_of_k_squares_by_trace

```python
sums_of_k_squares_by_trace(trace,k,length=True)
```
Computes the number of elements with trace `trace`, which can be represented as the sum of `k` squares.

##### Parameters
- `k,trace` -- positive int
- `length` -- boolean, optional, if True returns just the number of elements, if False returns the elements as well as their decompositions (in dictionary with the elements as the keys and list of decompositions as the value).

#### 🔷 compare_by_trace

```python
compare_by_trace(trace,min_k = 1, max_k = 7)
```
Computes `sums_of_k_squares_by_trace` for all _k_ between `min_k` and `max_k`, prints and returns the results (in array of ints).

##### Parameters
- `min_k, max_k, trace` -- positive int

#### 🔷 find_culprits_by_trace

```python
find_culprits_by_trace(trace,k, known_roots = {}):
```
Finds elements with trace `trace` which are sums of `k` squares but are not sums of `k-1` squares. 
Returns them in the list of lists with computed decompositions into the sums of squares (does not return the squares, but their square roots).

##### Parameters
- `k,trace` -- positive int
- `known_roots` -- dictionary, optional, can supply a dictionary of known square root of tuples for faster computing.


## Other functions (outside the class)

#### 🔷 isSquareFree

```python
isSquareFree(n)
```
Computes if positive integer `n` is squarefree.

##### Parameters
- `n` -- positive int

#### 🔷 safeFloor

```python
safeFloor(D,a,b,c)
```
Computes the floor of (a*sqrt(D)+b)/c. 

##### Parameters
- `D` -- non-negative int
- `a,b,c` -- int

#### 🔷 combinations

```python
combinations(n,k,l)
```
Generates all `k`-tuples (actually lists) with sum `n`, meaning (a1,a2,...,ak), a1 >= a2 >= a3... such that a1+a2+...+ak = `n`, such that a1<=`l`

##### Parameters
- `n,k,l` -- positive int

#### 🔷 take_one_from_each

```python
take_one_from_each(elements, indices,k)
```
Generates all the possible lists `res` of elements, such that `res[i]` is in `elements[indices[len(indices)-1-i]]`. Parameter `k` means that first `k` numbers in `indices` are ignored (i.e. the default parameter should be 0).

Warning: should be used as in the example below, due to the fact that the list from previous iteration is reused for the next iteration -- meaning if one wants to use the list later, e. g. example deep_copy has to be made
##### Parameters
- `elements` -- list of lists/sets
- indices -- list of ints, such that max(indices) < len(elements)
- 'k' -- positive int

##### Example use:


```python 
for item in take_one_from_each([[1,11],[2],[3,33,333]],[1,1,2],0):
    print(item)
```
Output:
```
[3,2,2]
[33,2,2]
[333,2,2]
```
#### 🔷 create_sums_of_2

```python
create_sums_of_2(elements1,elements2,max_ind)
```
Combines two list of sets of tuples into one using the rule that tuples in `results[i]` are all possible sums of one elements from `elements1[j]` and one elements from `elements2[k]` such that `j+k = i`. Counts only for `i <= max_ind`.

##### Parameters
- `elements1, elements2` -- list of sets of quadruples (same format as the attribute elements)
- `max_ind` --positive int


# Counting Permutation Patterns with Multidimensional Trees

This script accompanies the paper, and shows the following:
- Full rank is obtained for $\mathbb{S}_{\le 7}$ using pattern-trees with maximum vertex size $2$, as explained in Section 3.
- Full rank is obtained for $\mathbb{S}_{\le 5}$ using a combination of pattern-trees with maximum vertex size $1$ (corner-trees) and the special vertices corresponding to $\mathtt{3214}$ and $\mathtt{43215}$, as explained in Section 4.

## Technical Information
The script takes several minutes to complete and uses a few GBs of memory. It was tested with Python 3.10 and Sage 9.5 - please note that other configurations may have compatibility issues.

On a fresh install of Ubuntu 22.04, update and install Sage:
```
$ sudo apt update
$ sudo apt full-upgrade
$ sudo apt install sagemath
```
Verify versions:
```
$ python3 --version
Python 3.10.12
$ sage --version
SageMath version 9.5, Release Date: 2022-01-30
```
Then, run the script in this directory. Partial rank is displayed during the run to indicate progress.
```
$ python3 main.py
****** O(n^2 polylog n) for k <= 7 ******
Number of permutations: 5913 = 1! + 2! + 3! + 4! + 5! + 6! + 7!
Total 131072 rows. Rank: 5913

****** O(n^(7/4) polylog n) for k <= 5 ******
Number of permutations: 153 = 1! + 2! + 3! + 4! + 5!
Total 1024 rows. Rank: 153
```

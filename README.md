
# Cloning algorithm

This repository contains:

- c++ code used in reference [1] for the evaluation of the scaled cumulant generating functions (SCGFs) in:
    - Non-Markovian TASEP (Totally Asymmetric Simple Exclusion Process) `nonMarkovTASEP`,
    - Semi-Markov model for ion-channel gating with direction-time independence (DTI) `SemiMarkovDTI`,
    - Semi-Markov model for ion-channel gating without DTI `SemiMarkovNonDTI`.

- `cloning.py`, a Python implementation of the cloning algorithm for non-Markovian processes with heap sort, based on the examples at https://www.lpsm.paris//pageperso/lecomte/warwick-summer-school_2013.html. A basic Jupyter notebook demonstration (with Process-based parallelism) is [here](https://github.com/mcavallaro/cloning/blob/master/cloning.ipynb).

[1] Cavallaro, M., & Harris, R. J., [A framework for the direct evaluation of large deviations in non-Markovian processes](https://doi.org/10.1088/1751-8113/49/47/47LT02988), Journal of Physics A: Mathematical and Theoretical, 49, 47LT02 (2016) 171–176.

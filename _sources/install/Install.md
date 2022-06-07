# Installation Guide
This chapter walks you through the steps needed to install `dengo`.
Although it is not neccessary, it is however highly recommended users to install it a seperate environment, such that it does not mess with your current python packages and installations.

The repository of `dengo` locates at https://github.com/hisunnytang/dengo-merge.
```python
git clone https://github.com/hisunnytang/dengo-merge
cd dengo-merge
pip install -r requirements.txt
pip install -e .
```

[ChiantiPy](https://chiantipy.readthedocs.io/en/latest/) is a required package that contains the atomic database for astrophysical specroscopy data. Please read the their documentation or install-chiantipy.sh available in this repository.

[Sundials CVODE](https://computing.llnl.gov/projects/sundials/cvode) is a `LLNL` maintained solver for stiff and nonstiff ordinary differential equation (ODE) systems (initial value problem) given in explicit form $dydt = f(t,y)$. Current `dengo` templates work with CVODE with version greater thatn v5.3.0. If your reaction network has a sparse jacobian in nature, CVODE could take advantage of the sparsity, and use a sparse solver of Newton Iteration. [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) is a suite of sparse matrix algorithms that CVODE can take use to speed up the computation.

Note that CVODE and SuiteSparse are optional. While `Dengo` can function with the current installation, `Dengo` could make use of external ODE solver libraries to obtain additional performance gain. 

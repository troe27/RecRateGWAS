## known peculiarities


#### GMMAT fails with longitudinal data if the repeat-observations are grouped, and not interleaved.

- On my local setup, `glmmkin()` failed when repeated rows were grouped contiguously by individual, but the same logic worked when repeated observations were interleaved across individuals.
 the associate error message was:


```bash 
Duplicated id detected...
Assuming longitudinal data with repeated measures...
Error in asMethod(object) : the matrix is not triangular
```



#### GMMAT fails when the response variables are uniformely close to zero
- on my local setup the software failed when the phenotypes were too small. 
  Here, switching from crossover-per-basepair to crossover-per-Megabase helped. but z-score standardisation also worked, though not implemented here.
  The associated error message looked like this:

```bash
  Duplicated id detected...
  Assuming longitudinal data with repeated measures...
  Warning: CHOLMOD warning 'not positive definite' at file 'Cholesky/t_cholmod_rowfac_worker.c', line 449
  Status: FAIL
  Error: error in evaluating the argument 'x' in selecting a method for function 'chol2inv': error in evaluating the argument 'x' in selecting a method for function 't': leading principal minor
  of order 1 is not positive
``` 

OR

```bash
  Error in .local(x, ...) :
    internal_chm_factor: Cholesky factorization failed
``` 
OR

```bash
  Error: inv_sympd(): matrix is singular or not positive definite
```

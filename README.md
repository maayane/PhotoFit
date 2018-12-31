# PhotoFit
Calculating the evolution in time of the effective radius and effective temperature from multiple-band photometric light-curves (assuming blackbody SEDs).

[![PyPI](https://img.shields.io/pypi/v/SLAB-Diffusion.svg?style=flat-square)](https://pypi.python.org/pypi/SLAB-Diffusion)

```python
>>> import PhotoFit
>>> PhotoFit.calculate_T_and_R_in_time()
```

## Documenation

`PhotoFit` is a package to model the effective radius and temperature from multiple-bands photometry of a SN


## How to install the `PhotoFit` code?

### pip

`pip install PhotoFit`

### Python version
* `python 2`: higher than `2.7.10`
* `python 3`

### Required python packages
* `math`
* `numpy`
* `pylab`
* `math`
* `matplotlib`
* `pyphot`

## How to run the `PhotoFit` code?

### Edit the params.py file

After every modification of `params.py`, rerun

```python
>>> python setup.py install
```
from the `PhotoFit` directory

### Visualize the evolution of R and T for a given slab

```python
>>> import PhotoFit
>>> PhotoFit.calculate_T_and_R_in_time()
```

# PhotoFit
Calculating the evolution in time of the effective radius and effective temperature from multiple-band photometric light-curves (assuming blackbody SEDs).

THIS IS STILL UNDER CONSTRUCTION (Please wait for 7.01.2018 before installing anything)

[![PyPI](https://img.shields.io/pypi/v/SLAB-Diffusion.svg?style=flat-square)](https://pypi.python.org/pypi/SLAB-Diffusion)

```python
>>> import PhotoFit
>>> PhotoFit.calculate_T_and_R_in_time()
```

## Documenation

`PhotoFit` is a package to model the effective radius and temperature from multiple-bands photometry of a supernova.

### How does PhotoFit work?

1. Your data are usually taken at different dates for each band. The first thing to do is to define common epochs on which
the data will be interpolated, and calculate the interpolated flux and (more tricky) the interpolated errors. `PhotoFit` does this task using Monte Carlo Markov Chains simulations.
Although this task can be time-consuming (~0.5 min per band and per epoch), your can set the parameters (in the parameters file `params.py`) to only do this once for a given data set.
At the end of this first step, there is one SED per epoch.

2. `PhotoFit` then fits each SED with a blackbody model that has been corrected for:
    - the extinction: `PhotoFit` does this using Schlay & Finkbeiner (2011) and using the extinction curves of Cardelli et al. (1989).
    - the redshift
    - the effect of the filters transmission curves: `PhotoFit` does this using the `pyphot` package for synthetic photometry.

3. The fit itself can be done in two different ways (to choose from in the parameters file params.py):
    - Monte Carlo Markov Chain (with emcee). The advantage of this option is it gives you error bars on T and R. The disadvantage is that it is time-consuming
(~30 min per epoch for 100 walkers and 350 steps )
    - A linear fit with a grid of temperatures. The advantage is the speed. The disadvantage is you loose the error bars.

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
* `emcee`
* `matplotlib`
* `pyphot`

## How to run the `PhotoFit` code?

### Edit the params.py file

After every modification of `params.py`, rerun

```python
>>> python setup.py install
```
from your `PhotoFit` directory

### Calculate the evolution of R and T

The simplest way to run PhotoFit is
```python
>>> import PhotoFit
>>> Best=PhotoFit.calculate_T_and_R_in_time()
```
`Best` is a numpy array where the first column is the time (jd), the second column is the temperature (K) and the third column is the radius (cm).

By default, this will show you several plots by the time it finishes running:

If you do not want to see and save these plots, you can set `show_underlying_plots` to `False`, i.e. run

```python
>>> Best=PhotoFit.calculate_T_and_R_in_time(show_underlying_plots=False)
```
And if you want the code to tell you more about what it is doing at each step, you can set `verbose` to `True`:

```python
>>> Best=PhotoFit.calculate_T_and_R_in_time(verbose=True)
```
The result is stored in the `Results.txt` file of your output directory (defined in the `params.py` file). In addition to this file, the code
creates in this directory one sub-directory per epoch, with several files, the most important of which being `fit_result_FLux.pdf`, with a plot of
the data and the best fit:

show plot


### Visualize the evolution of R and T

### Visualize the evolution of L

```python
>>> import PhotoFit
>>> PhotoFit.()
```

### Visualize the spectral energy distributions (SEDs) at each epoch

```python
>>> import PhotoFit
>>> PhotoFit.calculate_T_and_R_in_time()
```

## The parameters file in details

## Give it a try with the test data!





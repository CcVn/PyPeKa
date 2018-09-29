# PyPeKa
An attempt to make a PK (then in a second step a PBPK) program in Python (this project is mainly for learning Python purpose)

My ultimate objective is to translate a huge script I made in a defunct Matlab competitor (acslX, RIP) to pilot PBPK models (the script was written in a version of the M langage, the Matlab langage, and the models were coded in the CSL langage used by this software which automatically translated it to C compiled code).

As a noob in Python, and to avoid reinventing the wheel, I have looked first for other examples made
- in Python : I found more or less nothing except this small script: https://github.com/alexrd/pk/blob/master/run_pk.py),
which may be useful to
- in R (found some packages that might be used as inspiration guide)
- and in Matlab : I found only one useful script, limited to PK models, but very appealing to me because 
  1) written in M langage (ie not using the proprietary Simbiology tool),
  2) written in an industrial context -and published-,
  3) it uses object oriente paradigm, which I wanted to use. 
  It can be found here in the Matlab community repository: https://fr.mathworks.com/matlabcentral/fileexchange/43521-preclinical-pkpd-modeling

So my first move will be to try to translate the Matlab package into Python with Numpy/Scipy modules.

If it works, I will try to extend its library with PBPK models (much more complex than 1-2-3 cpt PK models)

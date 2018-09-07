# PyPeKa
A foulish attempt to make a PK (then maybe PBPK) program in Python (for learning Python purpose)

This is public despite my very poor skills, because making it private requires a paid account ...

My primary goal is to translate a huge script I made in a defunct Matlab competitor software (acslX, RIP) to pilot PBPK models coded in the CSL langage used by this software (the script was written in a version of the M langage, the Matlab langage).

As a noob in Python, and to avoid reinventing the wheel, I have looked first for other examples made in Python (found more or less nothing), in R (found some packages that might be used as inspiration guide) and Matlab (found only one script, limited to PK models, but very appealing to me because 1) not coded using the proprietary Simbiology tool, 2) written in an industrial context -and published- and 3) it uses OOP). It can be found here in the Matlab community repository
https://fr.mathworks.com/matlabcentral/fileexchange/43521-preclinical-pkpd-modeling

So my first move will be to try to translate the Matlab package into Python/Numpy/Scipy.

If it works, I will try to extend its library with PBPK models (much more complex than 1-2-3 cpt PK models)

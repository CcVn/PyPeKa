#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 2018
From Lindhardt E, Gennemark P
https://fr.mathworks.com/matlabcentral/fileexchange/43521-preclinical-pkpd-modeling

All the main PK classes used
    SimFunction
    Model
    Experiment
    Dosing
    DataVariable

@author: Matlab code by Linhardt, converted to Python by CcVN
@date: 6/9/2018
"""
#import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats as stat
#from scipy.optimize import fmin as fminsearch
#from scipy.optimize import least_squares as lsq
#import LibNicePlots as nplt

class optimset(object):
    ''' Ersatz de optimset Matlab
     'Display','off', - displays no output.
     'Diagnostics','off', - does not display any diagnostics
     'MaxIter',2000, - Maximum number of iterations allowed is set to 2000.
     'TolFun',1e-10, - Termination tolerance on the function value.
     'TolX',1e-10 - Termination tolerance on x.
     '''

    def set_all_param(self, disp=0, diag=0, maxiter=2000, tolfun=1e-10, tolx=1e-10):
        self.disp = disp
        self.diag = diag
        self.maxiter = maxiter
        self.tolfun = tolfun
        self.tolx = tolx

    def __init__(self, disp=0, diag=0, maxiter=2000, tolfun=1e-10, tolx=1e-10):
        self.set_all_param(0, 0, 2000, 1e-10, 1e-10)
        #self.set_all_param(disp, diag, maxiter, tolfun, tolx)

    def get_val(self):
        return (self.disp, self.diag, self.maxiter, self.tolfun, self.tolx)

def non_linear_parameters_95_percent_confidence_interval(fvec, jac):
    ''' To replace a built-in Matlab function.  Returns the 95% confidence interval 
    on parameters from non-linear fit results.'''
    
    rss = np.sum(fvec**2) # residual sum of squares
    n, p = jac.shape # number of data points and parameters
    nmp = n - p # the statistical degrees of freedom
    ssq = rss / nmp # mean residual error
    J = np.matrix(jac) # the Jacobian
    c = np.linalg.inv(J.T*J) # covariance matrix
    pcov = c * ssq # variance-covariance matrix.
    # Diagonal terms provide error estimate based on uncorrelated parameters.
    # The sqrt convert from variance to std. dev. units.
    err = np.sqrt(np.diag(np.abs(pcov))) * 1.96  # std. dev. x 1.96 -> 95% conf
    # Here err is the full 95% area under the normal distribution curve. This
    # means that the plus-minus error is the half of this value
    return err

class BadInputException(Exception):
    ''' Base class for exceptions in this module. '''
    def __init__(self, identifier, message, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
        self.identifier = identifier
        print(message)
        self.message = message

class SimFunction(object):
    ''' SimFunction : Superclass for defining a dynamic mathematical model that can be simulated.
         No direct use of this class is recommended. 
         Examples of subclasses are PK_1comp_Ka_Cl_fun and PK_1comp_Ka_Cl_ode. 
         A SimFunction can either be a function or an ODE. 
         A Model is built of 1, 2 or 3 SimFunctions.

     Inputs:
      - parameterNames, vector of parameter names
      - parameterValues, vector of parameters
      - inputVariables, vector of input variable names
      - outputVariables, vector of output variable names

     Method: (dummy in the superclass)   TO BE DELETED?
      - Y = X.simulate(times, experiment, yInput, functionObj), simulates
         self - the object to simulate
         times - vector of time points to output
         experiment - the experiment to simulate (wrt dosing)
         yInput input variable values at time points in times, size|times|x|inputVariables|
         functionObj from which inputVariables can be calculated at arbitrary time-points

         Returns: Y matrix of output variable data, size |times|x|inputVariables|
    
      - outputMatchesInput
      - parametersAreCompatibleWith
    '''
    verbose = False

    def __init__(self, pNames={}, pValues=[], inputVar={}, outputVar={}):
        if pNames is not None: self.parameterNames = pNames
        if pValues is not None: self.parameterValues = pValues
        if inputVar is not None: self.inputVariables = inputVar
        if outputVar is not None: self.outputVariables = outputVar

    def initModelPD(self, pNames={}, pValues=[], inputVar={}, outputVar={}):
        ''' FIXME ''' #FIXME
        self.inputVariables = inputVar
        self.outputVariables = outputVar
        self.parameterValues = pValues
        if pValues is None:
            raise Exception('BadInput: At least one input argument is required.')
        if len(pValues) != len(self.parameterNames):
            raise Exception('BadInput: Number of parameters incorrect!')

    def __str__(self):
        message = 'Creating a ' + str(type(self).__name__) + ' type-like class with the following parameters/values\n'
        #ou self.__class__.__name__
        if self.calculateEffectCompartment > 0:
            message += '(effect compartment is activated)\n'
        for i, nom in enumerate(self.parameterNames):
            message += '\t-' + nom + ' : ' + str(self.parameterValues[i]) + '\n'
        return message

    def simulate(self, times, experiment, yInput, functionObj):
        '''  Dummy response in the super class '''
        return np.random.randn( (len(times), len(self.outputVariables)) )

    def outputMatchesInput(self, inputFunction):
        '''  TO DO '''
        yesno = True
        #nameidx = getnameidx(self.outputVariables, inputFunction.inputVariables) #renvoie les positions dans self.outputVariables ou apparaissent des elements de inputFunction.inputVariables
        for i, iV in enumerate(inputFunction.inputVariables):
            exists = False
            for j, oV in enumerate(self.outputVariables):
                if iV == oV:
                    exists = True

            if not(exists):
                yesno = False
                return yesno
        return yesno

    def parametersAreCompatibleWith(self, otherFunction):
        '''  TO DO '''
        yesno = True
        for i in enumerate(self.parameterNames):
            for j in enumerate(otherFunction.parameterNames):
                pName1 = self.parameterNames[i]
                pName2 = otherFunction.parameterNames[j]
                if pName1 == pName2: # strcmp(pName1, pName2) == 1
                    if self.parameterValues[i] != otherFunction.parameterValues[j]:
                        yesno = False
                        return yesno
        return yesno

class Model(object):
    ''' Model
     Defines a mathematical model

     Inputs:
     - FunctionPK, a PK model defined on closed form, or None
     - ODE, a PK model, a PD model, or a PKPD model defined on ODE form, or None
     - FunctionPD, a PD model defined on closed form, or None

     For an ODE that requires input variables, those variables must be contained
     in the output variables of FunctionPK.
     For a FunctionPD that requires input variables, those variables must be
     contained in the output variables of ODE, or if ODE is empty, in the output
     variables of FunctionPK.
     
     Methods:
     - fit(self, exps, lb, ub, ERROR_FUNCTION)
     - errorFunction(x,self,exps,ERROR_FUNCTION,lb)
     - fitLL
     - errorFunctionLL
     - simulate
     - utilitaries : getOutputVariables, verifyFunctionParameterOrder, 
                     verifyFunctionParameters, verifyInputChain, matrixSubset,
                     getParameterSubsetValues, getInputParametersFromFunctions,
                     getFreeParameterValues, getParameterValues
    '''

    verbose = False

    def __init__(self, FunctionPK = None, FunctionPKode = None, FunctionPD = None, Parameters = None, ISFREE = None):
        try:
            self.FunctionPK = FunctionPK if FunctionPK!=[] else None
            self.FunctionPKode = FunctionPKode if FunctionPKode!=[] else None
            self.FunctionPD = FunctionPD if FunctionPD!=[] else None
            self.Parameters = Parameters if Parameters!=[] else None
            self.ISFREE = ISFREE if ISFREE!=[] else None

            if not self.verifyInputChain():
                raise BadInputException('ResultChk:BadInput', 'Input chain is incorrect. Verify that output variables from one model part match inputs to the subsequent part.')
    
            self.Parameters = self.getInputParametersFromFunctions( (self.FunctionPK, self.FunctionPKode, self.FunctionPD) )
    
            print(self.verifyFunctionParameterOrder())
            if not self.verifyFunctionParameterOrder():
                raise BadInputException('ResultChk:BadInput', 'The function parameter values are not defined in order!')
    
            if not self.verifyFunctionParameters():
                raise BadInputException('ResultChk:BadInput', 'There is a conflict wrt parameter names!')
    
        except BadInputException:
            raise

    def fit(self, exps, lb, ub, ERROR_FUNCTION):
        '''
        Fit function
        '''
        if len(lb) != len(self.Parameters):
            raise BadInputException('ResultChk:BadInput', 'The length of lb does not match number of parameters.')
        if len(ub) != len(self.Parameters):
            raise BadInputException('ResultChk:BadInput', 'The length of ub does not match number of parameters.')
        if ERROR_FUNCTION < 0 or ERROR_FUNCTION > 4:
            raise BadInputException('ResultChk:BadInput', 'The value of ERROR_FUNCTION is illegal.')
        
        #=== Initialization ===
        # Check for free parameters
        self.ISFREE = np.zeros( (1, len(lb) ))
        for i in enumerate(lb):
            if lb[i]!=ub[i]:
                self.ISFREE[i]=1
        lbnew = lb[self.ISFREE>0]
        ubnew = ub[self.ISFREE>0]
        
        # Harmonize time points in each experiment
        for exp in exps:
            exp.harmonizeTimes()
        
        # 1. Are there output variables from this model for which we have experimental 
        # data to compare to? If not, throw an error.
        modelOutVar = self.getOutputVariables(self)
        ok = 0
        for exp in exps:
            for dv in exp.dataVariables:
                for outvar in modelOutVar:
                    if dv.name == outvar:
                        ok = 1
                        break
        if not ok:
            raise BadInputException('ResultChk:BadInput', 'There are no experimental data variables that match model outputs.')
        
        # 2. Fit data, error is calculated in errorFunction
        
        x0 = self.getFreeParameterValues
        options = optimset() 
        #Matlab code uses lsqnonlin
        [x, resnorm, residual, exitflag, output, lambda_, Jaco] = \
            np.least_squares(self.errorFunction, x0, lbnew, ubnew, options.get_val(), self, exps, ERROR_FUNCTION, lb)

# Matlab nlparci :  nonlinear regression parameter confidence intervalscollapse all in page
# Syntax: ci = nlparci(beta,resid,'jacobian',J) is an alternative syntax that also computes 95% confidence intervals. J is the Jacobian computed by nlinfit. If the 'robust' option is used with nlinfit, use the 'covar' input rather than the 'jacobian' input so that the required sigma parameter takes the robust fitting into account.
#        cifree = np.nlparci(x,residual,'jacobian',Jaco)
# no equivalent in numpy: found non_linear_parameters_95_percent_confidence_interval
            
        cifree = non_linear_parameters_95_percent_confidence_interval(x,Jaco)
        
        delta = (cifree[:,2] - cifree[:,1]) / 2
        #scipy t.ppf = Matlab tinv= Student's t inverse cumulative distribution function
        sefree = delta / stat.t.ppf(1 - 0.05/2, len(residual) - len(x))
        # Get nlparci's CI outputs, compute half their widths to get back to delta,
        # and divide by tinv(.975,v), where v is the degrees of freedom, equal to the
        # number of observations minus the number of parameters.
        # alpha=0.05, n-p degrees of freedom

        sdfree = np.sqrt(len(residual)) * sefree
        cvfree = 100*sdfree / x #FIXME? x' transposee # coefficient of variation
        self.setFreeParameterValues(x,lb)
        
        # Assign values to fixed parameters
        x = self.getParameterValues # return all parameters
        cv = np.nan*np.ones( (1,len(x)) )
        cv[self.ISFREE > 0] = cvfree
        se = np.nan*np.ones( (1,len(x)) )
        se[self.ISFREE > 0] = sefree
        ci = np.nan*np.ones( (len(x),2) )
        ci[self.ISFREE > 0,1] = cifree[:,1]
        ci[self.ISFREE > 0,2] = cifree[:,2]
     
        return [x, residual, ci, se, cv]
    

    def errorFunction(self,x,exps,ERROR_FUNCTION,lb):
        '''
        # Try some parameters and evaluate the error
        # ERROR_FUNCTION 0 unscaled
        #                1 scaled by yhat, the simulated value
        #                2 scaled by y, the experimental value
        #                3 y and yhat on log scale
        #                4 scaled by stdev
        '''
        self.setFreeParameterValues(x,lb)
        F = []
        for exp in exps:
            # simulate the model for exp e
            Y = self.simulate(self, exp.dataVariables[1].times, exp)
            modelOutVar = self.getOutputVariables(self)
            for dv in exp.dataVariables:
                # Find simulated timeseries from modelOutVar that matches dv.name with dv= dv 
                # find the COLUMN in Y to consider
                dataVariableName = dv.name
                COLUMN = -1
                for k, outvar in enumerate(modelOutVar):
                    if dataVariableName == outvar:
                        COLUMN = k
                        break
                if COLUMN < 0:
                    continue
                N = len(dv.times)
                F1 = np.zeros((N,1))
                for jj in range(N):
                    if ~np.isnan(dv.values[jj]):
                        if ERROR_FUNCTION == 0:
                            F1[jj] = Y[jj,COLUMN] - dv.values[jj]
                        elif ERROR_FUNCTION == 1 and Y[jj,COLUMN] != 0.0:
                            F1[jj] = (Y[jj,COLUMN] - dv.values[jj])/Y[jj,COLUMN]
                        elif ERROR_FUNCTION == 2 and dv.values[jj] != 0.0:
                            F1[jj] = (Y[jj,COLUMN] - dv.values[jj])/dv.values[jj]
                        elif ERROR_FUNCTION == 3 and Y[jj,COLUMN] != 0.0 and dv.values[jj] != 0.0:
                            F1[jj] = np.log(Y[jj,COLUMN]) - np.log(dv.values[jj])
                        elif ERROR_FUNCTION == 4 and ~np.isnan(dv.stdev[jj]) and dv.stdev[jj] != 0.0:
                            F1[jj] = (Y[jj,COLUMN] - dv.values[jj])/dv.stdev[jj]
                        else:
                            F1[jj] = 0
                    else:
                        F1[jj] = 0
                F.append(F1)
        return F
    

    def fitLL(self, exps, lb, ub, ERROR_FUNCTION, mapVariableNameToSigmaPara, MaxFunEvals=None, MaxIter=None, TolFun=None, TolX=None):
        if len(lb) != len(self.Parameters):
            raise BadInputException('ResultChk:BadInput', 'The length of lb does not match number of parameters.')
        if len(ub) != len(self.Parameters):
            raise BadInputException('ResultChk:BadInput', 'The length of ub does not match number of parameters.')
        if ERROR_FUNCTION < 0 or ERROR_FUNCTION > 5:
            raise BadInputException('ResultChk:BadInput', 'The value of ERROR_FUNCTION is illegal.')

        #--- Initialization

        # Check for free parameters
        self.ISFREE = np.zeros( (1,len(lb)) )
        print(lb)
        print(type(lb))
        for i, lb_ in enumerate(lb):
            if lb_ != ub[i]:
                self.ISFREE[i] = 1
        lbnew = lb[self.ISFREE>0]
#FIXME ? non utilise       ubnew = ub[self.ISFREE>0]

        # Set up a map from parameter name to index
        mapParaNameToIndex = {} #matlab : containers.Map() + keys() et values()
        count = 1
        for i in enumerate(lb):
            if self.ISFREE[i]:
                mapParaNameToIndex[str(self.Parameters[i])] = count
                count += 1
        #mapParaNameToIndex.keys()
        #mapParaNameToIndex.values)
        
        MAXFUNEVALS = MaxFunEvals if (MaxFunEvals is not None) else 200*len(lbnew)
        MAXITER     = MaxIter     if (MaxIter     is not None) else 200*len(lbnew)
        TOLFUN      = TolFun      if (TolFun      is not None) else 1e-4
        TOLX        = TolX        if (TolX        is not None) else 1e-4
        
        # Harmonize time points in each experiment
        for exp in exps:
            exp.harmonizeTimes()
        
        # 1. Are there output variables from this model
        # for which we have experimental data to compare to? If not,
        # throw an error.
        modelOutVar = self.getOutputVariables()
        ok = 0
        for exp in exps:
            for dv in exp.dataVariables:
                dataVariableName = dv.name
                for ov in modelOutVar:
                    if dataVariableName == ok:
                        ok = 1
                        break
        if not ok:
            raise BadInputException('ResultChk:BadInput', 'There are no experimental data variables that match model outputs.')
        
        # 2. Fit data, error is calculated in errorFunction
        options = optimset('Display','off',\
            'MaxFunEvals',MAXFUNEVALS,'Maxiter',MAXITER,'TolFun',TOLFUN,'TolX',TOLX)
        x0 = self.getFreeParameterValues
        # verifier suntaxe fmin
        [x, fval, exitflag, output] = np.fmin(self.errorFunctionLL,x0,\
            options,self,exps,ERROR_FUNCTION,lb,mapVariableNameToSigmaPara,mapParaNameToIndex)
# scipy.optimize.fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None, initial_simplex=None)[source]
#Minimize a function using the dowhn ill simplex algorithm. This algorithm only uses function values, not derivatives or second derivatives.

        self.setFreeParameterValues(x,lb) # Assign values to fixed parameters
        x = self.getParameterValues # return all parameters

        return [x, fval, exitflag, output]

    
    def errorFunctionLL(self,x,exps,ERROR_FUNCTION,lb,mapVariableNameToSigmaPara,mapParaNameToIndex):
        '''
        # Try some parameters and evaluate the error
        # ERROR_FUNCTION 0 unscaled
        #                1 scaled by yhat, the simulated value
        #                2 scaled by y, the experimental value
        #                3 y and yhat on log scale
        #                4 scaled by stdev
        '''    
        self.setFreeParameterValues(x,lb)
        F = 0
        for exp in exps:
            # simulate the model for exp 
            Y = self.simulate(self, exp.dataVariables[1].times, exp)
            modelOutVar = self.getOutputVariables(self)
            for dv in exp.dataVariables:
                # Find simulated timeseries from modelOutVar that matches dv.name
                # find the COLUMN in Y to consider
                dataVariableName = dv.name
                COLUMN = -1
                for k, outvar in enumerate(modelOutVar):
                    if dataVariableName == outvar:
                        COLUMN = k
                        break
                if COLUMN < 0:
                    continue

                try:
                    sigma_name = mapVariableNameToSigmaPara(dataVariableName)
                    sigma_index = mapParaNameToIndex(sigma_name)
                except:
                    raise BadInputException('ResultChk:BadInput', 'DataVariable ', dataVariableName, \
                        ' does not match a SIGMA parameter, or the SIGMA parameter has no proper index')
                    mapVariableNameToSigmaPara.keys()
                    mapVariableNameToSigmaPara.values()
                    mapParaNameToIndex.keys()
                    mapParaNameToIndex.values()

                sigma2 = x(sigma_index)
                N = len(dv.times)
                for jj in range(N):
                    if ~np.isnan(dv.values[jj]):  #tilde~ marche avec numpy
                        if ERROR_FUNCTION==0:
                            eps2 = (-Y[jj,COLUMN] + dv.values[jj])^2
                            LL = -0.5*(eps2/sigma2 + np.log(2*np.pi*sigma2))
                            F = F-LL
                        elif ERROR_FUNCTION==1:
                            if Y[jj,COLUMN] != 0.0 and dv.values[jj] != 0.0:
                                eps2 = (-np.log(Y[jj,COLUMN]) + np.log(dv.values[jj]))^2
                                LL = -0.5*(eps2/sigma2 + np.log(2*np.pi*sigma2))
                                #print([e.name ' ' dv.name  ': y = '\
                                #    num2str(dv.values[jj])\
                                #    ': sigma2 = ' num2str(sigma2)\
                                #    ': LL = ' num2str(LL) ': F = ' num2str(F)])
                                #http://www.statlect.com/normal_distribution_maximum_likelihood.htm
                                F = F-LL
                        #elif ERROR_FUNCTION==5
                        #    eps2 = (-Y[jj,COLUMN] + dv.values[jj])^2
                        #    sigma2 = Y[jj,COLUMN]^2 *sigmaA2 + sigmaB2  # sigmaA2 proportional error term, sigmaB2 additive error term
                        #    LL = -0.5*( eps2/sigma2 + np.log(2*np.pi*sigma2) )
                        #    F = F-LL
                        else:
                            raise BadInputException('ResultChk:BadInput', 'No valid ERROR_FUNCTION')
        if self.verbose:
            print(F)  #print(sigmaA2)   #print(sigmaB2)

        return F

    def simulate(self, times, experiment):
        ''' Sequentially simulate the PK model closed form (if any), the PK or PD model
        ODE form (if any) and the PD model closed form (if any)
        The function handle of PK models is transmitted to PD model (if required)
        '''
        Y, Y1, Y2 = [], [], []
        fh = []
        ySubset = []
        outputVars, inputVars = [], []
        #--- PK model, closed form
        if len(self.FunctionPK)>0:  #~isempty()
            f = self.FunctionPK
            Y1 = f.simulate(times, experiment, ySubset, fh)
            Y2 = Y1 # if FunctionPKode is empty, forward to functionPD
            Y = [Y, Y1]
            outputVars = f.outputVariables
            fh = f #input for the PD model ODE form
        #--- PK model or PD model, ODE form
        if len(self.FunctionPKode)>0:  #~isempty()
            f = self.FunctionPKode
            inputVars = f.inputVariables
            ySubset = self.matrixSubset(Y1, outputVars, inputVars)
            Y2 = f.simulate(times, experiment, ySubset, fh)
            Y = [Y, Y2]
            outputVars = f.outputVariables
            fh = f #input for the PD model closed form
        #--- PD model, closed form
        if len(self.FunctionPD)>0:  #~isempty()
            f = self.FunctionPD
            inputVars = f.inputVariables
            ySubset = self.matrixSubset(Y2, outputVars, inputVars)
            Y3 = f.simulate(times, experiment, ySubset, fh)
            Y = [Y, Y3]
        
        return Y
    
    def getOutputVariables(self):
        ''' Returns all output variables (from PK closed form, PD closed form,
        and PK or PD ODE form) in a list'''
        outputVars = []
        if len(self.FunctionPK) > 0:  #~isempty
            outputVars.append(self.FunctionPK.outputVariables)
        if len(self.FunctionPKode) > 0:   #~isempty
            outputVars.append(self.FunctionPKode.outputVariables)
        if len(self.FunctionPD) > 0:    #~isempty
            outputVars.append(outputVars, self.FunctionPD.outputVariables)
        return outputVars
    
    
    def verifyFunctionParameterOrder(self):
        ''' Return True if the function and the model have the same parameters
        in the same order, False otherwise '''
        functions = [self.FunctionPK, self.FunctionPKode, self.FunctionPD]
        for fun in functions:
            if fun: #~isempty
                selectedParameters = []
                for par in self.Parameters:
                    for parname in fun.parameterNames:
                        if par == parname:
                            selectedParameters.append(par)

                #Retun False if not the same number of parameters
                if len(selectedParameters) != len(fun.parameterNames):
                    return False
                #...or if not in the same order
                for p1, p2 in zip(selectedParameters, fun.parameterNames):
                    if p1 != p2:
                        return False

        return True
    
    def verifyFunctionParameters(self):
        ''' Return True if the PK model and the PD model have compatible parameters,
        False otherwise '''
        # IF PK model is closed form, PD model can be closed form or ODE
        if self.FunctionPK: #~isempty
            # If PD model uses ODE slot
            if self.FunctionPKode: #~isempty
                if not self.FunctionPK.parametersAreCompatibleWith(self.FunctionPKode):
                    return False
            # If PD model is closed form
            if self.FunctionPD:
                if not self.FunctionPK.parametersAreCompatibleWith(self.FunctionPD):
                    return False
        # IF PK model is ODE, PD model can only be closed form
        if self.FunctionPKode: #~isempty
            if self.FunctionPD:  #~isempty
                if not self.FunctionPKode.parametersAreCompatibleWith(self.FunctionPD):
                    return False
        return True
    
    def verifyInputChain(self):
        result = True
        # Make sure input/output from the functions match.
        if self.FunctionPK is not None:   #~isempty()
            print('PK function is not None...',repr(self.FunctionPKode))
            if self.FunctionPKode is not None:
                result = self.FunctionPK.outputMatchesInput(self.FunctionPKode)
            elif self.FunctionPD is not None:
                result = self.FunctionPK.outputMatchesInput(self.FunctionPD)
            if not (result): 
                return result
        
        elif self.FunctionPKode is not None:  #~isempty()
            print('PK function is None...',repr(self.FunctionPKode))
            if self.FunctionPD is not None:
                result = self.FunctionPKode.outputMatchesInput(self.FunctionPD)

        print('verify inputchain',result)
        return result
    
    def matrixSubset(self, Y, outputVars, inputVars):
        ''' Returns a subset of the output matrix from a function.
         Y = Output from function.
         outputVars = Output variables from fun1
         inputVars = Input variables to fun2
        '''
        ySubset = []
        for ivar in inputVars:
            for j, ovar in enumerate(outputVars):
                if ivar == ovar:
                    ySubset = [ ySubset, Y[:,j] ] #FIXME utiliser extend?
        return ySubset
    
    def getParameterSubsetValues(self, parameterNames, parameterValues):
        ''' Returns a subset of parameterValues where parameterNames == self.Parameters
        '''
        subset = []
        for i, par in enumerate(self.Parameters):
            for name in parameterNames:
                if par == name:
                    subset.append(parameterValues[i])
        return subset

    def getInputParametersFromFunctions(self, functions):
        parameters = []
        for fun in functions:
            if fun: # ~isempty(f)
                for pName in fun.parameterNames:
                    if not (self.parameterExists(parameters, pName)):
                        parameters.append(pName)
        return parameters
    
    def getFreeParameterValues(self):
        ''' ??? '''
        p = self.getParameterValues()
        pFreeValues = p[self.ISFREE>0]
        return pFreeValues
    
    def getParameterValues(self):
        ''' TO DO '''
        functions = [self.FunctionPK, self.FunctionPKode, self.FunctionPD]
        parameters = []
        pValues = []
        for fun in functions:
            if fun: # ~isempty(f):
                for j, pName in enumerate(fun.parameterNames):
                    if not (self.parameterExists(parameters, pName)):
                        parameters.append(pName)  #FIXME ; on met a jour localment mais on ne renvoie que les val
                        pValues.append(fun.parameterValues[j])
        return pValues
    
    def setFreeParameterValues(self, kfree, lb):
        '''
        # k vector of free parameters
        '''
        # Make sure len is the same
        if len(kfree) != sum(self.ISFREE):
            raise BadInputException('ResultChk:BadInput', 'Incorrect number of free parameter values!')
        k = lb
        j = 1
        for i, dummy in enumerate(k):
           if self.ISFREE[i]:
              k[i] = kfree[j]
              j += 1
        self.setParameterValues(k)
    
    def setParameterValues(self, k):
        '''
        # k vector of parameters
        '''
        # Make sure len is the same
        if len(self.getParameterValues) != len(k):
            raise BadInputException('ResultChk:BadInput', 'Incorrect number of input parameter values! Please check parameter k for correct number of input values.')
        if len(self.FunctionPK) > 0 : #~isempty(
            self.FunctionPK.parameterValues =    self.getParameterSubsetValues(self.FunctionPK.parameterNames, k)
        if len(self.FunctionPKode) > 0: #~isempty(
            self.FunctionPKode.parameterValues = self.getParameterSubsetValues(self.FunctionPKode.parameterNames, k)
        if len(self.FunctionPD) > 0: #~isempty(
            self.FunctionPD.parameterValues =    self.getParameterSubsetValues(self.FunctionPD.parameterNames, k)
    
    def parameterExists(self, parameters, param):
        exists = False
        for par in parameters:
            if par == param:
                return True
        return exists

class Experiment(object):
    '''
     Experiment
     Information about an experiment
    
     Inputs:
     - name, String
     - dosing, Dosing, the dosing schedule
     - dataVariables, vector of DataVariable, the measured data variables
     - plotColor, character [r=red, b=blue etc., see 'help plot' for details), or
                vector of three numbers (Red/Blue/Green in interval [0,1],
                e.g., [1 0 0]=red, [0 0 1]=blue). Default blue.
    
     Methods:
     - display(), displays information about this experiment
     - harmonizeTimes(), ensures that time vectors for all data variables are the
     same. Inserts NaN when required. See example below.
     - plot(figureNumber), plots data in a figure indexed figureNumber
    
    '''
    def __init__(self, name, dosing, dataVariables = [], plotColor = None):  # dataVariables vector of time series data
        try:
            self.name = name
            if isinstance(dosing,Dosing):
                self.dosing = dosing
            else:
                raise BadInputException('ResultChk:BadInput', 'Input dosing should be an object of the class Dosing')
    
            for dv in dataVariables:
                if not isinstance(dv,DataVariable):
                    raise BadInputException('ResultChk:BadInput', 'Input dataVariables should be vector of objects of the class DataVariable')
            self.dataVariables = dataVariables
    
            if plotColor is not None: self.plotColor = plotColor
            else:                     self.plotColor = 'b'
        except:
            raise

    def display(self):
        print('Experiment ',self.name)
        self.dosing.display()
        for dv in self.dataVariables:
            dv.display

    def harmonizeTimes(self):
        # Make sure that each dataVariable has the same time vector
        t = []
        for dv in self.dataVariables:
            t = np.union1d(t, dv.times)

        for dv in self.dataVariables:
            dvT = np.array(dv.times)
            dvV = np.array(dv.values)
            dvS = np.array(dv.stdev)
            newValues = []
            newStdDevs = []
            for tps in t:
                idx = np.where(dvT == tps)[-1] #idx = find(dvT == t[j], 1, 'last') # find the 1 last element satisfying the condition
                if idx.size>0:
                    newValues = [newValues,dvV[idx]]
                    newStdDevs = [newStdDevs,dvS[idx]]
                else: #if isempty[idx] #pour array: A.size == 0
                    newValues = [newValues, np.nan]
                    newStdDevs = [newStdDevs, np.nan]
            dv.times = t
            dv.values = newValues
            dv.stdev = newStdDevs

    def plot(self, ifigure):
        #plt.figure(ifigure)
        count = 1
        for dv in self.dataVariables:
            for j in [1,2]:
                plt.subplot(len(self.dataVariables),2,count)
                count += 1
                if j == 1:
                    plt.errorbar(dv.times, dv.values, dv.stdev, [self.plotColor, 'o-'])
                    plt.hold('on')
                else:
                    plt.semilogy(dv.times, dv.values, [self.plotColor, 'o-'])
                    plt.hold('on')

        count = 1
        for dv in self.dataVariables:
            for j in [1,2]:
                plt.subplot(len(self.dataVariables),2,count)
                count += 1
                plt.setp(plt.gca,'FontName','Helvetica')
                plt.setp(plt.gca,'FontSize',10)
                hXLabel = plt.xlabel('Time (' + str(self.dv.timeUnit) + ')')
                hYLabel = plt.ylabel(str(dv.name) + '(' + str(dv.valueUnit) + ')')
                hL = plt.legend(str(dv.name))
                plt.setp([hL,hXLabel,hYLabel])
# =============================================================================
#                 plt.setp([hL,hXLabel, hYLabel],
#                      'FontName','Helvetica',
#                      'Interpreter','tex',
#                      'FontSize',10)
#                 plt.setp(plt.gca,  \
#                      'Box'         , 'off'     ,
#                      'TickDir'     , 'out'     ,
#                      'YGrid'       , 'off'     ,
#                      'YGrid'       , 'off'     ,
#                      'XMinorTick'  , 'on'      ,
#                      'YMinorTick'  , 'on'      ,
#                      'XColor'      , [.3, .3, .3],
#                      'YColor'      , [.3, .3, .3],
#                      'LineWidth'   , 1)
#                 plt.legend('boxoff')
# =============================================================================

class Dosing(object):
    '''
     Dosing : Information about a dosing schedule
        (Note: for infusion, values are given as rates)

     Inputs:
      - admin, String, ev=extravascular, iv=intravenous, inf=infusion.
      - times, row vector of numbers, time points of doses.
      - values, row vector of numbers, dose values. NB: for infusion as rate.
      - duration, row vector of numbers, duration of infusion, empty for ev and iv.
      - timeUnit, String, unit of time points.
      - valueUnit, String, unit of values.

     Methods:
      - display(), displays information about this dosing
    '''
    def __init__(self, admin='ev', times=[], values=[], duration=[], timeUnit='min', valueUnit='umol'):
                    #admin: ev=extravascular, iv=intravenous, inf=infusion

        try:
            #Input controls
            if admin not in ['ev','iv','inf']:
                raise BadInputException('ResultChk:BadInput', 'admin should be ev, iv or inf')
    
            if times:
                times = np.array(times) #.reshape(1,len(times)) # mt, nt = times.shape
                if len(times.shape) != 1:
                    raise BadInputException('ResultChk:BadInput', 'Vector of times should be a row vector')
    
            if values:
                values = np.array(values) #.reshape(1,len(times))
                if len(values.shape) != 1:
                    raise BadInputException('ResultChk:BadInput', 'Vector of values should be a row vector')
    
            for i in range(len(times)-1):
                if times[i] >= times[i+1]:
                    raise BadInputException('ResultChk:BadInput', 'Elements in times must be strictly increasing.')
    
            if len(times) != len(values):
                raise BadInputException('ResultChk:BadInput', 'The len of values does not match the len of times.')
    
            if admin == 'inf':
                if duration:
                    duration = np.array(duration)
                    if len(duration.shape) != 1:
                        raise BadInputException('ResultChk:BadInput', 'Vector of duration should be a row vector')
                    if len(times) != len(duration):
                        raise BadInputException('ResultChk:BadInput', 'The len of duration does not match the len of times.')
                    for i in range(len(times)-1):
                        if times[i+1] < (times[i]+duration[i]):
                            raise BadInputException('ResultChk:BadInput', 'Infusions are overlapping.')
    
            #Creation of the Dosing object
            self.admin = admin
            self.times = times
            self.values = values
            self.duration = duration
            self.timeUnit = timeUnit
            self.valueUnit = valueUnit

        except BadInputException:
            raise

    def display(self):
        print(self.__str__())

    def __str__(self):
        try:
            text= (' Dosing\n')
            text += '  Admin: ' + str(self.admin) + '\n'
            text += '  Time: ' + str(self.times) + ' ' + str(self.timeUnit) + '\n'
            text += '  Values: ' + str(self.values) + ' ' + str(self.valueUnit) + '\n'
            if len(self.duration) > 0:
                text += '  Duration: ' + str(self.duration) + ' ' + str(self.timeUnit) + '\n'
        except:
            raise BadInputException('ResultChk:BadInput', 'Failed to convert some values')
        return text

class DataVariable(object):
    ''' DataVariable
    Information about measuered data of one variable

    Inputs:
    - name, String, name of data variable
    - times, row vector of numbers, sampling time points
    - values, row vector of numbers, sampled values
    - stdev, row vector of numbers, standard devation of sampled values
           if unknown, input [].
    - timeUnit, String, unit of time points
    - valueUnit, String, unit of values

    Methods:
    - display(): displays information about this data variable
    - y=AUC(t0,tend,scale,doplot): returns the AUC_t0_tend of the variable
           t0, start time, times[1] is default
           tend, end time, times(end) is default
           scale for interpolation, 'lin' or 'log', 'log' is default
                  For lin scale, negative values are assigned zero
           doplot, 0 = no plot, default
                   1 = plot data points and interpolating curve to
                       the current axes of the current figure
    '''

    def __init__(self, name, times=[], values=[], stdev=[], timeUnit='min', valueUnit='uM'):
        try:
            #Input checks
            #if times:
            times = np.array(times)# .reshape(1,times.size) #mt, nt = times.shape
            if len(times.shape) != 1:
                raise BadInputException('ResultChk:BadInput', 'Vector of times should be a row vector')
    
            #if values:
            values = np.array(values) #.reshape(1,times.size)
            if len(values.shape) != 1:
                raise BadInputException('ResultChk:BadInput', 'Vector of values should be a row vector')
    
            #if stdev:
            if not(stdev): #if len(stdev) == 0:
                stdev = np.nan * np.ones_like(times)
            stdev = np.array(stdev) #.reshape(1,times.size)
            if len(values.shape) != 1:
                raise BadInputException('ResultChk:BadInput', 'Vector of stdev should be a row vector')
    
#            if not (times.shape == values.shape and  values.shape == stdev.shape): #if len(times) != len(values) or len(times) != len(stdev):
            if not (times.shape == values.shape == stdev.shape): #if len(times) != len(values) or len(times) != len(stdev):
                print (times.shape,values.shape, stdev.shape)
                raise BadInputException('ResultChk:BadInput', 'Vectors of times, values and stdev must have the same shape/len.')
            for i in range(len(times)-1):
                if times[i] >= times[i+1]:
                    raise BadInputException('ResultChk:BadInput', 'ResultChk:BadInput', 'Elements in times must be increasing.')
    
            self.name = name
            self.times = times
            self.values = values
            self.stdev = stdev
            self.timeUnit = timeUnit
            self.valueUnit = valueUnit

        except:
            raise

    def display(self):
        print(self.__str__())

    def __str__(self):
        try:
            text= (' Data variable: ' + self.name + '\n')
            text += '  Time: ' + str(self.times) + ' ' + str(self.timeUnit) + '\n'
            text += '  Values: ' + str(self.values) + ' ' + str(self.valueUnit) + '\n'
            text += '  St devs: ' + str(self.stdev) + ' ' + str(self.valueUnit) + '\n'
        except:
            raise BadInputException('ResultChk:BadInput', 'ResultChk:BadInput', 'Failed to convert some values')
        return text

    def AUC(self, t0=None, tmax=None, scale_=None, doplot_=None):

        T0 = t0 if (t0 is not None) else self.times[0]
        TMAX = tmax if (tmax is not None) else self.times[-1]
        if scale_ is not None:
            scale = scale_
            if scale not in ['lin','log']:
                raise BadInputException('ResultChk:BadInput', 'ResultChk:BadInput', "Scale should be 'lin' or 'log'")
        else:
            scale = 'log'
        doplot = doplot_ if (doplot_ is not None) else 0

        N = 10000
        xq = np.linspace(T0,TMAX,N+1) # ou xq = np.arange(T0,TMAX+(TMAX-T0)/N,(TMAX-T0)/N)

        if scale == 'lin':
            Y = self.values
        elif scale == 'log':
            Y = np.log(self.values)
        elif scale == 'linlog':
            Y = -1  #FIXME ajouter linlog

        f = interpolate.interp1d(self.times, Y, kind='linear', fill_value='extrapolate')
        vq = f(xq)
        if scale == 'lin':
            vq[ vq<0 ] = 0
        elif scale == 'log':
            vq = np.exp(f(xq))
        elif scale == 'linlog':
            theauc = -1   #FIXME ajouter linlog

        #---AUC : calcul de la concentration moyenne puis on divise par l'intervalle de temps
        a = (vq[0]/2 + vq[-1]/2 + sum(vq[1:-1])) / (N-1)
        theauc = a*(TMAX-T0)

        if doplot > 0:
            settings = {}
            settings["fontsize"] = 12
            settings["lineWidth"] = 0.7
            settings["markerSize"] = 5
            if scale == 'lin':
                from nicePlot import nicePlot
                nicePlot(xq, vq, 'k', '', '-', 'interpoled', settings)
                plt.hold(True)
                nicePlot(self.times, Y, 'k', 'o', '', 'real', settings)
                plt.hold(False)
            else:
                pass
                #niceSemilogy(xq,vq,'k','','-',settings)
                #niceSemilogy(self.times,y,'k','o','',settings)

        print('AUC: '+str(round(theauc,1))+ ' ' + self.valueUnit + '.' + self.timeUnit + '\n')
        return theauc


if __name__ == "__main__": 
    #           EXAMPLES
    
    #======= Experiment
    # Example: Extravascular dosing at t=0h and t=24h with 10 and
    # 15umol/kg, respectively. Two measured variables, C (uM) and R (%),
    # the plot color is red (r or [1 0 0]).
    dose = Dosing('ev',[0, 24],[10, 15],[], 'h', 'umol/kg')
    t = [1, 10, 60]
    c = [0.2, 0.3, 0.04]
    c_std = [0.02, 0.06, 0.008]
    dv1 = DataVariable('C', t, c, c_std, 'min', 'uM')
    t = [1, 20, 60]
    r = [0.1, 0.8, 0.7]
    r_std = [0.01, 0.08, 0.07]
    dv2 = DataVariable('R', t, r, r_std, 'min', '#')
    exper = Experiment('test', dose, [dv1, dv2], 'r')
    exper.display()
    #---- 1st call to e.display() gives:
    # Experiment: test
    #  Data variable: C
    #   Time: 1  10  60 min
    #   Values: 0.2         0.3         0.04 uM
    #   St devs: 0.02        0.06         0.008 uM
    #  Data variable: R
    #   Time: 1  20  60 min
    #   Values: 0.1         0.8         0.7 %
    #   St devs: 0.01         0.08        0.07 %
    #
    
    #exper.plot(1) # Plot data in Figure 1
    exper.harmonizeTimes()
    exper.display()
    #---- 2nd call to e.display() gives:
    # Experiment: test
    #  Data variable: C
    #   Time: 1  10  20  60 min
    #   Values: 0.2         0.3         NaN        0.04 uM
    #   St devs: 0.02        0.06         NaN       0.008 uM
    #  Data variable: R
    #   Time: 1  10  20  60 min
    #   Values: 0.1         NaN         0.8         0.7 %
    #   St devs: 0.01         NaN        0.08        0.07 %
    
    #======= Dosing
    # Example 1: Extravascular dosing at t=0h and t=24h with 10 and 15umol/kg, respectively:
    d = Dosing('ev',[0, 24],[10, 15],[], 'h', 'umol/kg')
    print(d)
    
    # Example 2: Intravenous dosing at t=0min, t=10min and t=20min with 5, 10 and 15umol/kg, respectively:
    d = Dosing('iv',[0, 10, 20],[5, 10, 15],[], 'min', 'umol/kg')
    print(d)
    
    # Example 3: Infusion at t=0min with duration 10 minutes and at t=60min with duration
    # 15 minutes with 5 and 10umol/kg/min, respectively (infusion is given as rate):
    d = Dosing('inf',[0, 60],[5, 10],[10, 15], 'min', 'umol/kg/min')
    print(d)
    
    #======= DataVariables
    # Example 1:
    t = [1, 3, 8, 24]
    c = [0.2, 0.3, 0.1, 0.04]
    c_std = [0.02, 0.06, 0.03, 0.0008]
    dv = DataVariable('C', t, c, c_std, 'h', 'uM')
    dv.display()
    dv.AUC(0,24,'log',1)
    
    # Example 2, standard deviation unknown:
    t = [0, 1, 2, 3, 4, 6, 8, 10]
    c = [100, 71, 50, 35, 25, 12, 6.2, 3.1]
    c_std = []
    dv = DataVariable('C', t, c, c_std, 'min', 'uM')
    dv.display()
    dv.AUC(0,24,'lin')
    

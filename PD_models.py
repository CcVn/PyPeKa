# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 16:37:02 2018

@author: cvn_or
"""
import math
import PKObjects as pko

class PD_linear_fun(pko.SimFunction):
    ''' PD_linear_fun (Subclass of SimFunction)
        R = a *C(t)
        
     Inputs:
       - pValues, vector/scalar of the parameter [k]
       - inputVariableName, string (optional input, default 'C')
       - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
       - parameterNames, {'k'} 
       - parameterValues, assigned from input pValues
       - inputVariables, assigned from input inputVariableName
       - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
         times - vector of time points to output
         experiment - the experiment to simulate (wrt dosing)
         yInput input variable values at time points in times, size |times|x|inputVariables|
         defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)
         
     Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''
        
    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['a']
#FIXME: voir si on peut transferer dans SimFunction
#        super().__init__(pValues, inputVariableName, outputVariableName)
        self.inputVariables = inputVariableName
        self.outputVariables = outputVariableName
        if pValues is None:
            raise Exception('BadInput: At least one input argument is required.')
        self.parameterValues = pValues
        if len(pValues) != len(self.parameterNames):
            raise Exception('BadInput: Number of parameters incorrect!')
    
    def simulate(self, times, experiment, yInput, defObj):
        a, *others = self.parameterValues;
        Y = a * yInput;
        return Y


class PD_loglinear_fun(PD_linear_fun):
    ''' PD_loglinear_fun (Subclass of SimFunction)
        R = m * log (a + C(t))
        
     Inputs:
       - parameterNames, {'m','a'}
       - inputVariableName, string (optional input, default 'C')
       - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'k'} 
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

     Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['m','a']
        super().__init__(pValues, inputVariableName, outputVariableName)

    def simulate(self, times, experiment, yInput, defObj):
        m, a, *others = self.parameterValues
        Y = m * math.log(a + yInput)
        return Y


class PD_powerOf_fun(PD_linear_fun):
    ''' PD_powerOf_fun (Subclass of PD_linear_fun)
        R = a*C^b, a and b parameters
        
     Inputs:
      - parameterNames, {'a','b'}
      - inputVariableName, string (optional input, default 'C')
      - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'a', 'b'} 
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

      Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['m','a']
        super().__init__(pValues, inputVariableName, outputVariableName)        

    def simulate(self, times, experiment, yInput, defObj):
        a, b, *others = self.parameterValues
        Y = a * yInput**b
        return Y


class PD_Emax_fun(PD_linear_fun):
    ''' PD_Emax_fun (Subclass of PD_linear_fun)
        R = Emax *C(t) / (EC50 + C(t))
        
     Inputs:
      - pValues, vector of the two parameters [Emax, EC50]
      - inputVariableName, string (optional input, default 'C')
      - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'Emax','IC50'} 
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

     Return:s Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['Emax','EC50']
        super().__init__(pValues, inputVariableName, outputVariableName)        
    
    def simulate(self, times, experiment, yInput, defObj):
        Emax, EC50 = self.parameterValues
        Y = Emax*yInput / (EC50 + yInput)
        return Y


class PD_Emax_hill_fun(PD_linear_fun):
    ''' PD_Emax_hill_fun (Subclass of PD_linear_fun)
        R = Emax *C(t)^h / (EC50^h + C(t)^h)
        
     Inputs:
      - pValues, vector of the two parameters [Emax, EC50, h]
      - inputVariableName, string (optional input, default 'C')
      - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'Emax','EC50','h'} 
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, functionObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - functionObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

     Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['Emax','EC50','h']
        super().__init__(pValues, inputVariableName, outputVariableName)        
     
    def simulate(self, times, experiment, yInput, functionObj):
        Emax, EC50, h, *others = self.parameterValues
        Y = Emax * yInput**h / (EC50**h + yInput**h)
        return Y

        
class PD_Imax_fun(PD_linear_fun):
    ''' PD_Imax_fun (Subclass of PD_linear_fun)
        R = 1 - Imax *C(t) / (IC50 + C(t))
        
     Inputs:
      - pValues, vector of the two parameters [Imax, IC50]
      - inputVariableName, string (optional input, default 'C')
      - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'Imax','IC50'}
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

     Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['Imax','IC50']
        super().__init__(pValues, inputVariableName, outputVariableName)        

    def simulate(self, times, experiment, yInput, defObj):
        Imax, IC50 = self.parameterValues

        Y= 1 - Imax*yInput/(IC50 + yInput)
        return Y


class PD_Imax_hill_fun(PD_linear_fun):
    ''' PD_Imax_hill_fun (Subclass of PD_linear_fun)
        R = 1 - Imax *C(t)^h / (IC50^h + C(t)^h)
        
     Inputs:
      - pValues, vector of the two parameters [Imax, IC50, h]
      - inputVariableName, string (optional input, default 'C')
      - outputVariableName, string (optional input, default 'R')
    
     Fields defined in the superclass SimFunction
      - parameterNames, {'Imax','IC50','h'}
      - parameterValues, assigned from input pValues
      - inputVariables, assigned from input inputVariableName
      - outputVariables, assigned from input outputVariableName
    
     Method:
     Y = simulate(times, experiment, yInput, defObj), simulates
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput input variable values at time points in times, size |times|x|inputVariables|
      - defObj from which inputVariables can be calculated at arbitrary time-points (not implemented for this class)

     Returns: Y matrix of output variable data, size |times|x|outputVariables|
    '''

    def __init__(self, pValues= None, inputVariableName= ['C'], outputVariableName= ['R']):
        self.parameterNames = ['Imax','IC50']
        super().__init__(pValues, inputVariableName, outputVariableName)        
    
    def simulate(self, times, experiment, yInput, defObj):
        Imax, IC50, hill = self.parameterValues
        
        Y= Imax*yInput**hill / (IC50**hill + yInput**hill)   # 1 - ???
        return Y

#=====================================================================================
#  Examples 1 and 2 des differents fichiers
#----------------------
if __name__ == "__main__": 
    import numpy as np
    import PK_123comp_models_ODE as pkode
    import PK_123comp_models_fun as pkm
    import matplotlib as plt
    
    timeUnit='min'
    doseUnit='umol/kg'
    d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
    experiment = pko.Experiment('Ex', d, [])
    times = np.linspace(0,49,1)
    const = 0
    for model in ["Linear", "logLinear", "Emax", "Emax_hill", "Imax", "Imax_hill", "powerOf" ]:
        for effectComp in [0,1]: # Example 1: without PD; Ex 2 with PD
            if effectComp==0:
                sf1 = pkode.PK_1comp_Ka_Cl_ODE(const, effectComp, [0.86, 1.6, 0.84], ['G','Cu'])
                if model=="Linear":      sf2 = PD_linear_fun([5],'Cu','E')
                elif model=="logLinear": sf2 = PD_loglinear_fun([5, 1],'Cu','E')
                elif model=="Emax":      sf2 = PD_Emax_fun([1, 0.5],'Cu','E')
                elif model=="Emax_hill": sf2 = PD_Emax_hill_fun([1, 0.5, 3],'Cu','E')
                elif model=="Imax":      sf2 = PD_Imax_fun([1, 0.5],'Cu','E')
                elif model=="Imax_hill": sf2 = PD_Imax_hill_fun([1, 0.5, 3],'Cu','E')
                elif model=="powerOf":   sf2 = PD_powerOf_fun([1, 0.5])
                model = pko.Model(None, sf1, sf2)  # (sf1, None, sf2) dans powerOf
            else:
                sf1 = pkm.PK_1comp_Ka_Cl_fun(const, effectComp, [0.86, 1.6, 0.84, 0.1])
                if model=="Linear":      sf2 = PD_linear_fun([5], 'Ce')
                elif model=="logLinear": sf2 = PD_loglinear_fun([5, 1],'Ce')
                elif model=="Emax":      sf2 = PD_Emax_fun([1, 0.5], 'Ce')
                elif model=="Emax_hill": sf2 = PD_Emax_hill_fun([1, 0.5, 3], 'Ce')
                elif model=="Imax":      sf2 = PD_Imax_fun([1, 0.5],'Ce')
                elif model=="Imax_hill": sf2 = PD_Imax_hill_fun([1, 0.5, 3],'Ce')
                elif model=="powerOf":   sf2 = PD_powerOf_fun([1, 0.5],'Ce')
                model = pko.Model(sf1, None, sf2)

            Y = model.simulate(times,experiment)
            plt.plot(times,Y[:,0],'b')
            plt.hold(True)
            plt.plot(times,Y[:,1],'m')
            plt.plot(times,Y[:,2],'r')  #pas pour powerOf
            plt.legend(model.getOutputVariables())
            plt.xlabel(timeUnit)
            plt.hold(False)

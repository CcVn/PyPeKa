#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PK models with 1/2/3 compartment using ODE

@author: CcVN
@contact: cedric.vinson2@gmail.com
@since: 2018/09/06
@status: alpha
"""

import numpy as np
import matplotlib.pyplot as plt
import PKObjects as pko
import PK_123comp_models_fun as pkm
from scipy import integrate as integ

__author__ = "CcVN"
__contact__ = "cedric.vinson2@gmail.com"
__date__ = "2018"
__version__ = "0.0.1"

#[t,y,te,ye,ie] = ode15s(odefun,tspan,y0,options)
#[t,y] = ode15s(odefun,tspan,y0), where tspan = [t0 tf], integrates the system of differential equations
# y'=f(t,y) from t0 to tf with initial conditions y0. 
# Each row in the solution array y corresponds to a value returned in column vector t.
#t — Evaluation points (column vector): Evaluation points, returned as a column vector.
#If tspan contains two elements, [t0 tf], then t contains the internal evaluation points used to perform the integration.
#If tspan contains more than two elements, then t is the same as tspan.
#y — Solutions (rray)Solutions, returned as an array. Each row in y corresponds to the solution
# at the value returned in the corresponding row of t.


class PD_models_ode(pko.SimFunction):
    '''
    # Parent class for PD models with ODE
    #
    # Inputs:
    #  - pValues, vector of the PD parameters
            [R0, kout, Emax, EC50, n] (stimofLoss & stimofProd)
            [R0, kout, Imax, IC50, n] (inhofLoss & inhofProd)
            [KD, koff, Emax] (receptor occupancy)
    #  - inputVariableName, string (optional input, default 'C')
    #  - outputVariableName, string (optional input, default 'R')
    #
    # Fields defined in the superclass SimFunction
    #  - parameterNames:
            {'R0','kout','Emax','EC50','n'} (stimofLoss & stimofProd)
            {'R0','kout','Imax','IC50','n'} (inhofLoss & inhofProd)
            {'KD','koff','Emax'} (receptor occupancy)
    #  - parameterValues, assigned from input pValues
    #  - inputVariables, assigned from input inputVariableName
    #  - outputVariables, assigned from input outputVariableName
    #
    # Function:
    #  - Y=obj.simulate(times, experiment, yInput, functionObj), simulates object
    #  - times - vector of time points to output
    #  - experiment - the experiment to simulate (wrt dosing)
    #  - yInput input variable values at time points in times, size |times|x|inputVariables| (not implemented for this class)
    #  - functionObj from which inputVariables can be calculated at arbitrary time-points
    #  - returns Y matrix of output variable data, size |times|x|outputVariables|

    '''
    INDEX_FUNCTIONOBJ = 0

    def __init__(self, pValues=None, inputVariableName=['C'], outputVariableName=['R']):
        try:
            if pValues is None:
                raise pko.BadInputException(' At least one input argument is required.')
            else:
                self.parameterValues = pValues

            if inputVariableName:
                self.inputVariables = inputVariableName
            else:
                self.inputVariables = ['C']
            if outputVariableName:
                self.outputVariableName = outputVariableName
            else:
                self.outputVariables = ['R']

        except pko.BadInputException:
            print('Bad Input')

    def __str__(self):
        message = 'Creating a ' + str(type(self).__name__) + ' type-like class with the following parameters/values\n'
        for i, n in enumerate(self.parameterNames):
            message += '\t-' + n + ' : ' + str(self.parameterValues[i]) + '\n'
        return message
        
    def simulate(self, obj, times, experiment, yInput, functionObj):
        '''
        Assume times to be a vector
        '''
        times = np.array(times)
        if len(times.shape)==2 and times.shape[0]>1:
            times = times.T

        # Find the index of outputVariables of functionObj that
        # matches the inputVariable of this class, i.e., 'C' or 'Ce'
        obj.INDEX_FUNCTIONOBJ=0;
        for i in enumerate(functionObj.outputVariables):
            if functionObj.outputVariables[i] == self.inputVariables[0]:
                self.INDEX_FUNCTIONOBJ = i
        if self.INDEX_FUNCTIONOBJ == 0:
            raise pko.BadInputException('ResultChk:BadInput'+str('inputVariableName '+ \
                    self.inputVariables[0]+' not found in functionObj.outputVariables'))

        kn = self.parameterValues
        #Conversion in tuple for odeint; test if iterable in case of a single int
        try:
            iter(kn)
            knt = tuple(kn)
        except TypeError:
            kn = [kn]
            knt = tuple(kn)

        E = []
        N = 1
        init = np.zeros(1,N)
        init[0] = self.parameterValues(1) # assign R0
        # @@@@@@ PAS DANS MODELE RECEPTOR OCC???


        # Before the first dose
        #-----------------------
        if times[0] < experiment.dosing.times[0]:
            #select times before the first dose (and add the time of first dose if needed)
            t = times[ times < experiment.dosing.times[0] ]
            t = [ t, experiment.dosing.times[0]]
            memory = experiment.dosing.values[0] # Dose 1
            experiment.dosing.values[0] = 0.0 # temporarily remove the 1st dose

            if self.verbose: print('Before first dose, ', len(t), ' time points.')

            if experiment.dosing.admin in ['ev','iv']:
                #[Tdummy,Etem] = ode15s(@myodes,t,init,[],kn,experiment,functionObj,obj.INDEX_FUNCTIONOBJ)
                Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )

                #@@@@@@@@@@@@@@@@@@ GERER INDEX_FUNCTIONOBJ : ca sert a quoi?

            experiment.dosing.values[0] = memory # Dose 1 reintegrated

            E = Etem[0:len(t),:] #E.r_[ Etem[0:(len(t)-1),:] ] 
            init = Etem[-1,:]

        # Dose 1 to n_doses-1
        #---------------------
        for idx in range(len(experiment.dosing.times)-1):
            
            t = times[ times < experiment.dosing.times[idx+1] ]
            t = t[ t >= experiment.dosing.times[idx] ]

            #add the dosing time if it doesn't belong to the initial times vector
            addFront = 0
            if len(t)==0 or t[0] != experiment.dosing.times[idx]:
                np.insert(t, 0, experiment.dosing.times[idx])  #insert time at beginning of t
                addFront = 1

            #add the time of the next dose at the end of t, and deactivate it temporarily
            np.append(t, experiment.dosing.times[idx+1])
            memory = experiment.dosing.values[idx+1]
            experiment.dosing.values[idx+1] = 0.0 # temporarily remove

            if self.verbose: print('[Dose 1 to n_doses-1]: dose ', str(idx+1), ' at ', experiment.dosing.times[idx+1], 'h.', len(t), ' time points')

            #FIXME Etem=[]
            #           Etem = np.empty_like(Etem)  #FIXME Etem=[]

            #-- INTEGRATE: [[Tdummy,Etem] = ode15s(@myodes,t,init,[],kn,experiment,functionObj,obj.INDEX_FUNCTIONOBJ)
            Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )

                #@@@@@@@@@@@@@@@@@@ GERER INDEX_FUNCTIONOBJ : ca sert a quoi?

            experiment.dosing.values[idx+1] = memory # reinsert next dose
            if len(t)==2: Etem = Etem[[0,-1], :] # garder ligne 1 et finale
            Etem_slice = Etem[ addFront:len(t), : ] #tout sauf 1er si addFront et dernier (exclu en Python contrairement a Matlab)
            if len(E)>0:
                E = np.vstack( (E, Etem_slice) )
            else:
                E= Etem_slice  # si pas initialise

            init = Etem[-1,:] # last line, Matlab: Etem(length(t),:)

        # Dose n_doses
        #--------------
        t = times[times >= experiment.dosing.times[-1]] #-1 replaces idx = len(experiment.dosing.times)
        #Add the dosing time if it doesn't belong to the initial times vector
        addFront = 0
        if len(t) == 0 or t[0] != experiment.dosing.times[-1] :
            np.insert(t, 0, experiment.dosing.times[-1])  #t = [ experiment.dosing.times[-1], t]
            addFront = 1
        if self.verbose: print('Last dose at t=',experiment.dosing.times[-1],'. ', len(t), ' time points')

        if len(t) == 1:
            print('Error, t has len 1')   #FIXME mettre Exception
            Etem = init
        else:
            Etem = []
            #[Tdummy,Etem] = ode15s(@myodes,t,init,[],kn,experiment,functionObj,obj.INDEX_FUNCTIONOBJ)
            Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )
            if len(t) == 2:
                Etem = Etem[[0,-1], :] # garder ligne 1 et finale

        E_slice = Etem[ (addFront):len(t), : ] #tout sauf 1er si addFront
        E = np.vstack( (E, E_slice) ) # ou E = np.r_[ E, E_slice ]

        if self.verbose: print('output Y:\n',E)
        return E

    def myodes(t,y,k,experiment,functionObj,INDEX_FUNCTIONOBJ):
        '''Fake method for parent class'''
        pass

class PD_turnover_inhOfLoss_ode(PD_models_ode):
    ''' PD_turnover_inhOfLoss_ode (Subclass of SimFunction)
     A turnover model where the loss rate is inhibited.
     dR(t)/dt = kin - kout*(1-(Imax*C(t)^n)/(IC50^n+C(t)^n))*R(t),
     where kin=kout*R0,
     based on steady-state solution R'(t)=0=kin-kout*R(0) for C(0)=0

    pValues, vector of the five parameters [R0, kout, Imax, IC50, n]
    '''
    def __init__(self, pValues=None, inputVariableName=['C'], outputVariableName=['R']):
        try:
            super().__init__(pValues, inputVariableName, outputVariableName)
            #old syntax : super(PD_turnover_inhOfLoss_ode, self).__init__(pValues, inputVariableName, outputVariableName)

            self.parameterNames = ['R0','kout','Imax','IC50','n'] # inhofProd/inhOfLoss
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException(' Number of parameters incorrect!')
        except pko.BadInputException:
            print('Bad Input')

    def myodes(self,t,y,k,experiment,functionObj,INDEX_FUNCTIONOBJ):
        Csim = self.simulate(functionObj,t,experiment,[],[])
        C = Csim(INDEX_FUNCTIONOBJ)
        R0, kout, Imax, IC50, n = k
        kin = kout*R0
        dy = np.zeros(1,1)
        dy[0] = kin - kout*(1 - (Imax*C**n) / (IC50**n+C**n)) *y[0] # Inhibition of loss
        return dy

class PD_turnover_inhOfProd_ode(PD_turnover_inhOfLoss_ode):
    ''' PD_turnover_inhOfProd_ode
    A turnover model where the production rate is inhibited.
    dR(t)/dt = kin*(1-(Imax*C(t)^n)/(IC50^n+C(t)^n)) - kout*R(t),
    where kin=kout*R0,
    based on steady-state solution R'(t)=0=kin-kout*R(0) for C(0)=0

    - pValues, vector of the five parameters [R0, kout, Imax, IC50, n]
    - parameterNames, {'R0','kout','Imax','IC50','n'}
    '''
    def myodes(self,t,y,k,experiment,functionObj,INDEX_FUNCTIONOBJ):
        Csim = self.simulate(functionObj,t,experiment,[],[])
        C = Csim(INDEX_FUNCTIONOBJ)
        R0, kout, Imax, IC50, n = k
        kin = kout*R0
        dy = np.zeros(1,1)
        dy[0] = kin *(1 - (Imax*C**n)/(IC50**n+C**n)) - kout*y[0] #Inhibition of prod
        return dy

class PD_turnover_stimOfLoss_ode(PD_models_ode):
    '''PD_turnover_stimOfLoss_ode
    A turnover model where the loss rate is stimulated.
    dR(t)/dt = kin - kout*(1+(Emax*C(t)^n)/(EC50^n+C(t)^n))*R(t),
    where kin=kout*R0,
    based on steady-state solution R'(t)=0=kin-kout*R(0) for C(0)=0

    - pValues, vector of the five parameters [R0, kout, Emax, EC50, n]
    - parameterNames, {'R0','kout','Emax','EC50','n'}
    '''
    def __init__(self, pValues=None, inputVariableName=['C'], outputVariableName=['R']):
        try:
            super().__init__(pValues, inputVariableName, outputVariableName)
            #old syntax : super(PD_turnover_inhOfLoss_ode, self).__init__(pValues, inputVariableName, outputVariableName)

            self.parameterNames = ['R0','kout','Emax','EC50','n'] # inhofProd/inhOfLoss
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException(' Number of parameters incorrect!')
        except pko.BadInputException:
            print('Bad Input')

    def myodes(self, t, y, k, experiment, functionObj, INDEX_FUNCTIONOBJ):
        Csim = self.simulate(functionObj,t,experiment,[],[])
        C = Csim(INDEX_FUNCTIONOBJ)
        R0, kout, Emax, EC50, n = k
        kin = kout*R0
        dy = np.zeros(1,1)
        dy[0] = kin - kout*(1 + (Emax*C**n) / (EC50**n + C**n)) *y[0] # Stimulation of loss
        return dy

class PD_turnover_stimOfProd_ode(PD_turnover_stimOfLoss_ode):
    '''PPD_turnover_stimOfProd_ode
    A turnover model where the production rate is stimulated.
    dR(t)/dt = kin*(1+(Emax*C(t)^n)/(EC50^n+C(t)^n)) - kout*R(t)
    where kin=kout*R0,
    based on steady-state solution R'(t)=0=kin-kout*R(0) for C(0)=0

    - pValues, vector of the five parameters [R0, kout, Emax, EC50, n]
    - parameterNames, {'R0','kout','Emax','EC50','n'}
    '''
    def myodes(self, t, y, k, experiment, functionObj, INDEX_FUNCTIONOBJ):
        Csim = self.simulate(functionObj,t,experiment,[],[])
        C = Csim(INDEX_FUNCTIONOBJ)
        R0, kout, Emax, EC50, n = k
        kin = kout*R0
        dy = np.zeros(1,1)
        dy[0] = kin *(1 + (Emax*C**n)/(EC50**n+C**n)) - kout*y[0] #Stimulation of prod
        return dy

class PD_receptocc_ode(PD_models_ode):
    ''' PD_receptocc_ode 
    A receptor occupancy model; kon=koff/KD
    dR(t)/dt = kon*C(t)*(Emax-R(t)) - koff*R(t)
    '''
    def __init__(self, pValues=None, inputVariableName=['C'], outputVariableName=['R']):
        try:
            super().__init__(pValues, inputVariableName, outputVariableName)
            #old syntax : super(PD_turnover_inhOfLoss_ode, self).__init__(pValues, inputVariableName, outputVariableName)

            self.parameterNames = ['R0','kout','Emax','EC50','n'] # inhofProd/inhOfLoss
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException(' Number of parameters incorrect!')
        except pko.BadInputException:
            print('Bad Input')

    def myodes(self, t, y, k, experiment, functionObj, INDEX_FUNCTIONOBJ):
        Csim = self.simulate(functionObj,t,experiment,[],[])
        C = Csim(INDEX_FUNCTIONOBJ)
        KD, koff, Emax = k
        kon = koff/KD;
        dy = np.zeros(1,1)
        dy[0] = kon*C *(Emax - y[0]) - koff*y[0] # receptocc_ode
        return dy

#================================================================================================================
#           EXAMPLES WHICH WERE PUT INITIALLY IN THE DOCSTRINGS
#================================================================================================================
if __name__ == "__main__": 

    # PK_1comp_Ka_Cl_ode & PK_1comp_Ka_satCl_ode

    timeUnit = 'min'
    doseUnit = 'umol/kg'
    d = pko.Dosing('ev', [0, 24],[20, 10], [], timeUnit, doseUnit)
    experiment = pko.Experiment('Ex',d,[])
    times = np.arange(0,49,1) # [0:1:48]
    const = 0

    #For each model
    for m in ['inhOfLoss','inhOfProd', 'stimOfLoss', 'stimOfProd']:
        # With and without effectCompt activated
        for k in [0,1]:
            effectComp = k
            if k==0:
            # Example 1:
                if m=='inhOfLoss':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_inhOfLoss_ode([1, 0.1, 1, 1, 1])
                elif m=='inhOfProd':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_inhOfProd_ode([1, 0.1, 1, 1, 1])
                elif m=='stimOfLoss':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_inhOfProd_ode([1, 0.1, 1, 1, 1])
                elif m=='stimOfProd':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_stimOfProd_ode([1, 0.1, 1, 1, 1])  #inhofLoss dans matlab
            else:
            # Example 2:
                if m=='inhOfLoss':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86, 1.6, 0.84, 0.1],{'C','C2'})
                    sf2 = PD_turnover_inhOfLoss_ode([1, 0.1, 1, 1, 1],'C2','E')
                elif m=='inhOfProd':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86, 1.6, 0.84, 0.1],{'C','C2'})
                    sf2 = PD_turnover_inhOfProd_ode([1, 0.1, 1, 1, 1],'C2','E')
                elif m=='stimOfLoss':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_inhOfProd_ode([1, 0.1, 1, 1, 1])
                elif m=='stimOfProd':
                    sf1 = pkm.PK_1comp_Ka_Cl_fun(const,effectComp,[0.86,  1.6, 0.84])
                    sf2 = PD_turnover_stimOfProd_ode([1, 0.1, 1, 1, 1],'C2','E')
            model = pko.Model(sf1, sf2, [])
            Y = model.simulate(times,experiment)
            
            fig, ax = plt.subplots()
            if m in ['inhOfLoss','']:
                ax.plot(times, Y[:,0], 'b', times, Y[:,1], 'm')
                if k>1: ax.plot(times, Y[:,2], 'r')
            elif m in ['2comp_Ka']:
                pass
                #semilogy(times,[:,0],'g')
                #semilogy(times,[:,1],'b')
                #semilogy(times,Y[:,2],'r') 
                #if k>1: ax.plot(times,Y[:,3],'m')
            elif m in ['3comp_Ka']:
                pass
                #semilogy(times,Y[:,3],'k')
                #if k>1: ax.plot(times,Y[:,4],'k')

            ax.set_title('Sanity check')
            fig.legend(model.getOutputVariables())
            plt.xlabel(timeUnit)
            #plt.ylabel('uM')
            plt.show()
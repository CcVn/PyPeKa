#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PK models with one compartment

@author: CcVN
@contact: cedric.vinson2@gmail.com
@since: 2018/09/06
@status: alpha
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import PKObjects as pko

__author__ = "CcVN"
__contact__ = "cedric.vinson2@gmail.com"
__date__ = "2018"
__version__ = "0.0.1"

class PK_1comp_Ka_Cl_fun(pko.SimFunction):
    '''
     PK_1comp_Ka_Cl_fun :1 compartment, linear absorption, linear elimination, closed form 
      (Subclass of SimFunction)

     Inputs:
     - constants, integer 0, 1, 2 or 3
       0: microConstants, default. ka, V, k, (ke)
       1: volumeClearance. ka, Cl, V, (ke)
       2: microConstants and bioavailability. ka, V, k, F, (ke)
       3: volumeClearance and bioavailability. ka, Cl, V, F, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable C
       1: return the variables C and Ce (effect compartment)
     - pValues, vector of parameters (3 or 4 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','Ce'} 
                      Optional input argument. 
                      Defaults are {'G','C'} or {'G','C','Ce'} if effect comp.

     Fields defined in the superclass SimFunction
     - parameterNames, {'ka','V','k'} microConstants, no effect comp.
                       {'ka','V','k','ke'} microConstants, effect comp.
                       {'ka','Cl','V'} volumeClearance, no effect comp.
                       {'ka','Cl','V','ke'} volumeClearance, effect comp.
                       {'ka','V','k','F'} microConstants with F, no effect comp.
                       {'ka','V','k','F','ke'} microConstants with F, effect comp.
                       {'ka','Cl','V','F'} volumeClearance with F, no effect comp.
                       {'ka','Cl','V','F','ke'} volumeClearance  with F, effect comp.
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above

     Function:
     Y = obj.simulate(times, experiment, yInput, functionObj), simulates the object
        times - vector of time points to output
        experiment - the experiment to simulate (wrt dosing)
        yInput - input variable values at time points in times, size |times|x|inputVariables|
       obj from which inputVariables can be calculated at arbitrary time-points

     returns: Y matrix of output variable data, size |times|x|inputVariables|
    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['C']):
        try:
            self.nbComp = 1
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','k']   # ['k01','V','k10']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
            else:
                self.parameterValues = pValues

            # Modify parameters :
            if constants == 0:  pass #--> default
            elif constants == 1: self.parameterNames = ['ka','Cl','V']
            elif constants == 2: self.parameterNames = ['ka','V','k','F']
            elif constants == 3: self.parameterNames = ['ka','Cl','V','F']
            else: raise pko.BadInputException('ResultChk:BadInput', ' Constants not in [0,1,2,3]')

            if effectCompartment > 0:
                self.calculateEffectCompartment = 1
                if 'Ce' not in self.outputVariables: 
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(3)={'Ce'}
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0

            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            #Control the outputVariables argumetn
            if (len(outputVariables) == 1 and effectCompartment == 0) or \
                (len(outputVariables) == 2 and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+str(len(outputVariables)) + ': ' + str(outputVariables))

            self.infusion_mode = 0 

        except:
            print('ERROR')

#% %
    #============ Equations
    def simulate(self, times, experiment, yInput, functionObj):
        '''
        TAU is defined by the difference between times of dose 1 and 2
        Dose is defined by the dose 1
        '''
        # Reparametrize the model according to constants argument:
        k = self.parameterValues
        k01 = k[1]  # ka (toujours en position 1)
        if self.constants == 0: # 'ka','V','k'
            k01, V, k10, *autres = k
        elif self.constants == 1: # 'ka','Cl','V'
            V = k[3]
            k10 = k[1] / k[2]
        elif self.constants == 2: # 'ka','V','k','F'
            V = k[2]                 
            if experiment.dosing.admin in ['ev']:
                V = V / k[3] # divide by F gives V/F
            k10 = k[3]
        elif self.constants == 3: # 'ka','Cl','V','F'
            V = k[3] 
            if experiment.dosing.admin in ['ev']:
                V = V / k[3] # divide by F gives V/F
            k10 = k[1] / k[2] # Cl/F / V/F same as Cl/V

        C = np.zeros( (len(times),1) )
        alpha = k10
        #Dose1 = 0 # the first dose considered
        effet = (self.calculateEffectCompartment != 0)

        # EXTRAVASCULAR ADMINISTRATION
        #------------------------------
        if experiment.dosing.admin == 'ev':
            self.infusion_mode = 0

            A = (k01/V)/(k01-alpha)
            if effet:
                # for an effect compartment as well
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0*A / (ke0 - alpha)
                CE = -(AE*(k01 - alpha)) / (k01 - ke0)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    t_after_D = tps - experiment.dosing.times[idx]
                    if t_after_D >= 0:
                        AA = laDose *A
                        CC = AA
                        C[it] = C[it] + AA *np.exp(-alpha*t_after_D) \
                                      - CC *np.exp(  -k01*t_after_D)
                        if effet:
                            # for an effect compartment as well
                            AAE = laDose *AE
                            CCE = laDose *CE
                            DDE = AAE + CCE
                            C_E[it] = C_E[it] + AAE *np.exp(-alpha*t_after_D) \
                                              + CCE *np.exp(  -ke0*t_after_D) \
                                              - DDE *np.exp(  -k01*t_after_D)

        # INTRAVENOUS BOLUS
        #-------------------
        elif experiment.dosing.admin == 'iv':
            self.infusion_mode = 0

            if effet:
                # for an effect compartment as well
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0*A/(ke0-alpha)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    t_after_D = tps - experiment.dosing.times[idx]
                    if t_after_D >= 0.0:
                        tempA = A *np.exp(-alpha*t_after_D)
                        C[it] = C[it] + laDose*tempA
                        if effet:
                            tempA = AE *np.exp(-alpha*t_after_D)
                            tempC = AE *np.exp(  -ke0*t_after_D)
                            C_E[it] = C_E[it] + laDose*(tempA-tempC)

        # INTRAVENOUS INFUSION
        #----------------------
        elif experiment.dosing.admin == 'inf':
            self.infusion_mode = 1

            A = 1/V
            if effet:
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0 *A/(ke0 - alpha)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    inf_dur = experiment.dosing.duration[idx]
                    t_inf = tps - experiment.dosing.times[idx]
                    t_after_inf = t_inf - inf_dur

                    # After the infusion
                    if t_inf >= 0 and t_after_inf > 0:
                        temp = A/alpha*(1 - np.exp(-alpha*inf_dur)) *np.exp(-alpha*t_after_inf)
                        if effet:
                            tempE =  AE/alpha *(1 - np.exp(-alpha*inf_dur)) *np.exp(-alpha*t_after_inf)
                            tempE2 = AE/ke0   *(1 - np.exp(  -ke0*inf_dur)) *np.exp(  -ke0*t_after_inf)
                    # Current ongoing infusion
                    elif t_inf >= 0:
                        temp = A/alpha*(1 - np.exp(-alpha*(t_inf)))
                        if effet:
                            tempE  = AE/alpha *(1 - np.exp(-alpha*t_inf))
                            tempE2 = (AE/ke0) *(1 - np.exp(  -ke0*t_inf))

                    C[it] = C[it] + laDose*temp
                    if effet: C_E[it] = C_E[it] + laDose *(tempE - tempE2)

        if self.calculateEffectCompartment !=0:
            Y = (C, C_E)
        else:
            Y = C
        if self.verbose: print('output Y:\n',Y)
        return Y


class PK_2comp_Ka_Cl_fun(pko.SimFunction):
    '''
     PK_2comp_Ka_Cl_fun : 2 compartments, linear absorption, linear elimination, closed 
      (Subclass of SimFunction)

     Inputs:
     - constants, integer 0, 1, 2 or 3
        0: microConstants, default. ka, V, k, k12, k21, (ke)
        1: volumeClearance. ka, Cl, V1, Q, V2, (ke)
        2: microConstants and bioavailability, ka, V, k, k12, k21, F, (ke)
        3: volumeClearance and bioavailability. ka, Cl, V1, Q, V2, F, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable C
       1: return the variables C and Ce (effect compartment)
     - pValues, vector of parameters (3 or 4 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','Ce'} 
                      Optional input argument. 
                      Defaults are {'G','C'} or {'G','C','Ce'} if effect comp.

     Fields defined in the superclass SimFunction
     - parameterNames, {'ka','V','k','k12','k21'} microConstants, no effect comp.
                       {'ka','V','k','k12','k21','ke'} microConstants, effect comp.
                       {'ka','Cl','V1','Q','V2'} volumeClearance, no effect comp.
                       {'ka','Cl','V1','Q','V2','ke'} volumeClearance, effect comp.
                       {'ka','V','k','k12','k21','F'} microConstants with F, no effect comp.
                       {'ka','V','k','k12','k21','F','ke'} microConstants with F, effect comp.
                       {'ka','Cl','V1','Q','V2','F'} volumeClearance with F, no effect comp.
                       {'ka','Cl','V1','Q','V2','F','ke'} volumeClearance with F, effect comp.
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above

     Function:
     Y = obj.simulate(times, experiment, yInput, functionObj), simulates the object
        times - vector of time points to output
        experiment - the experiment to simulate (wrt dosing)
        yInput - input variable values at time points in times, size |times|x|inputVariables|
       obj from which inputVariables can be calculated at arbitrary time-points

     returns: Y matrix of output variable data, size |times|x|inputVariables|
    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['C']):
        try:
            self.nbComp = 2
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','k','k12','k21']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
            else:
                self.parameterValues = pValues
            
            # Modify parameters :
            if constants == 0:  pass #--> default
            elif constants == 1: self.parameterNames = ['ka','Cl','V1','Q','V2']
            elif constants == 2: self.parameterNames = ['ka','V','k','k12','k21','F']
            elif constants == 3: self.parameterNames = ['ka','Cl','V1','Q','V2','F']
            else: raise pko.BadInputException('ResultChk:BadInput', ' Constants not in [0,1,2,3]')
    
            if effectCompartment > 0:
                self.calculateEffectCompartment = 1
                if 'Ce' not in self.outputVariables: 
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(3)={'Ce'}
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0
    
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            #Control the outputVariables argumetn
            if (len(outputVariables) == 1 and effectCompartment == 0) or \
                (len(outputVariables) == 2 and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ str(len(outputVariables)) + ': ' + str(outputVariables))

            self.infusion_mode = 0 

        except:
            print('ERROR')

#% %
    #============ Equations
    def simulate(self, times, experiment, yInput, functionObj):
        '''
        TAU is defined by the difference between times of dose 1 and 2
        Dose is defined by the dose 1
        '''
        # Reparametrize the model according to constants argument:
        k = self.parameterValues
        k01 = k[0]  # ka (toujours en position 1)
        if self.constants == 0: # 'ka','V','k','k12','k21'
            k01, V, k10, k12, k21, *autres = k
        elif self.constants == 1: # 'ka','Cl','V1','Q','V2'
            V = k[2]
            k10 = k[1] / k[2]  # k = Cl / V1
            k12 = k[3] / k[2]  # k12 = Q / V1
            k21 = k[3] / k[4]  # k21 = Q / V2
        elif self.constants == 2: # 'ka','V','k','k12','k21','F'
            V = k[1]
            if experiment.dosing.admin in ['ev']:
                V = V / k[5] # divide by F gives V/F
            k10 = k[2]
            k12 = k[3]
            k21 = k[4]
        elif self.constants == 3: # 'ka','Cl','V1','Q','V2','F'
            V = k[2] 
            if experiment.dosing.admin in ['ev']:
                V = V / k[5] # divide by F gives V/F
            k10 = k[1] / k[2] # Cl/F / V/F same as Cl/V
            k12 = k[3] / k[2] # Q/ V1
            k21 = k[3] / k[4] # Q/ V2

        C = np.zeros( (len(times),1) )
        beta = 0.5 *(k12 + k21 + k10 - np.sqrt((k12 + k21 + k10)**2 - 4*k21*k10))
        alpha = k10 *k21/beta  # alpha = k10 en 1 comp
        effet = (self.calculateEffectCompartment != 0)

        # EXTRAVASCULAR ADMINISTRATION
        #------------------------------
        if experiment.dosing.admin == 'ev':
            self.infusion_mode = 0

            A = (k01/V)*(k21-alpha) / ((k01-alpha)*( beta - alpha))
            B = (k01/V)*(k21- beta) / ((k01- beta)*(alpha -  beta))
            if effet:
                # for an effect compartment as well
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0*A /(ke0 - alpha)
                BE = ke0*B /(ke0 - beta)
                CE = -(AE*(k01 - alpha) + BE*(k01 - beta)) / (k01 - ke0)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    t_after_D = tps - experiment.dosing.times[idx]
                    if t_after_D >= 0:
                        AA = laDose *A
                        BB = laDose *B
                        CC = AA + BB
                        C[it] =  C[it] + AA*np.exp(-alpha*t_after_D) \
                                       + BB*np.exp( -beta*t_after_D) \
                                       - CC*np.exp(  -k01*t_after_D)

                        if effet:
                            # for an effect compartment as well
                            AAE = laDose *AE
                            BBE = laDose *BE
                            CCE = laDose *CE
                            DDE = AAE + BBE + CCE
                            C_E[it] = C_E[it] + AAE *np.exp(-alpha*t_after_D) \
                                              + BBE *np.exp( -beta*t_after_D) \
                                              + CCE *np.exp(  -ke0*t_after_D) \
                                              - DDE *np.exp(  -k01*t_after_D)

        # INTRAVENOUS BOLUS
        #-------------------
        elif experiment.dosing.admin == 'iv':
            self.infusion_mode = 0

            A = (1/V) *(alpha - k21) / (alpha - beta)
            B = (1/V) *( beta - k21) / (beta - alpha)
            if effet:
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0*A / (ke0 - alpha)
                BE = ke0*B / (ke0 - beta)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    if (tps >= experiment.dosing.times[idx]):
                        tempA = A *np.exp(-alpha*(tps-laDose))
                        tempB = B *np.exp( -beta*(tps-laDose))
                        C[it] = C[it] + laDose*(tempA + tempB)

                        if self.calculateEffectCompartment !=0 :
                            tempA =      AE *np.exp(-alpha*(tps - laDose))
                            tempB =      BE *np.exp( -beta*(tps - laDose))
                            tempC = (AE+BE) *np.exp(  -ke0*(tps - laDose))
                            C_E[it] = C_E[it] + laDose *(tempA + tempB - tempC)


        # INTRAVENOUS INFUSION
        #----------------------
        elif experiment.dosing.admin == 'inf':
            self.infusion_mode = 1

            A = (1/V) *(alpha - k21) / (alpha - beta)
            B = (1/V) *( beta - k21) / (beta - alpha)
            if effet:
                C_E = np.zeros( (len(times),1) )
                ke0 = k[-1]
                AE = ke0*A / (ke0 - alpha)
                BE = ke0*B / (ke0 - beta)

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    inf_dur = experiment.dosing.duration[idx]
                    t_inf = tps - experiment.dosing.times[idx]
                    t_after_inf = t_inf - inf_dur

                    # After infusion 
                    if (t_inf >= 0 and t_after_inf > 0):
                        temp = A/alpha *(1 - np.exp(-alpha*inf_dur)) *np.exp(-alpha*t_after_inf) \
                             + B/beta  *(1 - np.exp( -beta*inf_dur)) *np.exp( -beta*t_after_inf)

                        if self.calculateEffectCompartment !=0 :
                            tempE =   AE/alpha   *(1 - np.exp(-alpha*inf_dur)) *np.exp(-alpha*t_after_inf) \
                                    + BE/beta    *(1 - np.exp( -beta*inf_dur)) *np.exp( -beta*t_after_inf)
                            tempE2 = (AE+BE)/ke0 *(1 - np.exp(  -ke0*inf_dur)) *np.exp(  -ke0*t_after_inf)

                    # Current ongoing infusion
                    elif t_inf >= 0:
                        temp = A/alpha *(1 - np.exp(-alpha*t_inf)) \
                             + B/beta  *(1 - np.exp( -beta*t_inf))

                        if self.calculateEffectCompartment !=0 :
                            tempE =       AE/alpha *(1 - np.exp(-alpha*t_inf)) \
                                       + BE/beta   *(1 - np.exp( -beta*t_inf))
                            tempE2 = (AE+BE)/ke0   *(1 - np.exp(  -ke0*t_inf))

                    C[it] = C[it] + laDose*temp
                    if effet: C_E[it] = C_E[it] + laDose *(tempE - tempE2)

        if effet:
            Y = (C, C_E)
        else:
            Y = C
        if self.verbose: print('output Y:\n',Y)
        return Y


class PK_3comp_Ka_Cl_fun(pko.SimFunction):
    '''
     PK_3comp_Ka_Cl_fun : 3 compartments, linear absorption, linear elimination, closed form
       (Subclass of SimFunction)

     Inputs:
     - constants, integer 0, 1, 2 or 3
        0: microConstants, default. ka, V, k, k12, k21, k13, k31
        1: volumeClearance. ka, Cl, V1, Q2, V2, Q3, V3
     - effectCompartment, integer 
        0: return the variable C [only allowed choice]
        1: effect-compartment [NOT IMPLEMENTED]
     - pValues, vector of parameters (3 or 4 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','C2','C3'} 
                      Optional input argument. Default is {'G','C','C2','C3'}

     Fields defined in the superclass SimFunction
     - parameterNames, {'ka','V','k','k12','k21','k13','k31'} microConstants
                       {'ka','Cl','V1','Q2','V2','Q3','V3'} volumeClearance
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above

     Function:
     Y = obj.simulate(times, experiment, yInput, functionObj), simulates the object
        times - vector of time points to output
        experiment - the experiment to simulate (wrt dosing)
        yInput - input variable values at time points in times, size |times|x|inputVariables|
       obj from which inputVariables can be calculated at arbitrary time-points

     returns: Y matrix of output variable data, size |times|x|inputVariables|
    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['C']):
        try:
            self.nbComp = 3
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','k','k12','k21','k13','k31']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
            else:
                self.parameterValues = pValues
            
            # Modify parameters :
            if constants == 0:  pass #--> default
            elif constants == 1: self.parameterNames = ['ka','Cl','V1','Q2','V2','Q3','V3']
            else: raise pko.BadInputException('ResultChk:BadInput', ' Constants not in [0,1]')
    
            if effectCompartment > 0:
                 raise pko.BadInputException('ResultChk:BadInput', 'Effect-compartment not implemented in PK_3comp_Ka_Cl_fun. Use PK_3comp_Ka_Cl_ode.')
            else:
                self.calculateEffectCompartment = 0
    
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            #Control the outputVariables argumetn
            if (len(outputVariables) == 1 and effectCompartment == 0):
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ \
                            str(len(outputVariables)) + ': ' + str(outputVariables))

            self.infusion_mode = 0

        except:
            print('ERROR')

#% %
    #============ Equations
    def simulate(self, times, experiment, yInput, functionObj):
        '''
        TAU is defined by the difference between times of dose 1 and 2
        Dose is defined by the dose 1
        '''
        # Reparametrize the model according to constants argument:
        #NB: Identique au code SS sauf definition de TAU en plus dans SS
        k = self.parameterValues
        k01 = k[0]  # ka (toujours en position 1)
        if self.constants == 0: # 'ka','V','k','k12','k21','k13','k31'
            k01, V, k10, k12, k21, k13, k31, *autres = k
        elif self.constants == 1: # 'ka','Cl','V1','Q2','V2','Q3','V3'
            V = k[2]
            k10 = k[1] / k[2]  # k = Cl / V1
            k12 = k[3] / k[2]  # k12 = Q / V1
            k21 = k[3] / k[4]  # k21 = Q / V2
            k13 = k[5] / k[2]  # k13 = Q3 / V1
            k31 = k[5] / k[6]  # k31 = Q3 / V3
        elif self.constants == 2: # 'ka','V','k','k12','k21','k13','k31','F'
            V = k[2]
            if experiment.dosing.admin in ['ev']:
                V = V / k[7] # divide by F gives V/F
            k10 = k[2]
            k12 = k[3]
            k21 = k[4]
            k13 = k[5]
            k31 = k[5]
        elif self.constants == 3: # 'ka','Cl','V1','Q','V2','Q3','V3','F'
            V = k[2] 
            if experiment.dosing.admin in ['ev']:
                V = V / k[7] # divide by F gives V/F
            k10 = k[1] / k[2] # Cl/F / V/F same as Cl/V
            k12 = k[3] / k[2] # Q/ V1
            k21 = k[3] / k[4] # Q/ V2
            k13 = k[5] / k[2]  # k13 = Q3 / V1
            k31 = k[5] / k[6]  # k31 = Q3 / V3

        a0 = k10*k21*k31
        a1 = k10*k31+k21*k31 + k21*k13 + k10*k21 + k31*k12
        a2 = k10 + k12 + k13 + k21 + k31
        p = a1 - a2*a2/3
        q = 2 *a2**3/27 - a1*a2/3 + a0
        r1 = np.sqrt(-p**3/27)
        r2 = 2 *(r1**(1/3))
        phi = math.acos(-q/2/r1)/3
        alpha =               -(math.cos(phi)*r2 - a2/3)
        beta =  -(math.cos(phi + 2*math.pi/3)*r2 - a2/3)
        gamma = -(math.cos(phi + 4*math.pi/3)*r2 - a2/3)
        C = np.zeros( (len(times),1) )
        #effet = (self.calculateEffectCompartment != 0)

        # EXTRAVASCULAR ADMINISTRATION
        #------------------------------
        if experiment.dosing.admin == 'ev':
            self.infusion_mode = 0

            A = (1/V) *(k01/(k01 - alpha))*((k21 - alpha) / (alpha -  beta))*((k31 - alpha)/(alpha - gamma))
            B = (1/V) *(k01/(k01 -  beta))*((k21 -  beta) / ( beta - alpha))*((k31 -  beta)/( beta - gamma))
            C_= (1/V) *(k01/(k01 - gamma))*((k21 - gamma) / (gamma -  beta))*((k31 - gamma)/(gamma - alpha))

                    
            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    t_after_D = tps - experiment.dosing.times[idx]
                    if t_after_D >= 0:
                        t1 =       A  *np.exp(-alpha*t_after_D)
                        t2 =       B  *np.exp( -beta*t_after_D)
                        t3 =       C_ *np.exp(-gamma*t_after_D)
                        t4 = (A+B+C_) *np.exp(  -k01*t_after_D)
                        C[it] = C[it] + laDose*(t1 + t2 + t3 - t4)

        # INTRAVENOUS BOLUS
        #-------------------
        elif experiment.dosing.admin == 'iv':
            self.infusion_mode = 0

            A = (1/V)*((k21-alpha)/(alpha-beta)) *((k31-alpha)/(alpha-gamma))
            B = (1/V)*((k21-beta )/(beta-alpha)) *((k31-beta )/(beta-gamma))
            C_= (1/V)*((k21-gamma)/(gamma-beta)) *((k31-gamma)/(gamma-alpha))

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    t_after_D = tps - experiment.dosing.times[idx]
                    if t_after_D >= 0:
                        tempA = A *np.exp(-alpha*t_after_D)
                        tempB = B *np.exp( -beta*t_after_D)
                        tempC = C_*np.exp(-gamma*t_after_D)
                        C[it] = C[it] + laDose *(tempA + tempB + tempC)

        # INTRAVENOUS INFUSION
        #----------------------
        elif experiment.dosing.admin == 'inf':
            self.infusion_mode = 1

            A = (1/V)*((k21-alpha)/(alpha-beta)) *((k31-alpha)/(alpha-gamma))
            B = (1/V)*((k21-beta )/(beta-alpha)) *((k31-beta )/( beta-gamma))
            C_= (1/V)*((k21-gamma)/(gamma-beta)) *((k31-gamma)/(gamma-alpha))

            for it, tps in enumerate(times):
                for idx, laDose in enumerate(experiment.dosing.values):
                    inf_dur = experiment.dosing.duration[idx]
                    t_inf = tps - experiment.dosing.times[idx]
                    t_after_inf = t_inf - inf_dur

                    # After infusion 
                    if (t_inf >= 0 and t_after_inf > 0):
                        temp = A/alpha  *(1 - np.exp(-alpha*inf_dur)) *np.exp(-alpha*t_after_inf) \
                             + B/beta   *(1 - np.exp( -beta*inf_dur)) *np.exp( -beta*t_after_inf) \
                             + C_/gamma *(1 - np.exp(-gamma*inf_dur)) *np.exp(-gamma*t_after_inf)
                    # Current ongoing infusion
                    elif t_inf >= 0:
                        temp =  A/alpha *(1 - np.exp(-alpha*t_inf)) \
                            +   B/beta  *(1 - np.exp( -beta*t_inf)) \
                            +  C_/gamma *(1 - np.exp(-gamma*t_inf))
                    C[it] = C[it] + laDose*temp

        Y = C
        if self.verbose: print('output Y:\n',Y)
        return Y



#================================================================================================================
#           EXAMPLES WHICH WERE PUT INITIALLY IN THE DOCSTRINGS
#================================================================================================================
if __name__ == "__main__": 

    # PK_1comp_Ka_Cl_ode & PK_1comp_Ka_satCl_ode

    timeUnit = 'min'
    doseUnit = 'umol'
    #Dosing varient selon les exemples
    d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
    experiment = pko.Experiment('Ex',d,[])
    times = np.arange(0,73,1) # [0:1:72]
    const = 0

    #For each model
    for m in ['1comp', '2comp', '3comp']:
        # With and without effectCompt activated
        for k in [0,1]:
            effectComp = k
            if k==0:
            # Example 1:
                if m=='1comp':
                    sf = PK_1comp_Ka_Cl_fun(const, effectComp, [0.86, 1.6,0.84])
                elif m=='2comp':
                    sf = PK_2comp_Ka_Cl_fun(const,effectComp,[0.86, 1.6, 0.84, 1.5, 0.1])
                elif m=='3comp':
                    sf = PK_3comp_Ka_Cl_fun(const,effectComp,[0.86, 1.6, 0.84, 1.5, 0.1, 1.4, 0.05])
            else:
            # Example 2: effectComp=1
                if m=='1comp':
                    sf = PK_1comp_Ka_Cl_fun(const, effectComp,[0.86, 1.6, 0.84, 1],['Cp','Ce'])
                elif m=='2comp':
                    sf = PK_2comp_Ka_Cl_fun( const, effectComp, [0.86, 1.6, 0.84, 1.5, 0.1, 1],['Cp','Ce'])
                elif m=='3comp':
                    print('PK 3comp model with effect compartment : Not implemented')
                    break
            Y = sf.simulate(times, experiment, [], [])
        
            fig, ax = plt.subplots()
            if m in ['1comp']:
                if k==0: ax.plot(times, Y[:,0], 'r')
                else: ax.plot(times, Y[0][:,0], 'b', Y[1][:,0], 'r')
            elif m in ['2comp']:
                if k==0: ax.semilogy(times, Y[:,0], 'r')
                else: ax.plot(times, Y[0][:,0], 'b', Y[1][:,0], 'r')
            elif m in ['3comp']:
                ax.semilogy(times,Y[:,0])

            ax.set_title('Sanity check')
            fig.legend(sf.outputVariables)
            plt.xlabel(timeUnit)
            plt.ylabel('uM')
            plt.show()
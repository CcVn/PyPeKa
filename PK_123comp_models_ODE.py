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

#scipy.integrate.ode(myodes).set_integrator('vode', method='bdf', order=15)
#If you want automatic switching between stiff/non-stiff then you can use,
#scipy.integrate.ode(myodes).set_integrator('lsoda')


class PK_1comp_Ka_Cl_ODE(pko.SimFunction):
    '''
     PK_1comp_Ka_Cl_ode (Subclass of SimFunction)
     1 compartment, linear absorption, linear elimination, ODE form

     Inputs:
     - constants, integer 0 or 1
       0: microConstants, default. ka, V, k, (ke)
       1: volumeClearance. ka, Cl, V, (ke)
       2: microConstants and bioavailability. ka, V, k, F, (ke)
       3: volumeClearance and bioavailability. ka, Cl, V, F, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable G and C
       1: return the variables G, C and Ce (effect compartment)
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
     Y = obj.simulate(times, experiment, yInput, functionObj), simulates the selfect
      - times - vector of time points to output
      - experiment - the experiment to simulate (wrt dosing)
      - yInput - input variable values at time points in times, size |times|x|inputVariables|
      - functionObj from which inputVariables can be calculated at arbitrary time-points

     Returns: Y matrix of output variable data, size |times|x|inputVariables|
    '''
    verbose = False

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['G','C']):
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
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(3)={'Ce'};
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0
    
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            #Control the outputVariables argumetn
            if (len(outputVariables) == (self.nbComp+1) and effectCompartment == 0) or \
                (len(outputVariables) == (self.nbComp+2) and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ \
                            str(len(outputVariables)) + ': ' + str(outputVariables))
                
            # Reparametrize the model according to constants argument:
            k = self.parameterValues
            k01 = k[1]  # ka (toujours en position 1)
            if self.constants == 0: # 'ka','V','k'
                k01, V, k10, *autres = k
            elif self.constants == 1: # 'ka','Cl','V'
                V = k[3]
                k10 = k[2] / k[3]
            elif self.constants == 2: # 'ka','V','k','F'
                V = k[2]                 
                if experiment.dosing.admin in ['ev']:
                    V = V / k[3] # divide by F gives V/F
                k10 = k[3]
            elif self.constants == 3: # 'ka','Cl','V','F'
                V = k[3] 
                if experiment.dosing.admin in ['ev']:
                    V = V / k[3] # divide by F gives V/F
                k10 = k[2] / k[3] # Cl/F / V/F same as Cl/V
            
            # Locally we use kn = [ka,V,k] or kn = [ka,V,k,ke]
            self.kn = [k01, V, k10]
            if self.calculateEffectCompartment != 0:
                self.kn.append(k[-1])
    
            self.F = 1  # was only present in PK_1comp_Ka_satCl_ODE for unknown reason
            self.infusion_mode = 0 
            #print(self)
            
        except:
            print('Bad Input')
            raise

#% %
    #============ Call to ODE solver
    def simulate(self, times, experiment, yInput, functionObj):
        '''
        Assume times to be a vector
        '''        
        times = np.array(times)
        if len(times.shape)==2 and times.shape[0]>1:
            times = times.T

        kn = self.kn
        #Conversion in tuple for odeint; test if iterable in case of a single int
        try:
            iter(kn)
            knt = tuple(kn)
        except TypeError:
            kn = [kn]
            knt = tuple(kn)

        if experiment.dosing.admin == 'ev':
            iDoseVariable = 0
            self.infusion_mode = 0
        elif experiment.dosing.admin == 'iv':
            iDoseVariable = 1
            iVolPara = 1 # Volume parameter
            self.infusion_mode = 0
        elif experiment.dosing.admin == 'inf':
            self.infusion_mode = 1

        E = [] #np.zeros((1,len(self.outputVariables)))
        init = np.zeros_like(self.outputVariables, dtype=float)
       
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
                #Etem = ode15s(    myodes, t, init, [], kn, experiment)
                Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )
            elif experiment.dosing.admin in ['inf']:
                Etem = integ.odeint(self.myodes_inf, init, t, args=(knt,), tfirst=True )

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

#           Etem = np.empty_like(Etem)  #FIXME Etem=[]
            if experiment.dosing.admin in ['ev','iv']:
                if experiment.dosing.admin == 'ev':
                    init[iDoseVariable] = init[iDoseVariable] + experiment.dosing.values[idx] * self.F
                elif experiment.dosing.admin == 'iv':
                    init[iDoseVariable] = init[iDoseVariable] + experiment.dosing.values[idx] / kn[iVolPara]
               #-- INTEGRATE: [Tdummy, Etem] = ode15s(@myodes, t, init, [], kn, experiment)
                Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )

            elif experiment.dosing.admin == 'inf':
                #[Tdummy, Etem] = ode15s(@myodes_inf, t, init, [], kn, experiment)
                Etem = integ.odeint(self.myodes_inf, init, t, args=(knt,), tfirst=True )

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
            print('Error, t has len 1')
            Etem = init
        else:
            Etem = []
            if experiment.dosing.admin in ['ev','iv']:
                if experiment.dosing.admin == 'ev':
                    init[iDoseVariable] = init[iDoseVariable] + experiment.dosing.values[-1] * self.F
                elif experiment.dosing.admin == 'iv':
                    init[iDoseVariable] = init[iDoseVariable] + experiment.dosing.values[-1] / kn[iVolPara]
                #Etem = ode15s(@myodes,t,init,[],kn,experiment)
                Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True ) 
            elif experiment.dosing.admin == 'inf':
                #Etem = ode15s(@myodes_inf,t,init,[],kn,experiment)
                Etem = integ.odeint(self.myodes, init, t, args=(knt,), tfirst=True )
            if len(t) == 2:
                Etem = Etem[[0,-1], :] # garder ligne 1 et finale

        E_slice = Etem[ (addFront):len(t), : ] #tout sauf 1er si addFront
        E = np.vstack( (E, E_slice) ) # ou E = np.r_[ E, E_slice ]

        if self.verbose: print('output Y:\n',E)
        return E

#% %
    #============ iv + ev, and infusion
    def myodes(self, t, y, k): #, experiment
        if len(y)==2: k01, V, k10 = k 
        elif len(y)==3: k01, V, k10, ke = k
        else: raise Exception
        
        if len(y)==2: dy = np.zeros((2,), dtype=float) # No effect compartment
        elif len(y)==3: dy = np.zeros((3,), dtype=float)

        # 'ka','V','k'
        dy[0] = -k01 * y[0]
        dy[1] =  k01 * y[0] / V - k10 * y[1]
        if len(y)==3: dy[2] = ke * (y[1] - y[2]) # An effect compartment, dy[2]

        return dy
 
    #============ Infusion
    def myodes_inf(self, t, y, k): #, experiment
        if len(y)==2: k01, V, k10 = k 
        elif len(y)==3: k01, V, k10, ke = k
        else: raise Exception
        
        inp = 0
        if self.infusion_mode > 0:  # preparation pour futur merge des myodes
            for idx, tps in enumerate(experiment.dosing.times):
                if  t >= tps and t < (tps + experiment.dosing.duration[idx]):
                    inp = experiment.dosing.values[idx]
                    break

        if len(y)==2: dy = np.zeros((2,), dtype=float) # No effect compartment
        else: dy = np.zeros((3,), dtype=float)

        # 'ka','V','k'
        dy[0] = 0.0  #-k01 * y[0]  # FIXME: calcule pour rien
        dy[1] = inp / V - k10 * y[1]
        if len(y)==3:  dy[2] = ke *(y[1] - y[2])   # An effect compartment, dy[2]
           
        return dy


class PK_1comp_Ka_satCl_ODE(PK_1comp_Ka_Cl_ODE):
    '''
     PK_1comp_Ka_satCl_ode (Subclass of Simdef)
     1 compartment, linear absorption, saturated elimination, ODE form

     Inputs:
     - constants, integer 0 or 1
       0: Without F: ka, V, Vmax, KM, (ke)
       1: With F: ka, V, Vmax, KM, F, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable C
       1: return both variables C and Ce (effect compartment)
     - pValues, vector of parameters (4 or 5 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','Ce'} 
                      Optional input argument. 
                      Defaults are {'G','C'} or {'G','C','Ce'} if effect comp.
    
     Fields defined in the superclass Simdef
     - Vm = maximum elimination rate (in amount per time unit)
     - Km = Michaelis-Menten constant (in concentration unit)
     - parameterNames, {'ka','V','Vmax','KM'}          without F, no effect comp.
                       {'ka','V','Vmax','KM','ke'}     without F, effect comp.
       parameterNames, {'ka','V','Vmax','KM''F',}      with F, no effect comp.
                       {'ka','V','Vmax','KM','F','ke'} with F, effect comp.
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above
    
     Function:
     Y=self.simulate(times, experiment, yInput, functionObj), simulates the object
     - times - vector of time points to output
     - experiment - the experiment to simulate (wrt dosing)
     - yInput input variable values at time points in times, size |times|x|inputVariables|
     - functionObj from which inputVariables can be calculated at arbitrary time-points
     
     returns Y matrix of output variable data, size |times|x|inputVariables|
    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['G','C']):
        try:
            self.nbComp = 1
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','Vmax','KM']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
    
            if constants==0:
                pass
                # default
            elif constants==1:
                self.parameterNames = ['ka','V','Vmax','KM','F']
            else:
                raise Exception('Error in PK_1comp_Ka_satCl_ode, constants not in [0,1]')
    
            if effectCompartment > 0:
                self.calculateEffectCompartment = 1
                if 'Ce' not in self.outputVariables: 
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(3)={'Ce'};
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0
            
            self.parameterValues = pValues
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            if (len(outputVariables)== (self.nbComp+1) and effectCompartment == 0) or \
                (len(outputVariables) == (self.nbComp+2) and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ \
                            str(len(outputVariables)) + ': ' + str(outputVariables))
    
            #Matlab code: was in function simulate
            self.k = self.parameterValues
            self.F = 1 if self.constants != 1 else self.parameterValues[4] # pas dans Ka_Cl_ODE, pourquoi?

            self.kn = self.k  # if applicable, ke is included here
            if self.constants==0:
                pass
            elif self.constants==1: # switch to microconstants
                # from 'ka','Cl','V1','Q','V2' to 'ka','V','k','k12','k21'
                self.kn[2] = k[3]
                self.kn[3] = self.k[2] / self.k[3]
                self.kn[4] = self.k[4] / self.k[3]
                self.kn[5] = self.k[4] / self.k[5]
            else:
                raise Exception('Error in simfunctionODE, constants not in [0,1]')

            self.infusion_mode = 0 

        except:
            raise

#    def simulate(self, times, experiment, ...) --> appel a classe parente

    def myodes(self, t, y, k):
        if len(y)==2:
            k01, V, Vmax, KM = k
        elif len(y)==3:
            k01, V, Vmax, KM, ke = k
        else: raise Exception

        if len(y)==2: dy = np.zeros((2,), dtype=float) # No effect compartment
        elif len(y)==3: dy = np.zeros((3,), dtype=float)
        else: raise Exception

        # 'ka','V','Vmax','KM'
        dy[0] = -k01 * y[0]
        admin = k01 * y[0]
        dy[1] = admin / V - (Vmax/V) *y[1] / (KM + y[1])
        if len(y)==3: dy[2] = ke * (y[1] - y[2]) # An effect compartment, dy[2]
    
        return dy

    def myodes_inf(self, t, y, k):
        #Unpacking identique myode
        if len(y)==2:
            k01, V, Vmax, KM = k
        elif len(y)==3:
            k01, V, Vmax, KM, ke = k
        else: raise Exception
        #Administration of the infusion : #FIXME to be merged with myodes
        if True:
            inp = 0
            for idx in range(len(experiment.dosing.times)):
                if  t >= experiment.dosing.times[idx] and t < (experiment.dosing.times[idx] + experiment.dosing.duration[idx]):
                    inp = experiment.dosing.values[idx]
                    break
        #Initialisation identique myode
        if len(y)==2: dy = np.zeros((2,), dtype=float) # No effect compartment # 'ka','V','k'
        elif len(y)==3: dy = np.zeros((3,), dtype=float)
        else: raise Exception

        # 'ka','V','Vmax','KM'
        dy[0] = -k01*y[0]
        admin = inp  # SEULE DIFF AVEC MODELE Ka_Cl avec unpacking de k
        dy[1] = admin/ V - (Vmax/V)*y[1] / (KM + y[1])
        if len(y)==3:  dy[2] = ke * (y[1] - y[2]) # An effect compartment, dy[2]

        return dy


class PK_2comp_Ka_Cl_ODE(PK_1comp_Ka_Cl_ODE):
    '''
     PK_2comp_Ka_Cl_ode : 2 compartments, linear absorption, linear elimination, ODE form
     (Subclass of SimFunction)

     Inputs:
     - constants, integer 0 or 1
       0: microConstants, default. ka, V, k, k12, k21, (ke)
       1: volumeClearance. ka, Cl, V1, Q, V2, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable C
       1: return both variables C and Ce (effect compartment)
     - pValues, vector of parameters (3 or 4 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','C2','Ce'}
                      Optional input argument. 
                      Defaults are {'G','C','C2'} or {'G','C','C2','Ce'} if effect comp.

     Fields defined in the superclass SimFunction
     - parameterNames, {'ka','V','k','k12','k21'} microConstants, no effect comp.
                       {'ka','V','k','k12','k21','ke'} microConstants, effect comp.
                       {'ka','Cl','V1','Q','V2'} volumeClearance, no effect comp.
                       {'ka','Cl','V1','Q','V2','ke'} volumeClearance, effect comp.
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above

     Function:
     Y=simulate(times, experiment, yInput, functionObj), simulates the object
        times - vector of time points to output
        experiment - the experiment to simulate (wrt dosing)
        yInput - input variable values at time points in times, size |times|x|inputVariables|
       obj from which inputVariables can be calculated at arbitrary time-points

     returns: Y matrix of output variable data, size |times|x|inputVariables|

    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['G','C','C2']):
        try:
            self.nbComp = 2
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','k','k12','k21']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
    
            if constants==0:
                pass # default
            elif constants==1:
                self.parameterNames = ['ka','Cl','V1','Q','V2']
            else:
                raise Exception('Error in PK_2comp_Ka_ODE, constants not in [0,1]')
    
            if effectCompartment > 0:
                self.calculateEffectCompartment = 1
                if 'Ce' not in self.outputVariables: 
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(4)={'Ce'};
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0
            
            self.parameterValues = pValues
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            if (len(outputVariables)== (self.nbComp+1) and effectCompartment == 0) or \
                (len(outputVariables) == (self.nbComp+2) and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ \
                            str(len(outputVariables)) + ': ' + str(outputVariables))

            self.k = self.parameterValues
            self.F = 1 #if self.constants != 1 else self.parameterValues[4] # pas dans Ka_Cl_ODE, pourquoi?

            #Matlab code: was in function simulate
            self.kn = self.k  # if applicable, ke is included here
            if self.constants==0:
                pass
            elif self.constants==1: # switch to microconstants
                # from 'ka','Cl','V1','Q','V2'
                # to 'ka','V','k','k12','k21'
                self.kn[1] = self.k[1]
                self.kn[2] = self.k[3]
                self.kn[3] = self.k[2] / self.k[3]
                self.kn[4] = self.k[4] / self.k[3]
                self.kn[5] = self.k[4] / self.k[5]
            else:
                raise Exception('Error in simfunctionODE, constants not in [0,1]')

        except:
            print('ERROR')  #FIXME

#    def simulate(self, times, experiment, ...) --> appel a classe parente

    def myodes(self, t, y, k):
        if len(y)==(self.nbComp+1):
            k01, V, k10, k12, k21 = k
        elif len(y)==(self.nbComp+2):
            k01, V, k10, k12, k21, ke = k  # remonte par rapport a Matlab
        else: raise Exception

        if len(y)==(self.nbComp+1): dy = np.zeros( (self.nbComp+1,), dtype=float) # No effect compartment
        elif len(y)==(self.nbComp+2): dy = np.zeros( (self.nbComp+2,), dtype=float)
        else: raise Exception

        # 'ka','V','k','k12','k21'
        dy[0] = -k01*y[0]
        dy[1] = k01*y[0]/V - k12*y[1] + k12*y[2] - k10*y[1]
        dy[2] = k21*y[1] - k21*y[2]
        if len(y)==(self.nbComp+2): # An effect compartment, dy[3]
            dy[3] = ke*(y[1] - y[3])

        return dy

    #============ Infusion
    def myodes_inf(self, t, y, k): #, experiment
        if len(y)==(self.nbComp+1):
            k01, V, k10, k12, k21 = k
        elif len(y)==(self.nbComp+2):
            k01, V, k10, k12, k21, ke = k  # remonte par rapport a Matlab
        else: raise Exception

        inp = 0
        if self.infusion_mode > 0:  # preparation pour futur merge des myodes
            for idx, tps in enumerate(experiment.dosing.times):
                if  t >= tps and t < (tps + experiment.dosing.duration[idx]):
                    inp = experiment.dosing.values[idx]
                    break

        if len(y)==(self.nbComp+1):   dy = np.zeros( (self.nbComp+1,), dtype=float) # No effect compartment
        elif len(y)==(self.nbComp+2): dy = np.zeros( (self.nbComp+2,), dtype=float)

        # 'ka','V','k'            
        dy[0] = -k01*y[0]  # FIXME: calcule pour rien
        dy[1] = inp / V - k12*y[1] + k12*y[2] - k10*y[1]
        dy[2] = k21*y[1] - k21*y[2]
        if len(y)==3:  dy[3] = ke*(y[1] - y[3])   # An effect compartment, dy[4]

        return dy


class PK_3comp_Ka_Cl_ODE(PK_1comp_Ka_Cl_ODE):
    '''
     PK_3comp_Ka_Cl_ode : 3 compartments, linear absorption, linear elimination, ODE form
     (Subclass of SimFunction)

     Inputs:
     - constants, integer 0 or 1
       0: microConstants, default. ka, V, k, k12, k21, k13, k31, (ke)
       1: volumeClearance. ka, Cl, V1, Q2, V2, Q3, V3, (ke)
     - effectCompartment, integer 0 or 1
       0: return the variable C
       1: return both variables C and Ce (effect compartment)
     - pValues, vector of parameters (3 or 4 parameters depending on the use 
              of an effect compartment)
     - outputVariables, vector of strings, e.g., {'G','C','C2','C3','Ce'}
                      Optional input argument. 
                      Defaults are {'G','C','C2','C3'} or {'G','C','C2','C3','Ce'} if effect comp.

     Fields defined in the superclass SimFunction
     - parameterNames, {'ka','V','k','k12','k21','k13','k31'} microConstants, no effect comp.
                       {'ka','V','k','k12','k21','k13','k31','ke'} microConstants, effect comp.
                       {'ka','Cl','V1','Q2','V2','Q3','V3'} volumeClearance, no effect comp.
                       {'ka','Cl','V1','Q2','V2','Q3','V3','ke'} volumeClearance, effect comp.
     - parameterValues, assigned from input pValues
     - inputVariables, []
     - outputVariables, as described above

     Function:
     Y=simulate(times, experiment, yInput, functionObj), simulates the object
        times - vector of time points to output
        experiment - the experiment to simulate (wrt dosing)
        yInput - input variable values at time points in times, size |times|x|inputVariables|
       obj from which inputVariables can be calculated at arbitrary time-points

     returns: Y matrix of output variable data, size |times|x|inputVariables|

    '''

    def __init__(self, constants, effectCompartment, pValues=None, outputVariables=['G','C','C2','C3']):
        try:
            self.nbComp = 3
            self.inputVariables = []
            self.outputVariables = outputVariables
            self.parameterNames = ['ka','V','k','k12','k21','k13','k31']
            self.constants = constants

            if pValues is None:
                raise pko.BadInputException('ResultChk:BadInput', ' At least three input arguments are required.')
    
            if constants==0:
                pass # default
            elif constants==1:
                self.parameterNames = ['ka','Cl','V1','Q2','V2','Q3','V3']
            else:
                raise Exception('Error in PK_2comp_Ka_ODE, constants not in [0,1]')
    
            if effectCompartment > 0:
                self.calculateEffectCompartment = 1
                if 'Ce' not in self.outputVariables: 
                    self.outputVariables.append('Ce')  #code origine force la position: obj.outputVariables(5)={'Ce'};
                self.parameterNames.append('ke')
            else:
                self.calculateEffectCompartment = 0
            
            self.parameterValues = pValues
            if len(pValues) != len(self.parameterNames):
                raise pko.BadInputException('ResultChk:BadInput', ' Number of parameters incorrect!')
            
            if (len(outputVariables)== (self.nbComp+1) and effectCompartment == 0) or \
                (len(outputVariables) == (self.nbComp+2) and effectCompartment == 1): 
                self.outputVariables = outputVariables
            else:
                raise pko.BadInputException('ResultChk:BadInput', ' Input argument outputVariables has incorrect len.'+ \
                            str(len(outputVariables)) + ': ' + str(outputVariables))

            self.k = self.parameterValues
            self.F = 1 #if self.constants != 1 else self.parameterValues[4] # pas dans Ka_Cl_ODE, pourquoi?

            #Matlab code: was in function simulate
            self.kn = self.k  # if applicable, ke is included here
            if self.constants==0:
                pass
            elif self.constants==1: # switch to microconstants
                # from 'ka','Cl','V1','Q2','V2','Q3','V3'
                # to 'ka','V','k','k12','k21','k13','k31'
                self.kn[1] = self.k[1]
                self.kn[2] = self.k[3]
                self.kn[3] = self.k[2] / self.k[3] # k12=Q2/V1
                self.kn[4] = self.k[4] / self.k[3] # k12=Q2/V1
                self.kn[5] = self.k[4] / self.k[5]
                self.kn[6] = self.k[6] / self.k[3] # k13=Q3/V1
                self.kn[7] = self.k[6] / self.k[7] # k31=Q3/V3
            else:
                raise Exception('Error in simfunctionODE, constants not in [0,1]')

        except:
            print('ERROR') #FIXME

#    def simulate(self, times, experiment, ...) --> appel a classe parente

    def myodes(self, t, y, k):
        if len(y)==(self.nbComp+1):
            k01, V, k10, k12, k21, k13, k31 = k
        elif len(y)==(self.nbComp+2):
            k01, V, k10, k12, k21, k13, k31, ke = k  # ke etait plus bas dans code Matlab
        else: raise Exception

        if len(y)==(self.nbComp+1): dy = np.zeros( (self.nbComp+1,), dtype=float) # No effect compartment
        elif len(y)==(self.nbComp+2): dy = np.zeros( (self.nbComp+2,), dtype=float)
        else: raise Exception

        # 'ka','V','k','k12','k21', 'k13', 'k31'
        dy[0] = -k01*y[0]
        dy[1] =  k01*y[0] / V - k12*y[1] + k21*y[2] - k13*y[1] + k31*y[3] - k10*y[1]  #corrige sur k12*y[2] et k13*y[3]
        dy[2] =  k12*y[1] - k21*y[2]  #etait k21*y[1]
        dy[3] =  k13*y[1] - k31*y[3]  #etait k31*y[1]
        if len(y)==(self.nbComp+2): # An effect compartment, dy[4]
            dy[4] = ke*(y[1] - y[4])

        return dy

    #============ Infusion
    def myodes_inf(self, t, y, k): #, experiment
        if len(y)==(self.nbComp+1):
            k01, V, k10, k12, k21, k13, k31 = k
        elif len(y)==(self.nbComp+2):
            k01, V, k10, k12, k21, k13, k31, ke = k  # remonte par rapport a Matlab
        else: raise Exception

        inp = 0
        if self.infusion_mode > 0:  # preparation pour futur merge des myodes
            for idx, tps in enumerate(experiment.dosing.times):
                if  t >= tps and t < (tps + experiment.dosing.duration[idx]):
                    inp = experiment.dosing.values[idx]
                    break

        if len(y)==(self.nbComp+1):   dy = np.zeros( (self.nbComp+1,), dtype=float) # No effect compartment
        elif len(y)==(self.nbComp+2): dy = np.zeros( (self.nbComp+2,), dtype=float)

        # 'ka','V','k','k12','k21', 'k13', 'k31'
        dy[0] = 0.0 #-k01*y[0]  # etait calcule pour rien
        dy[1] = inp / V - k12*y[1] + k21*y[2] - k13*y[1] + k31*y[3] - k10*y[1]
        dy[2] = k12*y[1] - k21*y[2]  # etait k21*y[1]
        dy[3] = k13*y[1] - k31*y[3]  # etait k31*y[1]
        if len(y)==3: dy[4] = ke*(y[1] - y[4])   # An effect compartment, dy[4]

        return dy

#================================================================================================================
#           EXAMPLES WHICH WERE PUT INITIALLY IN THE DOCSTRINGS
#================================================================================================================
if __name__ == "__main__": 

    # PK_1comp_Ka_Cl_ode & PK_1comp_Ka_satCl_ode

    timeUnit = 'min'
    doseUnit = 'umol'
    d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
    experiment = pko.Experiment('Ex',d,[])
    times = np.arange(0,73,1) # [0:1:72]
    const = 0

    #For each model
    for m in ['1comp', '1comp_satCl', '2comp', '3comp']:
        # With and without effectCompt activated
        for k in [0,1]:
            effectComp = k
            if k==0:
            # Example 1:
                if m=='1comp':
                    sf = PK_1comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84])
                elif m=='1comp_satCl':
                    sf = PK_1comp_Ka_satCl_ODE(const, effectComp, [0.86, 0.84, 1.6, 0.5])
                elif m=='2comp':
                    sf = PK_2comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84, 1.5, 0.1])
                elif m=='3comp':
                    sf = PK_3comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84, 1.5, 0.1, 1.4, 0.05])
            else:
            # Example 2: effectComp=1
                if m=='1comp':
                    sf = PK_1comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84, 1],     ['G','Cp','Ce'])
                elif m=='1comp_satCl':
                    sf = PK_1comp_Ka_satCl_ODE(const, effectComp, [0.86, 0.84, 1.6, 0.5, 1],['G','Cp','Ce'])
                elif m=='2comp':
                    sf = PK_2comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84, 1.5, 0.1, 1],['G','Cp','Cp2','Cpe'])
                elif m=='3comp':
                    sf = PK_3comp_Ka_Cl_ODE(   const, effectComp, [0.86, 1.6, 0.84, 1.5, 0.1, 1, 0.9, 1],['G','C','C2','C3','Ce'])
            Y = sf.simulate(times, experiment, [], [])
        
            fig, ax = plt.subplots()
            if m in ['1comp','1comp_satCl']:
                ax.plot(times, Y[:,0], 'r', times, Y[:,1], 'b')
                if k>1: ax.plot(times, Y[:,2], 'g')
            elif m in ['2comp']:
                pass
                #semilogy(times,[:,0],'g')
                #semilogy(times,[:,1],'b')
                #semilogy(times,Y[:,2],'r') 
                #if k>1: ax.plot(times,Y[:,3],'m')
            elif m in ['3comp']:
                pass
                #semilogy(times,Y[:,3],'k')
                #if k>1: ax.plot(times,Y[:,4],'k')

            ax.set_title('Sanity check')
            fig.legend(sf.outputVariables)
            plt.xlabel(timeUnit)
            plt.ylabel('uM')
            plt.show()
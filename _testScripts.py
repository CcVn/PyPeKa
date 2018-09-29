# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 17:28:04 2018
@author: Linhardt
"""
import matplotlib.pyplot as plt
import numpy as np
import PKObjects as pko
import PK_1comp_models as _pk
import PDmodels as _pd

class TestFailed(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, message, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
        print(message)

def testScripts(mytype=None):
    '''
    # testScripts
    # Test main functionality of scripts
    # Input:
    # mytype 0: all tests (default)
    #      1: Dosing (ok if no output)
    #      2: DataVariable (ok if no output)
    #      3: Experiment (ok if no output)
    #      4: Model (ok if no output)
    #      5: pko.SimFunction (9 plt.plots should be drawn)
    #      6: Parameter estimation iv (ok if no output except 1 warning)
    #      7: Parameter estimation ev
    #      8: Parameter estimation inf (ok if no output except 1 warning)
    #      9: Parameter estimation iv+ev (ok if no output)
    #      10: niceplt.plot, niceFigure, niceSemilogy (3 plt.plots should be drawn)
    #      11: _pk.PK_1/2comp_Ka_Cl_fun/ode  (ok if no output)
    #      12: PD models (ok if no output)
    #      13: PK closed form steady state equations (3 plt.plots should be drawn)
    #
    # Example: testScripts, to execute all tests
    #          testScripts(2), to execute test of the DataVariable class
    '''
    # INIT
    plt.close('all')
    timeUnit = 'min'
    valueUnit = 'uM'
    doseUnit = 'umol/kg'
    MYTYPE = mytype if mytype else 0


    # TEST DOSING
    if MYTYPE==0 or MYTYPE==1:
        print('== Test  1')
        d = pko.Dosing('ev', [0, 24], [10, 13], [], timeUnit, doseUnit)
        d = pko.Dosing('ev', [0], [10], [], timeUnit, doseUnit)
        d = pko.Dosing('ev', [0, 12, 24], [10, 11, 13], [], timeUnit, doseUnit)
        try:
            d = pko.Dosing('ev', [0, 24, 23], [10, 13, 1], [], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1a')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            raise TestFailed('Error, Dosing, test 1b')
        try:
            d = pko.Dosing('ev', [0, 24], [10, 13, 1], [], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1c')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1d')
        try:
            d = pko.Dosing('ev', [0, 24], [10, 13], [], timeUnit, doseUnit)  #[0.24] initialement transposed
            raise TestFailed('Error, expected an error message in Dosing, test 1e')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1f')
        try:
            d = pko.Dosing('ev', [0, 24], [10, 13].T, [], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1g')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1h')
        d = pko.Dosing('iv', [0, 24], [10, 13], [], timeUnit, doseUnit)
        d = pko.Dosing('iv', [0], [10], [], timeUnit, doseUnit)
        d = pko.Dosing('iv', [0, 12, 24], [10, 11, 13], [], timeUnit, doseUnit)
        d = pko.Dosing('inf', [0, 24], [10, 13], [4, 5], timeUnit, doseUnit)
        d = pko.Dosing('inf', [0], [10], [6], timeUnit, doseUnit)
        d = pko.Dosing('inf', [0, 12, 24], [10, 11, 13], [6, 8, 10], timeUnit, doseUnit)
        try:
            d = pko.Dosing('inf', [0, 10, 20], [10, 13, 1], [5, 11, 2], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1i')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1j')
        try:
            d = pko.Dosing('inf', [0, 24], [10, 13, 1], [2, 2], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1k')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1l')
        try:
            d = pko.Dosing('inf', [0, 24], [10, 13], [2, 2].T, timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1m')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1, test 1n')
        d = pko.Dosing('ev', [0, 24], [10, 15], [], 'h', 'umol/kg')
        d = pko.Dosing('iv', [0, 10, 20], [5, 10, 15], [], 'min', 'umol/kg')
        d = pko.Dosing('inf', [0, 60], [5, 10], [10, 15], 'min', 'umol/kg/min')
        try:
            d = pko.Dosing('aa', [0, 24], [10, 13, 1], [2, 2], timeUnit, doseUnit)
            raise TestFailed('Error, expected an error message in Dosing, test 1o')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Dosing, test 1p')


    # TEST DATA VARIABLE
    if MYTYPE==0 or MYTYPE==2:
        print('== Test  2')
        t = [1, 5, 10, 20, 30, 40, 50, 60]
        c = [0.2, 6, 8, 7, 5, 0.3, 0.09, 0.004]
        c01 = [i*0.1 for i in c]
        dv1 = pko.DataVariable('C', t, c, c01, timeUnit, valueUnit)
        t = [1, 20, 30, 40, 50, 55, 60]
        r = [0.1, 0.8, 0.7, 0.55, 0.3, 0.21, 0.1]
        r01 = [i*0.1 for i in r]
        dv2 = pko.DataVariable('R', t, r, r01, timeUnit, valueUnit)
        t = [1, 5, 10, 20, 30, 35, 50, 60]
        q = [0.15, 4, 5, 20, 55, 56, 60, 43]

#FIXME dv3 / dv4 /dv 5 inutilises???  dv5 deleted
        q01 = [i*0.1 for i in q]
        dv3 = pko.DataVariable('Q',t,q,q01, timeUnit, valueUnit)
        print('NB on origin test set: dv3 not used', dv3)
        t = [1, 10, 60]
        c = [0.2, 0.3, 0.04]
        c_std=[0.02, 0.06, 0.008]

        dv4 = pko.DataVariable('C', t, c, c_std, 'min', 'uM')
        print('NB on origin test set: dv4 not used', dv4)
        t = [1, 10, 60]
        c = [0.2, 0.3, 0.04]
        c_std=[]
        try:
            dv1 = pko.DataVariable('C', [1, 2], [1, 2, 3], [], timeUnit, valueUnit)
            raise TestFailed('Error, expected an error message in DataVariable, test 2a')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            print('Error, DataVariable, test 2b')
        try:
            dv1 = pko.DataVariable('C', [1, 2, 3], [1, 2, 3], [], timeUnit, valueUnit)  #FIXME : was 'C', [1, 2, 3]' in Matlab (transposed)
            raise TestFailed('Error, expected an error message in DataVariable, test 2c')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            print('Error, DataVariable, test 2d')
        try:
            dv1 = pko.DataVariable('C', [1, 2, 3], [1, 2, 3], [0.1, 0.2, 0.3], timeUnit, valueUnit)  #FIXME : was [1, 2, 3], [1, 2, 3]' in Matlab (transposed)
            raise TestFailed('Error, expected an error message in DataVariable, test 2e')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, DataVariable, test 2f')
        try:
            dv1 = pko.DataVariable('C', [1, 2, 3], [1, 2, 3], [0.1, 0.2, 0.3], timeUnit, valueUnit) #FIXME : was [0.1, 0.2, 0.3]' in Matlab (transposed)
            raise TestFailed('Error, expected an error message in DataVariable, test 2g')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, DataVariable, test 2h')
        try:
            dv1 = pko.DataVariable('C', [1, 3, 2], [1, 2, 3], [0.1, 0.2, 0.3], timeUnit, valueUnit)
            raise TestFailed('Error, expected an error message in DataVariable, test 2i')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, DataVariable, test 2j')


    # TEST EXPERIMENT
    if MYTYPE==0 or MYTYPE==3:
        print('== Test 3 : EXPERIMENT')
        dose = pko.Dosing('ev', [0, 24], [10, 15], [], 'h', 'umol/kg')
        t = [1, 10, 60]
        c = [0.2, 0.3, 0.04]
        c_std=[0.02, 0.06, 0.008]
        dv1 = pko.DataVariable('C', t, c, c_std, 'min', 'uM')
        t = [1, 20, 60]
        r=[0.1, 0.8, 0.7]
        r_std=[0.01, 0.08, 0.07]
        dv2 = pko.DataVariable('R', t, r, r_std, 'min', '%')
        e = pko.Experiment('test', dose, [dv1, dv2], 'r')
        e.harmonizeTimes()
        if sum(e.dataVariables[1].values==[0.2000, 0.3000 , np.nan, 0.0400])!=3:
            raise TestFailed('Error, Experiment')
        if sum(e.dataVariables[1].stdev==[0.02, 0.06, np.nan, 0.008])!=3:
            raise TestFailed('Error, Experiment')
        if sum(e.dataVariables[2].values==[0.1, np.nan, 0.8, 0.7])!=3:
            raise TestFailed('Error, Experiment')
        if sum(e.dataVariables[2].stdev==[0.01, np.nan, 0.08, 0.07])!=3:
            raise TestFailed('Error, Experiment')
        try:
            e = pko.Experiment('test', 'wrong', [dv1, dv2], 'r')
            raise TestFailed('Error, expected an error message in Experiment')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Experiment')
        try:
            e = pko.Experiment('test', dose, [dose, dose], 'r')
            raise TestFailed('Error, expected an error message in Experiment')
        except pko.BadInputException:
            pass # ok
        except Exception as e:
            print(e)
            raise TestFailed('Error, Experiment')


    # TEST MODEL
    if MYTYPE==0 or MYTYPE==4:
        print('== Test  4 : MODEL')

        # test -1
        f = _pk.PK_1comp_Ka_Cl_fun(0,0, [0.1, 1, 3])
        fode = _pd.PD_receptocc_ODE([1.0, 0.5, 100])
        m = pko.Model(f, fode, None)
        p = m.getParameterValues() # assigns p=[0.10 1.0 3.0 1.0 0.50 100.0]
        p[3] = 3.5
        m.setParameterValues(p) # assigns [0.10 1.0 3.5 1.0 0.50 100.0] to Model
        if not isSimilar(m.getParameterValues(), [0.10, 1.0, 3.5, 1.0, 0.50, 100.0],0.0001):
            raise TestFailed('Error, Model test -1')
        timeUnit = 'min'
        doseUnit = 'umol'
        valueUnit = 'uM'
        k_init = [1.5, 0.3, 0.17, 0.15]
        time_c = [0, 5, 10, 15, 20, 25, 30, 35, 40]
        print('time_c non used', time_c)
        c = [0,45.1427,16.9112,6.22335,2.28945,0.842243,0.309844,22.6853,8.49755]
        time_ce = [0, 5, 10, 15, 20, 30, 40]
        print('time_ce non used', time_ce)
        ce = [0,21.5009,23.9592,18.5689,12.7478,5.22205,13.9727]
        dosing = pko.Dosing('ev', [0, 30], [20, 10], [], timeUnit, doseUnit)
        print('dosing non used', dosing)
        model = pko.Model(_pk.PK_1comp_Ka_Cl_fun(0,1,k_init), None, None)

        # test 0
        times = [k for k in range(25)]
        experiment = []
        yInput = []
        functionObj=[]
        k = [0.1, 0.5]
        sf = pko.SimFunction(['k1','k11'],k,[],['C','Ce'])
        k = [0.6, 6]
        sf2 = pko.SimFunction(['k2','k22'],k,['E'],['A','B'])
        k = [0.7, 7]
        sf3 = pko.SimFunction(['k3','k33'],k,['B'],['E','F'])
        try:
            m1 = pko.Model(sf, sf2, sf3)
            raise TestFailed('Error, Model')
        except pko.BadInputException:
            pass # ok

        # TEST 1
        times = [k for k in range(25)]
        experiment = []
        yInput = []
        functionObj=[]
        k = [0.1, 0.5]
        sf = pko.SimFunction(['k1','k11'],k,[],['C','Ce'])
        Y = sf.simulate(times, experiment, yInput, functionObj)
        k = [0.6, 6]
        sf2 = pko.SimFunction(['k2','k22'],k,['Ce'],['A','B'])
        Y = sf2.simulate(sf, times, experiment, yInput, functionObj)
        k = [0.7, 7]
        sf3 = pko.SimFunction(['k3','k33'],k,['B'],['E','F'])
        Y = sf3.simulate(sf, times, experiment, yInput, functionObj)
        m1 = pko.Model(sf, sf2, sf3)
        multValue = 10
        #print('Get current k vector from Model')
        k = m1.getParameterValues()
        #print('Modify k vector: Multiply by 10.')
        knew = k * multValue
        #print('Set modified k vector to Model')
        m1.setParameterValues(knew)
        #print('Get current k vector from Model')
        kget = m1.getParameterValues()
        if knew == kget:
            pass # ok
        else:
            raise TestFailed('Error, Model test 1')

        # TEST 2
        times = [k for k in range(25)]
        experiment = []
        yInput = []
        functionObj=[]
        k = [0.1, 0.5]
        sf = pko.SimFunction(['k1','k11'],k,[],['C','Ce'])
        Y = sf.simulate(times, experiment, yInput, functionObj)
        k = [0.1, 6]
        sf2 = pko.SimFunction(['k1','k22'],k,['Ce'],['A','B'])
        Y = sf.simulate(times, experiment, yInput, functionObj)
        k = [0.5, 6]
        sf3 = pko.SimFunction(['k11','k22'],k,['B'],['E','F'])
        Y = sf.simulate(times, experiment, yInput, functionObj)
        m1 = pko.Model(sf, sf2, sf3)
        m1.Parameters
        #print('Get current k vector from Model')
        k = m1.getParameterValues()
        #print('Modify k vector')
        k=k*10
        #print('Set modified k vector to Model')
        #m1 = m1.setParameterValues(k)
        m1.setParameterValues(k)
        #print('Get current k vector from Model')
        k = m1.getParameterValues()
        m1.Function2.parameterNames
        x = m1.Function2.parameterValues
        if not isSimilar(x, [5, 60],0.0001):
            raise TestFailed('Error, Model test 2')
        m1.ODE.parameterNames
        x = m1.ODE.parameterValues
        if not isSimilar(x, [1, 60],0.0001):
            raise TestFailed('Error, Model test 2')

        # TEST 3
        times = [k for k in range(25)]
        experiment = []
        yInput = []
        functionObj=[]
        k = [0.1, 0.5]
        sf = pko.SimFunction(['k1','k11'],k,[],['C','Ce'])
        k = [0.6, 6]
        sf2 = pko.SimFunction(['k2','k22'],k,['Ce'],['A','B'])
        k = [0.5, 6]
        sf3 = pko.SimFunction(['k11','k22'],k,['B'],['E','F'])
        m1 = pko.Model(sf, sf2, sf3)
        k=np.ones(1,3)
        #print('Set modified k vector to Model')
# =============================================================================
#         success3 = False  # inutile???????
# =============================================================================
        try:
            m1 = m1.setParameterValues(k)
            raise TestFailed('Error, Model test 3')
        except pko.BadInputException:
            pass

        # TEST 4
        sf = _pk.PK_2comp_Ka_Cl_fun(0,0,np.ones(5,1))
        sfode = _pd.PD_receptocc_ODE([1, 1, 100])
        try:
            m1 = pko.Model(sf, sfode, None)
        except pko.BadInputException:
            raise TestFailed('Error, Model test 4')

        # TEST 5
        timeUnit = 'min'
        doseUnit = 'umol/kg/min'
        valueUnit = 'uM'
        d = pko.Dosing('inf', [0, 30], [5, 3], [4, 10], timeUnit, doseUnit)
        times = [0, 1, 2, 3, 4] +  [i for i in range(5,61,5)]
        Y1=[0, 14.2226, 20.0143, 24.0585, 27.4752, 16.9486, 12.2658,
            10.0389, 8.5605, 7.1696, 5.9074, 23.5884, 32.5343, 20.5914,
            17.1878, 14.2673, 11.8340]
        Y2=[0, 0.8288, 0.8613, 0.8887, 0.8956, 0.8382, 0.8152,
            0.7786, 0.7372, 0.7105, 0.6640, 0.8915, 0.9184, 0.8749,
            0.8509, 0.8264, 0.7933]
        dv1 = pko.DataVariable('Cu', times, Y1, [], timeUnit, valueUnit)
        dv2 = pko.DataVariable('E', times, Y2, [], timeUnit, valueUnit)
        e = pko.Experiment('testExp',d, [dv1, dv2])
        kpk = [1, 0.2, 0.2, 1.2, 0.3, 1] # true PK parameters
        pk = _pk.PK_2comp_Ka_Cl_fun(0,1,kpk,['Cu','Ce'])
        m = pko.Model(pk, None, None)
        lb = kpk*0.1
        ub = kpk*10
        lb[1]=1
        ub[1]=1
        [kpk,residual,ci,se,cv]=m.fit([e],lb,ub,0)
        # TEST 6
        doseTbid=[0, 7]
        for i in [1,2,3,4]:
            doseTbid = [doseTbid, i*24, i*24+7]
        doseT=doseTbid
        a=[  1.0934, 2.2054, 2.5035, 2.5071, 2.4038, 2.2622, 2.1080, 1.9507, 3.0503,
             3.3416, 3.3383, 3.2276, 3.0781, 2.9152, 2.7484, 2.5811, 2.4149, 2.2504,
             2.0879, 1.9278, 1.7704, 1.6161, 1.4651, 1.3180, 1.1753]
        doseV=np.ones(1,len(doseT))*1
        d1 = pko.Dosing('ev', [doseT], [doseV], [],'h','umol/kg')
        mestemps = [i for i in range(96,96+25,1)]
        dv11 = pko.DataVariable('Cp',mestemps,a.T, [],'h','uM')
        e1 = pko.Experiment('A',d1, [dv11],'k')
        exps=[e1]
        effectComp=0
        const = 0
        #kPk = [1 0.5 0.1 0.5]#'ka','V','Vmax','KM' CORRECT PARAS
        kPK = [1, 1 ,0.05, 0.7]#'ka','V','Vmax','KM'
        sf = _pk.PK_1comp_Ka_satCl_ODE(const, effectComp, kPK, ['G','Cp'])
        model = pko.Model(None, sf, None)
        #     Y=model.simulate(mestemps,exps[1])
        #     plt.figure(10000)
        #     plt.plot(mestemps,Y,'o-')
        #     Y[:,2]
        lb =  kPK*0.00001 # lower parameter bounds
        ub =  kPK*100000 # upper parameter bounds
        lb[1] = 1
        ub[1] = 1
        model.fit(exps,lb,ub,3)
        if not isSimilar(model.getParameterValues(), [1, 0.5, 0.1, 0.5],0.01):
            raise TestFailed('Error, Model test -4:6')
            print('Parameters: ', str(model.getParameterValues()))

    # TEST pko.SimFunction
    if MYTYPE==0 or MYTYPE==5:
        print('== Test  5 : pko.SimFunction')
        plt.close('all')
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [22, 8, 2, 5], timeUnit, doseUnit)
        experiment = pko.Experiment('testExp',d, [])
        times = [i*0.5 for i in range(72*2+1)]
        print('    Figure 1-3 tentatively successful if solid graphs are overwritten by')
        print('    dotted graphs in the plt.plots, and no error messages are reported.')
        # ONE COMP
        # 'ka','V','k'
        kmicro =  [ 0.85663, 1.5958, 0.8461 ]
        kmicro2 = [ 0.85663, 1.5958, 0.8461, 0.3 ]
        # 'ka','Cl','V'
        ka = kmicro[1]
        Cl = kmicro(2)*kmicro[3]
        V = kmicro(2)
        k = [ka,Cl,V]
        k2 = [ka, Cl, V, 1]
        adm = ['ev', 'iv', 'inf']
        plt.figure(1)
        for i in [1,2,3]:
            experiment.dosing.admin = str(adm[i])
            calculateEffectCompartment = 0
            for iconstant in [0,1]:
                plt.plot(3,4,(i-1)*4+iconstant+1)
                constants=iconstant
                sf1 = _pk.PK_1comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                sf2 = _pk.PK_1comp_Ka_Cl_ODE(constants,calculateEffectCompartment,kmicro)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro
                else:
                    sf1.parameterValues=k
                Y = sf1.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                del Y
                if iconstant==0:
                    sf2.parameterValues = kmicro
                else:
                    sf2.parameterValues=k
                Y = sf2.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Y[:,2],'g:')
                plt.setp([h1, h2],'LineWidth',2)
                plt.legend('C fun','C ode')
                plt.title('1 comp ', str(adm[i]), ' , no effectComp')
            calculateEffectCompartment = 1
            for iconstant in [0,1]:
                plt.plot(3,4,(i-1)*4+2+iconstant+1)
                constants=iconstant
                sf1 = _pk.PK_1comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro2)
                sf2 = _pk.PK_1comp_Ka_Cl_ODE(constants,calculateEffectCompartment,kmicro2)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro2
                else:
                    sf1.parameterValues=k2
                Y = sf1.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y[:,1],'b')
                plt.hold(True)
                h2 = plt.semilogy(times,Y[:,2],'k')
                del Y
                if iconstant==0:
                    sf2.parameterValues = kmicro2
                else:
                    sf2.parameterValues=k2
                Y = sf2.simulate(times,experiment, [], [])
                h3 = plt.semilogy(times,Y[:,2],'r:')
                h4 = plt.semilogy(times,Y[:,3],'g:')
                plt.setp([h1, h2, h3, h4],'LineWidth',2)
                plt.legend('C fun','Ce fun','C ode','Ce ode')
                plt.title('1 comp ', str(adm[i]), ' , effectComp')
        # TWO COMP
        # 'ka','V','k','k12','k21'
        kmicro =  [ 0.85663, 1.5958, 0.8461, 0.05013, 0.11506]
        kmicro2 = [ 0.85663, 1.5958, 0.8461, 0.05013, 0.11506, 0.2 ]
        # 'ka','Cl','V1','Q','V2'
        ka = kmicro[1]
        Cl = kmicro[2]*kmicro[3]
        V1 = kmicro[2]
        Q = kmicro[4]*kmicro[2]
        V2 = kmicro[4]*kmicro[2]/kmicro[5]
        k =  [ ka, Cl, V1, Q, V2]
        k2 = [ ka, Cl, V1, Q, V2, 0.2]
        plt.figure(2)
        for i in [1,2,3]:
            experiment.dosing.admin=str(adm[i])
            calculateEffectCompartment = 0
            for iconstant in [0,1]:
                plt.plot(3,4,(i-1)*4+iconstant+1)
                constants=iconstant
                sf1 = _pk.PK_2comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                sf2 = _pk.PK_2comp_Ka_Cl_ODE(constants,calculateEffectCompartment,kmicro)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro
                else:
                    sf1.parameterValues=k
                Y = sf1.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                del Y
                if iconstant==0:
                    sf2.parameterValues = kmicro
                else:
                    sf2.parameterValues=k
                Y = sf2.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Y[:,2],'g:')
                plt.setp([h1, h2],'LineWidth',2)
                plt.legend('C fun','C ode')
                plt.title('2 comp ', str(adm[i]), ' , no effectComp')
            calculateEffectCompartment = 1
            for iconstant in [0,1]:
                plt.plot(3,4,(i-1)*4+2+iconstant+1)
                constants=iconstant
                sf1 = _pk.PK_2comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro2)
                sf2 = _pk.PK_2comp_Ka_Cl_ODE(constants,calculateEffectCompartment,kmicro2)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro2
                else:
                    sf1.parameterValues=k2
                Y = sf1.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y[:,1],'b')
                plt.hold(True)
                h2 = plt.semilogy(times,Y[:,2],'k')
                del Y
                if iconstant==0:
                    sf2.parameterValues = kmicro2
                else:
                    sf2.parameterValues=k2
                Y = sf2.simulate(times,experiment, [], [])
                h3 = plt.semilogy(times,Y[:,2],'r:')
                h4 = plt.semilogy(times,Y[:,4],'g:')
                plt.setp([h1, h2, h3, h4],'LineWidth',2)
                plt.legend('C fun','Ce fun','C ode','Ce ode')
                plt.title('2 comp ', str(adm[i]), ' , effectComp')
        # 3 COMP
        # 'ka','V','k','k12','k21',k13,k31
        kmicro =  [ 0.85663, 1.5958, 0.8461, 0.05013, 0.11506, 0.07, 0.15     ]
        kmicro2 = [ 0.85663, 1.5958, 0.8461, 0.05013, 0.11506, 0.07, 0.15, 0.2]
        # 'ka','Cl','V1','Q','V2'
        ka = kmicro[1]
        Cl = kmicro[2]*kmicro[3]
        V1 = kmicro[2]
        Q = kmicro[4]*kmicro[2]
        V2 = kmicro[4]*kmicro[2]/kmicro[5]
        Q3 = kmicro[6]*kmicro[2]
        V3 = kmicro[6]*kmicro[2]/kmicro[7]
        k = [ka,Cl,V1,Q,V2,Q3,V3]
        k2 = [ka,Cl,V1,Q,V2,Q3,V3,0.2]
        plt.figure(3)
        for i in [1,2,3]:
            experiment.dosing.admin = str(adm[i])
            calculateEffectCompartment = 0
            for iconstant in [0,1]:
                plt.plot(3,2,(i-1)*2+iconstant+1)
                constants=iconstant
                sf1 = _pk.PK_3comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                sf2 = _pk.PK_3comp_Ka_Cl_ODE(constants,calculateEffectCompartment,kmicro)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro
                else:
                    sf1.parameterValues=k
                Y = sf1.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                del Y
                if iconstant==0:
                    sf2.parameterValues = kmicro
                else:
                    sf2.parameterValues=k
                Y = sf2.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Y[:,2],'g:')
                plt.setp([h1, h2],'LineWidth',2)
                plt.legend('C fun','C ode')
                plt.title('3 comp ', str(adm[i]), ' , no effectComp')
        print('    Figure 7-9 tentatively successful if two solid graphs are shown\n', \
              '    and no error messages are reported.')

        # F1 and ODE
        i=7
        #print('Test pko.SimFunction, produce figure ' int2str[i])
        plt.figure(i)
        experiment.dosing.admin='ev'
        sf = _pk.PK_2comp_Ka_Cl_fun(0,0,np.ones(5,1))
        sfode = _pd.PD_receptocc_ODE([1, 1, 100])
        m1 = pko.Model(sf, sfode, None)
        Y = m1.simulate(times, experiment)
        plt.plot(times,Y[:,1],'b')
        plt.hold(True)
        plt.plot(times,Y[:,2],'r')
        m1.getOutputVariables()

        # F1 and F2
        i=8
        #print('Test pko.SimFunction, produce figure ' int2str[i])
        plt.figure(i)
        experiment.dosing.admin='ev'
        sf = _pk.PK_2comp_Ka_Cl_fun(0,0,np.ones(5,1))
        sf2 = _pd.PD_Emax_fun([3, 2])
        m5 = pko.Model(sf, None, sf2)
        Y = m5.simulate(times, experiment)
        plt.plot(times,Y[:,1],'b')
        plt.hold(True)
        plt.plot(times,Y[:,2],'r')
        m5.getOutputVariables()

        # ODE and F2
        i=9
        #print('Test pko.SimFunction, produce figure ' int2str[i])
        plt.figure(i)
        experiment.dosing.admin='ev'
        sfo = _pk.PK_2comp_Ka_Cl_ODE(0,0,np.ones(5,1))
        sf2 = _pd.PD_Emax_fun([3, 2])
        m6 = pko.Model(None, sfo, sf2)
        Y = m6.simulate(times, experiment)
        plt.plot(times,Y[:,2],'b')
        plt.hold(True)
        plt.plot(times,Y[:,4],'r')
        m6.getOutputVariables()
        # ONE COMP NON-LINEAR
        # 'ka','V','k', 'KM'
        kmicro  = [ 0.85663, 1.5958, 0.8461, 2]
        kmicro2 = [ 0.85663, 1.5958, 0.8461, 2, 0.05 ]
        # 'ka','Cl','V'
        ka = kmicro[1]
        Cl = kmicro[2]*kmicro[3]
        V = kmicro[2]
        KM = kmicro[4]
        k = [ka, Cl, V, KM]
        k2 = [ka, Cl, V, KM, kmicro2[-1]]
        print('    Figure 11-13 tentatively successful if no error messages are reported.')
        for i in [11,12,13]:
            #print('Test pko.SimFunction, produce figure ', int2str[i])
            plt.figure(i)
            if i==13: experiment.dosing.admin='inf'
            if i==12: experiment.dosing.admin='iv'
            if i==11: experiment.dosing.admin='ev'
            plt.figure(i)
            calculateEffectCompartment = 0
            for iconstant in [0]:
                plt.plot(2,1,1)
                constants=iconstant
                sf1 = _pk.PK_1comp_Ka_satCl_ODE(constants,calculateEffectCompartment,kmicro)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro
                Y = sf1.simulate(times,experiment, [], [])
                plt.plot(times,Y[:,2],'b-')
                plt.hold(True)
                plt.legend('C ode')
            calculateEffectCompartment = 1
            for iconstant in [0]:
                plt.plot(2,1,2)
                constants=iconstant
                sf1 = _pk.PK_1comp_Ka_satCl_ODE(constants,calculateEffectCompartment,kmicro2)
                del Y
                if iconstant==0:
                    sf1.parameterValues = kmicro2
                Y = sf1.simulate(times,experiment, [], [])
                plt.plot(times,Y[:,2],'r:')
                plt.hold(True)
                plt.plot(times,Y[:,3],'m-')
                del Y
                plt.legend('C ode','Ce ode')

    # TEST PARAMETER ESTIMATION IV
    if MYTYPE==0 or MYTYPE==6:
        print('== Test  6')
        times = [i for i in range(0,73,1)]
        k = [1, 0.2, 0.2, 1.2, 0.3, 0.1]
        # SET UP EXP 1-2
        d1 = pko.Dosing('iv', [0, 30], [20, 10], [], timeUnit, doseUnit)
        e1 = pko.Experiment('testExp1',d1, [])
        sf1 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m1 = pko.Model(sf1, None, None)
        Y = m1.simulate(times, e1)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e1 = pko.Experiment('testExp1',d1, [dv1])
        e2 = pko.Experiment('testExp2',d1, [dv1, dv2])
        #e1.display
        #e1.plt.plot[1]
        # SET UP EXP 3
        del Y
        d3 = pko.Dosing('iv', [0], [50], [], timeUnit, doseUnit)
        e3 = pko.Experiment('testExp3',d3, [])
        sf3 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m3 = pko.Model(sf3, None, None)
        Y = m3.simulate(times, e3)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e3 = pko.Experiment('testExp3',d3, [dv1, dv2])
        #e3.display
        #e3.plt.plot(2)
        # Test 1
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        print('    Fit should warn below: Warning: Matrix is singular...')
        x = m1.fit([e1],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, np.nan], x, 0.001):
            raise TestFailed('Error, parameter estimation iv, test 1')
        # Test 2
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m1.fit([e1, e2],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation iv, test 2')
        # Test 3
        m3.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m1.fit([e3],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation iv, test 3')
        # Test 4
        del Y
        d4 = pko.Dosing('iv', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e4 = pko.Experiment('testExp4',d4, [])
        sf4_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5])
        sf4_2 = _pd.PD_Emax_fun([3, 1])
        m4 = pko.Model(sf4_1, None, sf4_2)
        Y = m4.simulate(times, e4)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e4 = pko.Experiment('testExp4',d4, [dv1, dv2])
        #e4.display
        #e4.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5, 6.5]
        m4.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m4.fit([e4],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 3.000, 1.000], x, 0.001):
            raise TestFailed('Error, parameter estimation iv, test 4')
        # Test 5
        del Y
        d5 = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e5 = pko.Experiment('testExp5',d4, [])
        sf5_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5],['C'])
        sf5_2 = _pd.PD_Emax_fun([3, 1],'C','E')
        m5 = pko.Model(sf5_1, None, sf5_2)
        Y = m6.simulate(times, e5)
        dv1 = pko.DataVariable('Cf', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ef', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e5 = pko.Experiment('testExp5',d5, [dv1, dv2])
        dv1 = pko.DataVariable('Cf', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ef', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e55 = pko.Experiment('testExp5',d5, [dv1, dv2])
        #e5.display
        #e5.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5,6.5]
        m5.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        try:
            m5.fit([e5, e55],lb,ub,0)
            raise TestFailed('Error, parameter estimation iv, test 5')
        except pko.BadInputException:
            pass # ok
        #         print('model output variables')
        #         m5.getOutputVariables
        #         print('data variable names')
        #         e5.dataVariables[1].name
        #         e5.dataVariables[2].name
        # SET UP EXP 6
        k = [1, 0.2, 0.2, 0.1]
        d6 = pko.Dosing('iv', [0, 30], [20, 10], [], timeUnit, doseUnit)
        e6 = pko.Experiment('testExp6',d6, [])
        sf6 = _pk.PK_1comp_Ka_Cl_fun(0,1,k)
        m6 = pko.Model(sf6, None, None)
        Y = m6.simulate(times, e6)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e6 = pko.Experiment('testExp6',d6, [dv1, dv2])
        #e6.display
        #e6.plt.plot[1]
        # Test 6
        m6.setParameterValues([1.1288, 0.2241, 0.1928, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m6.fit([e6],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation iv, test 6')
    # TEST PARAMETER ESTIMATION EV
    if MYTYPE==0 or MYTYPE==7:
        print('== Test  7')
        times = [i for i in range(73)]
        k = [1, 0.2, 0.2, 1.2, 0.3, 0.1]
        # SET UP EXP 1-2
        d1 = pko.Dosing('ev', [0, 30], [20, 10], [], timeUnit, doseUnit)
        e1 = pko.Experiment('testExp1',d1, [])
        sf1 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m1 = pko.Model(sf1, None, None)
        Y = m1.simulate(times, e1)
        #Y=Y+0.01*Y*randn(len(times),2)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e1 = pko.Experiment('testExp1',d1, [dv1])
        e2 = pko.Experiment('testExp2',d1, [dv1, dv2])
        #e2.display
        #e2.plt.plot[1]
        # SET UP EXP 3
        del Y
        d3 = pko.Dosing('ev', [0], [50], [], timeUnit, doseUnit)
        e3 = pko.Experiment('testExp3',d3, [])
        sf3 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m3 = pko.Model(sf3, None, None)
        Y = m3.simulate(times, e3)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e3 = pko.Experiment('testExp3',d3, [dv1, dv2])
        #e3.display
        #e3.plt.plot(2)
        # Test 1
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        [x,residual,ci,se,cv] =m1.fit([e1],lb,ub,0)
        if not isSimilar([1.0000, 0.2000, 0.2000, 1.2000, 0.3000, np.nan], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 1')
        # test 2
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x =m1.fit([e1, e2],lb,ub,0)
        if not isSimilar([1.0000, 0.2000, 0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 2')
        # test 3
        m3.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m1.fit([e3],lb,ub,3)
        if not isSimilar([1.0000, 0.2000,  0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 3')
        # test 4
        del Y
        d4 = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e4 = pko.Experiment('testExp4',d4, [])
        sf4_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5])
        sf4_2 = _pd.PD_Emax_fun([3, 1])
        m4 = pko.Model(sf4_1, None, sf4_2)
        Y = m4.simulate(times, e4)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e4 = pko.Experiment('testExp4',d4, [dv1, dv2])
        #e4.display
        #e4.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5, 6.5]
        m4.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m4.fit([e4],lb,ub,0)
        if not isSimilar([1.0000, 0.2000, 0.2000, 1.2000, 0.3000, 3.000, 1.000], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 4')
        # test 5
        del Y
        d5 = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e5 = pko.Experiment('testExp4',d4, [])
        sf5_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5])
        sf5_2 = _pd.PD_Emax_fun([3, 1])
        m5 = pko.Model(sf5_1, None, sf5_2)
        Y = m6.simulate(times, e5)
        dv1 = pko.DataVariable('Cf', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Rf', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e5 = pko.Experiment('testExp4',d5, [dv1, dv2])
        #e5.display
        #e5.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5, 6.5]
        m5.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        try:
            x = m5.fit([e5],lb,ub,3)
            raise TestFailed('Error, , parameter estimation ev, test 5')
        except pko.BadInputException:
            pass # ok
        timeUnit = 'min'
        doseUnit = 'umol'
        times = [i for i in range(0,51,5)] #[0:5:50]
        # SET UP EXP 6
        k = [1, 0.2, 0.2, 0.1]
        d6 = pko.Dosing('ev', [0, 30], [20, 10], [], timeUnit, doseUnit)
        e6 = pko.Experiment('testExp6',d6, [])
        sf6 = _pk.PK_1comp_Ka_Cl_fun(0,1,k)
        m6 = pko.Model(sf6, None, None)
        Y = m6.simulate(times, e6)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e6 = pko.Experiment('testExp6',d6, [dv1, dv2])
        #e6.display
        #e6.plt.plot[1]
        # test 6
        m6.setParameterValues([1.5, 0.3, 0.17, 0.15]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m6.fit([e6],lb,ub,0)
        if not isSimilar([1.0000, 0.2000, 0.2000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 6')
        # test 7
        timeUnit = 'min'
        doseUnit = 'umol'
        valueUnit = 'uM'
        k_init = [1.5, 0.3, 0.17, 0.15]
        times = [i for i in range(0,41,5)] #[0   5  10  15  20  25  30  35  40]
        c = [0,45.1427,16.9112,6.22335,2.28945,0.842243,0.309844,22.6853,8.49755]
        ce = [0,21.5009,23.9592,18.5689,12.7478,8.27831,5.22205,13.9918,13.9727]
        d7 = pko.Dosing('ev', [0, 30], [20, 10], [], timeUnit, doseUnit)
        model = pko.Model(_pk.PK_1comp_Ka_Cl_fun(0,1,k_init), None, None)
        dv1 = pko.DataVariable('C',times,c,c*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,ce,ce*0.1, timeUnit, valueUnit)
        ex2 = pko.Experiment('Ex2',d7, [dv1, dv2])
        lb = k*0.1
        ub = k*10
        x = model.fit([ex2],lb,ub,0)
        if not isSimilar([1.0000, 0.2000, 0.2000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 7')
        # test 8
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp8',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        sf2 = _pd.PD_Emax_hill_fun([3, 1, 0.5])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp8',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 3.5, 2, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 3, 1, 0.5], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 8')
        # test 9
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp9',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        sf2 = _pd.PD_linear_fun([0.5])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp9',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 0.5], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 9')
        # test 10
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp8',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        sf2 = _pd.PD_Imax_fun([1, 0.5])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp10',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 1.1, 0.2]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 1, 0.5], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 10')
        # test 11
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp11',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        sf2 = _pd.PD_Imax_hill_fun([1, 0.5, 2])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp11',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 1.1, 0.2, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 1, 0.5, 2], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 11')
        # test 12
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp12',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        sf2 = _pd.PD_loglinear_fun([5, 1])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp12',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 2, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 5, 0.5], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 12')
        # test 13
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp13',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        # R0, kout, Imax, IC50, n
        sf2 = _pd.PD_turnover_inhOfLoss_ODE([1, 0.1, 1, 0.5, 1])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp12',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.1, 0.2, 0.9, 0.6, 0.2, 0.9, 0.4, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 1, 0.1, 1, 0.5, 1], x, 0.01):
            raise TestFailed('Error, parameter estimation ev, test 13')
        # test 14
        del Y
        d = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e = pko.Experiment('testExp14',d, [])
        sf = _pk.PK_1comp_Ka_Cl_fun(0,0, [1, 0.5, 0.5])
        # R0, kout, Imax, IC50, n
        sf2 = _pd.PD_turnover_stimOfProd_ODE([1, 0.1, 1, 0.5, 1])
        m = pko.Model(sf, None, sf2)
        Y = m.simulate(times, e)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        ee = pko.Experiment('testExp14',d, [dv1, dv2])
        #     ee.display
        #     ee.plt.plot(2)
        kinit = [1.2, 0.25, 0.7, 0.6, 0.2, 0.8, 0.4, 1]
        m.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m.fit([ee],lb,ub,3)
        if not isSimilar([1, 0.5, 0.5, 1, 0.1, 1, 0.5, 1], x, 0.001):
            raise TestFailed('Error, parameter estimation ev, test 14')


    # TEST PARAMETER ESTIMATION INF
    if MYTYPE==0 or MYTYPE==8:
        print('== Test  8')
        times = [i for i in range(73)]
        k = [1, 0.2, 0.2, 1.2, 0.3, 0.1]
        # SET UP EXP 1-2
        d1 = pko.Dosing('inf', [0, 30], [20, 10], [5, 5], timeUnit, doseUnit)
        e1 = pko.Experiment('testExp1',d1, [])
        sf1 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m1 = pko.Model(sf1, None, None)
        Y = m1.simulate(times, e1)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e1 = pko.Experiment('testExp1',d1, [dv1])
        e2 = pko.Experiment('testExp2',d1, [dv1, dv2])
        #e1.display
        #e1.plt.plot[1]
        # SET UP EXP 3
        del Y
        d3 = pko.Dosing('inf', [0], [50], [20], timeUnit, doseUnit)
        e3 = pko.Experiment('testExp3',d3, [])
        sf3 = _pk.PK_2comp_Ka_Cl_fun(0,1,k)
        m3 = pko.Model(sf3, None, None)
        Y = m3.simulate(times, e3)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e3 = pko.Experiment('testExp3',d3, [dv1, dv2])
        #e3.display
        #e3.plt.plot(2)
        # Test 1
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        print('    Fit should warn below: Warning: Matrix is singular...')
        x = m1.fit([e1],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, np.nan], x, 0.001):
            raise TestFailed('Error, parameter estimation inf, test 1')
        # Test 2
        m1.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m1.fit([e1, e2],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation inf, test 2')
        # test 3
        m3.setParameterValues([1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m1.fit([e3],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation inf, test 3')
        # test 4
        del Y
        d4 = pko.Dosing('inf', [0, 30], [5, 3], [10, 10], timeUnit, doseUnit)
        e4 = pko.Experiment('testExp4',d4, [])
        sf4_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5])
        sf4_2 = _pd.PD_Emax_fun([3, 1])
        m4 = pko.Model(sf4_1, None, sf4_2)
        Y = m4.simulate(times, e4)
        dv1 = pko.DataVariable('C', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('R', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e4 = pko.Experiment('testExp4',d4, [dv1, dv2])
        #e4.display
        #e4.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5, 6.5]
        m4.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        x = m4.fit([e4],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 1.2000, 0.3000, 3.000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation inf, test 4')
        # test 5
        del Y
        d5 = pko.Dosing('ev', [0, 30], [5, 3], [], timeUnit, doseUnit)
        e5 = pko.Experiment('testExp4',d4, [])
        sf5_1 = _pk.PK_2comp_Ka_Cl_fun(0,0,k[0:5])
        sf5_2 = _pd.PD_Emax_fun([3, 1])
        m5 = pko.Model(sf5_1, None, sf5_2)
        Y = m6.simulate(times, e5)
        dv1 = pko.DataVariable('Cf', times, Y[:,1].T, Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Rf', times, Y[:,2].T, Y[:,2].T*0.1, timeUnit, valueUnit)
        e5 = pko.Experiment('testExp4',d5, [dv1, dv2])
        #e5.display
        #e5.plt.plot(2)
        kinit = [1.1288, 0.2241, 0.1928, 1.2128, 0.2753, 3.5, 6.5]
        m5.setParameterValues(kinit) # set false parameters
        lb = kinit*0.1
        ub = kinit*10
        try:
            m5.fit([e5],lb,ub,0)
            raise TestFailed('Error, parameter estimation inf, test 5')
        except pko.BadInputException:
            pass # ok
        # SET UP EXP 6
        k = [1, 0.2, 0.2, 0.1]
        d6 = pko.Dosing('inf', [0, 30], [20, 10], [5, 5], timeUnit, doseUnit)
        e6 = pko.Experiment('testExp6',d6, [])
        sf6 = _pk.PK_1comp_Ka_Cl_fun(0,1,k)
        m6 = pko.Model(sf6, None, None)
        Y = m6.simulate(times, e6)
        valueUnit = 'uM'
        dv1 = pko.DataVariable('C',times,Y[:,1].T,Y[:,1].T*0.1, timeUnit, valueUnit)
        dv2 = pko.DataVariable('Ce',times,Y[:,2].T,Y[:,2].T*0.1, timeUnit, valueUnit)
        e6 = pko.Experiment('testExp6',d6, [dv1, dv2])
        #e6.display
        #e6.plt.plot[1]
        # test 6
        m6.setParameterValues([1.1288, 0.2241, 0.1928, 0.1536]) # set false parameters
        lb = k*0.1
        ub = k*10
        x = m6.fit([e6],lb,ub,0)
        if not isSimilar([np.nan, 0.2000, 0.2000, 0.100], x, 0.001):
            raise TestFailed('Error, parameter estimation inf, test 6')


    # TEST PARAMETER ESTIMATION IV+EV
    if MYTYPE==0 or MYTYPE==9:
        print('== Test  9')
        valueUnit = 'uM'
        times = [i*0.5 for i in range(11)] #[0:0.5:5]
        # ----------------------------------------------------
        # TEST 9.1, 1-COMP MODEL
        # ----------------------------------------------------
        k = [1, 0.2, 0.2, 0.5]
        sf1 = _pk.PK_1comp_Ka_Cl_fun(2,0,k)
        m1 = pko.Model(sf1, None, None)
        sf2 = _pk.PK_1comp_Ka_Cl_ODE(2,0,k)
        m2 = pko.Model(sf2, None, None)
        # SET UP IV EXPERIMENT
        d1 = pko.Dosing('iv', [0], [10], [], timeUnit, doseUnit)
        #e1 = pko.Experiment('Exp1',d1, [])
        #Y = m1.simulate(times, e1)
        #Y=Y+0.01*Y*randn(len(times),1)
        Y = [   50.0000, 45.2419, 40.9365, 37.0409, 33.5160, 30.3265,
            27.4406, 24.8293, 22.4664, 20.3285, 18.3940]
        dv1 = pko.DataVariable('C',times,Y,Y*0.1, timeUnit, valueUnit)
        e1 = pko.Experiment('Exp1',d1, [dv1],'r')
        #e1.display
        #e1.plt.plot[1]
        # SET UP EV EXPERIMENT
        d2 = pko.Dosing('ev', [0], [30], [], timeUnit, doseUnit)
        #e2 = pko.Experiment('Exp2',d2, [])
        #Y = m1.simulate(times, e2)
        #Y=Y+0.01*Y*randn(len(times),1)
        Y = [0, 27.9663, 42.2673, 48.5333, 50.1548, 49.1668, 46.7836,
            43.7239, 40.4075, 37.0744, 33.8570]
        dv2 = pko.DataVariable('C',times,Y,Y*0.1, timeUnit, valueUnit)
        e2 = pko.Experiment('Exp2',d2, [dv2],'b')
        #e2.display
        #e2.plt.plot(2)
        # test 1a
        m1.setParameterValues([1.5, 0.4, 0.1,  0.7]) # set false parameters
        lb = k*0.1
        ub = k*10
        [x,residual,ci,se,cv] = m1.fit([e1, e2],lb,ub,3)
        if not isSimilar([1, 0.2, 0.2, 0.5], x, 0.01):
            raise TestFailed('Error, parameter estimation iv+ev, test 9.1a')
        # test 1b
        m2.setParameterValues([1.5, 0.4, 0.1, 0.7]) # set false parameters
        lb = k*0.1
        ub = k*10
        [x,residual,ci,se,cv] = m2.fit([e1, e2],lb,ub,3)
        if not isSimilar([1, 0.2, 0.2, 0.5], x, 0.01):
            raise TestFailed('Error, parameter estimation iv+ev, test 9.1b')
        # ----------------------------------------------------
        # TEST 2, 1-COMP MODEL, NON-LIN ELIMINATION
        # ----------------------------------------------------
        k = [2, 0.2, 0.3, 1, 0.5]
        sf = _pk.PK_1comp_Ka_satCl_ODE(1,0,k)
        m = pko.Model(None, sf1, None)
        # SET UP IV EXPERIMENT
        d1 = pko.Dosing('iv', [0], [1], [], timeUnit, doseUnit)
        e1 = pko.Experiment('Exp1',d1, [])
        Y = m.simulate(times, e1)
        np.random.seed(0)
        Y = Y + 0.01*Y*np.random.randn(len(times),2)
        dv1 = pko.DataVariable('C',times,Y[:,2].T, [], timeUnit, valueUnit)
        e1 = pko.Experiment('Exp1',d1, [dv1],'r')
        #e1.display
        #e1.plt.plot(10)
        # SET UP EV EXPERIMENT
        d2 = pko.Dosing('ev', [0], [3], [], timeUnit, doseUnit)
        e2 = pko.Experiment('Exp2',d2, [])
        Y = m.simulate(times, e2)
        Y = Y + 0.01*Y*np.random.randn(len(times),2)
        dv2 = pko.DataVariable('C',times,Y[:,2].T, [], timeUnit, valueUnit)
        e2 = pko.Experiment('Exp2',d2, [dv2],'b')
        #e2.display
        #e2.plt.plot(20)
        # test 2
        m.setParameterValues([10, 0.3, 2, 15, 0.7]) # set false parameters
        lb = k*0.1
        ub = k*10
        [x, residual, ci, se, cv] = m.fit([e1, e2],lb,ub,3)
        if not isSimilar(k,x,0.1):
            raise TestFailed('Error, parameter estimation iv+ev, test 9.2')


    # TEST NICEPLOT/NICEFIGURE/NICESEMILOGY
    if MYTYPE==0 or MYTYPE==10:
        print('== Test  10  :     DEACTIVATED')
        print('    Three plt.plots (figure 20-22 ) should be drawn')
#         settings.fontsize=16
#         settings.lineWidth=2
#         settings.markerSize=6
#         plt.figure(20)
#         monrange=[i for i in range (1,11,1)]
#         niceplt.plot(monrange,randn(10,1),'r','s',':',settings)
#         niceplt.plot(monrange,randn(10,1),'b','o','',settings)
#         niceplt.plot(monrange,randn(10,1),'g','','-',settings)
#         niceFigure(20,'This is an \Theta and \eta plt.plot','Time (h)','\muMol',settings)
#         plt.figure(21)
#         niceplt.plot(monrange,randn(10,1),'m','s',':',settings)
#         niceplt.plot(monrange,randn(10,1), [0.3, 0.3, 0.7],'>','-.',settings)
#         niceFigure(21,'','Time (h)','\muMol',settings)
#         plt.figure(22)
#         niceSemilogy(monrange,500+100*randn(10,1),'m','s',':',settings)
#         niceSemilogy(monrange,50+10*randn(10,1), [0.3, 0.3, 0.7],'>','-.',settings)
#         niceFigure(22,'','Time (h)','\muMol',settings)

    # TEST _pk.PK_1/2comp_Ka_Cl_fun/ode
    if MYTYPE==0 or MYTYPE==11:
        print('== Test  11')

        # _pk.PK_1comp_Ka_Cl_fun
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)]
        effectComp=0
        const = 0
        sf = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84])
        Y = sf.simulate(times,experiment, [], [])

        # _pk.PK_1comp_Ka_Cl_fun
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[i for i in range (73)] #[0:1:72]
        effectComp=1
        const = 0
        sf = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 1],['Cp','Cpe'])
        Y = sf.simulate(times,experiment, [], [])
        try:
            sf = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 1],['Cp'])
            raise TestFailed('Error, _pk.PK_1comp_Ka_Cl_fun')
        except pko.BadInputException:
            pass # ok
        try:
            effectComp=0
            sf = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6],['Cp'])
            raise TestFailed('Error, _pk.PK_1comp_Ka_Cl_fun')
        except pko.BadInputException:
            pass # ok

        # _pk.PK_1comp_Ka_Cl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[i for i in range (73)] #[0:1:72]
        effectComp=0
        const = 0
        sf = _pk.PK_1comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84])
        Y = sf.simulate(times,experiment, [], [])

        # _pk.PK_1comp_Ka_Cl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[i for i in range (73)] #[0:1:72]
        effectComp=1
        const = 0
        sf = _pk.PK_1comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84, 1],['G','Cp','Cpe'])
        Y = sf.simulate(times,experiment, [], [])
        try:
            sf = _pk.PK_1comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84, 1],['G','C'])
            raise TestFailed('Error, _pk.PK_1comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok
        try:
            effectComp=0
            sf = _pk.PK_1comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6],['G','Cp'])
            raise TestFailed('Error, _pk.PK_1comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok

        # _pk.PK_2comp_Ka_Cl_fun
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=0
        const = 0
        sf = _pk.PK_2comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 1.5, 0.1])
        Y = sf.simulate(times,experiment, [], [])

        # _pk.PK_2comp_Ka_Cl_fun
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=1
        const = 0
        sf = _pk.PK_2comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 1, 1, 1],['Cp','Cpe'])
        Y = sf.simulate(times,experiment, [], [])
        try:
            sf = _pk.PK_2comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 1, 1, 1],['Cp'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_fun')
        except pko.BadInputException:
            pass # ok
        try:
            effectComp=0
            sf = _pk.PK_2comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 1, 1],['Cp'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_fun')
        except pko.BadInputException:
            pass # ok

        # _pk.PK_2comp_Ka_Cl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=0
        const = 0
        sf = _pk.PK_2comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84, 1.5, 0.1])
        Y = sf.simulate(times,experiment, [], [])
        # _pk.PK_2comp_Ka_Cl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=1
        const = 0
        sf = _pk.PK_2comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84, 1.5, 0.1, 1],['G','Cp','Cp2','Cpe'])
        Y = sf.simulate(times,experiment, [], [])
        try:
            sf = _pk.PK_2comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84, 0.5, 0.5, 1],['Cp','a','b'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok
        try:
            effectComp=0
            sf = _pk.PK_2comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 1, 2],['G','C','C2'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok

        # _pk.PK_2comp_Ka_satCl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=0
        const = 0
        sf = _pk.PK_1comp_Ka_satCl_ODE(const,effectComp, [0.86, 0.84, 1.6, 0.5])
        Y = sf.simulate(times,experiment, [], [])
        # _pk.PK_2comp_Ka_satCl_ode
        timeUnit = 'min'
        doseUnit = 'umol'
        d = pko.Dosing('ev', [0, 24, 32, 48], [20, 10, 30, 20], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (73)] #[0:1:72]
        effectComp=1
        const = 0
        sf = _pk.PK_1comp_Ka_satCl_ODE(const,effectComp, [0.86, 0.84, 1.6, 0.5, 1],['G','Cp','Cpe'])
        Y = sf.simulate(times,experiment, [], [])
        try:
            sf = _pk.PK_2comp_Ka_satCl_ODE(const,effectComp, [0.86, 0.84, 1.6, 0.5, 1],['Cp','a'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok
        try:
            effectComp=0
            sf = _pk.PK_2comp_Ka_satCl_ODE(const,effectComp, [0.86, 0.84, 1.6, 0.5],['G','C'])
            raise TestFailed('Error, _pk.PK_2comp_Ka_Cl_ode')
        except pko.BadInputException:
            pass # ok
    
    
    # TEST PD models
    if MYTYPE==0 or MYTYPE==12:
        print('== Test  12')
        # _pd.PD_receptocc_ode
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=0
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84])
        Y1 = sf1.simulate(times,experiment)
        sf2 = _pd.PD_receptocc_ODE([0.2, 0.1, 1])
        Y2 = sf2.simulate(times,experiment, [],sf1)
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=1
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue'])
        sf2 = _pd.PD_receptocc_ODE([0.2, 0.1, 1],'Cue','E')
        model = pko.Model(sf1, sf2, None)
        Y = model.simulate(times,experiment)
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue'])
            sf2 = _pd.PD_receptocc_ODE([0.2, 0.1, 1],'Cue2','E')
            model = pko.Model(sf1, sf2, None)
            raise TestFailed('Error, _pd.PD_receptocc_ode 1')
        except pko.BadInputException:
            pass # ok
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue2'])
            sf2 = _pd.PD_receptocc_ODE([0.2, 0.1, 1],'Cue','E')
            model = pko.Model(sf1, sf2, None)
            raise TestFailed('Error, _pd.PD_receptocc_ode 2')
        except pko.BadInputException:
            pass # ok
        try:
            sf2 = _pd.PD_receptocc_ODE([0.2, 0.1],'Cue','E')
            raise TestFailed('Error, _pd.PD_receptocc_ode 3')
        except pko.BadInputException:
            pass # ok
        # _pd.PD_Emax_fun
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=0
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_ODE(const,effectComp, [0.86, 1.6, 0.84],['G','Cu'])
        sf2 = _pd.PD_Emax_fun([1, 0.5],'Cu','E')
        model = pko.Model(None, sf1, sf2)
        Y = model.simulate(times,experiment)
        # _pd.PD_Emax_fun
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=1
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1])
        sf2 = _pd.PD_Emax_fun([1, 0.5],'Ce')
        model = pko.Model(sf1, None, sf2)
        Y = model.simulate(times,experiment)
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue'])
            sf2 = _pd.PD_Emax_fun([0.2, 0.1],'Cue2','E')
            model = pko.Model(sf1, None, sf2)
            raise TestFailed('Error, _pd.PD_Emax_fun 1')
        except pko.BadInputException:
            pass # ok
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue2'])
            sf2 = _pd.PD_Emax_fun([0.2, 0.1],'Cue','E')
            model = pko.Model(sf1, None, sf2)
            raise TestFailed('Error, _pd.PD_Emax_fun 2')
        except pko.BadInputException:
            pass # ok
        try:
            sf2 = _pd.PD_Emax_fun([0.2],'Cue','E')
            raise TestFailed('Error, _pd.PD_Emax_fun 3')
        except pko.BadInputException:
            pass # ok
        # _pd.PD_powerOf_fun
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=0
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84])
        sf2 = _pd.PD_powerOf_fun([1, 0.5])
        model = pko.Model(sf1, None, sf2)
        Y = model.simulate(times,experiment)
        timeUnit = 'min'
        doseUnit = 'umol/kg'
        d = pko.Dosing('ev', [0, 24], [20, 10], [], timeUnit, doseUnit)
        experiment = pko.Experiment('Ex',d, [])
        times = [i for i in range (49)] #[0:1:48]
        effectComp=1
        const = 0
        sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1])
        sf2 = _pd.PD_powerOf_fun([1, 0.5],'Ce')
        model = pko.Model(sf1, None, sf2)
        Y = model.simulate(times,experiment)
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue'])
            sf2 = _pd.PD_powerOf_fun([0.2, 0.1],'Cue2','E')
            model = pko.Model(sf1, None, sf2)
            raise TestFailed('Error, _pd.PD_powerOf_fun 1')
        except pko.BadInputException:
            pass # ok
        try:
            sf1 = _pk.PK_1comp_Ka_Cl_fun(const,effectComp, [0.86, 1.6, 0.84, 0.1],['Cu','Cue2'])
            sf2 = _pd.PD_powerOf_fun([0.2, 0.1],'Cue','E')
            model = pko.Model(sf1, None, sf2)
            raise TestFailed('Error, _pd.PD_powerOf_fun 2')
        except pko.BadInputException:
            pass # ok
        try:
            sf2 = _pd.PD_powerOf_fun([0.2],'Cue','E')
            raise TestFailed('Error, _pd.PD_powerOf_fun 3')
        except pko.BadInputException:
            pass # ok


    # TEST CLOSED FORM STEADY STATE PK EQUATIONS
    if MYTYPE==0 or MYTYPE==13:
        print('== Test  13')
        print('    Figure 101-103 tentatively successful if solid graphs are overwritten by\n' \
              '    dotted graphs in the plt.plots, and no error messages are reported.')

        # ONE COMP
        kmicro = [0.5, 5, 0.05]
        plt.figure(101)
        adm=['iv','inf','ev']
        for j in [1,2,3]: #[1:3]        # for iv inf ev
            d = pko.Dosing(str(adm[j]), [i for i in range (0,121,12)],np.ones(1,11)*10,np.ones(1,11)*2,'','') #[0:12:12*10]
            experiment = pko.Experiment('testExp',d, [])
            for i in [1,2]: #without and withn effect compartment
                plt.plot(3,2,(j-1)*2+i)
                times = np.arange(0,132.1,0.1).tolist() #[0:0.1:132]
                #ou [i*0.1 for i in range (0,132*10+1,1)]
                constants=0
                calculateEffectCompartment = i-1
                if i==1:
                    sf0 = _pk.PK_1comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                    sf1 = _pk.PK_1comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment,kmicro)
                else:
                    sf0 = _pk.PK_1comp_Ka_Cl_fun(constants,calculateEffectCompartment, [kmicro, 0.2])
                    sf1 = _pk.PK_1comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment, [kmicro, 0.2])
                Y = sf0.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                Z = sf1.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Z,'g:')
                plt.setp([h1, h2],'LineWidth',2)
                if i==1:
                    plt.title('PK\_1comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', no effect-comp')
                else:
                    plt.title('PK\_1comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', effect-comp')
                plt.xlabel('Time')
                plt.ylabel('C')
                plt.ylim([1, 10])
        # TWO COMP
        kmicro = [0.5, 1.6, 0.3, 0.75, 0.15]#['ka','V','k','k12','k21','ke']
        plt.figure(102)
        adm=['iv','inf','ev']
        for j in range[1,2,3]: # for iv inf ev
            d = pko.Dosing(str(adm[j]), [i for i in range (0,121,12)] ,np.ones(1,11)*10,np.ones(1,11)*2,'','') #[0:12:12*10]
            experiment = pko.Experiment('testExp',d, [])
            for i in [1,2]: #without and withn effect compartment
                plt.plot(3,2,(j-1)*2+i)
                times = np.arange(0,132.1,0.1).tolist() #[0:0.1:132]
                constants=0
                calculateEffectCompartment = i-1
                if i==1:
                    sf0 = _pk.PK_2comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                    sf1 = _pk.PK_2comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment,kmicro)
                else:
                    sf0 = _pk.PK_2comp_Ka_Cl_fun(constants,calculateEffectCompartment, [kmicro, 0.2])
                    sf1 = _pk.PK_2comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment, [kmicro, 0.2])
                Y = sf0.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                Z = sf1.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Z,'g:')
                plt.setp([h1, h2],'LineWidth',2)
                if i==1:
                    plt.title('PK\_2comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', no effect-comp')
                else:
                    plt.title('PK\_2comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', effect-comp')
                    plt.xlabel('Time')
                plt.ylabel('C')
                plt.ylim([0.1, 10])
        # 3 COMP
        kmicro = [0.5, 1.6, 0.3, 0.75, 0.15, 0.6, 0.2]#['ka','V','k','k12','k21',k13,k31]
        plt.figure(103)
        adm=['iv','inf','ev']
        for j in [1,2,3]: # for iv inf ev
            d = pko.Dosing(str(adm[j]), [i for i in range (0,121,12)],np.ones(1,11)*10,np.ones(1,11)*2,'','')
            experiment = pko.Experiment('testExp',d, [])
            for i in [1]: #=1:1 #without and withn effect compartment
                plt.plot(2,2,j)
                times = np.arange(0,132,0.1).tolist()  #[i*0.1 for i in range(132*10+1)]
                constants=0
                calculateEffectCompartment = i-1
                if i==1:
                    sf0 = _pk.PK_3comp_Ka_Cl_fun(constants,calculateEffectCompartment,kmicro)
                    sf1 = _pk.PK_3comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment,kmicro)
                else:
                    sf0 = _pk.PK_3comp_Ka_Cl_fun(constants,calculateEffectCompartment, [kmicro, 0.2])
                    sf1 = _pk.PK_3comp_Ka_Cl_fun_SS(constants,calculateEffectCompartment, [kmicro, 0.2])
                Y = sf0.simulate(times,experiment, [], [])
                h1 = plt.semilogy(times,Y,'k')
                plt.hold(True)
                Z = sf1.simulate(times,experiment, [], [])
                h2 = plt.semilogy(times,Z,'g:')
                plt.setp([h1, h2],'LineWidth',2)
                if i==1:
                    plt.title('PK\_3comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', no effect-comp')
                else:
                    plt.title('PK\_3comp\_Ka\_Cl\_fun\_SS, ', str(adm[j]), ', effect-comp')
                plt.xlabel('Time')
                plt.ylabel('C')
                plt.ylim([0.1, 10])

def isSimilar(trueX,X,TOL):
# r=0, is not similar
# r=1, is similar
    if len(trueX)!=len(X):
        print('Error in function isSimilar. Vector lengths must be equal')
    r = False
    for i in range(len(X)):
        if not np.isnan(trueX[i]):
            if (abs(trueX[i])-abs(X[i]))>TOL:
                return
    r = True
    return r
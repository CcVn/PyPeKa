# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 16:46:22 2018

@author: cvn_or
"""
#%%
import matplotlib.pyplot as plt

class matplotlibSettings(object):
    def __init__(self, setting_dict):
        for key, val in setting_dict.items():    
            pass
            #eval("self."+str(key)) = eval("self."+str(val))
            
#%%
def nicePlot(x, y, col, markerType, lineStyle, label, settings):
    '''
    # nicePlot
    # Adjusts color, markers and lines to a plot.
    # Inputs:
    #   x, vector of x data
    #   y, vector of y data
    #   col, character [r=red, b=blue, see 'help plot' for details), or 
    #            vector of three numbers (Red/Blue/Green in interval [0,1],
    #            e.g., [1 0 0]=red, [0 0 1]=blue). Default blue.
    #   markerType, string, e.g. 'o', 's', see 'help plot'
    #   lineStyle, string, e.g. '-', ':', see 'help plot'
    #   settings, structure including at least the field markerSize and 
    #             lineWidth, e.g., settings.markerSize=7 settings.lineWidth=2 
    #
    # see also niceSemilogy
    '''
    h = plt.plot(x,y)
    #hold(True)
    plt.setp(h,'Color',col)
    
    if lineStyle:
        plt.setp(h, 'LineStyle',lineStyle, 'LineWidth',settings["lineWidth"])
    else:
        plt.setp(h, 'LineStyle','none')
    
    
    if markerType:
        plt.setp(h, 'Marker', markerType, 'MarkerEdgeColor', col, \
            'MarkerFaceColor', col, 'MarkerSize', settings["markerSize"])
    else:
        plt.setp(h, 'Marker','none')
    
    return h

#%%
def niceFigure(ifig, titleText, xLabel, yLabel, settings):
    '''
    # niceFigure
    # Adjusts fonts, colors and other settings to a figure.
    # Call after plotting has occured.
    # Inputs:
    #   ifig, integer for figure number
    #   titleText, string
    #   xLabel,string
    #   yLabel,string 
    #   settings, structure including at least the field fontsize (integer)
    #             e.g., settings.fontsize=14  
    '''
    plt.figure(ifig)
    #hold on
 #   plt.setp(plt.gcf,'Position',[50, 50, 500, 375],'Units','pixels')
    plt.setp(plt.gca,'FontName','Helvetica')
    plt.setp(plt.gca,'FontSize',settings["fontsize"])
    hTitle  = plt.title(titleText)
    hXLabel = plt.xlabel(xLabel)
    hYLabel = plt.ylabel(yLabel)
    plt.setp([hTitle, hXLabel, hYLabel], 
        'FontName','Helvetica', #'Interpreter','tex',
        'FontSize',settings["fontsize"])
    plt.setp(plt.gca,
        'Box'         , 'off'     ,
        'TickDir'     , 'out'     ,
        'YGrid'       , 'off'      ,
        'YGrid'       , 'off'      ,
        'XMinorTick'  , 'on'      ,
        'YMinorTick'  , 'on'      ,
        'XColor'      , [.1, .1, .1],
        'YColor'      , [.1, .1, .1],
        'LineWidth'   , 1)
    #plt.legend(boxoff) #fancybox=True?
    #get(plt.gca)

#%%
def niceSemilogx_y(x, y, col, markerType=None, lineStyle=None, settings=None, axis=0):
    ''' niceSemilogx_y : merge of niceSemilogx & niceSemilogy
     Adjusts color, markers and lines to a plot with logged x-axis or y-axis
     Inputs:
      - x, vector of x data
      - y, vector of y data
      - col, character [r=red, b=blue, see 'help plot' for details), or 
                vector of three numbers (Red/Blue/Green in interval [0,1],
                e.g., [1 0 0]=red, [0 0 1]=blue). Default blue.
      - markerType, string, e.g. 'o', 's', see 'help plot'
      - lineStyle, string, e.g. '-', ':', see 'help plot'
      - settings, structure including at least the field markerSize and 
                 lineWidth, e.g., settings.markerSize=7; settings.lineWidth=2; 
    
     see also nicePlot
    '''
    if axis==0:
        h = plt.semilogx(x, y)
    else:
        h = plt.semilogy(x, y)
    plt.hold(True)
    plt.setp(h, 'Color', col)
    
    if lineStyle is not None:
        plt.setp(h, 'LineStyle',settings["lineStyle"], \
                    'LineWidth',settings["lineWidth"])
    else:
        plt.setp(h, 'LineStyle','none')
    
    if markerType is not None:
        plt.setp(h, 'Marker', markerType,\
                   'MarkerEdgeColor', col,\
                   'MarkerFaceColor', col,\
                   'MarkerSize', settings["markerSize"])
    else:
        plt.setp(h, 'Marker','none')
    
    return h

def niceSemilogx(x, y, col, markerType=None, lineStyle=None, settings=None, axis=0):
    ''' niceSemilogx
     see niceSemilogx_y doc
    '''
    h = niceSemilogx_y(x, y, col, markerType, lineStyle, settings, axis=1)
    return h

def niceSemilogy(x, y, col, markerType=None, lineStyle=None, settings=None, axis=1):
    ''' niceSemilogy
     see niceSemilogx_y doc
    '''
    h = niceSemilogx_y(x, y, col, markerType, lineStyle, settings, axis=1)
    return h

#%%
#=======================================================================
if __name__ == "__main__":
    import numpy.random as nr
    # Example:
    plt.close('all')
    mysettings =  { 'fontsize': 14,
                  'fontsizeLegend': 10, 
                  'markerSize': 7,
#LineStyle: linestyle: [ '-' | '--' | '-.' | ':' | 'steps' | 'None' ] 
                  'lineStyle': '--',
                  'lineWidth': 2 }

#FIXME utiliser matplotlib.rcParams

    plt.figure(1)
    x=[i for i in range(10)]
    hA = nicePlot(x,     nr.randn(10), 'b', 'o', '-', 'A', mysettings)
    hB = nicePlot(x, 2 + nr.randn(10), 'r', 's', '--', 'B', mysettings)
    hC = nicePlot(x, 4 + nr.randn(10), 'm', 'x', ':', 'C', mysettings)
    hl = plt.legend(['A','B','C'],loc='upper left')

    #plt.setp([hl], 'FontName','Arial','Interpreter','tex', 'FontSize', settings['fontsizeLegend'])
    plt.xlim([0, 11])
    plt.ylim([-2, 7])
    #niceFigure(1,'The title','Time (h)','Conc. (uM)', settings)


    #Semilogy
    plt.figure(2)
    scale = [i for i in range(1,11,1)]
    niceSemilogy(scale, 1000+(500 * nr.randn(10,1)),'b','o','-',mysettings)
    plt.hold(True)
    niceSemilogy(scale,100+(50*nr.randn(10,1)),'r','s','--',mysettings)
    niceSemilogy(scale,10+(5*nr.randn(10,1)),'m','x',':',mysettings)
    hl = plt.legend(['A','B','C'],loc='upper right')   #, 'Location','NorthEast' in Matlab
    #plt.setp([hl], \
#         'FontName','AvantGarde',\
#         'FontSize',mysettings["fontsize"], \
#         'Interpreter','tex')
    plt.xlim([0, 11])
    plt.ylim([1, 10000])
    niceFigure(1,'The title','Time (h)','Conc. (uM)',mysettings)

    #NB: rc('lines', linewidth=2, color='r')  # a la place de settings        pass


#matplotlib.rc_context(rc=None, fname=None)[source]
#    Return a context manager for managing rc settings.
#matplotlib.rcParams
#with matplotlib.rc_context(fname='screen.rc'):
#    plt.plot(x, a)
#    with mpl.rc_context(fname='print.rc'):
#        plt.plot(x, b)
#    plt.plot(x, c)
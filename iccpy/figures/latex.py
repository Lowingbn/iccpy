import pylab as pl

# 240 pt figures (3.32 inches)
latexParams = {'figure.dpi':150,
               'figure.figsize':[3.32,2.49],
               'text.usetex':True,
               'text.fontsize':8,
               'font.size':8,
               'axes.labelsize':8,
               'figure.subplot.top':0.95,
               'figure.subplot.right':0.95,
               'figure.subplot.bottom':0.15,
               'figure.subplot.left':0.15,
               'lines.linewidth':0.5,
               'lines.markersize':3,
               'savefig.dpi':350}

if __name__=='__main__':
    for x in sorted(pl.rcParams):
        print x,pl.rcParams[x]
    pl.rcParams.update(latexParams)
    import numpy as np
    x = np.arange(100)/20.0
    y = np.sin(x)
    pl.plot(x,y, 'x')
    pl.xlabel(r'$\Delta x$')
    pl.ylabel(r'$T/T_0$')
    pl.show()

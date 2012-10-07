import pylab as pl


def panelplot(numRows, numCols):
    """ draw a panel of (numRows, numCols) plots with axes in the right places """
    subplots = []

    # no space between the panels
    pl.rcParams.update({'figure.subplot.wspace':0,'figure.subplot.hspace':0})
    
    for i in range(numRows*numCols):
        ax= pl.subplot(numRows,numCols,i+1)
        
        # only put left ticks on the left plots, right ticks on the right, lower on the lower
        left_tick = i%numCols==0
        right_tick = (i+1)%numCols==0
        
        
        for tick in ax.yaxis.get_major_ticks():
            tick.label1On=left_tick
            tick.label2On=right_tick

        lower_tick = i>=(numRows-1)*numCols
        for tick in ax.xaxis.get_major_ticks():
            tick.label1On=lower_tick
            
        
    
        subplots.append(ax)
        
    return subplots

if __name__=='__main__':
    subplots = panelplot(3,3)
    from numpy import arange
    x = arange(20)

    for ax in subplots:
        # set this to be the current axes
        pl.axes(ax)
        pl.plot(x,x)
        
        pl.ylim(0,19.9)
        pl.xlim(0,19.9)

pl.show()
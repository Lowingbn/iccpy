import iccpy.gadget
import matplotlib.pyplot as pl
import numpy as np
import iccpy.utils

sim_label = { 'aqa' : 'A', 'aqb' : 'B', 'aqc':'C', 'aqd':'D', 'aqe':'E' }
last_snapnum = { 'aqa2' : 1023, 'aqa3' : 511, 'aqa4' : 1023, 'aqb2' : 127, 'aqc2' : 127, 'aqd2' : 127, 'aqe2' : 127 }

r_200 = { 'aqa1' : 245.67, 'aqa2' : 245.88, 'aqa3' : 245.64, 'aqa4' : 245.70, 'aqa5' : 246.37, 'aqb2' : 187.70, 'aqb4' : 188.85,
          'aqc2' : 242.82, 'aqc4' : 243.68, 'aqd2' : 242.85, 'aqd4' : 243.60, 'aqe2' : 212.28, 'aqe4' : 213.63, 'aqf2' : 209.21,
          'aqf4' : 207.15 }
M_200 = { 'aqa1' : 183.9,  'aqa2' : 184.2,  'aqa3' : 183.6,  'aqa4' : 183.8,  'aqa5' : 185.3,  'aqb2' : 81.94,  'aqb4' : 83.45,
          'aqc2' : 177.4,  'aqc4' : 179.3,  'aqd2' : 177.4,  'aqd4' : 179.1,  'aqe2' : 118.5,  'aqe4' : 120.8,  'aqf2' : 113.5,
          'aqf4' : 110.1 }

merger_tree_filename = { 'aqa2' : '/gpfs/data/jch/Aquarius/Trees/Aq-A/2/trees/treedir_127/tree_127.0.hdf5',
                         'aqb2' : '/gpfs/data/d50wse/WMAP7_Trees/trees_Aq-B2/treedir_127/tree_127.0.hdf5',
                         'aqc2' : '/gpfs/data/d50wse/WMAP7_Trees/trees_Aq-C2/treedir_127/tree_127.0.hdf5',
                         'aqd2' : '/gpfs/data/d50wse/WMAP7_Trees/trees_Aq-D2/treedir_127/tree_127.0.hdf5',
                         'aqe2' : '/gpfs/data/d50wse/WMAP7_Trees/trees_Aq-E2/treedir_127/tree_127.0.hdf5'  }

def get_dir(sim_name):
    return "/gpfs/data/aquarius/halo_data/Aq-%s/%c/" % (sim_label[sim_name[0:3]], sim_name[3])

def load_last_snapshot(sim_name):
    return iccpy.gadget.load_snapshot(directory=get_dir(sim_name), snapnum=last_snapnum[sim_name])

def get_subhaloes(sim_name, snapnum=None):
    if snapnum==None:
        snapnum=last_snapnum[sim_name]

    catalogue = iccpy.gadget.SubfindCatalogue(get_dir(sim_name), snapnum)
    return catalogue.subhalo

def get_halo_centre(sim_name):
    return get_subhaloes(sim_name)[0].pot_min

def get_merger_tree(sim_name):
    return MergerTree(merger_tree_filename[sim_name])

def plot(plot_func, haloes=['A', 'B', 'C', 'D', 'E', 'F'], legend=None, tick_length=8, minor_tick_x_space=None, minor_tick_y_space=None):
    from matplotlib.ticker import MultipleLocator

    haloes = np.array(haloes)

    # no space between the panels
    pl.rcParams.update({'figure.subplot.wspace':0,'figure.subplot.hspace':0})

    all_haloes = np.array(['A', 'B', 'C', 'D', 'E', 'F'])

    plotIdxs = np.sort(iccpy.utils.match(haloes, all_haloes))

    numRows = 3
    numCols = 2

    for i in plotIdxs:
        ax = pl.subplot(numRows,numCols,i+1)
        plot_func(all_haloes[i], ax)

        #Tidy up plot
        if minor_tick_y_space is not None:
            ax.yaxis.set_minor_locator(MultipleLocator(minor_tick_y_space))    

        if minor_tick_x_space is not None:
            ax.xaxis.set_minor_locator(MultipleLocator(minor_tick_x_space))

        left_tick = i%numCols==0 or i-1 not in plotIdxs
        ax.yaxis.get_label().set_visible(left_tick)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1On=left_tick
            tick.tick1line.set_markersize(tick_length)
            tick.tick2line.set_markersize(tick_length)

        if left_tick and i-numCols in plotIdxs:
            lims = ax.get_ylim()
            ax.set_ylim(lims[0], 0.9999999999*lims[1])

        lower_tick = i>=(numRows-1)*numCols or i+numCols not in plotIdxs
        ax.xaxis.get_label().set_visible(lower_tick)
        for tick in ax.xaxis.get_major_ticks():
            tick.label1On=lower_tick
            tick.tick1line.set_markersize(tick_length)
            tick.tick2line.set_markersize(tick_length)

        for tick in ax.yaxis.get_minor_ticks() + ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(tick_length/2)
            tick.tick2line.set_markersize(tick_length/2)

        if lower_tick and i+1 in plotIdxs and (i+1)%numCols!=0:
            lims = ax.get_xlim()
            ax.set_xlim(lims[0], 0.9999999999*lims[1])

        if lower_tick and not left_tick and i<(numRows-1)*numCols:
            lims = ax.get_xlim()
            ax.set_xlim(lims[0]*1.0000000001, lims[1])

        if ax.get_legend() is not None:
            if legend is None:
                ax.get_legend().set_visible(False)
            elif legend is 'All' or legend is 'all' or all_haloes[i] in legend:
                ax.get_legend().draw_frame(False)
            else:
                ax.get_legend().set_visible(False)

def plot_test(halo, ax):
    ax.plot([1,2], [3,4])
    pl.legend(['FISH'])
    pl.ylabel('x')

if __name__=="__main__":
    #print load_last_snapshot("aqa4")
    #print get_halo_centre("aqa4")

    plot(plot_test, haloes=['A', 'B', 'C', 'D', 'E'], legend='A', minor_tick_x_space=0.025)
    #plot(plot_test, minor_tick_x_space=0.025)
    pl.show()

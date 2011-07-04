import os
from tables import openFile
from interp3d import interp3d
from numpy import float64, empty, array


elements = ('Calcium', 'Carbon', 'Iron', 'Magnesium', 'Neon', 'Nitrogen', 'Oxygen', 'Silicon', 'Sulphur')
helium_fracs = ('0.24', '0.25', '0.26', '0.27', '0.28', '0.29', '0.30')

class GimicCooling:
    def __init__(self, directory='/gpfs/data/rmdq85/wall_shock/Gadget/gadget3_cosma4/P-Gadget3-BG/BG_Tables/OldCoolingTables/',
                  helium_frac='0.24', Z=0, uv_background=False):
# directory='/data/rw13/dph3nlm/BG_Tables/LowRes/CoolingTables/'
        files = get_gimic_files(directory)
        redshifts = []
        
        name = os.path.join(directory, files[0])
        f = openFile(name, 'r')

        header = f.root.Header
        self.hydrogen_density_bins = header.Hydrogen_density_bins.read()
        self.temperature_bins = header.Temperature_bins.read()
        self.redshift = header.Redshift.read()

        f.close()

        solar_ratio = 10**Z
        cool_rates = empty((len(files), 
                            len(self.hydrogen_density_bins),
                            len(self.temperature_bins)), dtype=float64)

        for i,file in enumerate(files):
            name = os.path.join(directory, file)
            f = openFile(name, 'r')

            hhe = getattr(f.root.Metal_free, '_'.join(('Helium', helium_frac)))
            
            cooling = hhe.Cooling.read().astype(float64)
            heating = hhe.Heating.read().astype(float64)

            for element in elements:
                node = getattr(f.root, element)
                cooling += node.Cooling.read() * solar_ratio
                heating += node.Heating.read() * solar_ratio

            if uv_background:
                cool_rates[i] = cooling - heating
            else:
                cool_rates[i] = cooling

            redshifts.append(f.root.Header.Redshift.read())

            f.close()

        self.redshifts = array(redshifts, dtype=float64)
        self.cool_rates = cool_rates
        
    def cooling_rate(self, n_H, T, z):
        """ 
        interpolates the cooling rates at the given values
        n_H = hydrogen number density (/cm3)
        T = temperature (K)
        z = redshift
        
        all numpy broadcastable
        currently uses linear interpolation
        """
        

        cool = interp3d(self.redshifts, self.hydrogen_density_bins, self.temperature_bins, 
                        self.cool_rates, z, n_H, T)

        return cool
 

class ElementCooling:
    def __init__(self,  element, directory='/gpfs/data/rmdq85/wall_shock/Gadget/gadget3_cosma4/P-Gadget3-BG/BG_Tables/OldCoolingTables/',
              Z=0, uv_background=False):

                 #directory='/data/rw13/dph3nlm/BG_Tables/LowRes/CoolingTables/'
        if element not in elements:
            raise Exception('Element must be one of '+str(elements))

        files = get_gimic_files(directory)
        redshifts = []
        
        name = os.path.join(directory, files[0])
        f = openFile(name, 'r')

        header = f.root.Header
        self.hydrogen_density_bins = header.Hydrogen_density_bins.read()
        self.temperature_bins = header.Temperature_bins.read()
        self.redshift = header.Redshift.read()

        f.close()

        solar_ratio = 10**Z
        cool_rates = empty((len(files), 
                            len(self.hydrogen_density_bins),
                            len(self.temperature_bins)), dtype=float64)

        for i,file in enumerate(files):
            name = os.path.join(directory, file)
            f = openFile(name, 'r')

            node = getattr(f.root, element)
            cooling = node.Cooling.read().astype(float64) * solar_ratio
            heating = node.Heating.read().astype(float64) * solar_ratio

            if uv_background:
                cool_rates[i] = cooling - heating
            else:
                cool_rates[i] = cooling

            redshifts.append(f.root.Header.Redshift.read())

            f.close()

        self.redshifts = array(redshifts, dtype=float64)
        self.cool_rates = cool_rates
        
    def cooling_rate(self, n_H, T, z):
        """ 
        interpolates the cooling rates at the given values
        n_H = hydrogen number density (/cm3)
        T = temperature (K)
        z = redshift
        
        all numpy broadcastable
        currently uses linear interpolation
        """
        

        cool = interp3d(self.redshifts, self.hydrogen_density_bins, self.temperature_bins, 
                        self.cool_rates, z, n_H, T)

        return cool
 
def get_gimic_files(directory):
    """ find the files in the gimic directory """
    files = os.listdir(directory)
    files = sorted(f for f in files if f[:2] == 'z_')
    # ignore nocompton and photodis
    files =  files[:-2]
    return files



def test():
    """ plot a graph of the cooling rate """
    import pylab as pl
    from numpy import arange
    T = 1e4 * (1e8/1e4)**(arange(1000)/999.0)
    n_H = 0.1
    gimic = GimicCooling(helium_frac='0.24', Z=0, uv_background=False)
    lamb = gimic.cooling_rate(n_H, T, z=0.0)
    pl.loglog(T, lamb)
    pl.ylabel(r'$\Lambda/n_\mathrm{H}^2 \, (\mathrm{erg cm^3 /s})$')
    pl.xlabel(r'$T(\mathrm{K})$')
    pl.show()


if __name__=='__main__':

    test()

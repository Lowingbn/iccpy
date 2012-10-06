import unittest
from iccpy.particles.lattices import *
import numpy as np

class LatticeTestCase(unittest.TestCase):
        
    def testNearestNeighbourCCP(self):
        pos = ccp_lattice(5,10,15)
        idx = 21

        dist = np.sqrt(np.square(pos-pos[idx]).sum(1))
        
        dist = sorted(dist)
        self.assertEqual(dist[1], np.sqrt(0.5))
        self.assertEqual(pos.shape[0], 5*10*15*4)

    def testNearestNeighbourHCP(self):
        
        pos, box_size = hcp_lattice(5,10,15)
        idx = 21

        dist = np.sqrt(np.square(pos-pos[idx]).sum(1))
        
        dist = sorted(dist)
        self.assertAlmostEqual(dist[1], 1.0,10)
        self.assertEqual(pos.shape[0], 5*10*15*4)

ccpCase = LatticeTestCase('testNearestNeighbourCCP')
hcpCase = LatticeTestCase('testNearestNeighbourHCP')

if __name__=='__main__':
    unittest.main()

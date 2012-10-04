import unittest
import iccpy.gadget
import iccpy.gadget.binary_snapshot_io

class GadgetBinaryIOTestCase(unittest.TestCase):
    
    def test_load_snapshot(self):
        header, data = iccpy.gadget.binary_snapshot_io.read_snapshot_file('data/gadget_binary_snapshot_000')
        
        self.assertEqual(header['num_particles'][1], 2000)
        self.assertEqual(len(data['pos']), 2000)
        self.assertEqual(len(data['vel']), 2000)
        self.assertEqual(len(data['id']), 2000)
    
class GadgetSnapshotTestCase(unittest.TestCase):

    def test_load_snapshot(self):
        snapshot = iccpy.gadget.Snapshot(directory="data", filename="gadget_binary_snapshot", snapnum=0)

        self.assertEqual(snapshot.header.num_particles[1], 2000)
        self.assertEqual(len(snapshot.pos[1]), 2000)
        self.assertEqual(len(snapshot.vel[1]), 2000)
        self.assertEqual(len(snapshot.id[1]), 2000)

if __name__ == '__main__':
    unittest.main()

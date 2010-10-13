import unittest
import iccpy.gadget.snapshot

class GadgetSnapshotTestCase(unittest.TestCase):
    
    def test_load_snapshot(self):
        header, data = iccpy.gadget.snapshot.load_snapshot_file('data', 'gadget_binary_snapshot', 0)
        
        self.assertEqual(header['npart'][1], 2000)
        self.assertEqual(len(data['pos']), 2000)
        self.assertEqual(len(data['vel']), 2000)
        self.assertEqual(len(data['id']), 2000)
    
if __name__ == '__main__':
    unittest.main()

import unittest
import numpy as np
from atlas import display_group, ontology, lookup
from spatial import consensus, deduplicate, match_points


class SpatialChecks(unittest.TestCase):
    def test_no_consensus_chaining(self):
        points = {0: np.array([[0., 0.]]), 1: np.array([[0., 7.]]), 2: np.array([[0., 14.]])}
        self.assertEqual(len(consensus(points, 1, minimum=3)[0]), 0)

    def test_one_member_per_detector(self):
        points = {0: np.array([[0., 0.], [1., 0.]]), 1: np.array([[0., 1.]])}
        self.assertEqual(len(consensus(points, 1)[0]), 1)

    def test_maximum_matching_is_not_nearest_neighbor_greedy(self):
        a = np.array([[0., 0.], [0., 2.]])
        b = np.array([[0., 1.], [0., -1.]])
        self.assertEqual(int((match_points(a, b, 1, radius_um=1) >= 0).sum()), 2)

    def test_same_tile_cells_are_not_deduplicated(self):
        keep, _ = deduplicate(np.array([[0., 0.], [0., 1.]]), np.zeros((2, 2)), np.ones(2), 1)
        self.assertEqual(len(keep), 2)

    def test_empty_matching(self):
        self.assertEqual(match_points(np.empty((0, 2)), np.empty((0, 2)), 1).shape, (0,))


class AtlasChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nodes, cls.paths = ontology()

    def test_ca_fields_are_not_layers(self):
        for sid in [382, 423, 463]:
            self.assertEqual(display_group(sid, self.nodes, self.paths), sid)

    def test_named_singular_and_plural_layers(self):
        for sid, parent in {844: 985, 810: 48, 427: 895, 1140: 566, 1141: 566, 1142: 566}.items():
            self.assertEqual(display_group(sid, self.nodes, self.paths), parent)

    def test_yx_order_and_last_pixel_rounding(self):
        labels = np.array([[0, 12, 23], [31, 42, 53]], dtype=np.uint32)
        self.assertEqual(lookup([[1.2, 2.8]], labels).tolist(), [53])
        with self.assertRaises(ValueError):
            lookup([[2, 0]], labels)


if __name__ == '__main__':
    unittest.main()

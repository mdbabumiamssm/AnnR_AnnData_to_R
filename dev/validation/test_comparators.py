"""Independent negative controls. Run with the Python reference environment."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

read = load("compare_read")
rt = load("compare_roundtrip")

class NegativeControls(unittest.TestCase):
    def test_read_detects_tiny_value_change(self):
        c = read.Checker(Path(tempfile.gettempdir()))
        self.assertFalse(c.same_numeric("X", np.array([1.0 + 1e-7]), np.array([1.0])))
        self.assertTrue(c.issues)

    def test_read_detects_moved_index(self):
        c = read.Checker(Path(tempfile.gettempdir()))
        self.assertFalse(c.same_numeric("indices", np.array([0, 2]), np.array([0, 1])))
        self.assertTrue(c.issues)

    def test_roundtrip_detects_transpose(self):
        a = np.array([[1., 2.], [3., 4.]])
        self.assertFalse(rt.values_equal(a, a.T))
        self.assertTrue(rt.values_equal(a, a.copy()))

    def test_roundtrip_detects_tiny_value_change(self):
        self.assertFalse(rt.values_equal(np.array([1.]), np.array([1. + 1e-7])))

if __name__ == "__main__":
    unittest.main()

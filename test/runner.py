# Run tests. python3 test/runner.py

import unittest

# initialize the test suite
loader = unittest.TestLoader()
suite  = unittest.TestSuite()

import test_shex
import test_sparql

suite.addTests(loader.loadTestsFromModule(test_shex))
suite.addTests(loader.loadTestsFromModule(test_sparql))

# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)

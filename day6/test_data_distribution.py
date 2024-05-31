import unittest
import pandas as pd
import numpy as np
from data_distribution import load_data, calculate_statistics  # replace 'your_main_script' with the name of your script file

class TestStatistics(unittest.TestCase):

    def setUp(self):
        self.data = pd.DataFrame({
            'Dataset1': ['1', '3'],
            'Dataset2': ['4', '6']
        })

    def test_calculate_statistics(self):
        self.expected_means = {'Dataset1': 2, 'Dataset2': 5}
        self.expected_stds = {'Dataset1': 1, 'Dataset2': 1}

if __name__ == '__main__':
    unittest.main()

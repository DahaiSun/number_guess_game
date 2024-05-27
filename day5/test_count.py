import unittest
from library import count_lines, count_words, count_letters, process_file

class TestLibraryFunctions(unittest.TestCase):

    def test_count_lines(self):
        text = "Hello World\nHello World"
        self.assertEqual(count_lines(text), 2)

    def test_count_words(self):
        text = "Hello World\nHello World."
        self.assertEqual(count_words(text), 4)
    
    def test_count_letters(self):
        text = "Hello World"
        expected_counts = {chr(i): 0 for i in range(ord('a'), ord('z') + 1)}
        expected_counts['h'] = 1
        expected_counts['e'] = 1
        expected_counts['l'] = 3
        expected_counts['o'] = 2
        expected_counts['w'] = 1
        expected_counts['r'] = 1
        expected_counts['d'] = 1
        counts, total = count_letters(text)
        self.assertEqual(counts, expected_counts)
        self.assertEqual(total, 10)
    
if __name__ == '__main__':
    unittest.main()
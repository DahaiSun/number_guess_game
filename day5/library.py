
def count_lines(text):
    lines_count = len(text.split('\n'))
    return lines_count

def count_words(text):
    words_count = 0
    for line in text.split('\n'):
        words_count += len(line.split())
    return words_count

def count_letters(text):
    letter_counts = {chr(i): 0 for i in range(ord('a'), ord('z') + 1)}
    total_letters = 0
    for char in text.lower():
        if char in letter_counts:
            letter_counts[char] += 1
            total_letters += 1
    return letter_counts, total_letters

def process_file(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        text = file.read()
    total_lines = count_lines(text)
    total_words = count_words(text)
    letter_counts, total_letters = count_letters(text)
    return total_lines, total_words, letter_counts, total_letters
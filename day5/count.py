import sys
from library import process_file

def print_counts(filename):
    total_lines, total_words, letter_counts, total_letters = process_file(filename)
    for letter, count in sorted(letter_counts.items()):
        print(f"{letter}:{count}")
    print(f"Total letter count: {total_letters}")
    print(f"Total Lines: {total_lines}")
    print(f"Total words: {total_words}")
    print("done")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        print_counts(sys.argv[1])
    else:
        print("Usage: python count.py FILENAME.txt")

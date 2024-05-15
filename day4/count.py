import sys

def count(filename):
    letter_counts = {chr(i): 0 for i in range(ord('a'), ord('z') + 1)}
    total_letters = 0
    lines = 0
    words = 0

    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            for char in line.lower():
                if char in letter_counts:
                    letter_counts[char] += 1
                    total_letters += 1
            lines += 1
            words += len(line.split())

    for letter, count in  sorted(letter_counts.items()):
        print(f"{letter}:{count}")

    print (f"Totla letter count: {total_letters}")
    print (f"Total Lines: {lines}")
    print (f"Total words: {words}")
    print("done")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        count(sys.argv[1])
    else:
        print("Usage: python count.py FILENAME.txt")
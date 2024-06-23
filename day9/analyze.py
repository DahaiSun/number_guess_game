import sys
import re
from Bio import SeqIO
import argparse

def read_sequence(file_path):
    # 根据文件扩展名确定文件格式
    file_format = "fasta" if file_path.endswith((".fasta", ".fa")) else "genbank" if file_path.endswith(".gb") else None
    if file_format is None:
        raise ValueError("Unsupported file format. Please provide a FASTA or GeneBank file.")
    
    with open(file_path, "r") as file:
        record = SeqIO.read(file, format=file_format)
        return str(record.seq)

def duplicate(dna):
    length = 1
    result = ''
    while True:
        regex = r'([GATC]{' + str(length) + r'}).*\1'
        m = re.search(regex, dna)
        if m:
            result = m.group(1)
            length += 1
        else:
            break
    return result

def gc_content(sequence):
    return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

def high_GC_area(dna, window_size=50):
    max_gc_content = 0
    max_gc_sequence = ''
    for i in range(len(dna) - window_size + 1):
        window = dna[i:i + window_size]
        current_gc_content = gc_content(window)
        if current_gc_content > max_gc_content:
            max_gc_content = current_gc_content
            max_gc_sequence = window
    return max_gc_sequence, max_gc_content

def parse_arguments():
    parser = argparse.ArgumentParser(description="DNA Analyzer")
    parser.add_argument("file_path", type=str, help="Path to the FASTA or GeneBank file")
    parser.add_argument("function", type=str, choices=["duplicate", "high_GC_area"], help="Function to use: 'duplicate' or 'high_GC_area'")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    file_path = args.file_path
    function = args.function
    
    try:
        dna_sequence = read_sequence(file_path)
        
        if function == "duplicate":
            result = duplicate(dna_sequence)
            print(f"The longest repeated subsequence is: {result}")
        
        elif function == "high_GC_area":
            highest_gc_sequence, highest_gc_content = high_GC_area(dna_sequence)
            print(f"The 50bp region with the highest GC content is: {highest_gc_sequence} with a GC content of {highest_gc_content:.2f}%")
    
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

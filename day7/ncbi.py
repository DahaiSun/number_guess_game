import sys
import os
import csv
from datetime import datetime
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
from io import StringIO

def search_ncbi(term, number, organisms):
    Entrez.email = "sundahai63@gmail.com"
    ids = []

    query_term = f"{term} AND ({' OR '.join([f'{organism}[Organism]' for organism in organisms])})"
    handle = Entrez.esearch(db='protein', term=query_term, retmax=number)
    record = Entrez.read(handle)
    handle.close()
    ids.extend(record['IdList'])
    total = int(record['Count'])

    return ids[:number], total

def fetch_definitions(ids):
    definitions = []
    for id in ids:
        try:
            handle = Entrez.efetch(db='protein', id=id, rettype='gb', retmode="text")
            content = handle.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            handle.close()
            gb_record = SeqIO.read(StringIO(content), "genbank")
            definitions.append(gb_record.description)
        except HTTPError as e:
            print(f"Failed to fetch {id} from protein database: {e}")
        except Exception as e:
            print(f"An error occurred while fetching {id}: {e}")

    return definitions

def fetch_and_save(ids):
    filenames = []
    rettype = "gb"
    extension = "gb"

    for id in ids:
        try:
            handle = Entrez.efetch(db='protein', id=id, rettype=rettype, retmode="text")
            content = handle.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            filename = f"protein_{id}.{extension}"
            with open(filename, 'w') as file:
                file.write(content)
            handle.close()
            filenames.append(filename)
        except HTTPError as e:
            print(f"Failed to fetch {id} from protein database: {e}")
        except Exception as e:
            print(f"An error occurred while fetching {id}: {e}")

    return filenames

def log_results(log_file, date, term, number, total, organisms):
    log_exists = os.path.isfile(log_file)
    with open(log_file, 'a', newline='') as csvfile:
        fieldnames = ['date', 'term', 'max', 'total', 'database', 'organisms']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not log_exists:
            writer.writeheader()
        writer.writerow({'date': date, 'term': term, 'max': number, 'total': total, 'database': 'protein', 'organisms': ', '.join(organisms)})

def main():
    if len(sys.argv) < 4:
        print("What to search? Usage: python3 ncbi.py TERM NUMBER ORGANISMS")
        print("Example: python3 ncbi.py STAT5A 10 'Homo sapiens,Mus musculus'")
        sys.exit(1)
    
    term = sys.argv[1]
    number = int(sys.argv[2])
    organisms = sys.argv[3].split(',')

    ids, total = search_ncbi(term, number, organisms)
    print(f"Found {total} records for term '{term}' in protein database with organisms {', '.join(organisms)}")
    print("Wait for half a minute, Loading search results...")
    # Fetch definitions for displaying search results without saving files
    definitions = fetch_definitions(ids)

    print("\nSearch results:\n")
    for i, definition in enumerate(definitions):
        print(f"{i+1}. {definition}")

    # Ask user to choose which results to download
    choices = input("\nEnter the numbers of the items to download, separated by commas (e.g., 1,2,3,4): ")
    choices = [int(choice.strip()) for choice in choices.split(',')]
    chosen_ids = [ids[i-1] for i in choices]

    chosen_filenames = fetch_and_save(chosen_ids)
    
    for filename in chosen_filenames:
        print(f"Saved: {filename}")

    # Log search and download results
    log_results('search_log.csv', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), term, number, total, organisms)
    log_results('download_log.csv', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), term, len(chosen_ids), len(chosen_ids), organisms)

if __name__ == "__main__":
    main()

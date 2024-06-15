import sys
import os
import csv
from datetime import datetime
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
from io import StringIO
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext

def search_ncbi(term, number, organisms, db):
    Entrez.email = "sundahai63@gmail.com"
    ids = []

    query_term = f"{term} AND ({' OR '.join([f'{organism}[Organism]' for organism in organisms])})"
    handle = Entrez.esearch(db=db, term=query_term, retmax=number)
    record = Entrez.read(handle)
    handle.close()
    ids.extend(record['IdList'])
    total = int(record['Count'])

    return ids[:number], total

def fetch_definitions(ids, db):
    definitions = []
    for id in ids:
        try:
            handle = Entrez.efetch(db=db, id=id, rettype='gb', retmode="text")
            content = handle.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            handle.close()
            gb_record = SeqIO.read(StringIO(content), "genbank")
            definitions.append(gb_record.description)
        except HTTPError as e:
            print(f"Failed to fetch {id} from {db} database: {e}")
        except Exception as e:
            print(f"An error occurred while fetching {id}: {e}")

    return definitions

def fetch_and_save(ids, db, rettype):
    filenames = []
    extension = rettype.lower()

    for id in ids:
        try:
            handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode="text")
            content = handle.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            filename = f"{db}_{id}.{extension}"
            with open(filename, 'w') as file:
                file.write(content)
            handle.close()
            filenames.append(filename)
        except HTTPError as e:
            print(f"Failed to fetch {id} from {db} database: {e}")
        except Exception as e:
            print(f"An error occurred while fetching {id}: {e}")

    return filenames

def log_results(log_file, date, term, number, total, db, organisms):
    log_exists = os.path.isfile(log_file)
    with open(log_file, 'a', newline='') as csvfile:
        fieldnames = ['date', 'term', 'max', 'total', 'database', 'organisms']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not log_exists:
            writer.writeheader()
        writer.writerow({'date': date, 'term': term, 'max': number, 'total': total, 'database': db, 'organisms': ', '.join(organisms)})

def search():
    try:
        term = term_entry.get()
        number = int(number_entry.get())
        organisms = organisms_selector.get().split(',')
        db = db_selector.get()

        if not term or not number or not organisms:
            raise ValueError("Please fill in all fields before searching.")

        results_text.delete(1.0, tk.END)
        results_text.insert(tk.END, "Searching... Please wait.\n")
        results_text.see(tk.END)  # Scroll to the end
        app.update()

        ids, total = search_ncbi(term, number, organisms, db)
        results_text.delete(1.0, tk.END)
        results_text.insert(tk.END, f"Found {total} records for term '{term}' in {db} database with organisms {', '.join(organisms)}\n")
        results_text.insert(tk.END, "Wait for half a minute, Loading search results...\n")
        results_text.see(tk.END)  # Scroll to the end

        definitions = fetch_definitions(ids, db)

        results_text.insert(tk.END, "\nSearch results:\n")
        for i, definition in enumerate(definitions):
            results_text.insert(tk.END, f"{i+1}. {definition}\n")
        results_text.see(tk.END)  # Scroll to the end

        global search_ids
        search_ids = ids
    except ValueError as ve:
        messagebox.showerror("Input Error", "Check Input")
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")

def download():
    try:
        choices = [int(choice.strip()) for choice in download_entry.get().split(',')]
        if not choices:
            raise ValueError("Please enter the numbers of the items to download.")
        chosen_ids = [search_ids[i-1] for i in choices if i-1 < len(search_ids)]
        rettype = format_selector.get()

        chosen_filenames = fetch_and_save(chosen_ids, db_selector.get(), rettype)

        for filename in chosen_filenames:
            results_text.insert(tk.END, f"Saved: {filename}\n")
            results_text.see(tk.END)  # Scroll to the end

        term = term_entry.get()
        number = int(number_entry.get())
        organisms = organisms_selector.get().split(',')
        db = db_selector.get()
        log_results('search_log.csv', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), term, number, len(search_ids), db, organisms)
        log_results('download_log.csv', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), term, len(chosen_ids), len(chosen_ids), db, organisms)

        results_text.insert(tk.END, "\nDownload complete.\n")  # Provide feedback when download is finished
        results_text.see(tk.END)  # Scroll to the end
    except ValueError as ve:
        messagebox.showerror("Input Error", "Check Input")
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")

def exit_program():
    app.quit()

app = tk.Tk()
app.title("NCBI Search and Download")
app.geometry("600x800")

tk.Label(app, text="Search Term:").grid(row=0, column=0, padx=10, pady=5, sticky="w")
term_entry = tk.Entry(app)
term_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Label(app, text="How many results to show:").grid(row=1, column=0, padx=10, pady=5, sticky="w")
number_entry = tk.Entry(app)
number_entry.grid(row=1, column=1, padx=10, pady=5, sticky="ew")

tk.Label(app, text="Database:").grid(row=2, column=0, padx=10, pady=5, sticky="w")
db_selector = ttk.Combobox(app, values=["protein", "nucleotide"])
db_selector.grid(row=2, column=1, padx=10, pady=5, sticky="ew")
db_selector.current(0)

tk.Label(app, text="Organisms:").grid(row=3, column=0, padx=10, pady=5, sticky="w")
organisms_selector = ttk.Combobox(app, values=["Homo sapiens", "Mus musculus"])
organisms_selector.grid(row=3, column=1, padx=10, pady=5, sticky="ew")
organisms_selector.current(0)

search_button = tk.Button(app, text="Search", command=search)
search_button.grid(row=4, column=0, columnspan=2, pady=10, sticky="ew")

results_text = scrolledtext.ScrolledText(app, width=50, height=10)
results_text.grid(row=5, column=0, columnspan=2, padx=10, pady=5, sticky="nsew")

tk.Label(app, text="Enter the numbers of the items to download (comma-separated):").grid(row=6, column=0, columnspan=2, padx=10, pady=5, sticky="w")
download_entry = tk.Entry(app)
download_entry.grid(row=7, column=0, columnspan=2, padx=10, pady=5, sticky="ew")

tk.Label(app, text="File Format:").grid(row=8, column=0, padx=10, pady=5, sticky="w")
format_selector = ttk.Combobox(app, values=["gb", "fasta", "txt"])
format_selector.grid(row=8, column=1, padx=10, pady=5, sticky="ew")
format_selector.current(0)

download_button = tk.Button(app, text="Download", command=download)
download_button.grid(row=9, column=0, columnspan=2, pady=10, sticky="ew")

exit_button = tk.Button(app, text="Exit", command=exit_program)
exit_button.grid(row=10, column=0, columnspan=2, pady=10, sticky="ew")

# Configure grid resizing
app.grid_columnconfigure(1, weight=1)
app.grid_rowconfigure(5, weight=1)

app.mainloop()

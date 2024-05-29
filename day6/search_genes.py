import pandas as pd
import openpyxl
import tempfile
import os 

try:
    # Read data from the specified sheets
    data_df = pd.read_excel('normallized_reads.xlsx', sheet_name='Sheet1')
    genes_to_search_df = pd.read_excel('normallized_reads.xlsx', sheet_name='search_genes')
    print("Excel loaded successfully.")
except FileNotFoundError:
    print("failed")
    exit()
except Exception as e:
    print("error", e)
    exit()


matching_rows = pd.DataFrame()# Initialize an empty DataFrame for results

# search target genes listed in the "targets" column, search area is the "Gene" column in sheet,
# keep empty value for the none found genes
for gene_name in genes_to_search_df['targets'].values: 
   
    match = data_df[data_df['Gene'] == gene_name]
    if not match.empty:
        matching_rows = pd.concat([matching_rows, match], ignore_index=True)
    else: 
        no_match_row = pd.DataFrame({'Gene': [gene_name]})
        matching_rows = pd.concat([matching_rows, no_match_row], ignore_index=True)

# output search result in a temperate excel window
try:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
        temp_filename = tmp.name
    with pd.ExcelWriter(temp_filename, engine='openpyxl', mode='w') as writer:
        matching_rows.to_excel(writer, sheet_name='Matched_Genes', index=False)
    print("Completed successfully. Opening the temporary Excel file...")

 
    os.system(f'start excel "{temp_filename}"') # open the temorary excel file
except Exception as e:
    print("Failed to write to Excel file:", e)
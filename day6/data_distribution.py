import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os
import sys
import platform

def load_data(file_path, sheet_name):
    data = pd.read_excel(file_path, sheet_name=sheet_name)
    return data

def draw_histograms(data):
    plt.figure(figsize=(18, 8))

    plt.subplot(2, 1, 1)
    plt.hist(data['Dataset1'], bins=60, edgecolor='k', alpha=0.7)
    plt.title('Histogram of Dataset1')
    plt.xlabel('Value')
    plt.ylabel('Frequency')

    plt.subplot(2, 1, 2)
    plt.hist(data['Dataset2'], bins=60, edgecolor='k', alpha=0.7)
    plt.title('Histogram of Dataset2')
    plt.xlabel('Value')
    plt.ylabel('Frequency')

    plt.tight_layout()
    plt.show()

def calculate_statistics(data):
    mu_dataset1 = np.mean(data['Dataset1'])
    sigma_dataset1 = np.std(data['Dataset1'])

    mu_dataset2 = np.mean(data['Dataset2'])
    sigma_dataset2 = np.std(data['Dataset2'])

    results = pd.DataFrame({
        'Dataset': ['Dataset1', 'Dataset2'],
        'Mean': [mu_dataset1, mu_dataset2],
        'Standard Deviation': [sigma_dataset1, sigma_dataset2]
    })
    return results

def save_to_excel(results):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
            temp_filename = tmp.name
        with pd.ExcelWriter(temp_filename, engine='openpyxl', mode='w') as writer:
            results.to_excel(writer, sheet_name='statistics', index=False)

        if platform.system() == 'Windows':
            os.system(f'start excel "{temp_filename}"')
        else:
            print(f"Completed successfully. The temporary Excel file is saved at {temp_filename}")
            print("Please open the file manually.")
    except Exception as e:
        print("Failed to write to Excel file:", e)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 data_distribution.py <file_path> <sheet_name>")
        sys.exit(1)

    file_path = sys.argv[1]
    sheet_name = sys.argv[2]
    
    data = load_data(file_path, sheet_name)
    print("File loaded")

    draw_histograms(data)
    print("Figure generated")

    results = calculate_statistics(data)
    print("Calculation completed")

    save_to_excel(results)

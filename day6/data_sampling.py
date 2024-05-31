import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os

file_path = r'E:\github\wis_python_course_assignments\day6\StatisticsQ1.xlsx'
sheet_name = 'Sheet1'
data = pd.read_excel(file_path, sheet_name=sheet_name)
population1 = data['Dataset1']
population2 = data['Dataset2']

# sampling Distributions
def sampling_distribution(data, sample_size, num_samples, statistic):
    samples = np.random.choice(data, (num_samples, sample_size))
    if statistic == 'mean':
        return np.mean(samples, axis=1)
    elif statistic == 'median':
        return np.median(samples, axis=1)
    elif statistic == 'variance':
        return np.var(samples, axis=1)

# sampling from dataset 1
n1 = 10
num_samples = 1000

sampling_means1 = sampling_distribution(population1, n1, num_samples, 'mean')
sampling_medians1 = sampling_distribution(population1, n1, num_samples, 'median')
sampling_variances1 = sampling_distribution(population1, n1, num_samples, 'variance')

# sampling from dataset 2
n2 = 30

sampling_means2 = sampling_distribution(population2, n2, num_samples, 'mean')
sampling_medians2 = sampling_distribution(population2, n2, num_samples, 'median')
sampling_variances2 = sampling_distribution(population2, n2, num_samples, 'variance')

# plot sampling calculation results
plt.figure(figsize=(18, 12))

plt.subplot(3, 2, 1)
plt.hist(sampling_means1, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Mean (n=10) from Normal Population Dataset1')

plt.subplot(3, 2, 2)
plt.hist(sampling_means2, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Mean (n=30) from Skewed Population Dataset2')

plt.subplot(3, 2, 3)
plt.hist(sampling_medians1, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Median (n=10) from Normal Population Dataset1')

plt.subplot(3, 2, 4)
plt.hist(sampling_medians2, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Median (n=30) from Skewed Population Dataset2')

plt.subplot(3, 2, 5)
plt.hist(sampling_variances1, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Variance (n=10) from Normal Population Dataset1')

plt.subplot(3, 2, 6)
plt.hist(sampling_variances2, bins=60, edgecolor='k', alpha=0.7)
plt.title('Sampling Distribution of Variance (n=30) from Skewed Population Dataset2')

plt.tight_layout()
plt.show()
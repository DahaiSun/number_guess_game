
1. data_distribution.py

    Usage: python3 data_distribution.py <file_path> <sheet_name>"
    in my case the command should be：python3 data_distribution.py StatisticsQ1.xlsx Sheet1
    
    work-flow
    1) load excel file from command
    2) Draw histogram of 2 dataset 
    3) calcilate 𝜇 and 𝜎 for each dataset
    4) output result in a temporary window if using windows windows system, or print saving path of the output file

2. data_sampling.py

    1）From the normal population: Sampling distribution of the mean 𝑥̅, median and variance 𝑆
    
    2） From the skewed population: Sampling distribution of the mean 𝑥̅, median and variance 𝑆
2
for samples of size 𝑛 = 30
for samples of size 𝑛 = 10


3. search_gene.py
    1) load excel file
    2) search target genes listed in the "targets" column, search area is the "Gene" column in sheet
       <br> keep empty value for the none found genes
    3) output search results into temporary window

    how to use:
    1) match the file name with to the real file
    2) make a new sheet in the excel file named"search_genes"
      <br>make a column named"targets" 
       <br>list target genes under the column name:
       <br>for example:
        |targets|
        |-------|
        |STAT1|
        |STAT2|
        |STAT3|
        |STAT4|
    3) Data sheet name "sheet1"
        genes should under the column that named “Gene”
1. search_gene.py
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
2. data_distribution.py
    1) load excel file
    2) Draw histogram of 2 dataset 
    3) calcilate 𝜇 and 𝜎 for each dataset
    4) output result in a temporary window
    
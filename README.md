# The pipeline of DEG identification based on Oliver's data
## The download of raw data
Download the metadata about D. melanogaster using SRA Run selector of NCBI. Than run the python script 1run_download.py. 
python 1run_download.py /home/chuand/new_gene/data/dm_SraRunTable.txt /home/chuand/new_gene/data/dm_sra。
The first parameter /home/chuand/new_gene/data/dm_SraRunTable.txt is the metadata that downloaded by SRA Run selector. The second parameter is the path that the raw will be stored.

# The pipeline of DEG identification based on Oliver's data
## The download of raw data
Download the metadata about _D. melanogaster_ using SRA Run selector of NCBI. Than run the python script 1run_download.py. <br />
```
python 1run_download.py /home/chuand/new_gene/data/dm_SraRunTable.txt /home/chuand/new_gene/data/dm_sra
``` 
<br />
The first parameter /home/chuand/new_gene/data/dm_SraRunTable.txt is the metadata that downloaded by SRA Run selector. The second parameter is the path that the raw will be stored. <br />

## The identification of DEGs
Identification the sex-biased genes according to the parameter provided by the users <br />
```
nohup python3 deg.py "Drosophila melanogaster" gonad /home/chuand/new_gene/data/dm_sra /home/chuand/new_gene/data/dm_fastq /home/chuand/new_gene/data/dm_cleanfastq /home/chuand/new_gene/data/mapping /home/chuand/new_gene/data/bam /home/chuand/new_gene/data/mergebam  /home/chuand/new_gene/result/summary &
```
The first parameter "Drosophila melanogaster" is the species that we want to analysis, the script will retrieve information from the metadata (For example /home/chuand/new_gene/data/dm_SraRunTable.txt) that is downloaded from NCBI.
The second parameter gonad is the focal tissue, in which the DEGs will be analyzed.

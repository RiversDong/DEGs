# The pipeline of DEG identification based on Oliver's data
## The download of raw data
Download the metadata about _D. melanogaster_ using SRA Run selector of NCBI. Than run the python script 1run_download.py. <br />
```python
python 1run_download.py /home/chuand/new_gene/data/dm_SraRunTable.txt /home/chuand/new_gene/data/dm_sra
``` 
<br />
The first parameter /home/chuand/new_gene/data/dm_SraRunTable.txt is the metadata that downloaded by SRA Run selector. The second parameter is the path that the raw will be stored.

# The pipeline of TPM calculation and DEG identification based on GSE99574

## Acknowledgement
I am very appreciate to Chengchi Fang (Chinese Academy of Sciences) and Jianhai Chen (Sichuan University) for answering my questions and the friendly discussion about RNA-Seq sequencing technology.

## Information of softwares in my pipeline
<table class="MsoTableGrid" border="1" cellspacing="0" cellpadding="0" style="border-collapse:collapse;border:none;mso-border-top-alt:solid windowtext .5pt;
 mso-border-bottom-alt:solid windowtext .5pt;mso-yfti-tbllook:1184;mso-padding-alt:
 0cm 5.4pt 0cm 5.4pt;mso-border-insideh:none;mso-border-insidev:none">
 <tbody><tr style="mso-yfti-irow:0;mso-yfti-firstrow:yes">
  <td width="107" valign="top" style="width:80.05pt;border-top:solid windowtext 1.0pt;
  border-left:none;border-bottom:solid windowtext 1.0pt;border-right:none;
  mso-border-top-alt:solid windowtext .5pt;mso-border-bottom-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">Tool name<o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border-top:solid windowtext 1.0pt;
  border-left:none;border-bottom:solid windowtext 1.0pt;border-right:none;
  mso-border-top-alt:solid windowtext .5pt;mso-border-bottom-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">Version<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border-top:solid windowtext 1.0pt;
  border-left:none;border-bottom:solid windowtext 1.0pt;border-right:none;
  mso-border-top-alt:solid windowtext .5pt;mso-border-bottom-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">Availability<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:1">
  <td width="107" valign="top" style="width:80.05pt;border:none;mso-border-top-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK29"></a><a name="OLE_LINK28"><span style="mso-bookmark:OLE_LINK29"></span></a><span class="SpellE"><span style="mso-bookmark:OLE_LINK28"><span style="mso-bookmark:OLE_LINK29"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif;
  mso-font-kerning:0pt">fasterq</span></span></span></span><span style="mso-bookmark:OLE_LINK28"><span style="mso-bookmark:OLE_LINK29"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif;
  mso-font-kerning:0pt">-dump</span></span></span><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;mso-border-top-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">2.11.2<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;mso-border-top-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:2">
  <td width="107" valign="top" style="width:80.05pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK22"></a><a name="OLE_LINK21"><span style="mso-bookmark:OLE_LINK22"></span></a><span class="SpellE"><span style="mso-bookmark:OLE_LINK21"><span style="mso-bookmark:OLE_LINK22"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif">fastp</span></span></span></span><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif;
  mso-font-kerning:0pt"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">0.23.2<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">https://github.com/OpenGene/fastp<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:3">
  <td width="107" valign="top" style="width:80.05pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK25"></a><a name="OLE_LINK24"></a><a name="OLE_LINK23"><span style="mso-bookmark:OLE_LINK24"><span style="mso-bookmark:OLE_LINK25"><span lang="EN-US" style="mso-bidi-font-size:
  10.5pt;font-family:&quot;Times New Roman&quot;,serif">hisat2</span></span></span></a><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">2.2.1<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">http://daehwankimlab.github.io/hisat2<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:4">
  <td width="107" valign="top" style="width:80.05pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK33"></a><a name="OLE_LINK32"><span style="mso-bookmark:OLE_LINK33"></span></a><span class="SpellE"><span style="mso-bookmark:OLE_LINK32"><span style="mso-bookmark:OLE_LINK33"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif">samtools</span></span></span></span><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">1.14<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">http://www.htslib.org/download/<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:5">
  <td width="107" valign="top" style="width:80.05pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK35"></a><a name="OLE_LINK34"><span style="mso-bookmark:OLE_LINK35"></span></a><span class="SpellE"><span style="mso-bookmark:OLE_LINK34"><span style="mso-bookmark:OLE_LINK35"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif">featureCounts</span></span></span></span><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">2.0<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">https://sourceforge.net/projects/subread/files/subread-2.0.0<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:6;mso-yfti-lastrow:yes">
  <td width="107" valign="top" style="width:80.05pt;border:none;border-bottom:solid windowtext 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><a name="OLE_LINK36"></a><a name="OLE_LINK27"></a><a name="OLE_LINK26"><span style="mso-bookmark:OLE_LINK27"><span style="mso-bookmark:OLE_LINK36"><span lang="EN-US" style="mso-bidi-font-size:
  10.5pt;font-family:&quot;Times New Roman&quot;,serif">DESeq2</span></span></span></a><span lang="EN-US" style="mso-bidi-font-size:10.5pt;font-family:&quot;Times New Roman&quot;,serif"><o:p></o:p></span></p>
  </td>
  <td width="63" valign="top" style="width:47.6pt;border:none;border-bottom:solid windowtext 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">1.30.1<o:p></o:p></span></p>
  </td>
  <td width="384" valign="top" style="width:287.65pt;border:none;border-bottom:
  solid windowtext 1.0pt;mso-border-bottom-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt">
  <p class="MsoNormal"><span lang="EN-US" style="mso-bidi-font-size:10.5pt;
  font-family:&quot;Times New Roman&quot;,serif">https://bioconductor.org/packages/release/bioc/html/DESeq2.html<o:p></o:p></span></p>
  </td>
 </tr>
</tbody></table>

## 1. TPM calculation
Calculate the TPM using the wrapped script expression.py and the control file expression.ctl. All the needed parameters are provided in expression.ctl. If users provide “--download”, the raw data will be automatic download according to the configuration file
### 1.1 Usage
```
python expression.py expression.ctl
python expression.py expression.ctl --download
```
### 1.2 How to prepare the configuration file
All the data will be stored in the path specified by OUT parameter, for example fastq data were stored in "/home/chuand/new_gene/data/dm_gonad/fastq", and the other data follow the same rules for storing. The syntax of the configuration file follows the python syntax.
```
gtf = "/home/chuand/new_gene/data/BDGP6.22.97.gtf"
sra_run = "/home/chuand/new_gene/data/dm_SraRunTable.txt"
reference="/home/chuand/new_gene/data/ensembl_97/genome"
species="Drosophila melanogaster"
tissue = "gonad"
sra = "/home/chuand/new_gene/data/dm_sra"
fastq = "fastq"
clean_fastq ="clean_fastq"
mapping = "mapping"
bam = "bam"
merge_bam = "merge_bam"
OUT="/home/chuand/new_gene/data/dm_gonad"
countsummary="summary"
tpm_file = "/home/chuand/new_gene/data/dm_gonad/gonad_tpm"
```

## 2. The identification of DEGs

### 2.1 The download of raw data
Download the metadata about _D. melanogaster_ using SRA Run selector of NCBI. Than run the python script 1run_download.py. <br />
```
python 1run_download.py /home/chuand/new_gene/data/dm_SraRunTable.txt /home/chuand/new_gene/data/dm_sra
``` 
The first parameter /home/chuand/new_gene/data/dm_SraRunTable.txt is the metadata that downloaded by SRA Run selector. The second parameter is the path that the raw will be stored. <br />

### 2.2 Identification the sex-biased genes according to the parameters provided by the users <br />
```
nohup python3 deg.py "Drosophila melanogaster" gonad /home/chuand/new_gene/data/dm_sra /home/chuand/new_gene/data/dm_fastq /home/chuand/new_gene/data/dm_cleanfastq /home/chuand/new_gene/data/mapping /home/chuand/new_gene/data/bam /home/chuand/new_gene/data/mergebam  /home/chuand/new_gene/result/summary &
```
* The first parameter "Drosophila melanogaster" is the species that we want to analysis, the script will retrieve information from the metadata (For example /home/chuand/new_gene/data/dm_SraRunTable.txt) that is downloaded from NCBI. 
* The second parameter "gonad" is the focal tissue, in which the DEGs will be analyzed. the script will retrieve information from the metadata (For example /home/chuand/new_gene/data/dm_SraRunTable.txt) that is downloaded from NCBI. 
* The third parameter "/home/chuand/new_gene/data/dm_sra" is the path of raw data. 
* The fourth parameter "/home/chuand/new_gene/data/dm_fastq" is a path to store the data in fastq format. 
* The fifth parameter "/home/chuand/new_gene/data/dm_cleanfastq" is a path to store the clean data. 
* The sixth parameter "/home/chuand/new_gene/data/mapping" is a path to store the mapping data. 
* The seventh parameter "/home/chuand/new_gene/data/bam" is a path to store the mapping data with bam format. 
* The eighth path "/home/chuand/new_gene/data/mergebam" is a path to store the merged runs with one GSM ID. The ninth parameter "/home/chuand/new_gene/data/summary" is a file to store the read count.

### 2.3 The output of the script
* foldchange_gonad_w1118: The fold change of all genes in our analysis for tissue gonad of population w1118; 
* diff_gonad_w1118: The DEGs of all genes in our analysis for tissue gonad of population w1118; <br/>
* foldchange_gonad_orgR: The fold change of all genes in our analysis for tissue gonad of population orgR; <br/>
* diff_gonad_orgR: The DEGs of all genes in our analysis for tissue gonad of population orgR; <br/>

### 2.4 Fold change correlation of sex-biased genes based on testis and ovary
Figure a is my pipeline illustration for data analysis. b is sex-biased genes between w1118 and orgR. c is male-biased genes between w1118 and orgR. d is female-biased genes between w1118 and orgR. e, f and g are the correlation analysis for the overlapping parts of c, d and d, respectively. All the R values greater than 0.9 showing that the reliability of the sex-biased genes.<br/><br/>
 <img src="https://github.com/RiversDong/DEGs/blob/main/sexBiased.jpg" width=80% align="center">






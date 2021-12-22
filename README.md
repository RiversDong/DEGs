# The pipeline of DEG identification based on Oliver's data

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

## The download of raw data
Download the metadata about _D. melanogaster_ using SRA Run selector of NCBI. Than run the python script 1run_download.py. <br />
```
python 1run_download.py /home/chuand/new_gene/data/dm_SraRunTable.txt /home/chuand/new_gene/data/dm_sra
``` 
The first parameter /home/chuand/new_gene/data/dm_SraRunTable.txt is the metadata that downloaded by SRA Run selector. The second parameter is the path that the raw will be stored. <br />

## The identification of DEGs
Identification the sex-biased genes according to the parameters provided by the users <br />
```
nohup python3 deg.py "Drosophila melanogaster" gonad /home/chuand/new_gene/data/dm_sra /home/chuand/new_gene/data/dm_fastq /home/chuand/new_gene/data/dm_cleanfastq /home/chuand/new_gene/data/mapping /home/chuand/new_gene/data/bam /home/chuand/new_gene/data/mergebam  /home/chuand/new_gene/result/summary &
```
The first parameter "Drosophila melanogaster" is the species that we want to analysis, the script will retrieve information from the metadata (For example /home/chuand/new_gene/data/dm_SraRunTable.txt) that is downloaded from NCBI <br/>
The second parameter "gonad" is the focal tissue, in which the DEGs will be analyzed. the script will retrieve information from the metadata (For example /home/chuand/new_gene/data/dm_SraRunTable.txt) that is downloaded from NCBI.  <br/> 
The third parameter "/home/chuand/new_gene/data/dm_sra" is the path of raw data. <br/>
The fourth parameter "/home/chuand/new_gene/data/dm_fastq" is a path to store the data in fastq format. <br/>
The fifth parameter "/home/chuand/new_gene/data/dm_cleanfastq" is a path to store the clean data. <br/>
The sixth parameter "/home/chuand/new_gene/data/mapping" is a path to store the mapping data. <br/>
The seventh parameter "/home/chuand/new_gene/data/bam" is a path to store the mapping data with bam format <br/>
The eighth path "/home/chuand/new_gene/data/mergebam" is a path to store the merged runs with one GSM ID. <br/>
The ninth parameter "/home/chuand/new_gene/data/summary" is a file to store the read count. <br/>

## The output of the script
foldchange_gonad_w1118: The fold change of all genes in our analysis for tissue gonad of population w1118; <br/>
diff_gonad_w1118: The DEGs of all genes in our analysis for tissue gonad of population w1118; <br/>
foldchange_gonad_orgR: The fold change of all genes in our analysis for tissue gonad of population orgR; <br/>
diff_gonad_orgR: The DEGs of all genes in our analysis for tissue gonad of population orgR; <br/>

## Acknowledgement <br/>
I am very appreciate to Chengchi Fang in Institute of Hydrobiology, Chinese Academy of Sciences for answering my questions about RNA-Seq sequencing technology.





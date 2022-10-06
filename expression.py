import pandas as pd
import sys
import os
from multiprocessing import Pool
from download import *
from gtfparse import read_gtf

'''

'''

def pathCheck(OUT, sra, fastq, clean_fastq, mapping, bam, merge_bam, countsummary):
    fastq = os.path.join(OUT, "fastq")
    clean_fastq = os.path.join(OUT, "clean_fastq")
    mapping = os.path.join(OUT, "mapping")
    bam = os.path.join(OUT, "bam")
    merge_bam = os.path.join(OUT, "merge_bam")
    if not os.path.exists(fastq):
        os.mkdir(fastq)
    if not os.path.exists(clean_fastq):
        os.mkdir(clean_fastq)
    if not os.path.exists(mapping):
        os.mkdir(mapping)
    if not os.path.exists(bam):
        os.mkdir(bam)
    if not os.path.exists(merge_bam):
        os.mkdir(merge_bam)
    countsummary = os.path.join(OUT, countsummary)
    return sra, fastq, clean_fastq, mapping, bam, merge_bam, countsummary

def get_geo2run(species, tissue, sra_run):
    data = pd.read_csv(sra_run, sep = ",")
    targetSpTissue = data[(data["tissue"] == tissue) & (data["Organism"] == species)]
    row, column = targetSpTissue.shape
    geos = list(targetSpTissue["GEO_Accession (exp)"])
    runs = list(targetSpTissue["Run"])
    geo2run = {}
    for i in range(0, row):
        geo_Accession = geos[i]
        run = runs[i]
        if geo_Accession not in geo2run:
            geo2run[geo_Accession] = [run]
        else:
            geo2run[geo_Accession].append(run)
    return geo2run

def srr2fastq(species, tissue, sra, fastq, sra_run):
    data = pd.read_csv(sra_run, sep = ",")
    targetSpTissue = data[(data["tissue"] == tissue) & (data["Organism"] == species)]
    row, column = targetSpTissue.shape
    geos = list(targetSpTissue["GEO_Accession (exp)"])
    runs = list(targetSpTissue["Run"])
    geo2run = {}
    for i in range(0, row):
        geo_Accession = geos[i]
        run = runs[i]
        if geo_Accession not in geo2run:
            geo2run[geo_Accession] = [run]
        else:
            geo2run[geo_Accession].append(run)
    for i in geo2run:
        i_runs = geo2run[i]
        for j in i_runs:
            jfile = os.path.join(sra, j); jfile=os.path.join(jfile, j+".sra")
            cmd = "fastq-dump {0} --outdir {1} --split-files".format(jfile, fastq)
            os.system(cmd)

def fastqClean(fastq, clean_fastq, reportpath = "/home/chuand/new_gene/data/reports"):
    fastqs = os.listdir(fastq)
    for i in fastqs:
        iout = i.split("_")[0]+".fastq"
        ipath = os.path.join(fastq, i)
        ioutpath = os.path.join(clean_fastq, iout)
        ireport = os.path.join(reportpath, iout)
        cmd = "fastp --in1 {0} --out1 {1} -h {2}".format(ipath, ioutpath, ireport) #for single end data
        os.system(cmd)

def readmapping(clean_fastq, mapping, reference):
    fastqs = os.listdir(clean_fastq)
    cmds = []
    for i in fastqs:
        ipath = os.path.join(clean_fastq,i)
        iout = os.path.join(mapping, i.split(".")[0] + ".sam")
        cmd = "hisat2 -x {0} -U {1} -S {2}".format(reference, ipath, iout)
        cmds.append(cmd)
    pool = Pool(3)
    for i in cmds:
        pool.apply_async(os.system, args=(i,))
    pool.close()
    pool.join()

def sam2bam(mapping, bam):
    files = os.listdir(mapping)
    cmds = []
    for i in files:
        iin = os.path.join(mapping, i)
        iout = os.path.join(bam, i.split(".")[0]+".bam")
        cmd = "samtools view -bS {0} | samtools sort -@ 6 -o {1}".format(iin, iout)
        cmds.append(cmd)
    pool = Pool(6)
    for i in cmds:
        pool.apply_async(os.system, args=(i,))
    pool.close()
    pool.join()

def mergeBam(gsm2runs,bampath, mergebam):
    cmds = []
    for i in gsm2runs:
        runs = gsm2runs[i]
        tmp = []
        for j in runs:
            jpath = os.path.join(bampath, j+".bam")
            tmp.append(jpath)
        inbams = " ".join(tmp) 
        outbam = os.path.join(mergebam, i+".bam")
        cmd = "samtools merge -o {0} {1}".format(outbam, inbams)
        cmds.append(cmd)
    pool = Pool(6)
    for i in cmds:
        pool.apply_async(os.system, args=(i,))
    pool.close()
    pool.join()

def readcount(mergebam, countsummary, gtf):
    inbams = os.path.join(mergebam, "*.bam")
    # RNA-Seq测序中不包含内含子的信息
    cmd = "featureCounts -T 5 -a {0} -o {1} -t exon -g gene_id {2}".format(gtf, countsummary, inbams)
    os.system(cmd)

def builtmeta(tissue, species, sampleInfo = "/home/chuand/new_gene/data/GSE99574_All_samples_with_title.txt"):
    orgR = open("/home/chuand/new_gene/result/orgR.meta", "w")
    w1118 = open("/home/chuand/new_gene/result/w1118.meta", "w")
    orgR.write("sampleId\tgroup\ttitle\n")
    w1118.write("sampleId\tgroup\ttitle\n")
    sampleInfo = sampleInfo
    data = pd.read_csv(sampleInfo, sep="\t")
    targetSpTissue = data[(data["source name"] == tissue) & (data["organism"] == species)]
    row, column = targetSpTissue.shape
    titles = list(targetSpTissue["title"])
    geos = list(targetSpTissue["GEO accession number"])
    for i in range(0, row):
        i_title = titles[i]
        i_geo = geos[i].split(".")[0]
        if "orgR_" in i_title:
            if "_f_" in i_title:
                orgR.write(i_geo + "\t" + "female" + "\t" + i_title + "\n")
            else:
                orgR.write(i_geo + "\t" + "male" + "\t" + i_title + "\n")
        if "w1118_" in i_title:
            if "_f_" in i_title:
                w1118.write(i_geo + "\t" + "female" + "\t" + i_title + "\n")
            else:
                w1118.write(i_geo + "\t" + "male" + "\t" + i_title + "\n")
    w1118.close()
    orgR.close()
    outTissue = tissue.replace(" ","_")
    diff_w1118_meta = "/home/chuand/new_gene/result/w1118.meta"
    diff_w1118_out = "/home/chuand/new_gene/result/diff_"+outTissue+"_w1118"
    allfold_w1118_out = "/home/chuand/new_gene/result/foldchange_"+outTissue+"_w1118"
    diffcommand_w1118 = "{0} diffgene.r {1} {2} {3}".format(rcommand,\
            diff_w1118_meta,\
            diff_w1118_out,\
            allfold_w1118_out)

    diff_orgR_meta = "/home/chuand/new_gene/result/orgR.meta"
    diff_orgR_out = "/home/chuand/new_gene/result/diff_"+outTissue+"_orgR"
    allfold_orgR_out = "/home/chuand/new_gene/result/foldchange_"+outTissue+"_orgR"
    diffcommand_orgR = "{0} diffgene.r {1} {2} {3}".format(rcommand,\
            diff_orgR_meta, diff_orgR_out, allfold_orgR_out)
    os.system(diffcommand_orgR)
    os.system(diffcommand_w1118)

def count2tpm(countsummary, gtf, tpm_file):
    df = read_gtf(gtf); columns = ["seqname", "start", "end", "transcript_id", "gene_id"]
    exonInfo = df[df["feature"] == "exon"].loc[:,columns]
    genes = set(exonInfo["gene_id"]);gene2sumExonLen = {}
    for i in genes:
        tmp = exonInfo[exonInfo["gene_id"]==i]; row_num, col_num = tmp.shape
        exon_list = []
        for j in range(0, row_num):
            exon_list.extend(list(range(tmp.iloc[j]["start"], tmp.iloc[j]["end"]+1)))
        sumExonLen = len(set(exon_list))
        gene2sumExonLen[i] = sumExonLen
        
    read_data = pd.read_csv(countsummary, sep="\t", skiprows = 1); row_num, col_num = read_data.shape
    # calculate the tpm
    eff_gene_lens = []
    for i in range(0, row_num):
        gene = read_data.iloc[i]["Geneid"]
        eff_gene_len = gene2sumExonLen[gene]; eff_gene_lens.append(eff_gene_len)
    read_data["eff_gene_lens"] = eff_gene_lens
    count_pd = read_data.iloc[:,6:-1]; row_num, col_num = count_pd.shape
    count_pd_name = count_pd.columns.values
    tpm_pd=pd.DataFrame()
    tpm_pd["genes"]=read_data["Geneid"]
    for i in range(0, col_num):
        i_name=os.path.basename(count_pd_name[i])
        TPM = 10e6*((count_pd.iloc[:,i]/read_data["eff_gene_lens"])/sum(count_pd.iloc[:,i]/read_data["eff_gene_lens"]))
        TPM = TPM.round(3)
        tpm_pd[i_name] = TPM
    read_data.to_csv("del", sep="\t", index=False)
    tpm_pd.to_csv(tpm_file, sep="\t", index=False)

if __name__ == "__main__":
    argvs=sys.argv; ctrf=sys.argv[1]; ctr=open(ctrf).read(); exec(ctr);
    sra, fastq, clean_fastq, mapping, bam, merge_bam, countsummary = pathCheck(OUT, sra, fastq, clean_fastq, mapping, bam, merge_bam, countsummary)
    if "--download" in argvs:
        get_data(sra_run, sra)
    srr2fastq(species, tissue, sra, fastq,sra_run)
    fastqClean(fastq, clean_fastq)
    readmapping(clean_fastq, mapping, reference)
    sam2bam(mapping, bam)
    gsm2runs = get_geo2run(species, tissue, sra_run)
    mergeBam(gsm2runs, bam, merge_bam)
    readcount(merge_bam, countsummary, gtf)
    count2tpm(countsummary, gtf, tpm_file)
    #builtmeta(tissue, species)


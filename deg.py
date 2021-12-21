import pandas as pd
import sys
import os
from multiprocessing import Pool

species = sys.argv[1]
tissue = sys.argv[2]
srrpath = sys.argv[3]
fastqpath = sys.argv[4]
clearnfastqpath = sys.argv[5]
mappingpath = sys.argv[6]
bampath = sys.argv[7]
mergebam = sys.argv[8]
countsummary = sys.argv[9]
rcommand = "/home/chuand/.bin/miniconda3/envs/rbase/bin/Rscript"
#species = "Drosophila melanogaster"
#tissue = "abdomen without digestive or reproductive system"
#srrpath = "/home/chuand/new_gene/data/dm_sra"
#fastqpath = "/home/chuand/new_gene/data/dm_fastq"
def pathCheck():
    dm_fastq = "rm /home/chuand/new_gene/data/dm_fastq/*"
    dm_cleanfastq = "rm /home/chuand/new_gene/data/dm_cleanfastq/*"
    mapping = "rm /home/chuand/new_gene/data/mapping/*"
    bam = "rm /home/chuand/new_gene/data/bam/*"
    mergebam = "rm /home/chuand/new_gene/data/mergebam/*"
    os.system(dm_fastq)
    os.system(dm_cleanfastq)
    os.system(mapping)
    os.system(bam)
    os.system(mergebam)
    os.system("rm /home/chuand/new_gene/result/summary")
    os.system("rm /home/chuand/new_gene/result/summary.summary")
    os.system("rm /home/chuand/new_gene/result/orgR.meta")
    os.system("rm /home/chuand/new_gene/result/w1118.meta")

def get_geo2run(species, tissue, srrpath, fastqpath):
    data = pd.read_csv("/home/chuand/new_gene/data/dm_SraRunTable.txt", sep = ",")
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
def srr2fastq(species, tissue, srrpath, fastqpath):
    data = pd.read_csv("/home/chuand/new_gene/data/dm_SraRunTable.txt", sep = ",")
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

    srrs = os.listdir(srrpath)
    for i in geo2run:
        i_runs = geo2run[i]
        for j in i_runs:
            tmp = os.path.join(srrpath, j)
            jfile = os.path.join(tmp,j + ".sra")
            cmd = "fastq-dump {0} --outdir {1} --split-files".format(jfile, fastqpath)
            os.system(cmd)

def fastqClean(fastqpath, clearnfastqpath, reportpath = "/home/chuand/new_gene/data/reports"):
    fastq = os.listdir(fastqpath)
    for i in fastq:
        iout = i.split("_")[0]+".fastq"
        ipath = os.path.join(fastqpath, i)
        ioutpath = os.path.join(clearnfastqpath, iout)
        ireport = os.path.join(reportpath, iout)
        cmd = "fastp --in1 {0} --out1 {1} -h {2}".format(ipath, ioutpath, ireport) #for single end data
        os.system(cmd)
def readmapping(clearnfastqpath, mappingpath, reference = "/home/chuand/new_gene/data/ensembl_97/genome"):
    fastqs = os.listdir(clearnfastqpath)
    reference = reference
    cmds = []
    for i in fastqs:
        ipath = os.path.join(clearnfastqpath,i)
        iout = os.path.join(mappingpath, i.split(".")[0] + ".sam")
        cmd = "hisat2 -x {0} -U {1} -S {2}".format(reference, ipath, iout)
        cmds.append(cmd)
    pool = Pool(3)
    for i in cmds:
        pool.apply_async(os.system, args=(i,))
    pool.close()
    pool.join()
def sam2bam(mappingpath, bampath):
    files = os.listdir(mappingpath)
    cmds = []
    for i in files:
        iin = os.path.join(mappingpath, i)
        iout = os.path.join(bampath, i.split(".")[0]+".bam")
        cmd = "samtools view -bS {0} | samtools sort -@ 6 -o {1}".format(iin, iout)
        cmds.append(cmd)
    pool = Pool(6)
    for i in cmds:
        pool.apply_async(os.system, args=(i,))
    pool.close()
    pool.join()
def mergeBam(gsm2runs,bampath, mergebam):
    # samtools merge [options] -o out.bam [options] in1.bam ... inN.bam
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

def readcount(mergebam, countsummary, gtf="/home/chuand/new_gene/data/BDGP6.22.97.gtf"):
    gtf = gtf
    inbams = os.path.join(mergebam, "*.bam")
    # featureCounts -T 10 -a $gtf -o read.count -p -B -C -f -t exon -g gene_id  *.bam
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

if __name__ == "__main__":
    pathCheck()
    srr2fastq(species, tissue, srrpath, fastqpath) # transform sra to fastq format
    fastqClean(fastqpath, clearnfastqpath)
    readmapping(clearnfastqpath, mappingpath)
    sam2bam(mappingpath, bampath)
    gsm2runs = get_geo2run(species, tissue, srrpath, fastqpath)
    mergeBam(gsm2runs, bampath, mergebam)
    readcount(mergebam, countsummary, gtf="/home/chuand/new_gene/data/BDGP6.22.97.gtf")
    builtmeta(tissue, species)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ageData = pd.read_excel("F:\\New gene\\Supplementary Tables.xlsx", sheet_name = "S2")
w1118_diff = pd.read_csv("F:\\New gene\\DEGs\\diff_gonad_w1118", sep="\t")
w1118_diff_gene = set(w1118_diff.iloc[:,0])
orgR_diff = pd.read_csv("F:\\New gene\\DEGs\\diff_gonad_orgR", sep="\t")
orgR_diff_gene = set(orgR_diff.iloc[:,0])

w1118_all = pd.read_csv("F:\\New gene\\DEGs\\foldchange_gonad_w1118", sep="\t")
w1118_all_gene = w1118_all[w1118_all["log2FoldChange"].notna()].iloc[:,0]
w1118_all_gene = set(w1118_all_gene)

orgR_all = pd.read_csv("F:\\New gene\\DEGs\\foldchange_gonad_orgR", sep="\t")
orgR_all_gene = orgR_all[orgR_all["log2FoldChange"].notna()].iloc[:,0]
orgR_all_gene = set(orgR_all_gene)

def getGene(age, chromosome, population="w1118"):
    gene = ageData["g_id"][(ageData["branch"].isin(age)) & \
    (ageData["chromosome"].isin(chromosome)) & 
    (ageData["g_type"]=="protein_coding")]
    gene = set(gene)
    if population == "w1118":
        male_biased = w1118_diff.iloc[:,0][w1118_diff["log2FoldChange"]>1]
        male_biased = set(male_biased)
        male_biased = male_biased.intersection(gene)
        female_biased = w1118_diff.iloc[:,0][w1118_diff["log2FoldChange"]<-1]
        female_biased = set(female_biased)
        female_biased = female_biased.intersection(gene)
        unbiased = w1118_all_gene.difference(male_biased).difference(female_biased)
        unbiased = unbiased.intersection(gene)
        sum_all = sum([len(male_biased), len(female_biased), len(unbiased)])
        #print(population, age,chromosome, "male: ", len(male_biased)/sum_all)
        #print(population, age,chromosome, "female: ", len(female_biased)/sum_all)
        #print(population, age,chromosome, "unbiased: ", len(unbiased)/sum_all)
        #print("\n")
        ratio = [len(male_biased)/sum_all, len(female_biased)/sum_all, len(unbiased)/sum_all]
    if population == "orgR":
        male_biased = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]>1]
        male_biased = set(male_biased)
        male_biased = male_biased.intersection(gene)
        female_biased = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]<-1]
        female_biased = set(female_biased)
        female_biased = female_biased.intersection(gene)
        unbiased = orgR_all_gene.difference(male_biased).difference(female_biased)
        unbiased = unbiased.intersection(gene)
        sum_all = sum([len(male_biased), len(female_biased), len(unbiased)])
        #print(population, age,chromosome, "male: ", len(male_biased)/sum_all)
        #print(population, age,chromosome, "female: ", len(female_biased)/sum_all)
        #print(population, age,chromosome, "unbiased: ", len(unbiased)/sum_all)
        #print("\n")
        ratio = [len(male_biased)/sum_all, len(female_biased)/sum_all, len(unbiased)/sum_all]
    if population == "both":
        male_biased0 = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]>1]
        male_biased0 = set(male_biased0)
        female_biased0 = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]<-1]
        female_biased0 = set(female_biased0)
        unbiased0 = orgR_all_gene.difference(male_biased0).difference(female_biased0)

        male_biased1 = w1118_diff.iloc[:,0][w1118_diff["log2FoldChange"]>1]
        male_biased1 = set(male_biased1)
        female_biased1 = w1118_diff.iloc[:,0][w1118_diff["log2FoldChange"]<-1]
        female_biased1 = set(female_biased1)
        unbiased1 = w1118_all_gene.difference(male_biased1).difference(female_biased1)
        
        male_biased = male_biased0.intersection(male_biased1)
        male_biased=male_biased.intersection(gene)
        female_biased = female_biased0.intersection(female_biased1)
        female_biased = female_biased.intersection(gene)
        unbiased = unbiased0.intersection(unbiased1)
        unbiased = unbiased.intersection(gene)
        sum_all = sum([len(male_biased), len(female_biased), len(unbiased)])
        #print(population, age,chromosome, "male: ", len(male_biased)/sum_all)
        #print(population, age,chromosome, "female: ", len(female_biased)/sum_all)
        #print(population, age,chromosome, "unbiased: ", len(unbiased)/sum_all)
        #print("\n")
        ratio = [len(male_biased)/sum_all, len(female_biased)/sum_all, len(unbiased)/sum_all]
    return ratio

def draw(dataList):
    custom_lines = [Line2D([0], [0], color="#F5C767", lw=2),
                    Line2D([0], [0], color="#4198D7", lw=2),
                    Line2D([0], [0], color="#39A767", lw=2)]
    list_size = len(dataList)
    width = 0.33
    labels = ["Autosomes", "X chromosome"]
    for i in range(0, list_size):
        ax = fig.add_subplot(3,2,i+1)
        old = dataList[i][0]
        young = dataList[i][1]
        plt.bar([0, width, 2*width], old, width, color=["#F5C767", "#4198D7", "#39A767"])
        plt.bar([1.5, 1.5+width, 1.5+2*width], young, width, color=["#F5C767", "#4198D7", "#39A767"])
        plt.ylabel("Proportion of genes")
        plt.xticks([width, 1.5+width], labels)
        y_ticks = np.arange(0,1.25,0.25)
        plt.ylim([0,1])
        plt.yticks(y_ticks)
        if i==0:
            plt.legend(custom_lines, ['Male', 'Female', 'Unbiased'], loc="best")


if __name__ == "__main__":
    # w1118
    a_list = []
    a1 = getGene([-2, -1, 0], ["2L", "2R", "3L", "3R", 4], population="w1118")
    a2 = getGene([-2, -1, 0], ["X"], population="w1118")
    a_list.append(a1)
    a_list.append(a2)

    b_list = []
    b1 = getGene([1,2,3,4,5,6], ["2L", "2R", "3L", "3R", 4], population="w1118")
    b2 = getGene([1,2,3,4,5,6], ["X"], population="w1118")
    b_list.append(b1)
    b_list.append(b2)

    # orgR young
    c_list = []
    c1=getGene([-2, -1, 0], ["2L", "2R", "3L", "3R", 4], population="orgR")
    c2=getGene([-2, -1, 0], ["X"], population="orgR")
    c_list.append(c1)
    c_list.append(c2)

    d_list = []
    d1=getGene([1,2,3,4,5,6], ["2L", "2R", "3L", "3R", 4], population="orgR")
    d2=getGene([1,2,3,4,5,6], ["X"], population="orgR")
    d_list.append(d1)
    d_list.append(d2)
    
    # both
    e_list = []
    e1=getGene([-2, -1, 0], ["2L", "2R", "3L", "3R", 4], population="both")
    e2=getGene([-2, -1, 0], ["X"], population="both")
    e_list.append(e1)
    e_list.append(e2)
    
    f_list = []
    f1 = getGene([1,2,3,4,5,6], ["2L", "2R", "3L", "3R", 4], population="both")
    f2 = getGene([1,2,3,4,5,6], ["X"], population="both")
    f_list.append(f1)
    f_list.append(f2)
    fig = plt.figure()
    draw([a_list, b_list, c_list, d_list, e_list, f_list])
    fig.tight_layout()
    fig.savefig("F:\\New gene\\check_inMainTextSexBiased.pdf", dpi=600, transparent=True)


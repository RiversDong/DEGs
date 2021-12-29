import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit # used to fit the curve

def exponential_decay(x, a, b, d):
    return a*((1-d)*np.exp(b * x) + d)

parameters = {'axes.labelsize':5, 'ytick.labelsize':5, 'xtick.labelsize':5}
          
plt.rcParams.update(parameters)

def getChrPercentage(agefile="F:\\New gene\\Supplementary Tables.xlsx"):
    branch = [[-2], [-1], [0], [1], [2], [3], [4], [5,6]]
    age_data = pd.read_excel(agefile, sheet_name="S2")
    x_percentage = []
    autoChr_percentage = []
    for i in branch:
        i_gene_X = age_data[(age_data["branch"].isin(i)) &\
            (age_data["g_type"]=="protein_coding") &\
            (age_data["chromosome"].isin(["X"]))]["g_id"]
        #i_gene_X = age_data[(age_data["branch"].isin(i)) & (age_data["chromosome"].isin(["X"]))]["g_id"]
        i_gene_X = set(i_gene_X)
        i_gene_autoChrs = age_data[(age_data["branch"].isin(i)) &\
            (age_data["g_type"]=="protein_coding") &\
            (age_data["chromosome"].isin(["2L", "2R", "3L", "3R", 4]))]["g_id"]
        #i_gene_autoChrs = age_data[(age_data["branch"].isin(i)) & (age_data["chromosome"].isin(["2L", "2R", "3L", "3R", 4]))]["g_id"]
        i_gene_autoChrs = set(i_gene_autoChrs)
        sumXautoChr = len(i_gene_X) + len(i_gene_autoChrs)
        i_x_percentage = len(i_gene_X)/sumXautoChr
        i_autoChr_percentage = len(i_gene_autoChrs)/sumXautoChr
        x_percentage.append(i_x_percentage)
        autoChr_percentage.append(i_autoChr_percentage)
    return x_percentage, autoChr_percentage

ageData = pd.read_excel("F:\\New gene\\Supplementary Tables.xlsx", sheet_name = "S2")
w1118_diff = pd.read_csv("F:\\New gene\\DEGs\\diff_whole_body_w1118", sep="\t")
w1118_diff_gene = set(w1118_diff.iloc[:,0])
orgR_diff = pd.read_csv("F:\\New gene\\DEGs\\diff_whole_body_orgR", sep="\t")
orgR_diff_gene = set(orgR_diff.iloc[:,0])

w1118_all = pd.read_csv("F:\\New gene\\DEGs\\foldchange_whole_body_w1118", sep="\t")
w1118_all_gene = w1118_all[w1118_all["log2FoldChange"].notna()].iloc[:,0]
w1118_all_gene = set(w1118_all_gene)

orgR_all = pd.read_csv("F:\\New gene\\DEGs\\foldchange_whole_body_orgR", sep="\t")
orgR_all_gene = orgR_all[orgR_all["log2FoldChange"].notna()].iloc[:,0]
orgR_all_gene = set(orgR_all_gene)

def get_biase(age, chromosome, population="w1118"):
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
    if population == "orgR":
        male_biased = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]>1]
        male_biased = set(male_biased)
        male_biased = male_biased.intersection(gene)
        female_biased = orgR_diff.iloc[:,0][orgR_diff["log2FoldChange"]<-1]
        female_biased = set(female_biased)
        female_biased = female_biased.intersection(gene)
        unbiased = orgR_all_gene.difference(male_biased).difference(female_biased)
        unbiased = unbiased.intersection(gene)
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
    return male_biased, female_biased, unbiased

# xs,ys are the x and y values we pass in, ratio refers to our share (for example, 70% of men and 30% of women), sizes refer to the size of this point, and ax is the drawing function
def drawPieMarker(xs, ys, ratios, sizes, colors, ax):
    markers = []
    previous = 0

    for color, ratio in zip(colors, ratios):
        this = 2 * np.pi * ratio + previous
        x  = [0] + np.cos(np.linspace(previous, this, 30)).tolist() + [0]
        y  = [0] + np.sin(np.linspace(previous, this, 30)).tolist() + [0]
        xy = np.column_stack([x, y])
        previous = this
        markers.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes), 'facecolor':color, 'edgecolor':'black'})
    # scatter each of the pie pieces to create piesggg
    for marker in markers:
        ax.scatter(xs, ys, **marker, alpha=0.7)
 

if __name__ == "__main__":
    custom_lines = [Line2D([0], [0], color="#F5C767", lw=2),
                    Line2D([0], [0], color="#4198D7", lw=2),
                    Line2D([0], [0], color="#39A767", lw=2)]

    x_percentage, autoChr_percentage = getChrPercentage()
    branch_list = [[-2], [-1], [0], [1], [2], [3], [4], [5, 6]]
    chromosome_list = [["X"], ["2L", "2R", "3L", "3R", 4]]
    #chromosome_list = [["X"]]
    dataset = ["w1118", "orgR", "both"]
    fig = plt.figure()
    axindex = 1
    xpercentage = 0.157
    x_pos = [-75, -68, -64, -62.5, -58.5, -49.5, -28.5, -6.5]
    x_pos_fit = [-126, -65, -64, -62.5, -58.5, -49.5, -28.5, -6.5]
    myr = range(-75, 0, 5)
    for population in dataset:
        ax = fig.add_subplot(3,1,axindex)
        y_pos_x = x_percentage
        y_pos_autoChr = autoChr_percentage
        xy_index = 0
        onX = []
        onAutoChr = []
        maleBiased_x_percentage_list = []
        maleBiased_autoChr_percentage_list = []
        for branch in branch_list:
            male_biased_x, female_biased_x, unbiased_x = get_biase(branch, ["X"], population=population)
            male_biased_autoChr, female_biased_autoChr, unbiased_autoChr = get_biase(branch, \
                                                       ["2L", "2R", "3L", "3R", 4], population=population)
            
            non_biased_x = len(female_biased_x) + len(unbiased_x)
            non_biased_Chr = len(female_biased_autoChr) + len(unbiased_autoChr)
            non_biased_x_percentage = non_biased_x/(non_biased_x+non_biased_Chr)
            non_biased_x_percentage = round(non_biased_x_percentage, 3)
            non_biased_autoChr_percentage = 1 - non_biased_x_percentage
            onX.append(non_biased_x_percentage)
            onAutoChr.append(non_biased_autoChr_percentage)
            
            maleBiased_x_percentage = len(male_biased_x)/(len(male_biased_x) + len(male_biased_autoChr))
            maleBiased_x_percentage = round(maleBiased_x_percentage, 3)
            maleBiased_x_percentage_list.append(maleBiased_x_percentage)
            maleBiased_autoChr_percentage = 1 - maleBiased_x_percentage
            maleBiased_autoChr_percentage_list.append(maleBiased_autoChr_percentage)

            drawPieMarker([x_pos[xy_index]], [y_pos_x[xy_index]],\
                [maleBiased_x_percentage, 1-maleBiased_x_percentage], \
                [300], ['red', 'white'], ax)

            drawPieMarker([x_pos[xy_index]], [y_pos_autoChr[xy_index]],\
                [maleBiased_autoChr_percentage, 1-maleBiased_autoChr_percentage], \
                [300], ['blue', 'white'], ax)
                
            xy_index = xy_index + 1
        plt.axhline(y=xpercentage, color="grey", linestyle="--")
        ax.plot(x_pos, onX, "r.--", linewidth=1, alpha=0.7)
        ax.plot(x_pos, onAutoChr, "b.--", linewidth=1, alpha=0.7)
        plt.ylim([0,1])
        plt.ylabel("Proportion of genes")
        #plt.xlabel("Million years (branch)")
        y_ticks = np.arange(0,1.25,0.25)
        plt.yticks(y_ticks)
        axindex = axindex + 1
        plt.xticks(myr, myr, color = "green")
        popt, pcov = curve_fit(exponential_decay, x_pos_fit, maleBiased_x_percentage_list)
        NN = popt[0]
        rr = popt[1]
        dd = popt[2]
        print(NN, rr, dd)
        curve_x = np.linspace(-75, -3, 200)
        curve_y = exponential_decay(curve_x, NN, rr, dd)
        curve_y_auto = 1-curve_y
        ax.plot(curve_x, curve_y, "r", linewidth=1)
        ax.plot(curve_x, curve_y_auto, "b", linewidth=1)

    fig.tight_layout()
    fig.savefig("F:\\New gene\\newFigure4.pdf", dpi=600, transparent=True)



import os
import glob
import scipy.stats as stats
from scipy.stats import norm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import logging
from collections import defaultdict
from mageck_nest.gsea import *
import operator
import itertools
from mageck_nest.mleinstanceio import *
from mageck_nest.mageck_nest_PPI import *

def nest_output(self,mark):
    os.chdir(self.output_directory)
    gene_fdr_correction(self.allgenedict,self.adjust_method,self.negative_control,wald_only=True)
    genefile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.gene_summary.txt'
    sgrnafile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.sgrna_summary.txt'

    logging.info('Writing gene results to '+genefile)
    logging.info('Writing sgRNA results to '+sgrnafile)
    write_gene_to_file(self.allgenedict,genefile,betalabels=self.beta_labels,wald_only=True)
    write_sgrna_to_file(self.include_samples,self.allgenedict,self.invalid_gRNA_dict,sgrnafile,self.negative_control)

    gene_beta_values_list=betascore_plot(self.QC_path,genefile,negative_control=self.negative_control,wald_only=True)

    if self.QC_metric==True:
        PPI_diagnosis_output=open("%s/%s_QC_report.txt" %(self.QC_path,self.output_prefix),'a')
        return_NES=gsea_main(self.QC_path,genefile,gene_beta_values_list)

        PPI_diagnosis_insert=PPI_main(self,for_QC_output_only=True)
        PPI_diagnosis_output.write("mode: PPI-%s & outliers_removal-%s\n\n" %(mark[0],mark[1]))
        correlation_output=["rsquare:",PPI_diagnosis_insert[0][0][0],"; p_value:",PPI_diagnosis_insert[0][0][1]]
        regresssion_output=["slope:",PPI_diagnosis_insert[0][1][0],"; intercept:",PPI_diagnosis_insert[0][1][1]]
        PPI_diagnosis_output.write("\tcorrelation: %s\n" %"".join([str(i) for i in correlation_output]))
        PPI_diagnosis_output.write("\tregression: %s\n\n" %"".join([str(i) for i in regresssion_output]))

        for key,NES in return_NES.items():
            PPI_diagnosis_output.write("\tBeta: %s\n" %key)
            PPI_diagnosis_output.write("\tGSEA Normalized enrichment score: %s\n" %str(NES[0]))
            PPI_diagnosis_output.write("\tGSEA P-value: %s\n" %str(NES[1]))
            PPI_diagnosis_output.write("\tGSEA FDR:%s\n\n" %str(NES[2]))

def moffa_gene_set_file_preparation(gene_list):
    positive_control_genes_temp=['ARCN1', 'EIF3B', 'EIF3D', 'EIF5B', 'KPNB1', 'POLR2D', 'PRPF19',
    'PSMB2', 'PSMD1', 'RPL10A', 'RPL11', 'RPL23', 'RPL27', 'RPL32', 'RPL36', 'RPL4', 'RPL7', 'RPL9',
    'RPS11', 'RPS13', 'RPS14', 'RPS15A', 'RPS17', 'RPS26', 'RPSA', 'SNRNP200', 'SNRPD2', 'TUBA1B',
    'TUBA1C', 'U2AF1', 'COPB1', 'EEF2', 'EFTUD2', 'EIF3A', 'PSMA1', 'PSMA3', 'RAN', 'RPL12',
    'RPL18', 'RPL37', 'RPL37A', 'RPL5', 'RPL6', 'RPS18', 'RPS27A', 'RPS7', 'RPS8', 'RPS9', 'SF3B5',
    'SHFM1', 'CCT7', 'CDC5L', 'EIF4A3', 'HNRNPC', 'INTS4L2', 'KARS', 'NHP2L1', 'NUF2', 'PRPF31',
    'PSMB3', 'RPL14', 'RPL34', 'RPLP0', 'RPS3', 'RRM1', 'RUVBL2', 'SNRPD1', 'USP39', 'CCT3', 'CCT8',
    'EIF3G', 'NUP93', 'RPL7A', 'SF3A1', 'SFPQ', 'AQR', 'CLTC', 'EIF2S2', 'EIF3I', 'LSM6', 'PHB',
    'POLR2I', 'PSMD7', 'RPL19', 'RPS28', 'U2AF2', 'EIF3C', 'LSM4', 'RPL31', 'VCP', 'FTSJ3', 'NACA',
    'PSMA6', 'PSMD11', 'RILPL2', 'RPA2', 'RPL23A', 'RPL3', 'SUPT5H', 'COPZ1', 'FAM126A', 'PHF5A',
    'RPL24', 'RPLP2', 'RPS3A', 'SNRPE', 'TIMM10', 'CCT4', 'CDK11B', 'HAUS7', 'LSM3', 'LSM5', 'NXF1',
    'PABPN1', 'PSMA2', 'RPSAP9', 'COPS2', 'DYNC1H1', 'HSPA9', 'RPL13', 'RPL18A', 'RPS6', 'SF3B2',
    'SNRNP27', 'EIF3F', 'POLR2A', 'RPL10', 'SF3B3', 'SRSF1', 'TUBGCP2', 'ALYREF', 'RBM17', 'RPA1',
    'RPL13AP6', 'RPS24', 'RPS27', 'COPA', 'NUP133', 'PAPOLA', 'RFC4', 'RPS19', 'SF3B1', 'SMC3',
    'CCNB3', 'CDC40', 'CHD4', 'FUT9', 'HNRNPU', 'NAPA', 'PRPF3', 'RPS15', 'SDAD1', 'CSE1L', 'DDX18',
    'HSPE1', 'KCNJ5', 'NUP205', 'PRKAB2', 'PSMD6', 'RPL35A', 'SRCAP', 'UBA1', 'XPO1', 'EIF1AX', 'ETF1',
    'HINT3', 'HSPG2', 'INTS9', 'NOC4L', 'PELP1', 'PRPF38A', 'RPL17', 'RPL35', 'SKAP1', 'SMC2', 'SNRPG',
    'WDR12', 'ZNF207', 'CHMP2A', 'DNM2', 'GPR112', 'MS4A13', 'PHB2', 'POLA1', 'PRPF8', 'PSMC2', 'RPL13A',
    'RPS29', 'YY1', 'DDB1', 'EIF2B4', 'GNL3', 'MCM2', 'PRPF18', 'PSMC4', 'RPS20', 'WDR60', 'ADCYAP1R1',
    'COPS4', 'COPS6', 'CYP1A2', 'GAR1', 'GTDC2', 'HNRNPK', 'LIAS', 'MRPS31', 'NEDD8', 'PGLYRP4', 'PPP2R1A',
    'PRUNE', 'PSMC1', 'WDR61', 'XAB2', 'DYNC1I2', 'GUSBP4', 'RPS5', 'SF3B4', 'SRSF3', 'TLR9', 'TSTA3', 'BST2',
    'DAPK1', 'DDX46', 'DUOXA2', 'EIF4E3', 'POLR2F', 'RUVBL1', 'SPINT2', 'TRIP13', 'ZC3H13', 'BRIX1', 'CNOT3',
    'CPSF1', 'DMKN', 'FLT3LG', 'GTF3C4', 'NUP98', 'PTGDR2', 'RPL38', 'RPN2', 'SF3A2', 'SLC22A3', 'TFIP11',
    'THOP1', 'VPS28', 'ZAN', 'ALG10B', 'ERCC6L', 'HNRNPM', 'IGFBP2', 'NCOR2', 'PFDN2', 'QARS', 'RPL26',
    'RPL30', 'SRFBP1', 'SUPT6H', 'USP17L4', 'CCT6A', 'DDX49', 'EPHB4', 'EXOSC10', 'ISG20', 'MPL', 'NAPG',
    'NUDT21', 'PLK1S1', 'RPS4X', 'TEKT2', 'TEX11', 'TXLNG2P', 'WDR86', 'AURKB', 'C12orf66', 'CCNB2', 'CDK17',
    'COPS8', 'CTDNEP1', 'DDX51', 'ERH', 'EVX2', 'GREM1', 'GRIN1', 'HEATR1', 'MSANTD3', 'NID2', 'NUP54',
    'PAQR5', 'RBM47', 'SUPV3L1', 'ZBTB48']

    geneset_file=open("GSEA_analysis_on_essential_genes.gmt",'wt')
    positive_control_genes=[i.upper() for i in gene_list if i in positive_control_genes_temp]
    insert=["GSEA_analysis_on_essential_genes","GSEA_analysis_on_essential_genes"]+positive_control_genes
    geneset_file.write("%s\n" %"\t".join(insert))

def rnk_file_preparation(gene_beta_list,output_name,reverse=False):
    gene_beta_list.sort(key=operator.itemgetter(1),reverse=reverse)
    output=open("%s.rnk" %output_name,'wt')
    for i in gene_beta_list:
        try:
            float(i[1])
            if i[0]!="NA":
                output.write("%s\n" %"\t".join([str(k) for k in i]))
        except ValueError:
            print("wrong")

def gsea_main(target_path,file_name,gene_beta_values_list):
    return_NES=defaultdict(list)
    gene_list=[i[0] for i in list(itertools.chain.from_iterable(gene_beta_values_list.values()))]
    moffa_gene_set_file_preparation(gene_list)

    for short_column_name,gene_beta_list in gene_beta_values_list.items():
        output_name="%s_%s" %(file_name[:file_name.index(".txt")],short_column_name)
        rnk_file_preparation(gene_beta_list,output_name,reverse=False)
        ranked_file="%s.rnk" %output_name
        NES, pvalue,FDR=prerank(rnk="%s" %(ranked_file), gene_sets="GSEA_analysis_on_essential_genes.gmt",outdir=target_path)
        return_NES[short_column_name]=[NES, pvalue,FDR]

    remove_types=('*.gmt', '*.rnk')
    for remove_type in remove_types:
        for file in glob.glob(remove_type):
            os.remove(file)

    return return_NES

def betascore_plot(target_path,file_name,negative_control=None,wald_only=False):
    data=open(file_name,'rb')
    short_file_name=file_name[:file_name.index(".gene_summary.txt")]
    title=data.readline().decode().strip().split("\t")
    beta_value_columns=[title.index(i) for i in title if "|beta" in i]
    beta_value_names=[i[:i.index("|")] for i in [title[k] for k in beta_value_columns]]
    z_columns=[title.index(i) for i in title if "|z" in i]
    z_names=[i[:i.index("|")] for i in [title[k] for k in z_columns]]

    z_value_list={i:[] for i in z_names}
    beta_value_list={i:[] for i in beta_value_names}
    gene_beta_values_list={i:[] for i in beta_value_names}

    for line in data:
        elements=line.decode().strip().split("\t")
        if elements[0] not in negative_control:
            for i in range(len(z_names)):
                beta_value_list[beta_value_names[i]].append(float(elements[beta_value_columns[i]]))
                gene_beta_values_list[beta_value_names[i]].append([elements[0],float(elements[beta_value_columns[i]])])
                z_value_list[z_names[i]].append(float(elements[z_columns[i]]))

    # plot beta scores distribution
    for key,values in beta_value_list.items():
        beta_value_list[key]=[x for x in values if str(x) != 'nan' and abs(x)<10]
        fig, ax = plt.subplots()
        ax.set_xticks([0], minor=True)
        ax.xaxis.grid(True, which='minor',color='g',linestyle="-",linewidth=1)
        if len(beta_value_list[key])>3000:
            ax.hist(beta_value_list[key],bins=1000)
        else:
            ax.hist(beta_value_list[key],bins=100)

        plt.xlabel("Beta scores")
        plt.savefig("Hist of %s beta value %s .pdf" %(key,short_file_name))
        plt.close()
        os.rename("Hist of %s beta value %s .pdf" %(key,short_file_name), "/%s/Hist of %s beta value %s .pdf" %(target_path,key,short_file_name))

    # plot QQ-plot of z-scores
    for key,values in z_value_list.items():
        z_value_list[key]=[x for x in values if str(x) != 'nan' and abs(x)<10]
        fig=sm.qqplot(np.array(z_value_list[key]),stats.norm,fit=True, line='45')
        plt.savefig("QQplot of %s z score %s.pdf" %(key,short_file_name))
        plt.close()
        os.rename("QQplot of %s z score %s.pdf" %(key,short_file_name),"/%s/QQplot of %s z score %s.pdf" %(target_path,key,short_file_name))

    return gene_beta_values_list

from __future__ import print_function
import os
import sys
import math
import numpy as np
import pickle
from scipy import stats
import matplotlib.pyplot as plt
import operator
from collections import Counter,defaultdict
import logging
from scipy.stats import norm



def quantile_matching(value_list):
    mean_value=np.mean(value_list)
    normalized_value_list=[(i-mean_value) for i in value_list]
    abs_value=[abs(k) for k in normalized_value_list]
    abs_value.sort()
    var_value=float(np.percentile(abs_value,95))/norm.ppf(0.975)

    return((var_value))

def beta_non_PPI_prior_calculation(allgenedict,negative_control):
    '''
    For zero-centered beta prior, this function calculates the variance
    '''
    temp_beta=[(sk.nb_count.shape[1],sk.beta_estimate) for gene,sk in list(allgenedict.items()) if gene not in negative_control]
    temp_beta1=[i[1][i[0]:] for i in temp_beta]
    array_beta1=np.asarray(temp_beta1)
    beta1_prior_variance=[]
    for column in range(array_beta1.shape[1]):
        beta1=array_beta1[:,column]
        var_beta1=quantile_matching(beta1)
        beta1_prior_variance.append((var_beta1*2)**2)
    return(beta1_prior_variance)

def PPI_import_string_9(self):
    PPI_filename=os.path.join(os.path.dirname(__file__),'string_v9_network.p')
    self.PPI_data=pickle.load(open(PPI_filename,'rb'))

def PPI_weighting_rewiring(self):
    weighting_record=[]
    for gene, neighboring_gene_weighting in self.PPI_data.items():
        weighting_record+=list(neighboring_gene_weighting.values())
    weighting_correction_ratio=np.sum(weighting_record)/len(weighting_record)

    for gene, neighboring_gene_weighting in self.PPI_data.items():
        for neighboring_gene, weighting in neighboring_gene_weighting.items():
            neighboring_gene_weighting[neighboring_gene]=weighting/weighting_correction_ratio

def PPI_coverge_diagnosis(self):
    PPI_genes=self.PPI_data.keys()
    crispr_genes=self.allgenedict.keys()
    PPI_genes=[i.upper() for i in PPI_genes]
    crispr_genes=[i.upper() for i in crispr_genes]
    logging.info("Number of PPI genes: %s" %str(len(PPI_genes)))
    logging.info("Number of CRISPR genes: %s" %str(len(crispr_genes)))
    logging.info("Number of overlapped genes: %s" %str(len([i for i in PPI_genes if i in crispr_genes])))
    if len([i for i in PPI_genes if i in crispr_genes])<5000:
        logging.info("")
        return False
    else:
        return True

def PPI_network_plotting(self,for_QC_output_only=False,constant=0):
    os.chdir(self.file_directory)
    beta1_estimate=defaultdict(list)
    beta1_se_mat=defaultdict(list)
    for (tgid,tginst) in self.allgenedict.items():
        if tgid not in self.negative_control:
            tgid=tgid.upper()
            beta1_estimate[tgid]=tginst.beta_estimate.tolist()[tginst.nb_count.shape[1]:]
            beta1_se_mat[tgid]=tginst.beta_se_val

            for i in range(len(beta1_estimate[tgid])):
                if abs(beta1_estimate[tgid][i])>10:
                    beta1_estimate[tgid][i]=0

    true_central_beta=[]
    predicted_central_beta=[]
    gene_set=set(list(beta1_estimate.keys()))
    outliers_record=[]
    for gene,beta_value in beta1_estimate.items():
        #sgRNA_number=Counter(self.allgenedict[gene].eff_estimate)[1]
        if Counter([math.isnan(i) for i in beta_value])[True]==0:
            sgRNA_number=self.allgenedict[gene].nb_count.shape[1]
            original_surrounding_beta=[] # beta value of surrounding genes
            original_surrounding_beta_variance=[] # variance of beta of surrounding genes
            weighted_weighting=[] # the weighting given by the network data base
            for neighboring_gene, weighting in self.PPI_data[gene.upper()].items():
                if neighboring_gene in gene_set and Counter([math.isnan(i) for i in beta1_estimate[neighboring_gene]])[True]==0:
                    original_surrounding_beta.append(beta1_estimate[neighboring_gene])
                    weighted_weighting.append(weighting**1.5)
                    original_surrounding_beta_variance.append([1.0/(i**2) for i in beta1_se_mat[neighboring_gene]])

            if len(weighted_weighting)>0:
                original_surrounding_beta=np.matrix(np.vstack(original_surrounding_beta))
                original_surrounding_beta_variance=np.matrix(np.vstack(original_surrounding_beta_variance))
                weighted_weighting=np.matrix(np.vstack(weighted_weighting)).T
                weighted_weighting_sum=np.sum(weighted_weighting,axis=1)
                weighted_average=((weighted_weighting*original_surrounding_beta)/(weighted_weighting_sum+constant))
                weighted_average=weighted_average.tolist()[0]

                if for_QC_output_only==False:
                    beta_prior=[]
                    prior_SD=[i**(0.5) for i in self.non_PPI_beta_prior_variance]
                    for i in range(len(beta_value)):
                        predicted_beta=[k[i] for k in original_surrounding_beta]
                        predicted_beta_SD=quantile_matching(predicted_beta)
                        temp=((prior_SD[i]**2)*weighted_average[i])/((prior_SD[i])**2+(predicted_beta_SD)**2)
                        beta_prior.append(temp)
                    self.allgenedict[gene].prior_mean=beta_prior
                    self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance

                if Counter([math.isnan(i) for i in beta_value])[True]==0:
                    true_central_beta.append(beta_value)
                    predicted_central_beta.append(weighted_average)
            else:
                outliers_record.append(gene)
                if for_QC_output_only==False:
                    self.allgenedict[gene].prior_mean=[0]*len(beta_value)
                    self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance

        else:
            outliers_record.append(gene)
            if self_adjust==True:
                self.allgenedict[gene].prior_mean=[0]*len(beta_value)
                self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance

    true_central_beta=np.vstack(true_central_beta)
    predicted_central_beta=np.vstack(predicted_central_beta)

    if true_central_beta.shape!=predicted_central_beta.shape:
        sys.exit(1)

    if for_QC_output_only==True:
        beta_correlation_record=[]
        for k in range(true_central_beta.shape[1]):
            pearson_value=stats.pearsonr(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist())
            regression_value=np.polyfit(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist(),1)
            fit_fn = np.poly1d(regression_value)
            beta_correlation_record.append([pearson_value,regression_value])

            plt.scatter(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist(),color="r",s=0.1)
            plt.plot(predicted_central_beta[:,k].tolist(), fit_fn(predicted_central_beta[:,k].tolist()), '--k')
            plt.xlabel("Weighted average of neighboring genes")
            plt.ylabel("Center gene")
            file_name="{}_PPI_fitting_weighting_{}_beta1_{}_constant_{}.pdf".format(self.output_prefix,str(self.PPI_weighting),str(k),str(constant))
            plt.savefig(file_name)
            plt.close()
            os.rename("%s/%s" %(self.file_directory,file_name),"%s/%s" %(self.QC_path,file_name))

        return beta_correlation_record


def PPI_main(self,for_QC_output_only=False):
    constant=3
    if for_QC_output_only==True:
        beta_correlation_record=PPI_network_plotting(self,for_QC_output_only=True,constant=constant)
        return beta_correlation_record
    else:
        self.PPI_diagnosis=PPI_coverge_diagnosis(self)
        if self.PPI_diagnosis==True:
            PPI_network_plotting(self,for_QC_output_only=False,constant=constant)

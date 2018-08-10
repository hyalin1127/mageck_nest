import os
import pickle
import numpy as np
import matplotlib.pyplot as pylab
import operator
import glob
from collections import Counter,defaultdict
from scipy.stats import norm
import logging
from sklearn import linear_model

def total_known_outliers_number(self):
    outliers_number=0
    total_number=0
    for (tgid,tginst) in self.allgenedict.items():
        if tgid not in self.negative_control:
            eff_list=tginst.eff_estimate
            outliers_number+=len([i for i in eff_list if i==0])
            total_number+=len(eff_list)
    outliers_number_ratio=float(outliers_number)/total_number
    return outliers_number_ratio

def nest_selection_constnat_optimization(self):
    eff_likelihood_record=[]
    outliers_likelihood_record=[]
    for (tgid,tginst) in self.allgenedict.items():
        if tgid not in self.negative_control:
            likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
            likelihood_mean= np.matrix(np.mean(likelihood,axis=0))
            likelihood_list=likelihood_mean.tolist()[0]
            eff_list=tginst.eff_estimate
            for i in range(len(eff_list)):
                if eff_list[i]==0:
                    outliers_likelihood_record.append(likelihood_list[i])
                else:
                    eff_likelihood_record.append(likelihood_list[i])
    outliers_likelihood_record.sort()
    logging.info(outliers_likelihood_record)
    self.selection_constant=outliers_likelihood_record[int(len(outliers_likelihood_record)*0.9)]

def nest_selection(tginst,log_list,selection_constant):
    tginst.eff_estimate=[1]*len(tginst.sgrnaid)
    likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
    likelihood_mean= np.matrix(np.mean(likelihood,axis=0))
    likelihood_list=likelihood_mean.tolist()[0]
    for i in range(len(likelihood_list)):
        if likelihood_list[i]<selection_constant:
            tginst.eff_estimate[i]=0

def bayes_selection_constnat_optimization(self):
    removal_ratio=0.95

    likelihood_list=[]
    modified_likelihood_list=[]
    sgRNA_number=[] # total sgRNA_number

    for (tgid,tginst) in self.allgenedict.items():
        if tgid not in self.negative_control:
            #probability=tginst.sgrna_probability
            likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
            likelihood_mean= np.matrix(np.mean(likelihood,axis=0))
            #modified_likelihood_list+=[likelihood_mean.tolist()[0][i]-probability[i] for i in range(len(probability))]
            modified_likelihood_list+=likelihood_mean.tolist()[0]
            #probability_list+=probability
            likelihood_list+=(likelihood_mean.tolist()[0])
            sgRNA_number.append(len(tginst.sgrnaid))

    # selection_constant iteration start point
    sgRNA_number=[[sgRNA_number,count] for sgRNA_number,count in Counter(sgRNA_number).items()]
    sgRNA_number.sort(key=operator.itemgetter(1),reverse=True)
    common_sgRNA_number=sgRNA_number[0][0]
    modified_likelihood_list.sort(reverse=True)
    selection_constant=int(-modified_likelihood_list[int(len(modified_likelihood_list)*removal_ratio)]/(np.log(int(common_sgRNA_number*removal_ratio))-np.log(int(common_sgRNA_number*removal_ratio)-1)))

    # log_precaculation
    #log_list=defaultdict(list)
    for i in range(1,max([i[0] for i in sgRNA_number])+1):
        self.log_list[i]=[np.log(k+1) for k in range(i)]

    # selection constant searching
    constant_percentage={} # to record the percenge of on-target gRNAs given selection constant
    while selection_constant not in list(constant_percentage.keys()):
        total_sgRNA=0
        total_on_target_sgRNA=0
        for (tgid,tginst) in self.allgenedict.items():
            if tgid not in self.negative_control:
                likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
                likelihood_mean= np.matrix(np.mean(likelihood,axis=0)).tolist()
                temp=likelihood_mean[0]
                #probability=tginst.sgrna_probability
                #temp=[temp[i]-probability[i] for i in range(len(temp))]
                temp.sort(reverse=True)
                temp_accumulate=[np.sum([k for k in temp[:(i+1)]]) for i in range(len(temp))]
                modified_penalty=[selection_constant*i for i in self.log_list[len(tginst.sgrnaid)]]
                temp_accumulate_log_penalty=[sum(x) for x in zip(temp_accumulate,modified_penalty)]
                max_index=[i for i, j in enumerate(temp_accumulate_log_penalty) if j == max(temp_accumulate_log_penalty)]
                total_sgRNA+=len(tginst.sgrnaid)
                total_on_target_sgRNA+=(max_index[0]+1)
                if abs(total_sgRNA-3000)<10:
                    break
        constant_percentage[selection_constant]=float(total_on_target_sgRNA)/total_sgRNA
        if float(total_on_target_sgRNA)/total_sgRNA>removal_ratio:
            selection_constant-=1
        else:
            selection_constant+=1
        if [i>removal_ratio for i in list(constant_percentage.values())]==[True]*len(list(constant_percentage.values())):
            pass
        elif [i>removal_ratio for i in list(constant_percentage.values())]==[False]*len(list(constant_percentage.values())):
            pass
        else:
            break

    selection_percentage=[[constant,abs(percentage-removal_ratio)] for constant,percentage in constant_percentage.items()]
    selection_percentage.sort(key=operator.itemgetter(1))
    selection_constant=selection_percentage[0][0]
    self.selection_constant=selection_constant
    logging.info(selection_constant)

def bayes_selection(tginst,log_list,selection_constant):
    #likelihood_mean=np.mean(tginst.loglikelihood[i])
    likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
    likelihood_mean= np.matrix(np.mean(likelihood,axis=0)).tolist()[0]

    temp=[[tginst.sgrnaid[i],likelihood_mean[i]] for i in range(len(tginst.sgrnaid))]
    temp.sort(key=operator.itemgetter(1),reverse=True)


    temp_accumulate=[np.sum([k[1] for k in temp[:(i+1)]]) for i in range(len(temp))]

    modified_penalty=[selection_constant*i for i in log_list[len(tginst.sgrnaid)]]
    temp_accumulate_log_penalty=[sum(x) for x in zip(temp_accumulate,modified_penalty)]

    max_index=[i for i, j in enumerate(temp_accumulate_log_penalty) if j == max(temp_accumulate_log_penalty)]
    on_target_sgRNA=[temp[i][0] for i in range(max_index[-1]+1)]
    tginst.eff_estimate=[1]*len(tginst.sgrnaid)
    outliers_index=[orders for orders,sgRNA in enumerate(tginst.sgrnaid) if sgRNA not in on_target_sgRNA]
    for i in outliers_index:
        tginst.eff_estimate[i]=0

#!/usr/bin/env python3
'''
MAGeCK nest main entry
'''

import sys
import random
import math
import logging
import pickle
from collections import defaultdict

from mageck_nest.mleinstanceio import *
from mageck_nest.mleem import iteratenbem
from mageck_nest.mlemeanvar import MeanVarModel
from mageck_nest.mageckCount import normalizeCounts
from mageck_nest.bayes_selection import *
from mageck_nest.dispersion_characterization import *
from mageck_nest.mleargparse import *
from mageck_nest.mageck_nest_PPI import *
from mageck_nest.mageck_nest_output import *
from mageck_nest.mleclassdef import *
from mageck_nest.outliers_candidates import *

class Mageck_nest():
    def __init__(self, options):
        # Required
        self.count_table=options.count_table
        self.design_matrix=options.design_matrix

        # IO related
        self.beta_labels=options.beta_labels
        self.include_samples=options.include_samples
        self.output_prefix=options.output_prefix
        self.file_directory=os.getcwd()
        self.output_directory="%s/%s" %(os.getcwd(),self.output_prefix)
        if os.path.exists(self.output_directory)==False:
            os.makedirs(self.output_directory)

        # Normalization related
        self.adjust_method=options.adjust_method
        self.genes_varmodeling=options.genes_varmodeling
        self.negative_control=options.negative_control
        self.negative_control_gRNAs=None
        self.norm_method=options.norm_method

        # PPI and outliers removal
        self.outliers_removal=options.outliers_removal
        self.PPI_prior=options.PPI_prior
        self.QC_metric=options.QC_metric
        self.QC_path="/%s/QC_folder" %(self.output_directory)
        if os.path.exists(self.QC_path)==False:
            os.makedirs(self.QC_path)
        self.selection_constant=None
        self.PPI_diagnosis=None
        self.non_PPI_beta_prior_variance=None
        self.PPI_weighting=1
        self.suggested_remove_sgRNA=[]

        # Others
        self.allgenedict=None
        self.size_f=None
        self.mrm=None
        self.log_list=defaultdict(list)

    def nest_init(self):
        logging.info('Initiating ...')
        maxgene=np.inf

        if self.negative_control==None:
            self.allgenedict,self.invalid_gRNA_dict=read_gene_from_file(self.count_table,includesamples=self.include_samples,negative_control=self.negative_control)
        else:
            self.allgenedict,self.invalid_gRNA_dict,self.negative_control_gRNAs=read_gene_from_file(self.count_table,includesamples=self.include_samples,negative_control=self.negative_control)

        # calculate the size factor
        cttab_sel={}
        for (geneid,gk) in self.allgenedict.items():
            sgid=gk.sgrnaid
            sgreadmat=gk.nb_count.getT().tolist()
            for i in range(len(sgid)):
                cttab_sel[sgid[i]]=sgreadmat[i]
        if hasattr(self,'norm_method'):
            if self.norm_method!='none':
                self.size_f=normalizeCounts(cttab_sel,method=self.norm_method,returnfactor=True,reversefactor=True,negative_control_gRNAs=self.negative_control_gRNAs)
            else:
                self.size_f=None
        else:
            self.size_f=normalizeCounts(cttab_sel,returnfactor=True,reversefactor=True)
        logging.info('Size factor: '+','.join([str(x) for x in self.size_f]))
        desmat=self.design_matrix
        #------------------------------

        for (tgid,tginst) in self.allgenedict.items():
            if tgid not in self.negative_control:
                tginst.design_mat=desmat
                iteratenbem(tginst,debug=False,estimateeff=False,alpha_val=0.05,size_factor=self.size_f)
                tginst.w_estimate=[]

        self.non_PPI_beta_prior_variance=beta_non_PPI_prior_calculation(self.allgenedict,self.negative_control)
        PPI_import_string_9(self)
        PPI_weighting_rewiring(self)

    def nest_fitting(self):
        logging.info('Estimating dispersion factors ...')
        ngenes=0
        for (tgid,tginst) in self.allgenedict.items():
            if tgid not in self.negative_control:
                if ngenes<2000:
                    try:
                        sgrna_wide_dispersion_estimation_MAP_v2(tginst,self.design_matrix)
                        ngenes+=1
                    except:
                        pass

        logging.info('Modeling the mean and variance ...')

        self.mrm=MeanVarModel()
        self.mrm.model_mean_disp_by_glm(self.allgenedict,self.output_prefix,self.size_f)

    def nest_basic(self):
        logging.info('Calculating beta scores ...')

        self.suggested_remove_sgRNA=removal_suggestion(self)
        for (tgid,tginst) in self.allgenedict.items():
            if tgid not in self.negative_control:
                n_beta1=tginst.design_mat.shape[1]-1
                candidate_removed_tginst=[i for i in tginst.sgrnaid if i in self.suggested_remove_sgRNA]
                iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=True,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
                temp_non_PPI_beta_prior_variance=self.non_PPI_beta_prior_variance

                if abs(np.mean(tginst.beta_estimate[-n_beta1:]))>5:
                    temp_non_PPI_beta_prior_variance=[i/2 for i in temp_non_PPI_beta_prior_variance]
                    iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=True,size_factor=self.size_f,beta1_prior_var=temp_non_PPI_beta_prior_variance)
                if len(candidate_removed_tginst)>0 and len(tginst.sgrnaid)<30 and len(tginst.sgrnaid)>4 and abs(tginst.beta_estimate[-1])>5:
                    outliers_index=[orders for orders,sgRNA in enumerate(tginst.sgrnaid) if sgRNA in candidate_removed_tginst]
                    ratio_record=[]
                    tginst_record=[]
                    outliers_index_record=[]
                    tginst_2=copy.copy(tginst)
                    temp_2=tginst_2.beta_estimate[-n_beta1:]


                    for k in range(len(outliers_index)):
                        i=outliers_index[k]
                        self.outliers_removal=True
                        tginst.eff_estimate=[1]*len(tginst.sgrnaid)
                        tginst.eff_estimate[i]=0
                        iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=True,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
                        temp_non_PPI_beta_prior_variance=self.non_PPI_beta_prior_variance
                        if abs(np.mean(tginst.beta_estimate[-n_beta1:]))>5:
                            temp_non_PPI_beta_prior_variance=[i/2 for i in temp_non_PPI_beta_prior_variance]
                            iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=True,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=temp_non_PPI_beta_prior_variance)
                        tginst_1=copy.copy(tginst)
                        temp_1=tginst_1.beta_estimate[-n_beta1:]

                        ratio=np.log(abs(np.mean(temp_2)/np.mean(temp_1)))
                        ratio_record.append(abs(np.mean(temp_1)))
                        outliers_index_record.append(i)
                        if ratio>(5-0.2*len(tginst.sgrnaid)):
                            tginst_record.append(tginst_1)
                        else:
                            tginst_record.append(tginst_2)

                    ratio_record=[[i,j] for i,j in enumerate(ratio_record)]
                    ratio_record.sort(key=operator.itemgetter(1))

                    if ratio_record[0][1]>5 and len(ratio_record)>=2:
                        logging.info(ratio_record[0][1])
                        two_outliers=[outliers_index_record[ratio_record[0][0]],outliers_index_record[ratio_record[1][0]]]
                        tginst.eff_estimate=[1]*len(tginst.sgrnaid)
                        for i in two_outliers:
                            tginst.eff_estimate[i]=0
                        iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=True,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)

                    else:
                        self.allgenedict[tgid]=tginst_record[ratio_record[0][0]]

        for (tgid,tginst) in self.allgenedict.items():
            if tgid not in self.negative_control:
                if len(tginst.w_estimate)==0:
                    tginst.w_estimate=np.ones(len(tginst.sgrnaid))

        pickle.dump(self,open("/%s/%s_self_nest_major.p" %(self.output_directory,self.output_prefix),'wb'))

        #self=pickle.load(open("/%s/%s_self_nest_major.p" %(self.output_directory,self.output_prefix),'rb'))
        nest_output(self,["False","False"])

    def constant_optimization(self):
        self=pickle.load(open("/%s/%s_self_nest_major.p" %(self.output_directory,self.output_prefix),'rb'))

        if self.PPI_prior==True:
            logging.info('PPI validation...')
            PPI_main(self)
            if self.PPI_diagnosis==True:
                beta_prior_output(self)
            if self.PPI_diagnosis==False:
                logging.info("The overlapped number of input genes and PPI genes is less than <3000.")
                logging.info("PPI is not recommended.")

        if self.outliers_removal==True:
            logging.info('Estimate selection constant...')
            outliers_number_ratio=total_known_outliers_number(self)
            if outliers_number_ratio>0.05:
                nest_selection_constnat_optimization
                logging.info("Selection constant: %s" %self.selection_constant)
                for (tgid,tginst) in list(self.allgenedict.items()):
                    if tgid not in self.negative_control:
                        nest_selection(tginst,log_list=self.log_list,selection_constant=self.selection_constant)
            else:
                bayes_selection_constnat_optimization(self)
                logging.info("Selection constant: %s" %self.selection_constant)
                for (tgid,tginst) in list(self.allgenedict.items()):
                    if tgid not in self.negative_control:
                        bayes_selection(tginst,log_list=self.log_list,selection_constant=self.selection_constant)

        pickle.dump(self,open("/%s/%s_self_constant_optimization.p" %(self.output_directory,self.output_prefix),'wb'))

    def nest_iteration(self):
        logging.info('Final iteratioin for PPI or outliers removal...')
        os.chdir(self.output_directory)
        if self.PPI_prior==False and self.outliers_removal==True:
            marks=[[False,True]]
        elif self.PPI_prior==True and self.outliers_removal==True:
            marks=[[False,True],[True,False],[True,True]]
        elif self.PPI_prior==True and self.outliers_removal==False:
            marks=[[True,False]]

        for mark in marks:
            self=pickle.load(open("/%s/%s_self_constant_optimization.p" %(self.output_directory,self.output_prefix),'rb'))
            self.PPI_prior=mark[0]
            self.outliers_removal=mark[1]

            for (tgid,tginst) in self.allgenedict.items():
                if tgid not in self.negative_control:
                    iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=False,PPI_prior=self.PPI_prior,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
            #pickle.dump(self,open("/%s/%s_self_nest_iteration_PPI_%s_outliers_removal_%s.p" %(self.output_directory,self.output_prefix,self.PPI_prior,self.outliers_removal),'wb'))
            nest_output(self,[str(mark[0]),str(mark[1])])

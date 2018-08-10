'''
Reading and writing instances
'''


import sys
import numpy as np
from mageck_nest.mleclassdef import *
from collections import defaultdict
import itertools

import logging

def read_gene_from_file(filename,includesamples=None,negative_control=None):
    '''
    Reading gene models
    Parameters
    ----------
    filename
        file name of the read count table
    includesamples
        If not None, only samples in the includesampels are included
    '''
    # first, read count table
    allgenedict={}
    invalid_gRNA_dict=defaultdict(list)
    nline=0
    nsamples=0
    ngene=0
    negative_control_gRNAs=[]

    sampleids=[]
    sampleindex=[]
    sampleids_toindex={}
    hascsv=False
    if filename.upper().endswith('.CSV'):
        hascsv=True
        logging.info('Treating '+filename+' as csv file ...')
    for line in open(filename):
        nline+=1
        if hascsv==False:
            field=line.strip().split()
        else:
            field=line.strip().split(',')
        if nline==1:
            # The first line: check sample columns
            nsamples=len(field)-2
            sampleids=field[2:]
            for i in range(nsamples):
                sampleids_toindex[sampleids[i]]=i
            if includesamples != None:
                #logging.info('Loaded samples:'+';'.join(includesamples))
                for si in includesamples:
                    if si not in sampleids_toindex:
                        logging.error('Sample '+si+' cannot be found on the original read count table '+filename)
                        sys.exit(-1)
                sampleindex=[sampleids_toindex[si] for si in includesamples]
                #logging.info('Sample index: '+';'.join([str(x) for x in sampleindex]))
            else:
                sampleindex=[i for i in range(nsamples)]
            continue
        sgid=field[0].upper()
        gid=field[1].upper()
        if negative_control!=None:
            if gid in negative_control:
                negative_control_gRNAs.append(sgid)

        if gid not in allgenedict:
            sks=SimCaseSimple()
            sks.prefix=gid
            sks.nb_count=[]
            sks.sgrnaid=[]
            ngene+=1
            for i in sampleindex:
                sks.nb_count+=[[]]
            allgenedict[gid]=sks
        else:
            sks=allgenedict[gid]
        sks.sgrnaid+=[sgid]
        for i in range(len(sampleindex)):
            ni=sampleindex[i]
            try:
                #nrt=float(field[ni+2])+1 # add 1 pseudocount
                nrt=float(field[ni+2]) # no pseudocount
                sks.nb_count[i]+=[nrt]
            except ValueError:
                print('Error loading line '+str(nline))
    # end for loop
    logging.info('Loaded '+str(ngene)+' genes.')
    #
    # convert nb_count to matrix
    empty_gene=[]

    for (gene,ginst) in allgenedict.items():
        ginst.nb_count=np.matrix(ginst.nb_count)
        count_mean=np.mean(ginst.nb_count,axis=0)
        invalid_index=np.where(count_mean==0)[1]

        if len(invalid_index)!=0:
            invalid_gRNA_dict[gene]=[ginst.sgrnaid[i] for i in range(len(ginst.sgrnaid)) if i in invalid_index]

        temp=[ginst.sgrnaid[i] for i in range(ginst.nb_count.shape[1]) if i not in invalid_index]
        ginst.sgrnaid=temp

        invalid_cols = np.nonzero(ginst.nb_count.sum(axis=0) == 0)[1].tolist()
        temp=np.delete(ginst.nb_count,invalid_cols,axis=1)
        ginst.nb_count=temp+1

        if ginst.sgrnaid==[]:
            empty_gene.append(gene)

    allgenedict={key:values for key,values in list(allgenedict.items()) if key not in empty_gene}
    if negative_control==None:
        return allgenedict,invalid_gRNA_dict
    else:
        return allgenedict,invalid_gRNA_dict,negative_control_gRNAs

def write_gene_to_file(allgenedict,outfile,betalabels=None,wald_only=False):
    '''
    Write gene to file
    '''
    ofid=open(outfile,'w')
    tmpinst=allgenedict[list(allgenedict.keys())[0]]
    nbeta=len(tmpinst.beta_estimate)-(tmpinst.nb_count.shape[1])
    # print header
    # headerterms=['|beta','|z','|neg|p-value','|neg|fdr','|pos|p-value','|pos|fdr', '|neg|permutation', '|neg|permutation-fdr', '|pos|permutation','|pos|permutation-fdr' ] # one-sided
    # headerterms=['|beta','|z','|p-value','|fdr', '|permutation', '|permutation-fdr' ] # two-sided,using Wald test for p value
    headerterms=['|beta','|z','|p-value','|fdr','|wald-p-value','|wald-fdr' ] # two-sided,using permutation test for p value
    if wald_only==True:
        headerterms=['|beta','|z','|wald-p-value','|wald-fdr' ]
    if betalabels == None:
        # reportlabels='\t'.join(['\t'.join(['beta_'+str(i+1)+'|beta','beta_'+str(i+1)+'|z','beta_'+str(i+1)+'|p-value','beta_'+str(i+1)+'|permutation']) for i in range(nbeta)])  # two-sided
        reportlabels='\t'.join(['\t'.join(['beta_'+str(i+1)+ hst for hst in headerterms ]) for i in range(nbeta)]) # one-sided
    else:
        if len(betalabels)-1!=nbeta:
            raise ValueError('Beta labels do not match to columns of the design matrix.')
        # reportlabels='\t'.join(['\t'.join([sstt+'|beta',sstt+'|z',sstt+'|p-value',sstt+'|permutation']) for sstt in betalabels[1:]]) # two-sided
        reportlabels='\t'.join(['\t'.join([sstt+hst for hst in headerterms]) for sstt in betalabels[1:]])
    print('\t'.join(['Gene','sgRNA',reportlabels]),file=ofid)
    #ofid.write('\t'.join(['Gene','sgRNA',reportlabels]))
    # print for each gene
    for (tgid,tginst) in allgenedict.items():
        if wald_only==False:
            wfield=tginst.gene_to_printfield()
        else:
            wfield=tginst.gene_to_printfield(wald_only=True)
        #ofid.write('\t'.join(wfield))
        print('\t'.join(wfield),file=ofid)
    # end
    ofid.close()

def write_sgrna_to_file(include_samples,allgenedict,invalid_gRNA_dict,outfile,negative_control):
    '''
    Write gene to file
    '''
    ofid=open(outfile,'w')
    tmpinst=allgenedict[list(allgenedict.keys())[0]]
    mu_temp=["%s_mu" %i for i in include_samples]
    k_temp=["%s_k" %i for i in include_samples]
    sample_temp=mu_temp+k_temp
    #ofid.write('\t'.join(['Gene','sgRNA','eff']+sample_temp+["loglikelihood"]))
    print('\t'.join(['Gene','sgRNA','eff']+sample_temp+["loglikelihood"]),file=ofid)
    for (tgid,tginst) in allgenedict.items():
        if tgid not in negative_control:
            mu=tginst.mu_estimate.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
            k=tginst.sgrna_kvalue.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
            #fitting_alpha=np.matrix(tginst.dispersion_estimate).reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
            #MAP_alpha=np.matrix(tginst.MAP_sgrna_dispersion_estimate).reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
            likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1])).T
            likelihood_list= np.matrix(likelihood).tolist()
            for i in range(len(tginst.eff_estimate)):
                k_i=[str(g) for g in k[i,].tolist()[0]]
                mu_i=[str(g) for g in mu[i,].tolist()[0]]
                #fitting_a_i=[str(g) for g in fitting_alpha[i,].tolist()[0]]
                #insert=k_i+mu_i+[fitting_a_i[0]]
                #wfield=[tginst.prefix,tginst.sgrnaid[i],decformat(tginst.eff_estimate[i])]+insert+[str(j) for j in likelihood_list[i]]
                insert=k_i+mu_i
                wfield=[tginst.prefix,tginst.sgrnaid[i],decformat(tginst.eff_estimate[i])]+insert+[str(j) for j in likelihood_list[i]]
                #ofid.write('\t'.join(wfield))
                print('\t'.join(wfield),file=ofid)
            if tgid in list(invalid_gRNA_dict.keys()):
                for sgrna in invalid_gRNA_dict[tgid]:
                    wfield=[tgid,sgrna]+["NA"]*((tginst.nb_count.shape[0])*2+2)
                    #ofid.write('\t'.join(wfield))
                    print('\t'.join(wfield),file=ofid)
    ofid.close()

def beta_prior_output(self):
    output=open("/%s/%s_beta_prior.txt" %(self.output_directory,self.output_prefix),'wb')
    for (tgid,tginst) in list(self.allgenedict.items()):
        if tgid not in self.negative_control:
            insert=[tgid]+tginst.prior_mean+tginst.prior_variance
            output.write("{}\n".format("\t".join([str(k) for k in insert])).encode())

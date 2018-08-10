from __future__ import print_function, division, absolute_import
from functools import reduce

import os
import time
import numpy as np
import sys
from bs4 import BeautifulSoup
from numpy import in1d
from pandas import read_table, DataFrame
import sys
import warnings
import matplotlib.transforms as transforms
from matplotlib.colors import Normalize
from collections import OrderedDict
import pandas as pd
import matplotlib.pylab as plt

class _MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def gsea_plot(rank_metric, enrich_term, hit_ind, nes, pval, fdr, RES,
              phenoPos=None, phenoNeg=None, figsize =(6.5,6), **kwarg):
    """This is the main function for reproducing the gsea plot.

    :param rank_metric: rankings, rank_metric['rank'].values.
    :param enrich_term: gene_set name
    :param hit_ind: hit indexs of rank_metric['gene_name'] presented in gene set S.
    :param nes: Normalized enrichment scores.
    :param pval: nominal p-value.
    :param fdr: false discoveray rate.
    :param RES: ranking enrichment scores of all genes in rank_metric['gene_name'].
    :param phenoPos: phenotype lable, positive correlated.
    :param phenoNeg: phenotype lable, negative correlated.
    :param figsize: matplotlib figsize.
    :return: fig object of gsea plot.
    """
    #plt.style.use('classic')
    # center color map at midpoint = 0
    norm = _MidpointNormalize(midpoint=0)

    #dataFrame of ranked matrix scores
    x = rank_metric.index.values
    #figsize = (6,6)
    phenoP_label = phenoPos + ' (Positively Correlated)'
    phenoN_label = phenoNeg + ' (Negatively Correlated)'
    zero_score_ind = np.abs(rank_metric['rank']).argmin()
    z_score_label = 'Zero score at ' + str(zero_score_ind)
    nes_label = 'NES: '+ "{:.3f}".format(float(nes))
    pval_label = 'Pval: '+ "{:.3f}".format(float(pval))
    fdr_label = 'FDR: '+ "{:.3f}".format(float(fdr))
    im_matrix = rank_metric.ix[:,1:].T

    #in most case, we will have mangy plots, so do not display plots
    #It's also convinient to run this script on command line.
    plt.ioff()
    #GSEA Plots
    gs = plt.GridSpec(16,1)
    fig = plt.figure(figsize=figsize)
    #Ranked Metric Scores Plot
    ax1 =  fig.add_subplot(gs[11:])
    ax1.fill_between(x, y1= rank_metric['rank'], y2=0, color='#C9D3DB')
    ax1.set_ylabel("Ranked list metric",fontsize=14)
    ax1.text(.05, .9, phenoP_label, color='red', horizontalalignment='left', verticalalignment='top',
         transform=ax1.transAxes)
    ax1.text(.95, .05, phenoN_label, color='Blue', horizontalalignment='right', verticalalignment='bottom',
         transform=ax1.transAxes)

    # the x coords of this transformation are data, and the y coord are axes
    trans1 = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.vlines(zero_score_ind, 0, 1, linewidth=.5, transform=trans1, linestyles='--', color='grey')
    ax1.text(zero_score_ind, 0.5, z_score_label, horizontalalignment='center', verticalalignment='center',
             transform=trans1)
    ax1.set_xlabel("Rank in Ordered Dataset", fontsize=14)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(axis='both', which='both', top='off', right='off', left='off')
    ax1.locator_params(axis='y', nbins=5)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc) ))

    # use round method to control float number
    #ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  round(tick_loc, 1) ))

    #gene hits
    ax2 = fig.add_subplot(gs[8:10], sharex=ax1)

    # the x coords of this transformation are data, and the y coord are axes
    trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ax2.vlines(hit_ind, 0, 1,linewidth=.5,transform=trans2)
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off',labelleft='off')
    #colormap
    ax3 =  fig.add_subplot(gs[10],sharex=ax1)
    ax3.imshow(im_matrix, aspect='auto', norm=norm, cmap=plt.cm.seismic, interpolation='none') # cm.coolwarm
    ax3.spines['bottom'].set_visible(False)
    ax3.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off',labelleft='off')

    # Enrichment score plot
    ax4 = fig.add_subplot(gs[:8],sharex=ax1)
    ax4.plot(x,RES,linewidth=4,color ='#88C544')
    ax4.text(.1, .1, fdr_label, transform=ax4.transAxes)
    ax4.text(.1, .2, pval_label, transform=ax4.transAxes)
    ax4.text(.1, .3, nes_label, transform=ax4.transAxes)

    # the y coords of this transformation are data, and the x coord are axes
    trans4 = transforms.blended_transform_factory(ax4.transAxes, ax4.transData)
    ax4.hlines(0, 0, 1, linewidth=.5, transform=trans4, color='grey')
    ax4.set_ylabel("Enrichment score (ES)", fontsize=14)
    ax4.set_xlim(min(x), max(x))
    ax4.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off')
    ax4.locator_params(axis='y', nbins=5)
    # FuncFormatter need two argment, I don't know why. this lambda function used to format yaxis tick labels.
    ax4.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc)) )

    #fig adjustment
    fig.suptitle(enrich_term, fontsize=16)
    fig.subplots_adjust(hspace=0)
    #fig.tight_layout()
    plt.close(fig)

    return fig



def gsea_rank_metric(rnk):
    rank_metric = read_table(rnk,header=None)
    rank_metric.columns = ['gene_name','rank']
    rank_metric['rank2'] = rank_metric['rank']

    return rank_metric

def gsea_gmt_parser(gmt, min_size = 3, max_size = 1000, gene_list=None):
    with open(gmt) as genesets:
        genesets_dict = { line.strip("\n").split("\t")[0]: line.strip("\n").split("\t")[2:] for line in genesets.readlines()}

    #filtering dict
    if sys.version_info[0] == 3 :
        genesets_filter =  {k: v for k, v in genesets_dict.items() if len(v) >= min_size and len(v) <= max_size}
    elif sys.version_info[0] == 2:
        genesets_filter =  {k: v for k, v in genesets_dict.iteritems() if len(v) >= min_size and len(v) <= max_size}
    else:
        print("System failure. Please Provide correct input files")
        sys.exit(1)
    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())
        for subset in subsets:
            tag_indicator = in1d(gene_list, genesets_filter.get(subset), assume_unique=True)
            tag_len = sum(tag_indicator)
            if tag_len <= min_size or tag_len >= max_size:
                del genesets_filter[subset]
            else:
                continue

    filsets_num = len(genesets_dict) - len(genesets_filter)

    if filsets_num == len(genesets_dict):
        print("No gene sets passed throught filtering condition!!!, try new paramters again!\n" +\
              "Note: Gene names for gseapy is case sensitive." )
        sys.exit(1)
    else:
        return genesets_filter


def enrichment_score(gene_list, gene_set, weighted_score_type=1, correl_vector=None, esnull=None, rs=np.random.RandomState()):
    """This is the most important function of GSEAPY. It has the same algorithm with GSEA.

    :param gene_list:       The ordered gene list gene_name_list, rank_metric['gene_name']
    :param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
    :param weighted_score_type:  It's indentical to gsea's weighted_score method. weighting by the correlation
                            is a very reasonable choice that allows significant gene sets with less than perfect coherence.
                            options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
                            coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
                            might be appropriate. On the other hand, if one uses sets with largenumber of genes and only
                            a small subset of those is expected to be coherent, then one could consider using p > 1.
                            Our recommendation is to use p = 1 and use other settings only if you are very experienced
                            with the method and its behavior.

    :param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
                            the gene list. Or rankings, rank_metric['rank'].values
    :param esnull:          Only used this paramter when computing esnuall for statistial testing. set the esnull value
                            equal to the permutation number.
    :param rs:              Random state for initialize gene list shuffling. Default: np.random.RandomState(seed=None)

    :return:

     ES: Enrichment score (real number between -1 and +1)

     hit_index: index of a gene in gene_list, if gene included in gene_set.

     RES: Numerical vector containing the running enrichment score for all locations in the gene list .

    """

    axis = 0
    N = len(gene_list)

    #Test whether each element of a 1-D array is also present in a second array
    #It's more intuitived here than orginal enrichment_score source code.
    #use .astype to covert bool to intergers
    tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)

    if (weighted_score_type == 0 ):
        correl_vector = np.repeat(1, N)
    else:
        correl_vector = np.abs(correl_vector**weighted_score_type)


    #get indices of tag_indicator
    hit_ind = np.flatnonzero(tag_indicator).tolist()

    Nhint = np.sum(tag_indicator)
    sum_correl_tag = np.sum(correl_vector[tag_indicator.astype(bool)])

    # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
    if esnull:
        tag_indicator = tag_indicator.repeat(esnull).reshape(N, esnull).T
        correl_vector = correl_vector.repeat(esnull).reshape(N, esnull).T

        # gene list permutation
        for i in range(esnull):
            rs.shuffle(tag_indicator[i])

        # set axis to 1, because we have 2 dimentional array
        axis = 1
        Nhint = np.sum(tag_indicator, axis=axis).reshape(esnull, 1)
        sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis).reshape(esnull, 1)
    '''
    #similar results could be obtained when computing esnull using code below
    #
    if esnull:
        tag_null = np.empty((esnull, N))
        i=0
        while i < esnull:
            rs.shuffle(tag_indicator)
            tag_null[i] = tag_indicator
            i +=1
        axis = 1
        tag_indicator = tag_null
    '''


    #compute ES score, the code below is identical to gsea enrichment_score method.
    no_tag_indicator = 1 - tag_indicator
    Nmiss =  N - Nhint
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nmiss

    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)
    max_ES = np.max(RES, axis=axis)
    min_ES = np.min(RES, axis=axis)

    es = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)

    return es.tolist(), hit_ind, RES.tolist()


def ranking_metric(df, method, phenoPos, phenoNeg, classes, ascending):
    """The main function to rank an expression table.

   :param df:      gene_expression DataFrame.
   :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                   Others methods are:

                   1. 'signal_to_noise'

                      You must have at least three samples for each phenotype to use this metric.
                      The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
                      that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.”

                   2. 't_test'

                      Uses the difference of means scaled by the standard deviation and number of samples.
                      Note: You must have at least three samples for each phenotype to use this metric.
                      The larger the tTest ratio, the more distinct the gene expression is in each phenotype
                      and the more the gene acts as a “class marker.”

                   3. 'ratio_of_classes' (also referred to as fold change).

                      Uses the ratio of class means to calculate fold change for natural scale data.

                   4. 'diff_of_classes'

                      Uses the difference of class means to calculate fold change for log scale data

                   5. 'log2_ratio_of_classes'

                      Uses the log2 ratio of class means to calculate fold change for natural scale data.
                      This is the recommended statistic for calculating fold change for natural scale data.


   :param phenoPos: one of lables of phenotype's names.
   :param phenoNeg: one of lable of phenotype's names.
   :param classes:  a list of phenotype labels, to specify which column of dataframe belongs to what catogry of phenotype.
   :param ascending:  bool or list of bool. Sort ascending vs. descending.
   :return: returns correlation to class of each variable. same format with .rnk file. gene_name in first coloum,
            correlation in second column.

    visit here for more docs: http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
    """

    A = phenoPos
    B = phenoNeg
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T
    if method == 'signal_to_noise':
        sr = (df_mean[A] - df_mean[B])/(df_std[A] + df_std[B])
    elif method == 't_test':
        sr = (df_mean[A] - df_mean[B])/ np.sqrt(df_std[A]**2/len(df_std)+df_std[B]**2/len(df_std) )
    elif method == 'ratio_of_classes':
        sr = df_mean[A] / df_mean[B]
    elif method == 'diff_of_classes':
        sr  = df_mean[A] - df_mean[B]
    elif method == 'log2_ratio_of_classes':
        sr  =  np.log2(df_mean[A] / df_mean[B])
    else:
        print("Please provide correct method name!!!")
        sys.exit()
    sr.sort_values(ascending=ascending, inplace=True)
    df3 = sr.to_frame().reset_index()
    df3.columns = ['gene_name','rank']
    df3['rank2'] = df3['rank']

    return df3

def gsea_compute(data, gmt, n, weighted_score_type, permutation_type, method,
                 phenoPos, phenoNeg, classes, ascending, seed, prerank=False):
    """compute enrichment scores and enrichment nulls.

    :param data: prepreocessed expression dataframe or a pre-ranked file if prerank=True.
    :param gmt: all gene sets in .gmt file. need to call gsea_gmt_parser() to get results.
    :param n: permutation number. default: 1000.
    :param method: ranking_metric method. see above.
    :param phenoPos: one of lables of phenotype's names.
    :param phenoNeg: one of lable of phenotype's names.
    :param classes: a list of phenotype labels, to specify which column of dataframe belongs to what catogry of phenotype.
    :param weighted_score_type: default:1
    :param ascending: sorting order of rankings. Default: False.
    :param seed: random seed. Default: np.random.RandomState()
    :param prerank: if true, this function will compute using pre-ranked file passed by parameter data.

    :return:
      zipped results of es, nes, pval, fdr. Used for generating reportes and plotting.

      a nested list of hit indexs of input gene_list. Used for plotting.

      a nested list of ranked enrichment score of each input gene_sets. Used for plotting.

    """
    enrichment_scores = []
    w = weighted_score_type
    subsets = sorted(gmt.keys())
    dat = data.copy()
    if prerank:
        r2 =data.copy()
    else:
        r2 = ranking_metric(df=dat, method=method, phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending)
    ranking=r2['rank'].values
    gene_list=r2['gene_name']


    rank_ES = []
    hit_ind = []
    for subset in subsets:
        es,ind,RES = enrichment_score(gene_list=gene_list, gene_set=gmt.get(subset),
                              weighted_score_type=w, correl_vector=ranking)
        enrichment_scores.append(es)
        rank_ES.append(RES)
        hit_ind.append(ind)


    enrichment_nulls = [ [] for a in range(len(subsets)) ]
    rs = np.random.RandomState(seed)


    if permutation_type == "phenotype":

        dat2 = dat.T
        for i in range(n):
            dat2.apply(rs.shuffle, axis=0) #permutation classes
            r2 = ranking_metric(df=dat2.T, method=method, phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending )
            ranking2=r2['rank']
            gene_list2=r2['gene_name'].values


            for si,subset in enumerate(subsets):
                esn = enrichment_score(gene_list=gene_list2, gene_set=gmt.get(subset),
                                       weighted_score_type=w, correl_vector=ranking2)[0]
                enrichment_nulls[si].append(esn)
    else:
        for si,subset in enumerate(subsets):
            esn = enrichment_score(gene_list=gene_list, gene_set=gmt.get(subset), weighted_score_type=w,
                                   correl_vector=ranking, esnull=n, rs=rs)[0]
            enrichment_nulls[si] = esn # esn is a list, don't need to use append method.

    return gsea_significance(enrichment_scores, enrichment_nulls),hit_ind,rank_ES, subsets



def gsea_pval(es, esnull):
    """Compute nominal p-value.

    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign
    of the observed ES(S).
    """

    # to speed up, using numpy function to compute pval in parallel.
    es = np.array(es)
    esnull = np.array(esnull)
    try:
        condlist = [ es < 0, es >=0]
        choicelist = [np.sum(esnull < es.reshape(len(es),1), axis=1)/ np.sum(esnull < 0, axis=1) ,
                      np.sum(esnull >= es.reshape(len(es),1), axis=1)/ np.sum(esnull >= 0, axis=1)]
        pval = np.select(condlist, choicelist)

        return pval
    except:
        return np.repeat(1.0 ,len(es))



def normalize(es, enrNull):
    """normalize the ES(S,pi) and the observed ES(S), separetely rescaling
       the positive and negative scores by divident by the mean of the ES(S,pi).
    """

    try:
        if es == 0:
            return 0.0
        if es >= 0:
            meanPos = np.mean([a for a in enrNull if a >= 0])
            #print es, meanPos
            return es/meanPos
        else:
            meanNeg = np.mean([a for a in enrNull if a < 0])
            #print es, meanNeg
            return -es/meanNeg
    except:

        return 0.0 #return if according mean value is uncalculable

def gsea_significance(enrichment_scores, enrichment_nulls):
    """Compute nominal p-vals, normalized ES, and FDR q value.

        For a given NES(S) = NES* >= 0. The FDR is the ratio of the percantage of all (S,pi) with
        NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
        observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S) = NES* <= 0.
    """

    #print("Start to compute pvals..................................", time.ctime())

    #enrichmentPVals = []

    #compute pvals.
    enrichmentPVals = gsea_pval(enrichment_scores, enrichment_nulls).tolist()

    #new normalize enrichment score calculating method. this could speed up significantly.
    esnull_meanPos = []
    esnull_meanNeg = []

    es = np.array(enrichment_scores)
    esnull = np.array(enrichment_nulls)

    for i in range(len(enrichment_scores)):
        enrNull = esnull[i]
        meanPos = enrNull[enrNull >= 0].mean()
        esnull_meanPos.append(meanPos)


        meanNeg = enrNull[enrNull < 0 ].mean()
        esnull_meanNeg.append(meanNeg)


    pos = np.array(esnull_meanPos).reshape(len(es), 1)
    neg = np.array(esnull_meanNeg).reshape(len(es), 1)


    #compute normalized enrichment score and normalized esnull
    try:
        condlist1 = [ es >= 0, es < 0]
        choicelist1 = [ es/esnull_meanPos, -es/esnull_meanNeg ]
        nEnrichmentScores = np.select(condlist1, choicelist1).tolist()

        condlist2 = [ esnull >= 0, esnull < 0]
        choicelist2 = [ esnull/pos, -esnull/neg ]
        nEnrichmentNulls = np.select(condlist2, choicelist2).tolist()

    except:  #return if according nes, nesnull is uncalculable
        nEnrichmentScores = np.repeat(0.0, es.size).tolist()
        nEnrichmentNulls = np.repeat(0.0 , es.size).reshape(esnull.shape).tolist()

    vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])

    """
    def shorten(l, p=10000):

        #Take each len(l)/p element, if len(l)/p >= 2.

        e = len(l)/p
        if e <= 1:
            return l
        else:
            return [ l[i] for i in range(0, len(l), e) ]

    vals = shorten(vals) -> this can speed up second part. is it relevant TODO?

    """

    nvals = np.array(sorted(vals))
    nnes = np.array(sorted(nEnrichmentScores))
    fdrs = []

    for i in range(len(enrichment_scores)):
        nes = nEnrichmentScores[i]
        #this could be speed up twice
        if nes >= 0:
            allPos = int(len(vals) - np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(vals) - np.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(np.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(np.searchsorted(nvals, nes, side="right"))
            nesPos = int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(np.searchsorted(nnes, nes, side="right"))
        try:
            top = allHigherAndPos/float(allPos) #p value
            down = nesHigherAndPos/float(nesPos)

            fdrs.append(top/down)
        except:
            fdrs.append(1000000000.0)

    return zip(enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs)




def prerank(rnk, gene_sets, outdir='gseapy_out', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=1000, permutation_n=1000, weighted_score_type=1,
            ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None):

    #drop duplicates in ranking metrics.
    dat2 = gsea_rank_metric(rnk)
    dat2.drop_duplicates(subset='gene_name',inplace=True,keep='first')
    assert len(dat2) > 1

    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size=min_size, max_size=max_size, gene_list=dat2['gene_name'].values)

    #compute ES, NES, pval, FDR, RES
    warnings.filterwarnings("ignore")
    results,hit_ind,rank_ES, subsets = gsea_compute(data=dat2, n=permutation_n, gmt=gmt, weighted_score_type=weighted_score_type,
                                                    permutation_type='gene_set', method=None, phenoPos=pheno_pos, phenoNeg=pheno_neg,
                                                    classes=None, ascending=ascending, seed=seed, prerank=True)

    res = OrderedDict()
    for gs,gseale,ind,RES in zip(subsets, list(results), hit_ind, rank_ES):
        rdict = OrderedDict()
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat2.ix[ind,'gene_name'].tolist()
        rdict['hit_index'] = ind
        res[gs] = rdict


    res_df = pd.DataFrame.from_dict(res, orient='index')
    res_df.index.name = 'Enrich_terms'
    res_df.sort_values(by='fdr', inplace=True)
    res_final = res_df.head(graph_num)
    #res_df.to_csv('{a}/{b}.csv'.format(a=outdir, b='gseapy_reports'), float_format ='%.7f')

    for gs in res_final.index.values:
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=res.get(gs)['hit_index'],
                        nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'],
                        RES=res.get(gs)['rank_ES'], phenoPos=pheno_pos, phenoNeg=pheno_neg, figsize=figsize)
        fig.savefig('{a}/{b}.{c}'.format(a=outdir, b=gs, c=format), dpi=300,)

    for key,values in res.items():
        return(values["es"],values["pval"],values["fdr"])

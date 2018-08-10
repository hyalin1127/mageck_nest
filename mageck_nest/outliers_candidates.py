import os
from collections import Counter,defaultdict
import numpy as np
import operator
from scipy.stats import norm
import logging

def removal_suggestion(self):
    reverse_size_factor=[float(1)/i for i in self.size_f]

    transform_design_list=[tuple(i) for i in self.design_matrix.tolist()]
    categories=Counter(transform_design_list).keys()
    indexes=[]
    for category in categories:
        indexes.append([i for i in range(len(transform_design_list)) if transform_design_list[i]==category])
    include_samples_categories=[[self.include_samples[j] for j in i] for i in indexes]
    reverse_size_factor_categories=[[reverse_size_factor[j] for j in i] for i in indexes]

    all_combinations=[]
    for i in range(len(include_samples_categories)):
        for j in range(len(include_samples_categories)):
            if j>i:
                all_combinations.append((i,j))

    file=open(self.count_table)
    if (self.count_table).upper().endswith('.CSV'):
        field=file.readline().strip().split(",")
    else:
        field=file.readline().strip().split('\t')
    field_index=[[field.index(i) for i in j] for j in include_samples_categories]

    fc_record=defaultdict(list)
    fc_detailed_record=dict()
    for combination in all_combinations:
        fc_detailed_record[combination]=defaultdict(list)

    for line in file:
        if (self.count_table).upper().endswith('.CSV'):
            elements=file.readline().strip().split(",")
        else:
            elements=line.strip().split('\t')
        temp=[]
        for combination in all_combinations:
            values_1=np.multiply([float(elements[i])+1 for i in field_index[combination[0]]],reverse_size_factor_categories[combination[0]])
            values_2=np.multiply([float(elements[i])+1 for i in field_index[combination[1]]],reverse_size_factor_categories[combination[1]])
            fc_record[combination].append([elements[0].upper(),np.log(float(np.mean(values_1))/np.mean(values_2))])
            fc_detailed_record[combination][elements[1]].append([elements[0].upper(),np.log(float(np.mean(values_1))/np.mean(values_2))])

    suggested_remove_sgRNA=[]
    for combination,fc in fc_record.items():
        fc_value=[i[1] for i in fc]
        fc_value_median=np.median(fc_value)
        abs_fc=[abs(i-fc_value_median) for i in fc_value]
        abs_fc.sort()
        var_fc=float(np.percentile(abs_fc,68))/norm.ppf(0.84)
        for gene,detailed_fc in fc_detailed_record[combination].items():
            gene_median=np.median([i[1] for i in detailed_fc])
            for k in detailed_fc:
                if abs(k[1]/var_fc)>1.5 and abs(k[1]-gene_median)>1.5:
                    suggested_remove_sgRNA.append(k[0])
                    #if gene=="PCNA":
                    #    logging.info(combination)
                    #    logging.info(k[0])
        #suggested_remove_sgRNA+=[i[0] for i in fc if abs(i[1]/var_fc)>1.5]
        '''
        suggested_remove_sgRNA+=[i[0] for i in fc[:int((len(fc)*0.05))]]
        suggested_remove_sgRNA+=[i[0] for i in fc[int(len(fc)*0.95):]]
        '''
    return suggested_remove_sgRNA

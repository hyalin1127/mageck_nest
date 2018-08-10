
import argparse,textwrap

def crisprseq_parseargs():
    """
    Parsing mageck arguments.
    """
    parser=argparse.ArgumentParser(description='MAGeCK-NEST: performs sgRNA, gene and pathway analysis on CRISPR-Cas9 screening data.')
    # definition of sub commands
    parser.add_argument('-v', '--version',action='version',version='%(prog)s 0.0.1')
    subparser=parser.add_subparsers(help='commands to run mageck',dest='subcmd')
    subm_nest=subparser.add_parser('nest',help='Perform estimation of gene essentiality.')

    reqgroup=subm_nest.add_argument_group(title='Required arguments',
    description='')
    reqgroup.add_argument('-k','--count_table',required=True,help='Provide a tab-separated count table. Each line in the table should include sgRNA name (1st column), target gene (2nd column) and read counts in each sample.')
    reqgroup.add_argument('-d','--design_matrix',required=True,help='Provide a design matrix, either a quoted string of the design matrix or a file name. If using quoted string, for instance, "1,0;1,1" can be used for 2-sample conditions like 0_day and 4_weeks, and --include-samples should be specified as "0_day,4_weeks". If a file is given, the row of the design matrix must match the order of the samples in the count table, or the order of the samples by the --include-samples option.')

    iogroup=subm_nest.add_argument_group(title='Optional arguments for input and output',
    description='')
    iogroup.add_argument('-n','--output-prefix',default=None,help='The prefix of the output file(s).')
    iogroup.add_argument('-i', '--include-samples', help='Specify the sample labels if the design matrix is not given by file in the --design-matrix option. Sample labels are separated by ",", and must match the labels in the count table.')
    iogroup.add_argument('-b', '--beta-labels', help='Specify the labels of the variables (i.e., beta), if the design matrix is not given by file in the --design-matrix option. Should be separated by ",", and the number of labels must equal to (# columns of design matrix), including baseline labels. Default value: "bata_0,beta_1,beta_2,...".')

    genegroup=subm_nest.add_argument_group(title='Optional arguments for normalization.',
    description='"--norm-method control -e non_essential_genes_control" is recommended, which you have to specify negative control names, such as AAVS1')
    genegroup.add_argument("--norm-method",choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (nonormalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option). Defautl is median.')
    genegroup.add_argument("-e","--negative_control",help="The name of negative controls, such as AAVS1.",default=[])
    genegroup.add_argument('--genes-varmodeling',default="2000",help='The number of genes for mean-variance modeling. Default is 2000.')
    genegroup.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds). Defautl is FDR.')

    advancegroup=subm_nest.add_argument_group(title='Optional arguments for PPI incorporation and outliers removal',
    description='')
    advancegroup.add_argument("-o","--outliers_removal",action='store_true',help="Speicify whehter you want to remove outliers and recalculate..")
    advancegroup.add_argument("-p","--PPI_prior",action='store_true',help="Specify whether you want to incorporate PPI as prior")
    advancegroup.add_argument("-q","--QC_metric",action='store_true',help="Specify whether you want to derive quality control metrics")

    advancegroup=subm_nest.add_argument_group(title='Example',
    description='python3 mageck_nest.py nest -n mageck_nest_cell_line_A -i day_0,week_4 -k readcount_table.txt --norm-method control -e AAVS1 -d "1,0;1,1" -q')
    #--------------------------------------------------------------------------------
    args=parser.parse_args()

    if args.subcmd == None:
        parser.print_help()
        sys.exit(0)
    if args.output_prefix==None:
        args.output_prefix="mageck_nest_{}".format(args.count_table)
    if args.negative_control!=[]:
        args.negative_control=args.negative_control.split(',')
        args.negative_control=[i.upper() for i in args.negative_control]

    return args

def postargs(args):
    '''
    post-processing of argument parsing
    '''
    # configure logging information
    logging.basicConfig(level=10,
        format='%(levelname)-5s @ %(asctime)s.%(msecs)03d: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        # stream=sys.stderr,
        filename=args.output_prefix+'.log',
        filemode='w'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s.%(msecs)03d: %(message)s ','%a, %d %b %Y %H:%M:%S')
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    logging.info('Parameters: '+' '.join(sys.argv))

    from mledesignmat import parse_designmat

    try:
        import scipy
        from scipy.stats import nbinom
    except ImportError:
        logging.error('Cannot find scipy (required for mle approach). Please check your scipy installation.')
        sys.exit(-1)
    try:
        import numpy as np
        import numpy.linalg as linalg
    except ImportError:
        logging.error('Cannot find numpy (required for mle approach). Please check your numpy installation.')
        sys.exit(-1)
    # parsing design matrix
    (desmat,sampleid,betalabel)=parse_designmat(args.design_matrix)
    args.design_matrix=desmat

    # parsing sample label
    if sampleid ==None:
        # design matrix is provided as a string
        if args.include_samples !=None:
            args.include_samples=args.include_samples.split(',')
            if len(args.include_samples) != desmat.shape[0]:
                logging.error('The number of samples in the --include-samples option do not match rows in design matrix.')
                sys.exit(-1)
        if args.beta_labels!=None:
            args.beta_labels=args.beta_labels.split(',')
            if len(args.beta_labels) != desmat.shape[1]:
                logging.error('The number of labels in the --beta-labels option do not match columns in design matrix.')
                sys.exit(-1)
    else:
        # design matrix is provided as file
        if args.include_samples !=None:
            logging.error('Sample labels are included in the design matrix file '+args.design_matrix+'. The --include-samples option should not be used.')
            sys.exit(0)
        if args.beta_labels!=None:
            logging.error('Beta labels are included in the design matrix file '+args.design_matrix+'. The --beta-labels option should not be used.')
            sys.exit(0)
        args.include_samples=sampleid
        args.beta_labels=betalabel
        if len(args.include_samples) != desmat.shape[0]:
            logging.error('The number of samples in the --include-samples option do not match rows in design matrix.')
            sys.exit(-1)
        if len(args.beta_labels) != desmat.shape[1]:
            logging.error('The number of labels in the --beta-labels option do not match columns in design matrix.')
            sys.exit(-1)
    # log design matrix and column, row labels
    logging.info('Design matrix:')
    for desmat_1line in str(desmat).split('\n'):
        logging.info(desmat_1line)
    if args.beta_labels != None:
        logging.info('Beta labels:'+','.join(args.beta_labels))
    if args.include_samples != None:
        logging.info('Included samples:'+','.join(args.include_samples))

    return args

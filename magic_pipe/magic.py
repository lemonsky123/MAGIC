import argparse
from .utils import *


def require_arg(args_list):
    # warn the user to input the required args
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values not in args_list:
                msg = "Args of {} must be chosen from {}".format(f"-{self.dest}", ",".join(args_list))
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequireArgs


def require_argl():
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if type(values) != int:
                msg = "Args of {} must be integer".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            elif values < 1:
                msg = "Args of {} must be greater than or equal to 1".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            else:
                cpu_count = multiprocessing.cpu_count()
                if values > cpu_count:
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Local cpu count: {cpu_count}, "
                          f"processes user provide: {values} Adjusting to {cpu_count}")
                    # values = cpu_count
            setattr(args, self.dest, values)

    return RequireArgs


def _vcf_to_matrix_(args):
    from .vcf2mat import parse_vcf
    if args.e:
        if float(args.frq_low) > float(args.frq_up):
            raise ValueError(f"Lower limit for {args.e} can not be greater than upper limit!")
    parse_vcf(path_vcf=args.path_vcf, path_out=args.path_out, variant_type=args.variant, e=args.e, cadd=args.cadd,
              LoF=args.lof, frq_low=args.frq_low, frq_up=args.frq_up, path_corr=args.path_corr)


def _fisher_perm_skat_acat_(args):
    # fisher test and permutation
    from .fisher_test_perm import parse_snp_matrix, get_case_control_samples
    path_out_fisher = f"{args.o}.fisher"
    case_samples, control_samples = get_case_control_samples(args.c)
    parse_snp_matrix(path_matrix=args.p, case_samples=case_samples, control_samples=control_samples, method=args.m,
                     path_out=path_out_fisher, num_process=args.l, progress_bar=args.progress, hypothesis=args.s,
                     path_pathway=args.w)
    # SKAT and ACAT analysis
    from .skat_acat import get_case_control_samples_num, parse_covariate, read_in_snp_matrix, calculate_skat, read_in_snp_matrix_with_pathway
    path_out_skat = f"{args.o}.skat"
    case_samples_num, control_samples_num = get_case_control_samples_num(args.c)
    if args.X:
        covariate = parse_covariate(args.X, case_samples_num, control_samples_num, args.c)
    else:
        covariate = None
    if not args.w:
        genes_snp_genotypes = read_in_snp_matrix(args.p)
    else:
        # if pathway is provided, the keys of genes_snp_genotypes are pathway names, not genes
        genes_snp_genotypes = read_in_snp_matrix_with_pathway(path_matrix=args.p, path_pathway=args.w)
    # calculate_skat(genes_snp_genotypes, case_samples_num, control_samples_num, args.o, covariate)
    calculate_skat(genes_snp_genotypes=genes_snp_genotypes, case_samples_num=case_samples_num,
                   control_samples_num=control_samples_num, path_out=path_out_skat, X=covariate)
    # merge results of fisher test and SKAT into one file
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_out_fisher}"))
    fisher_header = []
    gene_fisher_skat_results_items = defaultdict(list)
    with open(path_out_fisher) as fisher:
        for line in fisher:
            items = line.strip().split("\t")
            if line.startswith("gene"):
                fisher_header = items[1:]
            else:
                gene_fisher_skat_results_items[items[0]] = items[1:]
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_out_skat}"))
    skat_header = []
    with open(path_out_skat) as skat:
        for line in skat:
            items = line.strip().split("\t")
            if line.startswith("gene"):
                skat_header = items[1:]
            else:
                gene_fisher_skat_results_items[items[0]].extend(items[1:])
    column_nums = len(fisher_header) + len(skat_header)
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {args.o}"))
    with open(args.o, "w+") as out:
        out.write("gene\t" + "\t".join(fisher_header + skat_header) + "\n")
        for each_gene in gene_fisher_skat_results_items:
            this_gene_fisher_skat_result = gene_fisher_skat_results_items[each_gene]
            if len(this_gene_fisher_skat_result) != column_nums:
                print(yellow(f"[WARNING]: {each_gene} has no fisher or skat results, please "
                             f"check! It will not be written to {args.o}"))
            else:
                out.write(each_gene + "\t" + "\t".join(gene_fisher_skat_results_items[each_gene]) + "\n")


def _logistic_(args):
    from .logistic import parse_matrix, parse_matrix_without_gene_set
    if args.path_group and args.gene_type:
        parse_matrix(path_matrix=args.p, path_out=args.o, path_pathway=args.w, path_group=args.path_group,
                     gene_type=args.gene_type)
    else:
        parse_matrix_without_gene_set(path_matrix=args.p, path_out=args.o)


def _scPBS_(args):
    from .scPBS import parse_sparse_matrix, parse_vcf2mat, logistic
    genes_with_exp_gt_02, data_divide = parse_sparse_matrix(path_h5ad=args.path_h5ad)
    samples_variants_count_in_cell_mat = parse_vcf2mat(path_vcf2mat=args.path_vcf2mat, path_group=args.path_group,
                                                       genes_with_exp_gt_02=genes_with_exp_gt_02)
    logistic(samples_variants_count_in_cell_mat=samples_variants_count_in_cell_mat, progress=args.progress,
             data_divide=data_divide, path_out=args.path_out)


def _seurat_to_h5ad_(args):
    seurat_to_h5ad(path_rds=args.path_rds, path_h5ad=args.path_h5ad, path_rscript=args.path_rscript)


def _seurat_plot_cluster_(args):
    seurat_plot_cluster(path_rds=args.path_rds, path_cor=args.path_cor, path_pdf=args.path_pdf,
                        path_rscript=args.path_rscript)


def main():
    parser = argparse.ArgumentParser(prog='magic',
                                     description="A toolkit for WES/WGS data analysis\n"
                                                 "magic toolkits have five main tools: \n"
                                                 "\t%(prog)s vcf2mat: Extract biological features and information"
                                                 " from the vcf file generated by GATK\n"
                                                 "\t%(prog)s RVAS: Perform fisher test, permutation test, SKAT "
                                                 "and ACAT analysis for genes extracted from the vcf2mat command\n"
                                                 "\t%(prog)s logistic: Run logistic on gene sample matrix\n"
                                     "\t%(prog)s scPBS: single-cell Polygenic Burden Score\n"
                                     "\t%(prog)s rds2h5ad: Convert Seurat rds to h5ad\n"
                                     "\t%(prog)s seurat_plot: Feature plot for correlation results\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(title="modules", help='magic modules, use -h/--help for help')
    sub_vcf_to_matrix = subparsers.add_parser("vcf2mat", description="Convert the vcf file to a matrix, "
                                                                     "each row stands for a snp")
    sub_fisher_skat = subparsers.add_parser("RVAS", description="Perform fisher test, permutation test and SKAT "
                                                                "and ACAT analysisfor each gene generated by vcf2mat")
    sub_logistic = subparsers.add_parser("logistic", description="Run logistic on gene sample matrix")
    sub_scPBS = subparsers.add_parser("scPBS",
                                      description="single-cell Polygenic Burden Score")
    sub_seurat_to_h5ad = subparsers.add_parser("rds2h5ad", description="Convert Seurat rds to h5ad")
    sub_seurat_plot_cluster = subparsers.add_parser("seurat_plot", description="Feature plot for correlation results")
    # vcf2mat
    sub_vcf_to_matrix.add_argument("-p", help="Path to reference vcf file", action="store", dest="path_vcf", type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-o", help="Path to output file", action="store", dest="path_out", type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-e", help="Variable to filter, must be chosen from MAX_AF, MAX_AF_POPS, AFR_AF, "
                                              "AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, gnomAD_AF, gnomAD_AFR_AF, "
                                              "gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, "
                                              "gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, default is not given. "
                                              "If not given, won't use the database to filter variants", action="store",
                                   dest="e", type=str, default="",
                                   choices=["MAX_AF", "MAX_AF_POPS", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
                                            "AA_AF",
                                            "EA_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "",
                                            "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF",
                                            "gnomAD_SAS_AF"])
    sub_vcf_to_matrix.add_argument("-f", help="Lower limit for -e, a float number to filter the variants by database, "
                                              "default is negative infinity. eg. when -f 0, it means that if gnomAD_AF "
                                              "> 0, then the transcript is output", action="store", dest="frq_low",
                                   type=str, default="-inf")
    sub_vcf_to_matrix.add_argument("-F", help="Upper limit for -e, a float number to filter the variants by database, "
                                              "default is positive infinity. eg. when -F 0.05, it means that if "
                                              "gnomAD_AF <= 0.05, then the transcript is output", action="store",
                                   dest="frq_up", type=str, default="inf")
    sub_vcf_to_matrix.add_argument("-v", help="Variant types to filter the vcf, must be chosen from lof/"
                                              "missense_damage/missense_benign/synonymous/3UTR/5UTR, default is lof",
                                   action="store", dest="variant", type=str, default="lof",
                                   choices=["lof", "missense_damage", "missense_benign", "synonymous", "3UTR", "5UTR"])
    sub_vcf_to_matrix.add_argument("-l", help="Chosen from HC/LC, default is not given, if given, HC or LC WON'T be "
                                              "used to filter for lof variant", action="store", dest="lof", type=str,
                                   choices=["HC", "LC", ""], default="")
    sub_vcf_to_matrix.add_argument("-s", help="Threshold for CADD_PHRED, used in missense_damage/missense_benign"
                                              "/synonymous, default is 15", action="store", dest="cadd", type=str,
                                   default=15)
    sub_vcf_to_matrix.add_argument("-c", help="Path to put corrected variants(AF > 0.05, 0/0 to 1/1 and 1/1 to 0/0)",
                                   action="store", dest="path_corr", type=str, required=True)
    sub_vcf_to_matrix.set_defaults(func=_vcf_to_matrix_)
    # fisher test and permutation test
    sub_fisher_skat.add_argument("-p", help="Path to input matrix with columns representing samples and rows "
                                            "representing snps by vcf2mat", action="store", dest="p", type=str,
                                 required=True)
    sub_fisher_skat.add_argument("-c", help="Path to input file with case and control samples, must be two columns, "
                                            "first col for samples, the other for case/control", action="store",
                                 dest="c", type=str, required=True)
    sub_fisher_skat.add_argument("-m", help="Method to use in the adjustment of p values, must be chosen from "
                                            "bonferroni/sidak/holm-sidak/holm/simes-hochberg/hommel/fdr_bh/fdr_by/"
                                            "fdr_tsbh/fdr_tsbky, default is fdr_bh",
                                 action=require_arg(["bonferroni", "sidak", "holm-sidak", "holm",
                                                     "simes-hochberg", "hommel", "fdr_bh", "fdr_by", "fdr_tsbh",
                                                     "fdr_tsbky"]),
                                 type=str, default="fdr_bh")
    sub_fisher_skat.add_argument("-l", help="Processes to use in the fisher test and permutation step, default is 10",
                                 action=require_argl(), type=int, default=10)
    sub_fisher_skat.add_argument("-progress", help="Whether to show the progress bar, default is off",
                                 action="store_true", default=False)
    sub_fisher_skat.add_argument("-s", help="Alternative hypothesis for fisher exact test, must be chosen from "
                                            "two-sided/less/greater, less and greater is one-sided",
                                 action=require_arg(["two-sided", "less", "greater"]), type=str, default="two-sided")
    sub_fisher_skat.add_argument("-w",
                                 help="Path to pathway table, each row stands for a pathway, the 1st and 2nd col are "
                                      "pathway name while the other cols are genes in this pathway. e.g. "
                                      "pathway\tpathway\tgene1\tgene2...",
                                 action="store", default=None)
    # SKAT and ACAT
    sub_fisher_skat.add_argument("-X", help="Path to covariates file, default id not provided", action="store",
                                 dest="X", type=str, default=None)
    sub_fisher_skat.add_argument("-o", help="Path to put output matrix with the following information: "
                                            "gene\tcase\tcase_wild\tcontrol\tcontrol_wild\tp_value\tp_adj"
                                            "\tp_permutation\tburden_p1\tburden_p2\tskat_p1\tskat_p2\tacatv_p1"
                                            "\tacatv_p2\tskato_p1\tskato_p2\tacato_p", action="store", dest="o",
                                 type=str, required=True)
    sub_fisher_skat.set_defaults(func=_fisher_perm_skat_acat_)
    # logistic
    sub_logistic.add_argument("-p",
                              help="Path to input matrix with columns representing samples and rows representing snps "
                                   "by vcf2mat", action="store", dest="p", type=str, required=True)
    sub_logistic.add_argument("-o", help="Path to put output file with two columns, one for sample and the other for "
                                         "variant nums", action="store", dest="o", type=str, required=True)
    sub_logistic.add_argument("-w", help="Path to a pathway file in gmt format", action="store", dest="w", type=str,
                              default=None)
    sub_logistic.add_argument("-g", help="Path to group file with two columns, one for sample name, another for sample "
                                         "type, case/control", action="store", dest="path_group", type=str, default="")
    sub_logistic.add_argument("-t", help="Type of gene name, ensemble or symbol", action="store", dest="gene_type",
                              type=str, choices=["ensemble", "symbol", ""], default="")
    sub_logistic.set_defaults(func=_logistic_)
    # scPBS
    sub_scPBS.add_argument("-p", help="Path to input matrix with columns representing samples and rows representing "
                                      "snps by vcf2mat", action="store", dest="path_vcf2mat", type=str, required=True)
    sub_scPBS.add_argument("-g", help="Path to group file with two columns, one for sample name, another for sample "
                                      "type, case/control", action="store", dest="path_group", type=str, required=True)
    sub_scPBS.add_argument("-progress", help="Whether to show the progress bar, default is off",
                           action="store_true", dest="progress", default=False)
    sub_scPBS.add_argument("-b", help="Path to h5ad format file with gene expression in each cell, can be converted "
                                      "from Seurat", action="store", dest="path_h5ad", default=False)
    sub_scPBS.add_argument("-o", help="Path to put correlation result of odds ratio and gene expression",
                           action="store", dest="path_out", default=False)
    sub_scPBS.set_defaults(func=_scPBS_)
    # seurat to h5ad
    sub_seurat_to_h5ad.add_argument("-r", help="Path to rds file saved from Seurat", action="store", dest="path_rds",
                                    type=str, required=True)
    sub_seurat_to_h5ad.add_argument("-o", help="Path to put converted h5ad file", action="store", dest="path_h5ad",
                                    type=str, required=True)
    sub_seurat_to_h5ad.add_argument("-rscript", help="Path to Rscript executable", action="store", dest="path_rscript",
                                    type=str, default="Rscript")
    sub_seurat_to_h5ad.set_defaults(func=_seurat_to_h5ad_)
    # seurat plot cluster
    sub_seurat_plot_cluster.add_argument("-r", help="Path to rds file saved from Seurat", action="store",
                                         dest="path_rds", type=str, required=True)
    sub_seurat_plot_cluster.add_argument("-c", help="Path to correlation result of odds ratio and gene expression",
                                         action="store", dest="path_cor", required=True)
    sub_seurat_plot_cluster.add_argument("-o", help="Path to put pdf of feature plot", action="store", dest="path_pdf",
                                         required=True)
    sub_seurat_plot_cluster.add_argument("-rscript", help="Path to Rscript executable", action="store",
                                         dest="path_rscript", type=str, default="Rscript")
    sub_seurat_plot_cluster.set_defaults(func=_seurat_plot_cluster_)
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()

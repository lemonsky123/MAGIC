import logging
import rpy2.robjects.numpy2ri as rpyn
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import IntVector, Formula, DataFrame, NULL, pandas2ri, default_converter, conversion
from .utils import *


class FilterWarningMessage(logging.Filter):
    def filter(self, record):
        return not re.search(r"Warning messages:", record.getMessage()) or len(record.getMessage().strip()) != 0


logging.basicConfig(
    level=logging.WARNING,
    format="[%(levelname)s FROM R]: %(message)s",
)
rpy2_logger.addFilter(FilterWarningMessage())


base = importr("base")
skat = importr("SKAT")
acat = importr("ACAT")


def get_case_control_samples_num(path_samples):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_samples}"))
    case_samples = []
    control_samples = []
    for each_item in csv.DictReader(open(path_samples), delimiter="\t"):
        if each_item["type"] == "case":
            case_samples.append(each_item["sample"])
        elif each_item["type"] == "control":
            control_samples.append(each_item["sample"])
        else:
            raise InvalidSampleTypeError(f"Sample types must be case/control, got {each_item['type']}")
    return len(case_samples), len(control_samples)


def read_in_snp_matrix(path_matrix):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_matrix}"))
    genes_snp_genotypes = defaultdict(dict)
    with open(path_matrix) as file:
        for line in file:
            if line.startswith("variant"):
                continue
            items = line.strip().split("\t")
            snp = items[0]
            gene = items[1]
            # replace . with 0
            # genotype = IntVector(list(map(lambda k: 0 if k == "." else int(k), items[7:])))
            # replace . with 9
            genotype = IntVector(list(map(lambda k: 9 if k == "." else int(k), items[7:])))
            # old style step1, samples starts from the 6th
            # genotype = IntVector(list(map(lambda k: 0 if k == "." else int(k), items[5:])))
            genes_snp_genotypes[gene][snp] = genotype
    return genes_snp_genotypes


def read_in_snp_matrix_with_pathway(path_matrix, path_pathway):
    pathway_snp_genotypes = defaultdict(dict)
    genes_snp_genotypes = read_in_snp_matrix(path_matrix=path_matrix)
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_pathway}"))
    with open(path_pathway, encoding='utf-8') as file:
        for line in file:
            items = line.strip().split("\t")
            pathway = items[0]
            this_pathway_genes = items[2:]
            for each_gene in this_pathway_genes:
                if each_gene in genes_snp_genotypes:
                    this_gene_genotypes = genes_snp_genotypes[each_gene]
                    # Each snp in the vcf file is unique
                    pathway_snp_genotypes[pathway].update(this_gene_genotypes)
    return pathway_snp_genotypes


def parse_covariate(path_covariate, case_samples_num, control_samples_num, path_samples):
    print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Covariate file is provided, reading from {path_covariate}"))
    cov_data = pd.read_csv(path_covariate, sep="\t", header=None, index_col=0, engine="c", low_memory=False)
    with conversion.localconverter(default_converter + pandas2ri.converter):
        cov_data_r = conversion.py2rpy(cov_data)
    cov_data_mat = base.as_matrix(cov_data_r)
    if cov_data_mat.nrow != case_samples_num + control_samples_num:
        raise SamplesNotEqualError(f"samples in {path_covariate} not equal to those in {path_samples}")
    return cov_data_mat


def calculate_skat(genes_snp_genotypes, case_samples_num, control_samples_num, path_out, X):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start performing SKAT and ACAT analysis"))
    if not os.path.exists(os.path.dirname(path_out)):
        print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Creating directory {os.path.dirname(path_out)}"))
        if os.path.dirname(path_out):
            os.makedirs(os.path.dirname(path_out))
    trait = IntVector([1] * case_samples_num + [0] * control_samples_num)
    if X:
        print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Constructing NULL model with covariates"))
        fmla = Formula("trait ~ X")
        env = fmla.environment
        env['trait'] = trait
        env['X'] = X
        obj = skat.SKAT_Null_Model(fmla, out_type="D")
        obj_acat = acat.NULL_Model(trait, Z=X)
    else:
        print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Constructing NULL model without covariates"))
        fmla = Formula("trait ~ 1")
        env = fmla.environment
        env['trait'] = trait
        env['1'] = 1
        obj = skat.SKAT_Null_Model(fmla, out_type="D")
        obj_acat = acat.NULL_Model(trait, Z=NULL)
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {path_out}"))
    with open(path_out, "w+") as out:
        out.write("gene\tburden_p1\tburden_p2\tskat_p1\tskat_p2\tacatv_p1\tacatv_p2\tskato_p1\tskato_p2\tacato_p\n")
        for each_gene in genes_snp_genotypes:
            this_gene_snp_genotypes = DataFrame(genes_snp_genotypes[each_gene])
            this_gene_snp_genotypes_mat = base.as_matrix(this_gene_snp_genotypes)
            burden_pvalue = base.c(
                skat.SKAT(this_gene_snp_genotypes_mat, obj, r_corr=1, weights_beta=IntVector([1, 1])).rx2("p.value"),
                skat.SKAT(this_gene_snp_genotypes_mat, obj, r_corr=1, weights_beta=IntVector([1, 25])).rx2("p.value"))
            skat_pvalue = base.c(
                skat.SKAT(this_gene_snp_genotypes_mat, obj, r_corr=0, weights_beta=IntVector([1, 1])).rx2("p.value"),
                skat.SKAT(this_gene_snp_genotypes_mat, obj, r_corr=0, weights_beta=IntVector([1, 25])).rx2("p.value"))
            acatv_pvalue = base.c(acat.ACAT_V(this_gene_snp_genotypes_mat, obj_acat, weights_beta=IntVector([1, 1])),
                                  acat.ACAT_V(this_gene_snp_genotypes_mat, obj_acat, weights_beta=IntVector([1, 25])))
            skato_pvalue = base.c(
                skat.SKAT(this_gene_snp_genotypes_mat,
                          obj, method="optimal.adj", weights_beta=IntVector([1, 1])).rx2("p.value"),
                skat.SKAT(this_gene_snp_genotypes_mat,
                          obj, method="optimal.adj", weights_beta=IntVector([1, 25])).rx2("p.value"))
            acato_pvalue = acat.ACAT(base.c(burden_pvalue, skat_pvalue, acatv_pvalue))
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                each_gene,
                burden_pvalue[0],
                burden_pvalue[1],
                skat_pvalue[0],
                skat_pvalue[1],
                acatv_pvalue[0],
                acatv_pvalue[1],
                skato_pvalue[0],
                skato_pvalue[1],
                rpyn.rpy2py(acato_pvalue)[0]
            ))
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Finish writing to {path_out}"))

from .utils import *


warnings.filterwarnings("ignore")


def parse_group_file(path_group):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_group}"))
    total_sample_type = OrderedDict()
    for each_item in csv.DictReader(open(path_group), delimiter="\t"):
        if each_item["type"] == "case":
            total_sample_type[each_item["sample"]] = 1
        elif each_item["type"] == "control":
            total_sample_type[each_item["sample"]] = 0
    return total_sample_type


def parse_sparse_matrix(path_h5ad):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from  {path_h5ad}"))
    # sparse_matrix = io.mmread(path_matrix)
    # mat_dense = sparse_matrix.toarray()
    # var_names = np.genfromtxt(path_rownames, dtype=str)
    # col_names = np.genfromtxt(path_colnames, dtype=str)
    data_raw = read_h5ad(path_h5ad)
    # the matrix in h5ad from Seurat is cell x gene, we must convert it to gene x cell
    data_mat = data_raw.to_df().T
    print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Drop rows with all zeros"))
    data_mat = data_mat.loc[(data_mat != 0).any(axis=1)]
    # data_mat = pd.DataFrame(mat_dense, columns=col_names, index=var_names)
    data_exp = dict()
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Divide the exp of each cell by the max "
                f"exp in each gene"))
    for each_gene in data_mat.index:
        this_gene_exp_values = data_mat.loc[[each_gene]].values
        # if not this_gene_exp_values.any():
        #    print(yellow(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Expression of {each_gene} is 0 in"
        #                 f" all cells, drop it"))
        #    continue
        this_gene_exp_values_divide_max = this_gene_exp_values / np.max(this_gene_exp_values, axis=1)[:, None]
        data_exp[each_gene] = this_gene_exp_values_divide_max[0]
    data_divide = pd.DataFrame.from_dict(data_exp, orient='index')
    data_divide.columns = data_mat.columns
    # genes with (exp / max) > 0.2
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Get genes with (exp / max exp) > 0.2"))
    genes_with_exp_gt_02 = defaultdict(set)
    for each_cell in data_divide.columns:
        this_cell_genes_gt_02 = data_divide[each_cell].index[data_divide[each_cell].values > 0.2]
        for each_gene in this_cell_genes_gt_02:
            genes_with_exp_gt_02[each_gene].add(each_cell)
    return genes_with_exp_gt_02, data_divide


def parse_vcf2mat(path_vcf2mat, path_group, genes_with_exp_gt_02):
    total_sample_type = parse_group_file(path_group=path_group)
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_vcf2mat}"))
    samples = []
    total_snp = defaultdict(lambda: 0)
    samples_variants_count_in_cell = defaultdict(lambda: defaultdict(lambda: 0))
    with open(path_vcf2mat) as file:
        for line in file:
            items = line.strip().split("\t")
            snp = items[0]
            # if one variant occurs more than once, keep only the first occurred one
            if total_snp[snp] != 0:
                continue
            else:
                total_snp[snp] += 1
            if line.startswith("variant"):
                samples = items[7:]
            elif not samples:
                raise ValueError(f"Sample names not found in {path_vcf2mat}")
            else:
                this_gene_cells = genes_with_exp_gt_02[items[1]]
                if this_gene_cells:
                    for each_sample_genotype in zip(samples, items[7:]):
                        if each_sample_genotype[1] == "1" or each_sample_genotype[1] == "2":
                            for each_cell in this_gene_cells:
                                samples_variants_count_in_cell[each_sample_genotype[0]][each_cell] += 1
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Add group label to each sample"))
    # samples_with_variant = list(samples_variants_count_in_cell.keys())
    for each_sample in samples:
        this_sample_group = total_sample_type[each_sample]
        samples_variants_count_in_cell[each_sample]["group"] = this_sample_group
    # convert the sample_cell_variant count dict to pandas dataframe, and fill NaN with 0
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Convert sample cell variants count matrix "
                f"to pandas dataframe"))
    samples_variants_count_in_cell_mat = pd.DataFrame.from_dict(samples_variants_count_in_cell, orient='index').fillna(0)
    return samples_variants_count_in_cell_mat


def logistic(samples_variants_count_in_cell_mat: pd.DataFrame, progress: bool, data_divide: pd.DataFrame,
             path_out: str):
    cell_names = samples_variants_count_in_cell_mat.columns.tolist()
    # if there are "-" in cell names, replace them with "_" for the logistic analysis. smf.logit will throw an exception
    # if there are "-" in cell names.
    cell_names_sub = list(map(lambda k: re.sub("-", "_", k), cell_names))
    cell_odds = dict()
    samples_variants_count_in_cell_mat.columns = cell_names_sub

    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Run logistic on all cells"))
    if progress:
        pbar = tqdm(total=len(cell_names_sub))
    for each_cell in cell_names_sub:
        if each_cell == "group":
            continue
        try:
            model = smf.logit(f"group ~ {each_cell}", data=samples_variants_count_in_cell_mat).fit(disp=0)
        except Exception:
            model = smf.logit(f"group ~ {each_cell}", data=samples_variants_count_in_cell_mat).fit_regularized(
                method='l1', alpha=0.1, disp=0)
        model_odds = pd.DataFrame(np.exp(model.params), columns=['OR'])
        model_odds['z-value'] = model.pvalues
        model_odds[['2.5%', '97.5%']] = np.exp(model.conf_int())
        cell_odds[each_cell] = model_odds.loc[each_cell, "OR"]
        if progress:
            pbar.update()
    if progress:
        pbar.close()

    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Calculate Correlation between odds ratio "
                f"and gene expression"))
    cells = data_divide.columns
    cells_sub = list(map(lambda k: re.sub("-", "_", k), cells))
    if len(cells) + 1 != len(cell_names):
        raise SamplesNotEqualError(f"Cell nums should be the same after stating variants for each sample")
    gene_corr_pvalue = defaultdict(dict)
    odds_ratio_per_cell = [cell_odds[each_cell] for each_cell in cells_sub]
    genes = data_divide.index
    if progress:
        pbar = tqdm(total=len(genes))
    for each_gene in genes:
        this_gene_dict = data_divide.loc[each_gene].to_dict()
        this_gene_items = [this_gene_dict[each_cell] for each_cell in cells]
        rho = pearsonr(this_gene_items, odds_ratio_per_cell)
        gene_corr_pvalue[each_gene]["correlation"] = rho.correlation
        gene_corr_pvalue[each_gene]["pvalue"] = rho.pvalue
        if progress:
            pbar.update()
    if progress:
        pbar.close()

    gene_corr_pvalue_sorted = sorted(gene_corr_pvalue.items(), key=lambda k: k[1]["correlation"], reverse=True)
    with open(path_out, "w+") as out:
        out.write("gene\todds ratio\tpvalue\n")
        for each_item in gene_corr_pvalue_sorted:
            out.write("{}\t{}\t{}\n".format(each_item[0], each_item[1]["correlation"], each_item[1]["pvalue"]))

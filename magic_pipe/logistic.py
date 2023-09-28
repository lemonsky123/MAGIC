from .utils import *


def parse_group_file(path_group):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_group}"))
    total_sample_type = OrderedDict()
    for each_item in csv.DictReader(open(path_group), delimiter="\t"):
        if each_item["type"] == "case":
            total_sample_type[each_item["sample"]] = 1
        elif each_item["type"] == "control":
            total_sample_type[each_item["sample"]] = 0
        else:
            raise ValueError(f"Sample type must be case or control, got {each_item['type']}")
    return total_sample_type


def parse_pathway_file(path_pathway):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_pathway}"))
    total_set_names = []
    gene_sets = defaultdict(set)
    with open(path_pathway) as file:
        for line in file:
            items = line.strip().split("\t")
            pathway = items[0]
            genes = items[2:]
            total_set_names.append(pathway)
            for each_gene in genes:
                gene_sets[each_gene].add(pathway)
    total_set_names = list(dict.fromkeys(total_set_names))
    return gene_sets, total_set_names


def parse_matrix(path_matrix, path_out, path_pathway, path_group, gene_type):
    total_sample_type = parse_group_file(path_group=path_group)
    gene_sets, total_set_names = parse_pathway_file(path_pathway=path_pathway)
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_matrix}"))
    samples = []
    total_snp = defaultdict(lambda: 0)
    samples_variants_count = defaultdict(lambda: defaultdict(lambda: 0))
    with open(path_matrix) as file:
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
                raise ValueError(f"Sample names not found in {path_matrix}")
            else:
                # this_gene_sets = gene_sets[items[1]]
                if gene_type == "symbol":
                    this_gene_sets = gene_sets[items[1]]
                elif gene_type == "ensemble":
                    this_gene_sets = gene_sets[items[2]]
                if this_gene_sets:
                    for each_sample_genotype in zip(samples, items[7:]):
                        if each_sample_genotype[1] == "1" or each_sample_genotype[1] == "2":
                            for each_set in this_gene_sets:
                                samples_variants_count[each_sample_genotype[0]][each_set] += 1
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {path_out}"))
    with open(path_out, "w+") as out:
        out.write("sample\t{}\tstatus\n".format("\t".join(total_set_names)))
        for each_sample in total_sample_type:
            this_sample_variant_nums = [str(samples_variants_count[each_sample][each_set]) for each_set in total_set_names]
            out.write("{}\t{}\t{}\n".format(each_sample, "\t".join(this_sample_variant_nums), total_sample_type[each_sample]))


def parse_matrix_without_gene_set(path_matrix, path_out):
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_matrix} without gene sets and group names"))
    samples = []
    total_snp = defaultdict(lambda: 0)
    samples_variants_count = defaultdict(lambda: 0)
    with open(path_matrix) as file:
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
                raise ValueError(f"Sample names not found in {path_matrix}")
            else:
                for each_sample_genotype in zip(samples, items[7:]):
                    if each_sample_genotype[1] == "1" or each_sample_genotype[1] == "2":
                        samples_variants_count[each_sample_genotype[0]] += 1
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {path_out}"))
    with open(path_out, "w+") as out:
        for each_sample in samples:
            out.write("{}\t{}\n".format(
                each_sample,
                samples_variants_count[each_sample]
            ))

from .utils import *


def format_genotype(genotype, reverse=False):
    genotype_items = genotype.split(":")
    if not reverse:
        if genotype_items[0] == "0/0":
            return "0"
        elif genotype_items[0] == "1/1":
            return "2"
        elif genotype_items[0] == "0/1" or genotype_items[0] == "1/0":
            return "1"
        else:
            return "."
    else:
        if genotype_items[0] == "0/0":
            # replace 0/0 to 1/1
            return "2"
        elif genotype_items[0] == "1/1":
            # replace 1/1 to 0/0
            return "0"
        elif genotype_items[0] == "0/1" or genotype_items[0] == "1/0":
            return "1"
        else:
            return "."


def parse_vcf(path_vcf, path_out, variant_type, e, cadd, LoF, frq_low, frq_up, path_corr):
    formats = []
    sample_names = []
    snp_out_format = ["0/0", "0/1", "1/1", "1/0", "./."]
    pattern_CSQ = re.compile(r";CSQ=(.*?)$")
    pattern_AC = re.compile(r"AC=(.*?);AF=(.*?);AN=(.*?);")
    missense_types = ["inframe_deletion", "inframe_insertion", "missense_variant", "stop_lost"]
    lof_types = ["stop_gained", "frameshift_variant", "start_lost", "splice_acceptor_variant", "splice_donor_variant"]
    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_vcf}"))
    out = open(path_out, "w+")
    corr = open(path_corr, "w+")
    with open(path_vcf) as file:
        for line in file:
            if line.startswith("#"):
                if re.search(r'Format: (.*?)">', line):
                    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start parsing format line"))
                    formats = re.findall(r'Format: (.*?)">', line)[0].split("|")
                elif line.startswith("#CHROM"):
                    print(green(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start collecting sample names"))
                    sample_names = line.strip().split("\t")[9:]
                    out.write("variant\tsymbol\tENSG\tENST\tsample_stat\tvariant_stat\tinfo\t" + "\t".join(sample_names) + "\n")
            elif not formats:
                raise ValueError(f"Format line not found in {path_vcf}")
            elif not sample_names:
                raise ValueError(f"Sample name not found in {path_vcf}")
            else:
                items = line.strip().split("\t")
                # print(pattern_AC.findall(items[7]))
                ac, af, an = list(map(float, pattern_AC.findall(items[7])[0]))
                samples_info_items = list(map(lambda k: k.split(":")[0], items[9:]))
                # count the frequency of "0/0", "1/1", "1/0", "0/1" and "./."
                snp_type_counts = Counter(samples_info_items)
                snp_type_out = []
                for each_snp in snp_out_format:
                    if each_snp in snp_type_counts:
                        snp_type_out.append(str(snp_type_counts[each_snp]))
                    else:
                        snp_type_out.append("0")
                # whether to output the corrected variants, if true, the variant will be output
                whether_output_corrected = False
                if af <= 0.05:
                    genotypes = list(map(format_genotype, items[9:]))
                else:
                    whether_output_corrected = True
                    genotypes = list(map(lambda k: format_genotype(k, reverse=True), items[9:]))
                transcripts = pattern_CSQ.findall(items[7])[0].split(",")
                if len(transcripts) > 1:
                    raise ValueError(f"More than one transcript detected in {items[7]}")
                this_snp_trans_type = ""
                each_transcript = transcripts[0]
                this_transcript_items = each_transcript.split("|")
                # replace the empty string in this_transcript_items with 0
                this_transcript_format_dict_rep0 = dict(zip(formats, list(map(lambda k: k if k else 0, this_transcript_items))))
                genom_af = this_transcript_format_dict_rep0["gnomAD_AF"]
                if any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, lof_types))):
                    if LoF and this_transcript_format_dict_rep0["LoF"] == LoF:
                        this_snp_trans_type = "lof"
                    else:
                        this_snp_trans_type = "lof"
                elif any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, missense_types))):
                    if (re.search(r"^deleterious\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^probably_damaging\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"]))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) > cadd:
                        this_snp_trans_type = "missense_damage"
                    if (re.search(r"^tolerated\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^benign\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"]))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) < cadd:
                        this_snp_trans_type = "missense_benign"
                if this_transcript_format_dict_rep0["Consequence"] == "synonymous_variant" and float(this_transcript_format_dict_rep0["CADD_PHRED"]) < cadd:
                    this_snp_trans_type = "synonymous"
                # if not float(this_transcript_format_dict_rep0[e]) < frequency:
                #    continue
                # variant_stat AC/AF/AN/gnomad_AF
                this_snp_variant_stat = f"{ac}/{af}/{an}/{genom_af}"
                # whether to output this variant
                whether_output = False
                if not e:
                    whether_output = True
                elif float(frq_low) < float(this_transcript_format_dict_rep0[e]) <= float(frq_up):
                    whether_output = True
                if whether_output:
                    if variant_type == this_snp_trans_type:
                        if whether_output_corrected:
                            corr.write(line)
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            f"{items[0]}_{items[1]}_{items[3]}_{items[4]}",
                            this_transcript_format_dict_rep0['SYMBOL'],
                            this_transcript_format_dict_rep0['Gene'],
                            this_transcript_format_dict_rep0['Feature'],
                            "/".join(snp_type_out),
                            this_snp_variant_stat,
                            each_transcript,
                            "\t".join(genotypes)
                        ))
    out.close()
    corr.close()

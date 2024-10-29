import polars
import numpy as np
import sys
from pathlib import Path

# Paths
project_dir = Path("/data/teamgdansk/mwaleron/carmen-analysis")

data_dir = project_dir.joinpath("data")
temp_dir = project_dir.joinpath("temp")

data_subs_dir = data_dir.joinpath("subsidiary-files")
data_pub_dir = data_dir.joinpath("to-be-published")

gencode_data_dir = data_subs_dir.joinpath("GENCODE")
gencodev39_stripped = gencode_data_dir.joinpath("gencode.v39.annotation.stripped.versions.parquet")

contig_scaffold_list_file = data_pub_dir.joinpath("scaff_all_expanded.tsv")
pogo_annotations_file = data_dir.joinpath("carmen-mapped-protein-annotations-pogo.parquet")

genestats_file = temp_dir.joinpath(f"stats/{sys.argv[1]}_stats.parquet")
exonstats_file = temp_dir.joinpath(f"exons/{sys.argv[1]}_exons.parquet")
nonexonstats_file = temp_dir.joinpath(f"nonexons/{sys.argv[1]}_nonexons.parquet")

def find_contiguous_true_ranges(seq):
    contiguous_true_ranges = []
    start_idx = None

    for i in range(len(seq)):
        if seq[i] == True:
            if start_idx is None:
                start_idx = i
        elif start_idx is not None:
            contiguous_true_ranges.append((start_idx, i))
            start_idx = None

    # Check if there's an open range at the end
    if start_idx is not None:
        contiguous_true_ranges.append([start_idx, len(seq)])

    return contiguous_true_ranges

pogo_annotations = polars.scan_parquet(pogo_annotations_file).with_columns(Chromosome = polars.lit("chr")+polars.col("Chromosome"))
gencodev39 = polars.scan_parquet(gencodev39_stripped)
contig_scaffold_list = polars.scan_csv(contig_scaffold_list_file, separator="\t")

target_gene = gencodev39.filter(
    polars.col("gene_id") == sys.argv[1],
    polars.col("feature") == "gene"
).collect()

target_exons = gencodev39.filter(
    polars.col("gene_id") == sys.argv[1],
    polars.col("feature") == "exon"
).collect()

target_peptides = pogo_annotations.filter(
    polars.col("Chromosome") == target_gene["seqname"],
    ( polars.col("Gene_start").is_between(target_gene["start"], target_gene["end"]) ) |
    ( polars.col("Gene_end").is_between(target_gene["start"], target_gene["end"]))
).collect()

print(target_gene)
print(target_exons)
print(target_peptides)

target_gene_length = target_gene["end"][0]-target_gene["start"][0]

target_exon_coverage = np.full(target_gene_length, False)
for exon in target_exons.iter_rows(named=True):
    exon_start = exon["start"]-target_gene["start"]
    exon_end = exon["end"]-target_gene["start"]
    target_exon_coverage[exon_start[0]:exon_end[0]]=True

target_peptide_coverage = np.full(target_gene_length, False)
for peptide in target_peptides.iter_rows(named=True):
    peptide_start = peptide["Gene_start"]-target_gene["start"]
    peptide_end = peptide["Gene_end"]-target_gene["start"]
    target_peptide_coverage[peptide_start[0]:peptide_end[0]]=True

exon_bases_covered_by_peptides = np.bitwise_and(target_exon_coverage,target_peptide_coverage)
non_exon_bases_covered_by_peptides = np.bitwise_and(target_peptide_coverage, ~target_exon_coverage)

print(np.sum(target_exon_coverage))
print(np.sum(exon_bases_covered_by_peptides))
result = polars.DataFrame(
    {
        "gene_id" : target_gene["gene_id"][0],
        "exon_percentage" : (np.sum(target_exon_coverage) / target_gene_length) * 100,
        "exon_coverage" : (np.sum(exon_bases_covered_by_peptides) / np.sum(target_exon_coverage)) * 100,
        "non_exon_coverage" : (np.sum(non_exon_bases_covered_by_peptides) / np.sum(~target_exon_coverage)) * 100,
    },
    schema=
    {
        "gene_id":polars.Utf8, 
        "exon_percentage":polars.Float64, 
        "exon_coverage":polars.Float64, 
        "non_exon_coverage":polars.Float64, 
    }
)
result_exons = polars.DataFrame(
    {
        "gene_id" : target_gene["gene_id"][0],
        "exon_ranges" : find_contiguous_true_ranges(target_exon_coverage),
    }
)
result_non_exons = polars.DataFrame(
    {
        "gene_id" : target_gene["gene_id"][0],
        "exon_ranges" : find_contiguous_true_ranges(non_exon_bases_covered_by_peptides),
    }
)
print(result_exons)


result.write_parquet(genestats_file)
result_exons.write_parquet(exonstats_file)
result_non_exons.write_parquet(nonexonstats_file)
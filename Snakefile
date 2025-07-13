import scripts.getGCForGlobalGenome

def getGenomeSize(wildcards):
    return scripts.getGCForGlobalGenome.getGCAndGenomeSize("data_" + wildcards.species + "/genome.fa")["only_called_nieb"]

def getGCgenome(wildcards):
    return scripts.getGCForGlobalGenome.getGCAndGenomeSize("data_" + wildcards.species + "/genome.fa")["gc"]

rule blockAndGaps:
    input:
        "data_{species}/chain_{rank}.gz"
    output:
        "data_{species}/gaps_{rank}.bed",
        "data_{species}/blocks_{rank}.bed"
    shell:
        "python3 scripts/getGapsPos.py {input} {output}"

rule get_microsat:
    input:
        "data_{species}/microsat.msdb"
    output:
        "data_{species}/microsat.bed"
    shell:
        "cut -f1,2,3,4,5,6 {input} > {output}"

rule sort_microsat:
    input:
        "data_{species}/microsat.bed"
    output:
        "data_{species}/microsat_annotated.bed"
    shell:
        "python3 scripts/sort_microsat.py {input} {output}"

rule createGenomeFile:
    input:
        "data_{species}/genome.fa"
    output:
        "data_{species}/genome.genome"
    shell:
        "python3 scripts/createGenomeFile.py {input} {output}"

rule cleanGaps:
    input:
        gaps = "data_{species}/gaps_{rank}.bed",
        blocks = "data_{species}/blocks_{rank}.bed"
    output:
        "data_{species}/cleaned_gaps_{rank}.bed"
    shadow: "shallow"
    shell:
        """
        bedtools subtract -a {input.gaps} -b {input.blocks} > tmp.bed
        bedtools sort -i tmp.bed > tmp.sort.bed
        bedtools merge -i tmp.sort.bed > tmp.merge.bed
        Rscript scripts/filter_per_size.R tmp.merge.bed {output}
        """

rule getSeqs:
    input:
        bed = "data_{species}/cleaned_{type}.bed", # either uniq or gaps_1
        genome = "data_{species}/genome.fa"
    output:
        "data_{species}/cleaned_{type}.fa"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output}"

rule checkNs:
    input:
        "data_{species}/cleaned_{type}.fa"
    output:
        "data_{species}/Ns_{type}.csv"
    shell:
        "python3 scripts/countNPerSeq.py {input} {output}"

rule filter_Ns:
    input:
        "data_{species}/Ns_{type}.csv"
    output:
        "data_{species}/not_too_many_Ns_{type}.csv"
    shell:
        "awk -F'\t' '$2 > 10 {{print $1}}' {input} > {output}"

rule filter_gaps:
    input:
        gaps = "data_{species}/cleaned_{type}.bed",
        Ns = "data_{species}/not_too_many_Ns_{type}.csv"
    output:
        "data_{species}/filtered_{type}.bed"
    shell:
        "grep -f <(sed 's/:/\t/g' {input.Ns} | sed 's/-/\t/g') -v {input.gaps} > {output}"

rule mergeBlocks:
    input:
        "data_{species}/blocks_1.bed",
        "data_{species}/blocks_2.bed"
    output:
        "data_{species}/blocks_total.bed"
    shell:
        "bedtools merge -i <(cat {input} | bedtools sort -i -) > {output}"

rule getUniqGaps:
    input:
        gaps1 = "data_{species}/filtered_gaps_1.bed",
        gaps2 = "data_{species}/filtered_gaps_2.bed"
    output:
        "data_{species}/cleaned_uniq.bed"
    shell:
        "bedtools intersect -a {input.gaps1} -b {input.gaps2} > {output}"

rule getCovNIEBxMicrosat:
    input:
        microsat = "data_{species}/microsat_annotated.bed",
        niebs = "data_{species}/niebs.bed"
    output:
        "data_{species}/nieb_microsat.tsv"
    shell:
        "scripts/intersectKeepingNames.bin +a {input.niebs} +b {input.microsat} +m hit +o {output}"


rule categorizeNIEBs:
    input:
        niebs = "data_{species}/niebs.bed",
        gaps = "data_{species}/cleaned_uniq.bed",
        blocks = "data_{species}/blocks_total.bed"
    output:
        "data_{species}/niebs_annotated.bed"
    shadow: "shallow"
    shell:
        """
        bedtools intersect -f 1 -a {input.niebs} -b {input.gaps} > tmp.gaps.bed
        bedtools intersect -f 1 -a {input.niebs} -b {input.blocks} > tmp.blocs.bed
        bedtools intersect -wa -a {input.niebs} -b {input.gaps} | bedtools subtract -a - -b tmp.gaps.bed > tmp.overlapping.bed
        Rscript scripts/annotateNIEBs.R tmp.gaps.bed tmp.blocs.bed tmp.overlapping.bed {output}
        """

rule getUniqueFromOverlap:
    input:
        niebs = "data_{species}/niebs_annotated.bed",
        gaps = "data_{species}/cleaned_uniq.bed"
    output:
        "data_{species}/uniq_part_of_overlap_niebs.bed"
    shell:
        "bedtools intersect -a <(grep overlap {input.niebs}) -b {input.gaps} > {output}"

rule uniqueAssociationToSat:
    input:
        niebs = "data_{species}/uniq_part_of_overlap_niebs.bed",
        microsat = "data_{species}/microsat_annotated.bed"
    output:
        "data_{species}/nieb_microsat_uniq.bed"
    shell:
        "bedtools intersect -b {input.niebs} -a {input.microsat} > {output}"

rule getCovNIEBoverlap:
    input:
        "data_{species}/uniq_part_of_overlap_niebs.bed",
        "data_{species}/nieb_microsat_uniq.bed",
        "data_{species}/microsat_gap_niebs.bed",
        "data_{species}/niebs_annotated.bed"
    output:
        "results_{species}/cov_niebs_overlap_microsat.csv"
    shell:
        "Rscript scripts/get_cov_niebs_overlap_microsat.R {input} {output}"

rule getCovGapXMicrosat:
    input:
        microsat = "data_{species}/microsat_annotated.bed",
        annotated = "data_{species}/niebs_annotated.bed"
    output:
        "data_{species}/microsat_gap_niebs.bed"
    shell:
        "scripts/intersectKeepingNames.bin +a {input.annotated} +b {input.microsat} +o {output} +m both"

rule show_cov:
    input:
       gaps = "data_{species}/cleaned_uniq.bed",
       blocks = "data_{species}/blocks_total.bed",
       genome = "data_{species}/genome.genome",
       niebs_annotated = "data_{species}/niebs_annotated.bed"
    output:
        touch("data_{species}/.cov_genome_done")
    shell:
        "Rscript scripts/check_cov.R {input}"

rule getCovMicrosat:
    input:
        niebs = "data_{species}/niebs.bed",
        microsat = "data_{species}/microsat_annotated.bed",
        genome = "data_{species}/genome.genome",
        niebs_annotated = "data_{species}/niebs_annotated.bed",
        niebXmicrosat = "data_{species}/nieb_microsat.tsv",
        microsatXgapXnieb = "data_{species}/microsat_gap_niebs.bed"
    output:
        "data_{species}/stat_coverage_microsat.csv"
    shell:
        "Rscript scripts/get_stat_coverage.R {input} {output}"

rule getStatCategories:
    input:
        "data_{species}/niebs_annotated.bed",
    output:
        "data_{species}/stat_niebs_categories.csv"
    shell:
        "Rscript scripts/stat_niebs_per_category.R {input} {output}"

rule plot_GC_microsat:
    input:
        microsat = "data_{species}/microsat.bed"
    output:
        "results_{species}/gc_histogram_microsat.png"
    shell:
        "Rscript scripts/plot_GC_microsat.R {input} {output}"

rule plotSizeBlocksGaps:
    input:
        blocks = "data_{species}/blocks_total.bed",
        gaps = "data_{species}/cleaned_uniq.bed"
    output:
        "results_{species}/stats_gaps_blocks.csv",
        "results_{species}/plot_blocks_gaps.png"
    shell:
        "Rscript scripts/plot_distribution_sizes.R {input} {output}"

rule plotCoverageBarplot:
    input:
        niebs = "data_{species}/niebs.bed",
        microsat = "data_{species}/microsat_annotated.bed",
        genome = "data_{species}/genome.genome",
        niebs_annotated = "data_{species}/niebs_annotated.bed",
        microsatXgapXnieb = "data_{species}/microsat_gap_niebs.bed",
        uniq_overlap = "data_{species}/uniq_part_of_overlap_niebs.bed",
        microsat_uniq = "data_{species}/nieb_microsat_uniq.bed"
    output:
        "data_{species}/savestate_cov.csv",
        "results_{species}/plot_covs.png"
    shell:
        "Rscript scripts/plot_coverage_microsat_V2.R {input} {output}"


rule getCoverage:
    input:
        aoe = "data_{species}/cleaned_gaps_{type}.aoe",
        bed = "data_{species}/niebs.bed"
    output:
        niebs = "data_{species}/niebs_counts_{type}.csv",
        ref = "data_{species}/ref_counts_{type}.csv"
    shell:
        """
        scripts/countFeatures.bin +a {input.aoe} +b {input.bed} +o {output.niebs} +d +k source +p whole
        scripts/countFeatures.bin +a {input.aoe} +b {input.aoe} +o {output.ref} +d +k source +p whole
        """

rule sumBases:
    input:
        niebs = "data_{species}/niebs_counts_{type}.csv",
        total = "data_{species}/ref_counts_{type}.csv"
    output:
        niebs = "data_{species}/niebs_total_{type}.csv",
        total = "data_{species}/ref_total_{type}.csv"
    shell:
        """
        python3 scripts/sum_sizes.py {input.niebs} {output.niebs}
        python3 scripts/sum_sizes.py {input.total} {output.total}
        """
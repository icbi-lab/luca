
process SCQC_MERGE_STATS {
    input:
        path(stats_tsv)

    output:
        path("qc_stats_all.tsv")

    script:
    """
    mkdir out
    head -n1 ${stats_tsv[0]} > out/qc_stats_all.tsv
    for f in *.tsv; do
        tail -n+2 \$f >> out/qc_stats_all.tsv
    done
    mv out/qc_stats_all.tsv .
    """
}

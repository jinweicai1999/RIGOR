# RIGOR: four-dimension benchmark results

This directory provides gene-level scores and paired comparisons for the 26-gene retrospective benchmark. The original Full report for each gene was used to select reference entries previously judged full or partial matches. Each selected entry was scored as 1, 0.5 or 0. Scores are computed within M1–M4 and combined using the selected-entry counts (38, 11, 28, 15 of 92) as weights. Weights are renormalized when a gene has no selected entries in a dimension. These are post-hoc concordance scores on a Full-selected subset, not an estimate of complete-report accuracy. Reference entries may overlap, and semantic assessments were made with AI assistance without independent blinded adjudication.

## Files

- `per_gene_condition_scores.csv`: the overall score for each gene and condition. No claim-by-claim scoring records are included.
- `condition_summary.csv`: number of genes, mean and sample standard deviation.
- `paired_tests_and_BH.csv`: observed paired differences, raw and BH-adjusted P values, test definitions and correction families.
- `standardized_26_gene_reports.zip`: the 156 reports from the standardized Full and five single-layer ablation conditions (26 per condition).
- `standardized_report_sha256.csv`: identifiers and SHA-256 hashes for the 156 archived reports.
- `verify_summaries.py`: recompute the descriptive statistics, standardized exact tests, and both BH adjustments using only Python's standard library.

The 26-gene output records from the original baseline campaign, including the original Full, annotation-only, literature-only and deterministic no-LLM outputs, are in the [earlier public benchmark archive](https://github.com/jinweicai1999/RIGOR/tree/9bc5ae70eaeece44ba91a4b4a039a16b4f72fcc7/supplementary_data/revision_20261018). That archive also retains historical five-dimension results, which must not be confused with the four-dimension results here.

## Analysis sets and comparisons

The **baseline** set contains four conditions: the original Full report, annotation/GO/InterPro only, pre-cutoff literature only, and deterministic no-LLM. Six paired comparisons were assessed by 199,999 within-gene label swaps per comparison (seed 20261027), reselecting entries and recomputing weights after each swap. BH correction covered these six comparisons. The tests are post-hoc, depend on exchangeability and concern the selected-entry statistic.

The **standardized ablation** set contains a different standardized Full report and five leave-one-layer-out reports per gene. Each was scored against the same fixed selected entries. Comparisons use two-sided exact random-sign signed-rank tests with zero differences omitted and averaged ranks for tied absolute differences. BH correction covers the five standardized Full-versus-ablation tests. The standardized Full and original Full are distinct outputs and must not be combined as one condition.

The reviewer-requested gene-level scores and complete report outputs are provided here and in the earlier archive. Individual reference-entry credits and the detailed adjudication records are not included; consequently these files permit recalculation of descriptive summaries and gene-level paired statistics, but not independent reconstruction or adjudication of every entry-level score.

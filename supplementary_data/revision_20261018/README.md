# RIGOR revision evidence package — 2026-10-18

[Download the complete archive](RIGOR_revision_data_20261018.zip) · [Archive SHA-256](RIGOR_revision_data_20261018.zip.sha256)

This package publishes existing archived results and the associated audit materials.
It does not represent new model calls, new wet-lab experiments, or completion of all
outstanding manuscript revisions. The existing application code is unchanged.

## Contents

| Directory inside the archive | Material |
|---|---|
| `01_baselines` | Four conditions × 26 genes; paired comparisons; reference texts; worked example; original scoring script; portable scoring replay; archived normalized outputs and 338 raw benchmark API responses |
| `02_ablation` | Existing standardized-output ablation summaries and paired/sensitivity comparisons; these use their separate full-evidence reference mean of 79.85, rather than 80.82 |
| `03_same_evidence_comparison` | C2DP40 structured-versus-synthesized display and its source/interpretation limits |
| `04_stability` | Original 200-run archive, 800-call manifest, run-level citations, and 19,900 cited-PMID-set comparisons |
| `05_free_recall` | 26 evidence-free recall responses, adjudications, and retained API records |
| `06_route_audit` | Existing network-free P1 source-path audit with fixtures and stubbed PubMed results |
| `07_figures` | Response Figures R1/R2 and their plotting script; Figure R2 is the single Jaccard-distribution plot |
| `08_api_manifest` | Per-call manifest for the original 338-call benchmark campaign |
| `09_provenance` | Local resource sources, available version labels, original-date status, per-file hashes, and trait-specific GWAS backgrounds |

## Readable provenance lookups

- [Evidence-layer lookup](provenance/Evidence_Layer_Provenance_1018.csv)
- [244-file provenance and SHA-256 inventory](provenance/Evidence_File_Provenance_SHA256_1018.csv)
- [Trait-specific GWAS background lookup](provenance/GWAS_Trait_Background_1018.csv)
- [Epigenomic accession lookup](provenance/Epigenomic_Accession_Lookup_1018.csv)
- [Official/local RGAP expression-file identity check](provenance/Expression_RGAP_Official_File_Match_1018.json)
- [Provenance audit](provenance/Provenance_Audit_1018.json)

The 244 retained files include supporting and derived artifacts; they are not 244
independent evidence inputs. The author reports source download dates of 2026-04-10 for annotation, expression,
co-expression and promoter-epigenomic resources, and 2026-04-16 for GWAS association
files. No independent download log was recovered; this reporting basis is explicit in
the lookup. Creation/copy timestamps, index-build dates and verification dates remain separate.
The three core RAP-DB GFFs exactly match the official 2025-03-19 release. The ten GWAS traits map to the RiceVarMap 529-accession panel of Xie et al. (2015),
supported by 40/40 trait-SNP numerical reconstructions matching NMISS, BETA, STAT
and P at the stored precision. The source check uses four biallelic SNPs, not a full
GWAS rerun, and does not recover the historical PLINK command/version.
RiceENCODE is cited as the 2021 resource; its retained 2026-04-10 local snapshot is
identified by the file/accession/hash records, without inventing a provider release number.
Some exact upstream support-file/processed-BED revisions remain unverified. Embedded workbook dates and local index-build times are not download
dates. A current file hash alone does not establish its identity with every historical
benchmark input.

The retained expression workbook is byte-identical to the official RGAP/MSU Release 7
download, as documented in the linked identity-check record. That match does not
reconstruct the original project download date.

## Checks and replay

The complete archive is 11,476,057 bytes. SHA-256:

```text
72443aaa82272ff9a120439278f29dcf340ebd8b90233eab7844dcf57e353890
```

After extracting the archive, run:

```bash
python3 RIGOR_revision_data_20261018/01_baselines/reproduce_alignment.py
```

This standard-library-only replay compares all 104 baseline scores and their component
scores with the included table. The retained check reports zero mismatches. It does not
contact a model provider. `MANIFEST_SHA256.csv` inside the archive covers all other 47
packaged files. See the [packaging-check record](packaging_audit_20261018.json).

The packaging audit records the state before publication; its upload-status field is
not a live repository-status indicator. Archive CRC checks and a credential-pattern
scan passed, but such checks are not a guarantee of scientific validity or absolute
security. Scoring results measure retrospective text/evidence alignment, not genome-wide
prediction accuracy. The deferred 26-gene prior-literature audit and wet-lab documentation
are not presented here as completed.

## Source-date correction

[Official release/download ordering checks](provenance/Resource_Release_and_Download_Check_1018.csv) · [Author-reported dates and scope](provenance/Download_Date_Author_Confirmation_1018.json) · [RAP-DB archive identity checks](provenance/Official_Resource_Identity_Checks_1018.json) · [GWAS portal/local spot checks](provenance/GWAS_Portal_Local_Spot_Check_1018.csv)

This correction changes provenance metadata only; analysis outputs, scoring scripts, and raw API records are unchanged. Resource-publication dates are not silently substituted for exact processed-file build dates.

## GWAS source-panel resolution — 2026-10-19

[40 numerical source checks](provenance/GWAS_Source_Reconstruction_40_checks_1019.csv) · [Methods and scope](provenance/GWAS_Source_Reconstruction_1019.json) · [Portable verification script](provenance/verify_gwas_source.py) · [Public phenotype input](provenance/RiceVarMap_phenos.csv) · [Public genotype input](provenance/RiceVarMap_genotypes_1.csv)

The source panel has 529 accessions after QC of the original 533-accession panel
([Xie et al., 2015](https://doi.org/10.1073/pnas.1515919112)); the phenotype table has
521–529 nonmissing values depending on the trait. Per-SNP NMISS can be lower.
The prior portal/local P-value spot check was not a like-for-like model comparison
and is not used to reject the population assignment. One of five queried markers
had deletion/multi-allelic portal display calls and was not recoded. The four
biallelic markers give 40/40 matching trait-SNP checks. No historical scientific
outputs, scoring results, or provider API records were changed.

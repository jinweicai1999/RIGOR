# RIGOR — Rice Integrative Gene functiOn Reasoner

RIGOR is an interpretable, evidence-prioritized AI framework that converts fragmented public
rice evidence into traceable functional and breeding-oriented inference. It supports single-gene
analysis, batch prioritization of candidate gene sets, and prediction for regulatory ncRNAs.

- Web resource: <http://research.gzsys.org.cn/rigor/>
- Corresponding authors: Wenlong Zhao (zhaowlong7@mail.sysu.edu.cn), Yuchan Zhang (zhyuchan@mail.sysu.edu.cn)

## Repository layout

```
.
├── src/                     RIGOR source code
│   ├── prompts/             runtime prompt templates (system_prompts.py)
│   └── utils/               scoring, ranking and post-processing utilities
├── apps/                    web dashboard and command-line entry points
├── data/                    packaged local evidence resources (see note below)
└── supplementary_data/      benchmark registry and reproducibility data
```

> **Note on `data/`.** The packaged annotation, co-expression, promoter-epigenomic and GWAS
> resources used at inference time are large reference files assembled from third-party sources
> (RAP-DB/RGAP, RGAP expression, eRice/RiceENCODE chromatin data, published GWAS summary
> statistics). They are not redistributed here; see the manuscript and its Supplemental
> Materials and Methods for the exact source, version and access route of each layer.

## Supplementary data

| File | Contents |
|---|---|
| `supplementary_data/Supplementary_Data_S1.csv` | 200-run stability audit (run status, generation time, output length, provider-returned model, response IDs, cited and unsupported PMIDs, output SHA-256) |
| `supplementary_data/Supplementary_Data_S2.zip` | Drought and plasma-membrane batch-validation statistics: `primary_tests.csv`, `scenario_score_profile.csv`, `plasma_membrane_gene_set.csv`, `mapping_audit.json`, `data_quality_audit.csv` |
| `supplementary_data/Supplementary_Data_S3.zip` | 26-gene evidence-layer ablation: `condition_summary.csv`, `paired_ablation_comparisons.csv`, `paired_ablation_comparisons_no_M5.csv`, `integrity_audit.csv` |
| `supplementary_data/26_Gene_Benchmark_Registry.xlsx` | 26-gene time-delayed benchmark registry: overview, selected cohort, ground-truth references, screening log, candidate pool, source versions, leakage audit |

The frozen evidence snapshot, the populated prompts, all 200 complete model outputs, the
800-response provider manifest and per-file SHA-256 hashes are distributed with the revision's
Supplementary Data and are also attached to the [Releases](../../releases) page of this repository
(`Data_archive_200runs.zip`).

## Model configuration

Production runs used **deepseek-flash** (temperature 0.1, thinking mode disabled,
`max_tokens = 4000`) through an OpenAI-compatible endpoint. The provider updated this model's API
alias during the revision period, so the August 2026 stability audit files record the identifier
returned at that time (`deepseek-v4-flash`) while later benchmark runs returned `deepseek-flash`;
both refer to the same production model.

## Citation

Jin, W., Feng, J., Zhao, R., Qin, Y., Chen, Z., Cao, Z., Zheng, B., Dong, Y., Zhao, W., and
Zhang, Y. RIGOR: an evidence-prioritized AI framework for systematic inference of rice gene
function and breeding relevance. *(citation to be completed on publication)*

## License

MIT — see [LICENSE](LICENSE).

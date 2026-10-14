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
│   └── utils/               scoring, ranking and post-processing utilities
├── apps/                    web dashboard and command-line entry points
├── prompts/                 line-anchored export of runtime prompt expressions
├── scripts/                 data-index and figure-generation utilities
└── supplementary_data/      benchmark registry and reproducibility data
```

The source and runtime prompts are included in this repository. A historical dated dashboard
backup is omitted; `apps/web_dashboard.py` is the current
dashboard entry point and reads its credential from the server environment.
Prompt expressions are embedded in the current `apps/` and `src/agents/` code. Their
auditable text and source locations are exported in `prompts/runtime_prompt_templates.md`;
`scripts/export_runtime_prompts.py` regenerates that file. The empty placeholder
`src/prompts/system_prompts.py` is not a runtime prompt source.

> **Note on local `data/`.** The packaged annotation, co-expression, promoter-epigenomic and GWAS
> resources used at inference time are large reference files assembled from third-party sources
> (RAP-DB/RGAP, RGAP expression, eRice/RiceENCODE chromatin data, published GWAS summary
> statistics). They are not redistributed here. The supplementary methods describe their
> scientific sources; a file-level provenance and checksum table is being prepared separately.

## Local setup

Install the packages in `requirements.txt`, mount the locally licensed evidence files under
`data/` (or configure the corresponding paths), and set `DEEPSEEK_API_KEY` in the server
environment before starting `streamlit run apps/web_dashboard.py`. No API credential is stored
in the repository. Without the local evidence bundle, the source alone does not reproduce
the complete inference outputs.

## Supplementary data

| File | Contents |
|---|---|
| `supplementary_data/Supplementary_Data_S1.csv` | 200-run stability audit (run status, generation time, output length, provider-returned model, response IDs, cited and unsupported PMIDs, output SHA-256) |
| `supplementary_data/Supplementary_Data_S2.zip` | Drought and plasma-membrane batch-validation statistics: `primary_tests.csv`, `scenario_score_profile.csv`, `plasma_membrane_gene_set.csv`, `mapping_audit.json`, `data_quality_audit.csv` |
| `supplementary_data/Supplementary_Data_S3.zip` | 26-gene evidence-layer ablation: `condition_summary.csv`, `paired_ablation_comparisons.csv`, `paired_ablation_comparisons_no_M5.csv`, `integrity_audit.csv` |
| `supplementary_data/26_Gene_Benchmark_Registry.xlsx` | 26-gene time-delayed benchmark registry, three sheets: **Selected 26** (the locked cohort, with RAP ID, class, earliest public date, PMID/DOI, experimental validation, gold-standard status and selection basis), **Screening Log** (all 225 screened PubMed records with dates, detected genes, experimental terms, screen score, inclusion/exclusion decision and reason) and **Leakage Audit** (target-specific pre-cutoff exposure check and confidence for each benchmark gene) |

The frozen evidence snapshot, populated prompts, all 200 complete model outputs, the
800-response provider manifest, and per-file SHA-256 hashes are retained in the separate
`Data_archive_200runs.zip` revision archive. This archive is not part of the source tree above.

## Model configuration

The reported runs requested a DeepSeek flash-class model through an OpenAI-compatible endpoint
with temperature 0.1, thinking mode disabled, and `max_tokens = 4000`. The August 2026
stability responses returned `deepseek-v4-flash`; later benchmark responses returned
`deepseek-flash`. These are the raw API-returned identifiers. Backend equivalence across
provider aliases, the model knowledge cutoff, and context-window size cannot be verified
from the retained responses.

## Citation

Jin, W., Feng, J., Zhao, R., Qin, Y., Chen, Z., Cao, Z., Zheng, B., Dong, Y., Zhao, W., and
Zhang, Y. RIGOR: an evidence-prioritized AI framework for systematic inference of rice gene
function and breeding relevance. *(citation to be completed on publication)*

## License

MIT — see [LICENSE](LICENSE).

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
ANNOTATION_DIR = DATA_DIR / "annotations"
EXPRESSION_DIR = DATA_DIR / "expression"
PROCESSED_DIR = DATA_DIR / "processed"
NCRNA_DIR = DATA_DIR / "ncrna"
EPIGENOMICS_DIR = DATA_DIR / "epigenomics"
CHROMATIN_DIR = EPIGENOMICS_DIR / "chromatin_accessibility"
HISTONE_DIR = EPIGENOMICS_DIR / "histone_modification"
GWAS_INFO_DIR = DATA_DIR / "gwas"
NETWORK_DIR = DATA_DIR / "network"

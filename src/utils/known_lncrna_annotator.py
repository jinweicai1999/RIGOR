"""
Has lncRNA Annotationer

Parse lncRNA.fa, extract Gene group positions and LOC ID
For position-based FunctionInference
"""

from pathlib import Path 
from src .project_paths import NCRNA_DIR 


class KnownLncRnaAnnotator :
    def __init__ (self ,fasta_file =None ):
        self .fasta_file =Path (fasta_file )if fasta_file else NCRNA_DIR /"lncRNA.fa"
        self .lncrna_db ={}

    def build_database (self ):
        """Parse lncRNA.fa, build Database"""
        if not self .fasta_file .exists ():
            return 

        with open (self .fasta_file )as f :
            for line in f :
                if line .startswith ('>'):
                #Format: >XR_001539667.2#OSA-LNC9267644#LOC9267644
                    parts =line .strip ()[1 :].split ('#')
                    if len (parts )>=3 :
                        refseq_id =parts [0 ]
                        lnc_id =parts [1 ]
                        loc_id =parts [2 ]

                        self .lncrna_db [lnc_id .lower ()]={
                        'RefSeq_ID':refseq_id ,
                        'LNC_ID':lnc_id ,
                        'LOC_ID':loc_id 
                        }

    def query (self ,lnc_id ):
        """Query lncRNA Information"""
        return self .lncrna_db .get (lnc_id .lower ())

    def has_lncrna (self ,lnc_id ):
        """Check if exists"""
        return lnc_id .lower ()in self .lncrna_db 

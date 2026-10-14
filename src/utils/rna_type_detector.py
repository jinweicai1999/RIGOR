"""
Intelligent RNA Type Recognition System

Based on the following strategies:
1. Sequence length judgment
2. ID pattern matching
3. Sequence database alignment
"""

import os 
from pathlib import Path 
from src .project_paths import NCRNA_DIR 


class RnaTypeDetector :
    def __init__ (self ):
        self .base_dir =NCRNA_DIR 

        #Sequence length threshold
        self .length_rules ={
        'miRNA':(18 ,25 ),#Mature miRNA: 18–25 nt
        'miRNA_precursor':(60 ,300 ),#miRNA precursor: 60–300 nt
        'phasiRNA':(21 ,21 ),#phasiRNA: fixed 21 nt
        'lncRNA':(200 ,100000 )# lncRNA: >200nt
        }

        #Sequence database cache
        self ._seq_db ={}

    def _load_fasta (self ,fasta_file ):
        """Load FASTA file into dictionary"""
        if fasta_file in self ._seq_db :
            return self ._seq_db [fasta_file ]

        seq_dict ={}
        filepath =self .base_dir /fasta_file 

        if not filepath .exists ():
            return seq_dict 

        with open (filepath )as f :
            current_id =None 
            current_seq =[]

            for line in f :
                line =line .strip ()
                if line .startswith ('>'):
                    if current_id :
                        seq_dict [current_id .lower ()]=''.join (current_seq ).upper ()
                    current_id =line [1 :].split ()[0 ]#Take the first field as the ID
                    current_seq =[]
                else :
                    current_seq .append (line )

            if current_id :
                seq_dict [current_id .lower ()]=''.join (current_seq ).upper ()

        self ._seq_db [fasta_file ]=seq_dict 
        return seq_dict 

    def detect_by_id (self ,rna_id ):
        """Identify RNA type based on ID pattern"""
        rna_id_lower =rna_id .lower ()

        #miRNA pattern
        if any (pattern in rna_id_lower for pattern in ['osa-mir','osa-let','osa-mirn']):
            return 'miRNA'

            #phasiRNA pattern
        if 'phasi'in rna_id_lower :
            return 'phasiRNA'

            #lncRNA pattern
        if 'lncn'in rna_id_lower or rna_id_lower .startswith ('tcons'):
            return 'lncRNA'

        return None 

    def detect_by_length (self ,sequence ):
        """Determine RNA type based on sequence length"""
        seq_len =len (sequence )

        candidates =[]
        for rna_type ,(min_len ,max_len )in self .length_rules .items ():
            if min_len <=seq_len <=max_len :
                candidates .append (rna_type )

                #Priority: Exact match > Range match
        if 'phasiRNA'in candidates :
            return 'phasiRNA'
        elif 'miRNA'in candidates :
            return 'miRNA'
        elif 'miRNA_precursor'in candidates :
            return 'miRNA_precursor'
        elif 'lncRNA'in candidates :
            return 'lncRNA'

        return None 

    def detect_by_database (self ,sequence ):
        """Search for matches in SequenceDatabase"""
        sequence =sequence .upper ().strip ()

        #Check each database (Order determines Priority; mature miRNA library takes precedence over precursor library)
        databases =[
        ('novel_mature_miRNA.fa','miRNA'),
        ('mirbase_osa_mature.fa','miRNA'),
        ('phasiRNA.fa','phasiRNA'),
        ('novel_lncRNA.fa','lncRNA'),
        ('novel_precursor.fa','miRNA_precursor')
        ]

        for fasta_file ,rna_type in databases :
            seq_dict =self ._load_fasta (fasta_file )
            if sequence in seq_dict .values ():
                return rna_type 

        return None 

    def detect (self ,rna_id =None ,sequence =None ):
        """
        Comprehensive RNA type identification

        Priority:
        1. Exact match in SequenceDatabase (most reliable)
        2. ID pattern matching
        3. Sequence length determination
        """
        results =[]

        #Strategy 1: SequenceDatabase matching
        if sequence :
            db_type =self .detect_by_database (sequence )
            if db_type :
                return db_type ,'database_match'

                #Strategy 3: Sequence length
            len_type =self .detect_by_length (sequence )
            if len_type :
                results .append ((len_type ,'length'))

                #Strategy 2: ID pattern
        if rna_id :
            id_type =self .detect_by_id (rna_id )
            if id_type :
                results .append ((id_type ,'id_pattern'))

                #Return the first matching result
        if results :
            return results [0 ]

        return None ,'unknown'

    def get_sequence_from_id (self ,rna_id ):
        """Retrieve sequence from database using ID"""
        rna_id_lower =rna_id .lower ()

        #Mature miRNA library takes priority (exact match + prefix fuzzy match, but not across libraries)
        #Priority: miRBase known > novel mature > phasiRNA > lncRNA > precursor (last resort)
        mature_databases =[
        'mirbase_osa_mature.fa',
        'novel_mature_miRNA.fa',
        'phasiRNA.fa',
        'novel_lncRNA.fa',
        ]
        for fasta_file in mature_databases :
            seq_dict =self ._load_fasta (fasta_file )
            #Exact match
            if rna_id_lower in seq_dict :
                return seq_dict [rna_id_lower ]
                #Prefix fuzzy match (handles two cases):
                #1. Library ID with suffix: Osa-miRN01_RiceNCexp → Matches input Osa-miRN01
                #2. Library ID with arm identifier: osa-mir168a-5p → Matches input osa-mir168a (returns 5p strand)
            for seq_id ,sequence in seq_dict .items ():
                if seq_id .startswith (f"{rna_id_lower }_")or seq_id .startswith (f"{rna_id_lower }-"):
                    return sequence 

                    #Precursor library is last resort—only checked if no mature library matches above
        precursor_dict =self ._load_fasta ('novel_precursor.fa')
        if rna_id_lower in precursor_dict :
            return precursor_dict [rna_id_lower ]
        for seq_id ,sequence in precursor_dict .items ():
            if seq_id .startswith (f"{rna_id_lower }_")or rna_id_lower in seq_id :
                return sequence 

        return None 

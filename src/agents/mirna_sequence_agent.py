import pandas as pd 
from collections import defaultdict 

class MiRnaSequenceAgent :
    """sRNA candidate target Gene prediction engine based on sequence complementarity (plant standard)"""

    #Plant miRNA scoring threshold
    PENALTY_THRESHOLD =4.0 #Only keep if total penalty ≤ 4.0
    CORE_PENALTY_MAX =2.0 #Core region penalty cap (positions 2–12)
    CLEAVAGE_POS =(9 ,10 )#1-based cleavage site (positions 10–11, 0-based 9–10)

    def __init__ (self ,annotator ):
        self .annotator =annotator 
        self .complement ={'A':'U','U':'A','G':'C','C':'G','T':'A'}

    def _reverse_complement (self ,seq ):
        """Reverse complement"""
        seq =seq .upper ().replace ('T','U')
        return ''.join (self .complement .get (b ,'N')for b in reversed (seq ))

    def _score_pairing (self ,mirna_seq ,target_seq ):
        """
        Plant miRNA scoring rules (referencing psRNATarget / plant AGO1 cleavage mechanism):

        Penalty system (lower = better):
          Perfect pairing (A-U, G-C)    : +0.0
          G-U wobble                 : +0.5
          Mismatch                   : +1.0
          Gap/Bulge                  : +2.0

        Elimination rules (veto):
          1. Cleavage site (1-based positions 10, 11; i.e., 0-based 9, 10) has mismatch or gap → eliminate
          2. ≥3 consecutive mismatches → eliminate
          3. Core region (1-based positions 2–12; i.e., 0-based 1–11) penalty > 2.0 → eliminate
          4. Total penalty > 4.0 → eliminate

        Returns (penalty, is_valid): lower penalty is better; is_valid=False means eliminated
        """
        mirna_seq =mirna_seq .upper ().replace ('T','U')
        target_seq =target_seq .upper ().replace ('T','U')

        #miRNA 5' end corresponds to mRNA target: compare target reverse complement with miRNA position-by-position
        #Actually, miRNA pairs with reverse complement of mRNA, equivalent to mirna_rc compared with target order
        mirna_rc =self ._reverse_complement (mirna_seq )

        min_len =min (len (mirna_rc ),len (target_seq ))

        total_penalty =0.0 
        core_penalty =0.0 
        consec_mismatch =0 
        max_consec =0 

        for i in range (min_len ):
            m =mirna_rc [i ]
            t =target_seq [i ]
            pos_1based =i +1 #Convert to 1-based indexing

            #Cleavage site check (positions 10, 11)
            at_cleavage =pos_1based in (10 ,11 )

            if m =='N'or t =='N':
                pen =2.0 #Treat as gap
                consec_mismatch +=1 
            elif m ==t :
                pen =0.0 
                consec_mismatch =0 
            elif (m =='G'and t =='U')or (m =='U'and t =='G'):
                pen =0.5 
                consec_mismatch =0 
                #Wobble at cleavage site allowed but penalty continues
                if at_cleavage :
                    pass #Wobble at cleavage site is acceptable
            else :
                pen =1.0 
                consec_mismatch +=1 
                #Mismatch at cleavage site → automatic rejection
                if at_cleavage :
                    return 999.0 ,False 

            max_consec =max (max_consec ,consec_mismatch )

            total_penalty +=pen 
            if 2 <=pos_1based <=12 :#Core region
                core_penalty +=pen 

                #≥3 consecutive mismatches → Eliminated
        if max_consec >=3 :
            return 999.0 ,False 

            #Core region penalty exceeds limit →Eliminated
        if core_penalty >self .CORE_PENALTY_MAX :
            return total_penalty ,False 

            #Total penalty exceeds limit → Eliminated
        if total_penalty >self .PENALTY_THRESHOLD :
            return total_penalty ,False 

        return total_penalty ,True 

    def predict_targets (self ,mirna_seq ,top_n =10 ,rna_label ="miRNA"):
        """
        Predict target genes (plant standard: Scan CDS + 3'UTR full length)

        Parameters:
            mirna_seq : small RNA sequence (mature form, ≤25 nt)
            top_n     : return top N candidates
            rna_label : label (miRNA / phasiRNA)

        Returns:
            DataFrame with columns: Gene_ID, Symbol, Penalty, Scoring_Mode, PMIDs
        """
        mirna_seq =mirna_seq .upper ().replace ('T','U')
        if not all (b in 'AUCG'for b in mirna_seq ):
            raise ValueError ("Sequence contains invalid characters; only A/U/C/G/T are supported")
        if rna_label =="phasiRNA":
            if not (21 <=len (mirna_seq )<=24 ):
                raise ValueError (f"phasiRNA abnormal length：{len (mirna_seq )}nt（should be 21-24nt）")
        elif not (18 <=len (mirna_seq )<=25 ):
            raise ValueError (f"miRNA abnormal length：{len (mirna_seq )}nt（should be 18-25nt）")

        print (f"  🔬 SequenceValidationthrough：{len (mirna_seq )}nt")
        print (f"  🧬 {rna_label } Sequence：{mirna_seq }")
        print (f"  📡 CurrentlyScanGenegroup（CDS + 3'UTR）...")

        results =[]
        gene_count =0 
        window_size =len (mirna_seq )

        for gene_id in self .annotator .db .keys ():
            gene_info =self .annotator .query (gene_id )
            if not gene_info :
                continue 

            gene_count +=1 

            #Scan CDS + 3'UTR full length (mock implementation, to be replaced with real genome sequence)
            full_seq =self ._generate_mock_mrna (gene_id )

            best_penalty =999.0 
            for i in range (len (full_seq )-window_size +1 ):
                window =full_seq [i :i +window_size ]
                penalty ,valid =self ._score_pairing (mirna_seq ,window )
                if valid and penalty <best_penalty :
                    best_penalty =penalty 

            if best_penalty <=self .PENALTY_THRESHOLD :
                symbol =gene_info .get ('Primary_Symbol')or (
                gene_info .get ('Synonyms',[''])[0 ]if gene_info .get ('Synonyms')else 'No'
                )
                results .append ({
                'Gene_ID':gene_id ,
                'Symbol':symbol ,
                'Penalty':best_penalty ,
                'PMIDs':gene_info .get ('PMIDs',[]),
                'Scoring_Mode':'plant_psrnatarget',
                })

            if gene_count %5000 ==0 :
                print (f"    HasScan {gene_count } eachGene...")

        print (f"  ✅ ScanCompleted！totalScan {gene_count } eachGene")

        df =pd .DataFrame (results )
        if df .empty :
            return df 

            #Lower penalty is better; sorted in ascending order
        df =df .sort_values (by ='Penalty',ascending =True )
        return df .head (top_n )

    def _generate_mock_mrna (self ,gene_id ):
        """Generate simulated full-length mRNA sequence (CDS + 3'UTR, for demonstration)"""
        import hashlib 
        hash_val =int (hashlib .md5 (gene_id .encode ()).hexdigest (),16 )
        bases =['A','U','C','G']
        seq =[]
        for i in range (500 ):#500nt simulated full-length mRNA
            seq .append (bases [(hash_val >>(i %32 ))%4 ])
        return ''.join (seq )

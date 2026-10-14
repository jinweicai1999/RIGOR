import pandas as pd 

class LncRnaAnnotator :
    def __init__ (self ,gtf_file ,fa_file =None ):
        self .gtf_file =gtf_file 
        self .fa_file =fa_file 
        self .db ={}

    def build_database (self ):
        """Parse GTF file to build lncRNA position database"""
        print ("⏳ Currently loading lncRNA gene group annotation...")

        with open (self .gtf_file ,'r')as f :
            for line in f :
                if line .startswith ('#'):
                    continue 

                parts =line .strip ().split ('\t')
                if len (parts )<9 :
                    continue 

                chrom =parts [0 ]
                feature =parts [2 ]
                start =int (parts [3 ])
                end =int (parts [4 ])
                strand =parts [6 ]
                attrs =parts [8 ]

                #Only process transcript lines
                if feature !='transcript':
                    continue 

                    #Parse gene_id
                gene_id =None 
                for attr in attrs .split (';'):
                    attr =attr .strip ()
                    if attr .startswith ('gene_id'):
                        gene_id =attr .split ('"')[1 ]
                        break 

                if not gene_id :
                    continue 

                    #Store information
                if gene_id not in self .db :
                    self .db [gene_id ]={
                    'ID':gene_id ,
                    'Chromosome':chrom ,
                    'Start':start ,
                    'End':end ,
                    'Strand':strand ,
                    'Length':end -start +1 
                    }

        print (f"✅ Successfulload {len (self .db )} eachlncRNAGenegroup positionInformation")

    def query (self ,lnc_id ):
        """Query lncRNA information"""
        #Remove possible transcript suffix
        base_id =lnc_id .split ('.')[0 ]
        return self .db .get (base_id )

    def query_by_coordinate (self ,chrom ,start ,end ,min_overlap =0.8 ):
        """
        Find similar lncRNAs by gene group coordinates.

        Similarity defined as: overlap_len / min(query_len, lnc_len),
        i.e., the proportion of the overlap between user interval and lncRNA interval relative to the shorter one.

        Return a list sorted by similarity in descending order, each item contains:
          {'id': str, 'similarity': float, 'chromosome': str,
           'start': int, 'end': int, 'strand': str}
        Only return entries with similarity >= min_overlap.
        """
        #Normalize chromosome names: unify to lowercase chrN format
        def _norm_chrom (c ):
            c =c .strip ().lower ()
            if not c .startswith ('chr'):
                c ='chr'+c 
            return c 

        query_chrom =_norm_chrom (chrom )
        query_len =end -start +1 
        hits =[]

        for lnc_id ,info in self .db .items ():
            lnc_chrom =_norm_chrom (info ['Chromosome'])
            if lnc_chrom !=query_chrom :
                continue 

            lnc_start =info ['Start']
            lnc_end =info ['End']
            lnc_len =lnc_end -lnc_start +1 

            #Calculate overlap
            overlap_start =max (start ,lnc_start )
            overlap_end =min (end ,lnc_end )
            overlap_len =max (0 ,overlap_end -overlap_start +1 )

            if overlap_len ==0 :
                continue 

            similarity =overlap_len /min (query_len ,lnc_len )
            if similarity >=min_overlap :
                hits .append ({
                'id':lnc_id ,
                'similarity':round (similarity ,4 ),
                'chromosome':info ['Chromosome'],
                'start':lnc_start ,
                'end':lnc_end ,
                'strand':info ['Strand'],
                })

        hits .sort (key =lambda x :x ['similarity'],reverse =True )
        return hits 

    def find_nearby_genes (self ,lnc_id ,protein_annotator ,distance =10000 ):
        """Find protein-coding genes near lncRNA"""
        lnc_info =self .query (lnc_id )
        if not lnc_info :
            return []

        chrom =lnc_info ['Chromosome']
        lnc_start =lnc_info ['Start']
        lnc_end =lnc_info ['End']

        nearby_genes =[]
        checked_count =0 
        same_chrom_count =0 
        sample_printed =0 

        #Iterate over all protein genes
        for gene_id ,gene_info in protein_annotator .db .items ():
            gene_loc =gene_info .get ('Location','')
            if not gene_loc :
                continue 

            checked_count +=1 

            #Print first 5 samples for debugging
            if sample_printed <5 :
                print (f"[example{sample_printed +1 }] {gene_id }: Location='{gene_loc }'")
                sample_printed +=1 

                #Parse location: chr01:1000-2000 (+) or chr01:1000-2000
            try :
            #Remove possible strand direction information
                gene_loc_clean =gene_loc .split ('(')[0 ].strip ()
                gene_chrom ,pos =gene_loc_clean .split (':')
                gene_start ,gene_end =map (int ,pos .split ('-'))
            except :
                continue 

                #Check if on the same chromosome (remove chr prefix and leading zeros)
            lnc_chrom_norm =chrom .lower ().replace ('chr','').lstrip ('0')or '0'
            gene_chrom_norm =gene_chrom .lower ().replace ('chr','').lstrip ('0')or '0'

            if sample_printed <=5 and checked_count <=10 :
                print (f"  lnc_chrom_norm='{lnc_chrom_norm }', gene_chrom_norm='{gene_chrom_norm }'")

            if lnc_chrom_norm !=gene_chrom_norm :
                continue 

            same_chrom_count +=1 

            #Calculate distance
            if gene_end <lnc_start :
                dist =lnc_start -gene_end 
            elif gene_start >lnc_end :
                dist =gene_start -lnc_end 
            else :
                dist =0 #Overlap

            if dist <=distance :
                nearby_genes .append ({
                'Gene_ID':gene_id ,
                'Distance':dist ,
                'Info':gene_info 
                })

        print (f"[debug] lncRNAchromosome: '{chrom }' (after standardization: '{chrom .lower ().replace ('chr','')}')")
        print (f"[debug] Checkdone {checked_count } eachHasLocation'sGene")
        print (f"[debug] among which {same_chrom_count } one on same chromosome")
        print (f"[debug] found {len (nearby_genes )} in the±{distance }bpin the range")

        #Sort by distance
        nearby_genes .sort (key =lambda x :x ['Distance'])
        return nearby_genes 

import urllib .parse 
import re 
import os 

class RiceGeneAnnotator :
    def __init__ (self ,gff_path ,mapping_path ):
        self .gff_path =gff_path 
        self .mapping_path =mapping_path 
        self .db ={}
        self .alias_index ={}

    def build_database (self ):
        lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"
        print ("Loading transcripts.gff annotation database (with literature and GO terms)..."if lang =="en"else "Currently loading transcripts.gff annotation library (including literature and GO annotations)...")

        #1. Parsing RAP-MSU mapping table
        rap_to_msu ={}
        with open (self .mapping_path ,'r')as f :
            for line in f :
                parts =line .strip ().split ('\t')
                if len (parts )>=2 and parts [1 ].strip ()!='None':
                    rap_id =parts [0 ].strip ()
                    msu_ids =[]
                    for model_id in [m .strip ()for m in parts [1 ].strip ().split (',')if m .strip ()]:
                        if model_id not in msu_ids :
                            msu_ids .append (model_id )
                    rap_to_msu [rap_id ]=msu_ids 
                    for msu in msu_ids :
                        self .alias_index [msu .lower ()]=rap_id 

                        #2. Parse transcripts.gff (focusing specifically on mRNA level)
        with open (self .gff_path ,'r')as f :
            for line in f :
                if line .startswith ('#'):continue 
                parts =line .strip ().split ('\t')

                #We only look at mRNA because the richest annotations are here
                if len (parts )<9 or parts [2 ]!='mRNA':continue 

                attrs ={}
                for item in parts [8 ].split (';'):
                    if '='in item :
                        k ,v =item .split ('=',1 )
                        attrs [k ]=urllib .parse .unquote (v )

                rap_id =attrs .get ('Locus_id')
                if not rap_id :continue 

                #💡 Core: Extract Literature PMID
                raw_pmid =attrs .get ('Literature_PMID','')
                pmids =[p .strip ()for p in raw_pmid .split (',')]if raw_pmid else []

                #💡 Optimization details: Split GO terms into a list for easier multi-line printing later
                raw_go =attrs .get ('GO','')
                go_terms =[g .strip ()for g in raw_go .split (',')]if raw_go else []

                raw_orf =attrs .get ('ORF_evidence','')
                uniprot_id =self ._extract_uniprot_id (raw_orf )

                #Extract various former names
                #💡 Precisely extract official Gene Symbol, sorted by authority priority
                primary_symbol =""
                for preferred_key in ['RAP-DB Gene Symbol Synonym(s)','CGSNL Gene Symbol','Oryzabase Gene Symbol Synonym(s)']:
                    if preferred_key in attrs :
                        primary_symbol =attrs [preferred_key ].split (',')[0 ].strip ()
                        break 

                        #Extract all former names as a PubMed search pool
                symbols =[]
                for k ,v in attrs .items ():
                    if 'Synonym'in k or 'CGSNL'in k :
                        symbols .extend ([s .strip ()for s in v .split (',')])

                record ={
                "RAP_ID":rap_id ,
                "Transcript_ID":attrs .get ('ID',''),
                "Transcript_IDs":[attrs .get ('ID','')]if attrs .get ('ID','')else [],
                "Primary_Symbol":primary_symbol ,#🌟 New: Save the most authoritative name separately
                "MSU_IDs":rap_to_msu .get (rap_id ,[]),
                "UniProt_ID":uniprot_id ,
                "Chromosome":parts [0 ],
                "Start":int (parts [3 ]),
                "End":int (parts [4 ]),
                "Strand":parts [6 ],
                "TSS":int (parts [3 ])if parts [6 ]=="+"else int (parts [4 ]),
                "Location":f"{parts [0 ]}:{parts [3 ]}-{parts [4 ]} ({parts [6 ]})",
                "Description":attrs .get ('Note','Unknown'),
                "Synonyms":list (set (symbols )),
                "PMIDs":pmids ,
                "GO":go_terms ,
                "InterPro":attrs .get ('InterPro','N/A'),
                "Subcellular":attrs .get ('Subcellular location','N/A')
                }

                if rap_id in self .db :
                    existing =self .db [rap_id ]
                    transcript_id =record .get ("Transcript_ID")
                    if transcript_id and transcript_id not in existing .setdefault ("Transcript_IDs",[]):
                        existing ["Transcript_IDs"].append (transcript_id )
                    if not existing .get ("Transcript_ID")and record .get ("Transcript_ID"):
                        existing ["Transcript_ID"]=record ["Transcript_ID"]
                    if not existing .get ("UniProt_ID")and record .get ("UniProt_ID"):
                        existing ["UniProt_ID"]=record ["UniProt_ID"]
                    if existing .get ("InterPro")in ("","N/A")and record .get ("InterPro")not in ("","N/A"):
                        existing ["InterPro"]=record ["InterPro"]
                    if existing .get ("Subcellular")in ("","N/A")and record .get ("Subcellular")not in ("","N/A"):
                        existing ["Subcellular"]=record ["Subcellular"]
                    if not existing .get ("PMIDs")and record .get ("PMIDs"):
                        existing ["PMIDs"]=record ["PMIDs"]
                    continue 

                self .db [rap_id ]=record 
                self .alias_index [rap_id .lower ()]=rap_id 
                for sym in symbols :
                    self .alias_index [sym .lower ()]=rap_id 

        print (f"Gene annotation database loaded: {len (self .db )} rice genes indexed."if lang =="en"else f"GeneAnnotationDatabaseHasload：total included {len (self .db )} riceGene。")

    def query (self ,search_term ):
        if not search_term :
            return None 

        raw =search_term .strip ()
        candidates =[raw ]

        #Common variants: Remove Transcript suffix (e.g., .1) / Remove LOC_ prefix / Restore LOC_ prefix
        if "."in raw :
            candidates .append (raw .split (".")[0 ])

        no_loc =raw [4 :]if raw .upper ().startswith ("LOC_")else raw 
        candidates .append (no_loc )
        candidates .append (f"LOC_{no_loc }")

        #Additional completion: For each candidate, remove Transcript suffix again
        #and restore .1 variant (alias_index stores MSU ID with Transcript suffix)
        expanded =[]
        for item in candidates :
            expanded .append (item )
            if "."in item :
                expanded .append (item .split (".")[0 ])
            else :
            #Fallback to .1: LOC_Os01g01010 → LOC_Os01g01010.1
                expanded .append (item +".1")

                #Deduplicate and query in order
        seen =set ()
        for term in expanded :
            key =term .lower ().strip ()
            if not key or key in seen :
                continue 
            seen .add (key )
            target_rap =self .alias_index .get (key )
            if target_rap and target_rap in self .db :
                return self .db [target_rap ]
        return None 

    @staticmethod 
    def _extract_uniprot_id (raw_orf ):
        if not raw_orf :
            return ""
        match =re .search (r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b',raw_orf )
        if match :
            return match .group (1 )
        return ""

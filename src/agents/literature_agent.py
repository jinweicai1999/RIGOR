import os 
import re 
import ssl 
from urllib .error import HTTPError 
from datetime import datetime 
from Bio import Entrez ,Medline 
from dotenv import load_dotenv 
from src .utils .llm_client import chat_completion_text ,make_openai_client 

load_dotenv ()

class LiteratureAgent :
    def __init__ (self ):
        Entrez .email =os .getenv ("ENTREZ_EMAIL","jinweicai1999@gmail.com")
        if os .getenv ("BIO_ENTREZ_INSECURE_SSL","1")=="1":
            ssl ._create_default_https_context =ssl ._create_unverified_context 

        self .client =make_openai_client (
        api_key =os .getenv ("LLM_API_KEY"),
        base_url =os .getenv ("LLM_BASE_URL","https://api.longcat.chat/openai"),
        )
        self .model_name =os .getenv ("LLM_MODEL_NAME","LongCat-Flash-Lite")
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"
        self .pubmed_max_date =os .getenv ("PUBMED_MAX_DATE","").strip ()
        self ._pubmed_search_cache ={}
        self ._pubmed_fetch_cache ={}
        self ._pubmed_rate_limit_callback =None 

    def _txt (self ,zh ,en ):
        return en if self .lang =="en"else zh 

    def _get_esearch_kwargs (self ,term ,retmax ,sort ="relevance"):
        kwargs ={
        "db":"pubmed",
        "term":term ,
        "retmax":retmax ,
        "sort":sort ,
        }
        if self .pubmed_max_date :
            kwargs ["datetype"]="pdat"
            kwargs ["maxdate"]=self .pubmed_max_date 
        return kwargs 

    def search_pubmed_ids (self ,term ,retmax =5 ,sort ="relevance"):
        cache_key =(term ,retmax ,sort ,self .pubmed_max_date )
        if cache_key in self ._pubmed_search_cache :
            return list (self ._pubmed_search_cache [cache_key ])
        try :
            handle =Entrez .esearch (**self ._get_esearch_kwargs (term ,retmax ,sort ))
            record =Entrez .read (handle )
            handle .close ()
            ids =record .get ("IdList",[])
            self ._pubmed_search_cache [cache_key ]=list (ids )
            return ids 
        except HTTPError as e :
            if getattr (e ,"code",None )==429 :
                if self ._pubmed_rate_limit_callback :
                    self ._pubmed_rate_limit_callback ("search",term )
                print (self ._txt ("PubMed请求过于频繁，请耐心等待或稍后重试。","PubMed rate limit detected. Please wait patiently or retry later."))
                return []
            raise 

    def _parse_record_date (self ,rec ):
        for key in ("DP","DEP","EDAT","CRDT"):
            value =rec .get (key )
            if not value :
                continue 
            text =str (value )
            m =re .search (r"(\d{4})(?:[ /-](\d{1,2}))?(?:[ /-](\d{1,2}))?",text )
            if not m :
                continue 
            year =int (m .group (1 ))
            month =int (m .group (2 )or 1 )
            day =int (m .group (3 )or 1 )
            try :
                return datetime (year ,month ,day ).date ()
            except Exception :
                try :
                    return datetime (year ,1 ,1 ).date ()
                except Exception :
                    return None 
        return None 

    def _parse_date_text (self ,text ):
        if not text :
            return None 
        text =str (text ).strip ()
        if not text :
            return None 
        m =re .search (r"(\d{4})(?:[-/ ](\d{1,2}))?(?:[-/ ](\d{1,2}))?",text )
        if not m :
            return None 
        year =int (m .group (1 ))
        month =int (m .group (2 )or 1 )
        day =int (m .group (3 )or 1 )
        try :
            return datetime (year ,month ,day ).date ()
        except Exception :
            try :
                return datetime (year ,1 ,1 ).date ()
            except Exception :
                return None 

    def _record_within_cutoff (self ,rec ):
        if not self .pubmed_max_date :
            return True 
        try :
            cutoff =datetime .strptime (self .pubmed_max_date ,"%Y/%m/%d").date ()
        except Exception :
            return True 
        rec_date =self ._parse_record_date (rec )
        if rec_date is None :
            return True 
        return rec_date <=cutoff 

    def _convert_to_pmids (self ,raw_ids ):
        valid_pmids =[]
        for item in raw_ids :
            item =item .strip ()
            if not item :continue 

            if item .isdigit ():
                valid_pmids .append (item )
            elif item .upper ().startswith ("DOI:"):
                doi =item [4 :].strip ()
                print (self ._txt (f"detected DOI ({doi })，Currentlyconvert to PMID...",f"DOI detected ({doi }), converting to PMID..."))
                try :
                    ids =self .search_pubmed_ids (f"{doi }[AID]",retmax =1 ,sort ="relevance")
                    if ids :
                        valid_pmids .append (ids [0 ])
                except :
                    pass 
        return list (set (valid_pmids ))

    def _fetch_medline_by_ids (self ,id_list ):
        if not id_list :return ""
        cache_key =tuple (str (pid )for pid in id_list )
        if cache_key in self ._pubmed_fetch_cache :
            return self ._pubmed_fetch_cache [cache_key ]
        try :
            fetch_handle =Entrez .efetch (db ="pubmed",id =id_list ,rettype ="medline",retmode ="text")
            records =Medline .parse (fetch_handle )

            paper_text =""
            for rec in records :
                if not self ._record_within_cutoff (rec ):
                    continue 
                pmid =rec .get ("PMID","Unknown")
                title =rec .get ("TI","No title available")
                journal =rec .get ("JT","Unknown Journal")
                abstract =rec .get ("AB","No abstract available")
                doi ="N/A"
                for lid in rec .get ("LID",[]):
                    if "[doi]"in lid :doi =lid .replace (" [doi]","")

                paper_text +=f"- PMID: {pmid }\n  Title: {title }\n  Journal: {journal } (DOI: {doi })\n  Abstract: {abstract }\n\n"
            self ._pubmed_fetch_cache [cache_key ]=paper_text 
            return paper_text 
        except HTTPError as e :
            if getattr (e ,"code",None )==429 :
                if self ._pubmed_rate_limit_callback :
                    self ._pubmed_rate_limit_callback ("fetch",",".join (str (pid )for pid in id_list [:5 ]))
                print (self ._txt ("PMID FetchFailed: PubMed请求过于频繁，请耐心等待或稍后重试。","PMID fetch failed: PubMed rate limit detected. Please wait patiently or retry later."))
                return ""
            print (self ._txt (f"PMID FetchFailed: {e }",f"PMID fetch failed: {e }"))
            return ""
        except Exception as e :
            print (self ._txt (f"PMID FetchFailed: {e }",f"PMID fetch failed: {e }"))
            return ""

    def _build_smart_query (self ,gene_dict ):
        """🌟 Ultimate robust version: Smartly construct PubMed search query, specialized in handling various quirky naming conventions 🌟"""
        raw_terms =[]

        # 1. Extract all possible names (including Primary Symbol and Synonyms)
        primary_sym =gene_dict .get ('Primary_Symbol','').strip ()
        if primary_sym :raw_terms .append (primary_sym )

        for syn in gene_dict .get ('Synonyms',[]):
            syn =syn .strip ()
            # Filter out purely numeric or extremely short meaningless names
            if len (syn )>2 and not syn .isdigit ():
                raw_terms .append (syn )

        if not raw_terms :
            return ""

        final_terms =[]

        for term in raw_terms :
        # --- Extremely strict string cleaning (to prevent quote errors from breaking the search) ---
        # Remove leading/trailing quotes and extra spaces, convert multiple internal spaces to single spaces
            clean_term =term .replace ('"','').replace ("'","").strip ()
            clean_term =re .sub (r'\s+',' ',clean_term )

            if len (clean_term )<3 :continue 

            # A. Add original form
            final_terms .append (f'"{clean_term }"[Title/Abstract]')

            # B. If not containing Os, force-add Os (be careful not to add extra spaces or quotes)
            if not clean_term .lower ().startswith ('os'):
                final_terms .append (f'"Os{clean_term }"[Title/Abstract]')

                # C. 🌟 Ultimate family name extractor (specialized in handling PP2C 52 and RRM1-30) 🌟
                # Logic: Whether it's connected (PP2C52), hyphenated (RRM-30), or spaced (PP2C 52), always remove the trailing numeric part!
                # Regex explanation: \s*[-_]?\s*\d+$ matches trailing whitespace, hyphen/underscore, and a sequence of digits
            family_name =re .sub (r'\s*[-_]?\s*\d+$','',clean_term ).strip ()

            if family_name and family_name !=clean_term and len (family_name )>=3 :
                final_terms .append (f'"{family_name }"[Title/Abstract]')
                if not family_name .lower ().startswith ('os'):
                    final_terms .append (f'"Os{family_name }"[Title/Abstract]')

                    # 2. Deduplicate and assemble ultimate query
        unique_terms =list (set (final_terms ))

        # To prevent the query from becoming excessively large (e.g., some genes have dozens of aliases), we take at most the first 8 high-quality search terms
        term_str =" OR ".join (unique_terms [:8 ])
        query =f"({term_str }) AND (rice[Title/Abstract] OR Oryza sativa[Title/Abstract])"

        return query 

    def fetch_abstracts_route_a (self ,raw_pmids ,gene_info =None ):
        """Route A: Extract official PMID literature + expand search for potential functional literature"""
        all_blocks =[]

        # 1. Extract officially documented landmark literature
        print (self ._txt ("[Route A] Extracting official PMID literature...","[Route A] Extracting official PMID literature..."))
        clean_pmids =self ._convert_to_pmids (raw_pmids )
        if clean_pmids :
            official_text =self ._fetch_medline_by_ids (clean_pmids [:10 ])
            if official_text :
                block =f"\n\n### 【officialLiteraturearea：Hasrecorded coreFunction】 ###\n{official_text }"
                all_blocks .append (block )

                # 2. Expand search: Find potential other functions based on gene description
        if gene_info :
            print (self ._txt ("[Route A+] Expanding search for potential functional literature...","[Route A+] Expanding search for additional functional literature..."))
            additional_pmids ,has_cross_species =self ._search_target_gene_literature (gene_info ,use_symbol_search =True )
            # Remove already-has official PMIDs to avoid duplication
            additional_pmids =[pid for pid in additional_pmids if pid not in clean_pmids ]
            if additional_pmids :
                additional_text =self ._fetch_medline_by_ids (additional_pmids [:10 ])
                if additional_text :
                    species_note ="(Includes cross-species homologous protein literature)"if has_cross_species else ""
                    block =f"\n\n### 【expandLiteraturearea：potential additionalFunctionclue{species_note }】 ###\n{additional_text }"
                    all_blocks .append (block )

        return "".join (all_blocks )

    def _extract_protein_keywords_with_llm (self ,note ):
        """Use LLM to extract core protein name, common abbreviations, and function keywords from Note"""
        if not note or len (note .strip ())<5 :
            return []

        if self .lang =="en":
            prompt =f"""You are a bioinformatics expert. This is a functional description of a rice gene:

"{note }"

Identify:
1) core protein domains/modules (e.g., RING finger, C6HC, HSP70, G3BP)
2) major functional keywords (e.g., ubiquitin ligase, stress granule, chaperone)

Output only 2-3 core English search keywords, comma-separated, no explanation.
Example: RING finger, ubiquitin ligase
Example: G3BP, stress granule"""
        else :
            prompt =f"""you are biologicalInformationexpert in。This is a riceGene'sFunctionDescription：

"{note }"

please identify：
1. The core protein corresponding to this descriptionDomainorFunctionmodule（if RING finger、C6HC、HSP70、G3BP such as）
2. this protein'sPrimaryFunctionkeywords（if ubiquitin ligase、stress granule、chaperone such as）

onlyOutput 2-3 The top 5 core English search keywords，separated by commas，do not explain。
for example：RING finger, ubiquitin ligase
for example：G3BP, stress granule"""

        try :
            keywords_str =chat_completion_text (
            self .client ,
            model =self .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
            #Compatible with Chinese/Japanese comma, semicolon, line break, dunhao, etc., to prevent tokenization failure after translation
            tokens =re .split (r"[,\n;；，、|]+",keywords_str )
            keywords =[]
            for token in tokens :
                clean =token .strip ().strip ("\"'`")
                if not clean :
                    continue 
                    #Remove common numbering prefixes like "1) RING finger"
                clean =re .sub (r"^\d+[\)\.\-:\s]+","",clean ).strip ()
                #Discard overly long sentences, keep phrases that resemble retrieval terms
                if len (clean .split ())>8 :
                    continue 
                keywords .append (clean )
                #Deduplicate while preserving order
            seen =set ()
            dedup =[]
            for kw in keywords :
                k =kw .lower ()
                if k in seen :
                    continue 
                seen .add (k )
                dedup .append (kw )
            return dedup [:3 ]#Return at most 3 keywords
        except Exception as e :
            print (self ._txt (f"LLM keyword extractionFailed: {e }",f"LLM keyword extraction failed: {e }"))
            return []

    def _fallback_terms_from_symbol_and_note (self ,gene_info ):
        """When LLM fails to extract enough terms, fall back to symbol/synonym and lightweight keyword extraction from note."""
        terms =[]
        primary_sym =gene_info .get ('Primary_Symbol','').strip ()
        if primary_sym and len (primary_sym )>2 :
            terms .append (primary_sym )

        for syn in gene_info .get ('Synonyms',[])[:5 ]:
            syn =str (syn ).strip ()
            if len (syn )>2 and not syn .isdigit ():
                terms .append (syn )

                #Extract suspected protein keywords from description (e.g., G3BP, PP2C, RING finger)
        note =str (gene_info .get ('Description','')or '')
        if note :
            note_tokens =re .findall (r"[A-Za-z][A-Za-z0-9\-/]{2,}",note )
            for tk in note_tokens [:12 ]:
                low =tk .lower ()
                if low in {"similar","putative","protein","domain","family","unknown"}:
                    continue 
                terms .append (tk )

                #Deduplicate while preserving order
        seen =set ()
        out =[]
        for t in terms :
            key =t .lower ()
            if key in seen :
                continue 
            seen .add (key )
            out .append (t )
        return out 

    def _search_target_gene_literature (self ,gene_info ,use_symbol_search =False ):
        """Search Literature for the target Gene: prioritize rice, and if not found, search cross-species homologous protein Literature

        Args:
            gene_info: GeneInformation dictionary
            use_symbol_search: Whether to use Gene symbol directly for search (used in Has GeneRoute A)
        """
        import re 
        search_terms =[]

        if use_symbol_search :
        # Route A: Directly use Gene symbol and aliases
            primary_sym =gene_info .get ('Primary_Symbol','').strip ()
            if primary_sym and len (primary_sym )>2 :
                search_terms .append (primary_sym )
                # Take the first 2 aliases
            for syn in gene_info .get ('Synonyms',[])[:2 ]:
                syn =syn .strip ()
                if len (syn )>2 and not syn .isdigit ():
                    search_terms .append (syn )
        else :
        # Route B: Extract core keywords (e.g., G3BP) from Note using LLM
            note =gene_info .get ('Description','')
            if note :
                llm_keywords =self ._extract_protein_keywords_with_llm (note )
                if llm_keywords :
                    print (self ._txt (f"LLM identify keywords: {', '.join (llm_keywords )}",f"LLM keywords: {', '.join (llm_keywords )}"))
                    search_terms .extend (llm_keywords )
                    # Route B forces robust fallback: tokenize symbol/synonym/note
            if len (search_terms )<2 :
                fallback_terms =self ._fallback_terms_from_symbol_and_note (gene_info )
                if fallback_terms :
                    print (self ._txt (
                    f"LLMinsufficient keywords，fallback use symbol/synonym keywords: {', '.join (fallback_terms [:4 ])}",
                    f"LLM keywords insufficient, fallback to symbol/synonym terms: {', '.join (fallback_terms [:4 ])}"
                    ))
                    search_terms .extend (fallback_terms )

                    # Perform deduplication while preserving order again
        dedup_terms =[]
        seen_terms =set ()
        for term in search_terms :
            clean =str (term ).strip ()
            if len (clean )<3 :
                continue 
            key =clean .lower ()
            if key in seen_terms :
                continue 
            seen_terms .add (key )
            dedup_terms .append (clean )
        search_terms =dedup_terms 

        if not search_terms :
            return [],False 

        all_pmids =[]
        cross_species_pmids =[]

        # 2. First search rice Literature
        for term in search_terms [:2 ]:
            query =f'"{term }"[Title/Abstract] AND (rice[Title/Abstract] OR Oryza[Title/Abstract])'
            print (self ._txt (f"SearchriceLiterature: {term } ...",f"Searching rice literature: {term } ..."))
            try :
                all_pmids .extend (self .search_pubmed_ids (query ,retmax =8 ,sort ="relevance"))
            except Exception as e :
                print (self ._txt (f"SearchFailed: {e }",f"Search failed: {e }"))

                # 3. If rice Literature is insufficient, search other plant Literature
        if len (all_pmids )<3 and search_terms :
            main_keyword =search_terms [0 ]
            query =f'"{main_keyword }"[Title/Abstract] AND (plant[Title/Abstract] OR Arabidopsis[Title/Abstract] OR maize[Title/Abstract] OR wheat[Title/Abstract])'
            print (self ._txt (f"riceLiteratureInadequate，expandSearchplant homologous proteins: {main_keyword } ...",f"Rice literature insufficient. Expanding to plant homologs: {main_keyword } ..."))
            try :
                plant_pmids =self .search_pubmed_ids (query ,retmax =8 ,sort ="relevance")
                all_pmids .extend (plant_pmids )
                if plant_pmids :
                    cross_species_pmids .extend (plant_pmids )
            except Exception as e :
                print (self ._txt (f"plantLiteratureSearchFailed: {e }",f"Plant literature search failed: {e }"))

                # 4. If plant Literature is still insufficient, finally search human/mouse Literature
        if len (all_pmids )<3 and search_terms :
            main_keyword =search_terms [0 ]
            query =f'"{main_keyword }"[Title/Abstract] AND (human[Title/Abstract] OR mouse[Title/Abstract])'
            print (self ._txt (f"plantLiteraturestill insufficient，expandSearchmammalian homologous protein: {main_keyword } ...",f"Plant literature still insufficient. Expanding to mammalian homologs: {main_keyword } ..."))
            try :
                mammal_pmids =self .search_pubmed_ids (query ,retmax =8 ,sort ="relevance")
                all_pmids .extend (mammal_pmids )
                if mammal_pmids :
                    cross_species_pmids .extend (mammal_pmids )
            except Exception as e :
                print (self ._txt (f"mammalianLiteratureSearchFailed: {e }",f"Mammalian literature search failed: {e }"))

        return list (set (all_pmids )),len (cross_species_pmids )>0 

    def fetch_abstracts_route_b (self ,network_genes ,target_gene_info =None ,progress_callback =None ):
        all_text_blocks =[]
        print (self ._txt ("[Route B] Start Literature Search based on Symbol...","[Route B] Start symbol-based related literature retrieval..."))

        # 🌟 New: Prioritize searching Literature for the target Gene itself (based on Note and key InterPro)
        has_cross_species_lit =False 
        if target_gene_info :
            target_pmids ,has_cross_species_lit =self ._search_target_gene_literature (target_gene_info )
            if target_pmids :
                text =self ._fetch_medline_by_ids (target_pmids [:10 ])
                if text :
                    target_id =target_gene_info .get ('RAP_ID','')
                    target_sym =target_gene_info .get ('Primary_Symbol','')
                    species_note ="(Includes cross-species homologous protein literature)"if has_cross_species_lit else ""
                    block =f"\n\n### 【exclusiveLiteraturearea：targetGene {target_id } ({target_sym }) homologous proteinsFunctionclue{species_note }】 ###\n{text }"
                    all_text_blocks .append (block )

        total_genes =len (network_genes )
        max_literature_records =max (1 ,total_genes *10 )
        cumulative_literature_records =0 
        if progress_callback :
            progress_callback (
            done =0 ,
            total =total_genes ,
            gene_data =None ,
            retrieved_for_gene =0 ,
            cumulative =0 ,
            maximum =max_literature_records ,
            )

        def _rate_limit_callback (stage ,detail ):
            if progress_callback :
                progress_callback (
                done =idx -1 ,
                total =total_genes ,
                gene_data =gene_data ,
                retrieved_for_gene =0 ,
                cumulative =cumulative_literature_records ,
                maximum =max_literature_records ,
                rate_limited =True ,
                rate_limit_stage =stage ,
                rate_limit_detail =detail ,
                )

        previous_rate_limit_callback =self ._pubmed_rate_limit_callback 
        self ._pubmed_rate_limit_callback =_rate_limit_callback 
        try :
            iterable_network_genes =list (enumerate (network_genes ,start =1 ))
            for idx ,gene_data in iterable_network_genes :
                rap_id =gene_data .get ("RAP_ID","")
                primary_sym =gene_data .get ("Primary_Symbol","")
                raw_pmids =gene_data .get ("PMIDs",[])
                query =None 

                # Prioritize extracting its officially recorded PMIDs
                clean_pmids =self ._convert_to_pmids (raw_pmids )
                gene_fetch_ids =list (clean_pmids )
                target_count =10 -len (gene_fetch_ids )

                # If less than 10 papers, aggressively expand family Literature via PubMed
                if target_count >0 :
                    query =self ._build_smart_query (gene_data )
                    if query :
                        print (self ._txt (f"IntelligentSearch {primary_sym or rap_id }: {query [:50 ]}...",f"Smart search {primary_sym or rap_id }: {query [:50 ]}..."))
                        try :
                            for pid in self .search_pubmed_ids (query ,retmax =target_count ,sort ="relevance"):
                                if pid not in gene_fetch_ids :# Deduplicate
                                    gene_fetch_ids .append (pid )
                        except Exception as e :
                            print (self ._txt (f"SearchFailed: {e }",f"Search failed: {e }"))

                            # 🌟 Key breakthrough: Group Fetch and physically isolate!
                if gene_fetch_ids :
                    gene_fetch_ids =gene_fetch_ids [:10 ]
                    text =self ._fetch_medline_by_ids (gene_fetch_ids )
                    if text :
                    # Forcefully label this Literature with a "brand" so the large model cannot process sequentially!
                        block =f"\n\n### 【exclusiveLiteraturearea：comradeGene {rap_id } ({primary_sym })】 ###\n{text }"
                        all_text_blocks .append (block )
                        cumulative_literature_records +=len (gene_fetch_ids )

                if progress_callback :
                    progress_callback (
                    done =idx ,
                    total =total_genes ,
                    gene_data =gene_data ,
                    retrieved_for_gene =len (gene_fetch_ids )if gene_fetch_ids else 0 ,
                    cumulative =cumulative_literature_records ,
                    maximum =max_literature_records ,
                    )
        finally :
            self ._pubmed_rate_limit_callback =previous_rate_limit_callback 

        print (self ._txt (f"final acquisition {len (all_text_blocks )} eachLiteratureblock。",f"Retrieved {len (all_text_blocks )} literature block(s)."))
        return "".join (all_text_blocks )

    def summarize (self ,gene_info ,is_route_a ,abstracts_text ):
        target_gene =gene_info ['RAP_ID']
        note =gene_info .get ('Description','No')
        interpro =gene_info .get ('InterPro','No')
        go_terms =", ".join (gene_info .get ('GO',[]))if gene_info .get ('GO')else "No"
        subcell =gene_info .get ('Subcellular','No')

        # 🌟 Orphan Gene emergency kit: If no Literature exists and it's not a celebrity gene
        if not abstracts_text .strip ():
            if is_route_a :
                return "> ⚠️ Target Gene has PMID records but failed to successfully fetch Literature content."
            else :
            # Enter pure hardware blind guess mode
                orphan_prompt =f"""
                You are an extremely rigorous expert in rice molecular biology。
                targetGene '{target_gene }' a complete"orphan dark matter"。neitherHasexclusiveLiterature，neitherHasin 8456 SampletotalExpressionNetworkfind any comrades。
                
                please combine it only remains【hardware clues】：
                - FunctionDescription (Note)：{note }
                - InterPro Domain：{interpro }
                - GO Annotation：{go_terms }
                - subcellular localization：{subcell }
                
                perform pure bioinformatics logicInference。Must strictly use the following template：
                "although thisGenelacks coExpressionNetworkwith directLiteraturesupport，but based on its havingHas's [fill inDomainorGOName] feature，we reasonably infer targetGene {target_gene } possible..."
                （if evenDomainnoneHas，please truthfully state thisGenecurrently completelyNomethodInference，must rely on wetExperimentValidation。avoid rigid application。）
                """
                return chat_completion_text (
                self .client ,
                model =self .model_name ,
                messages =[{"role":"user","content":orphan_prompt }],
                temperature =0.1 
                )

                # ==========================================
                # Original Map-Reduce batch reading engine (keep unchanged below)
                # ==========================================
        raw_papers =abstracts_text .split ("- PMID:")

        target_gene =gene_info ['RAP_ID']
        interpro =gene_info .get ('InterPro','No')
        go_terms =", ".join (gene_info .get ('GO',[]))if gene_info .get ('GO')else "No"
        subcell =gene_info .get ('Subcellular','No')

        if is_route_a :
            desc =f"""CurrentlyanalyzeHasknown starGene '{target_gene }'。
【GeneitselfAnnotation】：
- FunctionDescription（Description）：{note }
- InterPro Domain：{interpro }
- GO Annotation：{go_terms }

Literaturedivided into two parts：
1. 【officialLiteraturearea】：Hasrecorded coreFunction（officialPMID）
2. 【expandLiteraturearea】：Potential additional insights from homologous protein searchFunction"""

            task =f"""Please strictly follow the followingStructureOutputreport：

**1. Hasrecorded coreFunction**
based on officialLiteraturearea，Summary '{target_gene }' 'sHasbiological knowledgeFunctionwith mechanism。

**2. potential additionalFunctionexplore**
extension-basedLiteraturearea（Especially cross-species homologous protein studies），InferencethisGenemay haveHasotherFunction。
if extendedLiteraturearea is empty，then indicates"temporaryNoadditionalFunctionclue"。"""
        else :
        # Use LLM to extract core protein keywords
            llm_keywords =[]
            if note and note !='No'and 'similar to'in note .lower ():
                llm_keywords =self ._extract_protein_keywords_with_llm (note )

                # Construct homology hints
            homology_hint =""
            if llm_keywords :
                homology_hint =f"\n  → core homologous protein：{llm_keywords [0 ]}（{', '.join (llm_keywords [1 :])}）"

            desc =f"""targetGene '{target_gene }' Functionunknown。
【targetGeneitselfAnnotation - highestPriority】：
- FunctionDescription（Description）：{note if note and note !='No'else ' No clear description'}{homology_hint }
- InterPro Domain：{interpro }
- GO Annotation：{go_terms }

following are co-ExpressioncomradeGene'sLiterature（used forAuxiliaryValidation）。"""

            task =f"""Please strictly follow the followingStructureOutputreport：

**1. comradeGeneFunctioninventory**
List each comrade individuallyGene'sFunction（format：**[GeneID]**：Literatureindicates...）

**2. targetGene {target_gene } potential mechanismInference**
【⚠️ Inferencerule】：
1. **highestPriority**：FunctionDescription（Description）field。if this field contains"Similar to XXX"，this clearly indicates thatGeneisXXXhomologous proteinsGene，Inferencemust first state"thisGene'sFunctiondescription indicates it isXXXhomologous proteins"
2. **timePriority**：InterPro DomainandGOAnnotationprovided moleculesFunctionEvidence
3. **AuxiliaryValidation**：totalExpressioncomradeGeneand cross species homologous proteinLiterature

ifSearchreached cross-speciesLiterature，please cite specificPMIDexplain species and（if"based on Arabidopsis/human homologous protein research（PMID: xxx）"）。

please use3-4give rigorous in one sentenceInference，must prioritize based onFunctionhomology describedInformation。"""

        prompt =f"""
        You are an extremely rigorous expert in rice molecular biology。you haveHasan exceptionally wide reading scope，hate academic fraud andNoend brain complementation。
        {desc }

        【LiteratureDatapool】：
        {abstracts_text [:500000 ]}

        【⚠️ InformationPriorityhierarchy - must strictly comply with】：
        1. **Geneitself Note Description**（if "Similar to XXX"）- this is the most directFunctionpoint to，Priorityhighest！
        2. **key InterPro Domain**（particularlyFunctionnatureDomain）- provide molecularFunctiondirect evidence
        3. **totalExpressioncomradeLiterature** - used forValidationand supplement the aboveInference，rather than dominantInference

        【extremely strict anti hallucinationWarning】：
        1. species isolation：clearly study wheat、Homologs in barley and other species，Please directly ignore thisConclusion！
        2. prioritize selfAnnotation：ifGene Note middleHasExplicitFunctionDescription（if "Similar to Ras-GTPase-activating protein"），must treat it asInferencecore，cannot be co-ExpressionLiteratureSubmerged！
        3. ExplicitTarget：Inferencepart must and only targetGene '{target_gene }' as subject！
        4. dehydration removal：not relevant，directly write in inventory"NoHaseffective comradesLiteraturesupport"。

        【task requirements】：
        {task }
        
        【⚠️ mandatory complete citation rules（never discounted）】：
        As a rigorous researcher，you must at the bottom of the report，create a【Appendix：all supporting this reportHasHasefficiencyLiterature】chapter。
        you must in Part I（comradeGeneFunctioninventory）mentioned in【each entryFunction】corresponding toLiterature，List all in complete format！
        strictly prohibit selecting only"core"list！prohibit omitting titles！
        
        must use the following format（EachLiteratureoccupy one line）：
        * 《Literaturefull title》 - journal name
        """

        try :
            print (self ._txt (" Model currently performing global Literature reasoning...","Model is performing global literature reasoning..."))
            return chat_completion_text (
            self .client ,
            model =self .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            return f"> ❌ large model reasoningInterfacedisconnect orReport Error: {e }"

    def fetch_abstracts_route_c (self ,deg_infos :list [dict ])->str :
        """ Route C: DEG cluster Literature Fetch, limited to rice species, no cross-species fallback allowed."""
        if not deg_infos :
            return ""

        all_pmids =[]
        for gene_data in deg_infos :
            rap_id =gene_data .get ("RAP_ID","")
            primary_sym =gene_data .get ("Primary_Symbol","")or rap_id 
            raw_pmids =gene_data .get ("PMIDs",[])
            gene_pmids =list (self ._convert_to_pmids (raw_pmids ))

            remaining =max (0 ,5 -len (gene_pmids ))
            if remaining >0 :
                query =self ._build_smart_query (gene_data )
                if query :
                    print (self ._txt (
                    f"[pathway C] Search DEG Gene {primary_sym } riceLiterature...",
                    f"[Route C] Searching rice literature for DEG gene {primary_sym }..."
                    ))
                    try :
                        for pid in self .search_pubmed_ids (query ,retmax =remaining ,sort ="relevance"):
                            if pid not in gene_pmids :
                                gene_pmids .append (pid )
                    except Exception as e :
                        print (self ._txt (f"SearchFailed: {e }",f"Search failed: {e }"))

            all_pmids .extend (gene_pmids [:5 ])

        merged_pmids =list (dict .fromkeys (all_pmids ))
        print (self ._txt (
        f"[pathway C] merged to obtain {len (merged_pmids )} unique chain PMID。",
        f"[Route C] Collected {len (merged_pmids )} unique PMIDs after merging."
        ))
        return self ._fetch_medline_by_ids (merged_pmids )

    def summarize_deg_cluster (self ,cluster_name :str ,deg_infos :list [dict ],abstracts_text :str )->str :
        gene_names =[info .get ("Primary_Symbol")or info .get ("RAP_ID")for info in deg_infos ]
        gene_names =[x for x in gene_names if x ]
        gene_list_text =", ".join (gene_names [:60 ])if gene_names else " No available gene name"

        gene_context =[]
        for info in deg_infos [:30 ]:
            rap_id =info .get ("RAP_ID","")
            symbol =info .get ("Primary_Symbol","")or rap_id 
            desc =info .get ("Description","")or " No clear description"
            interpro =info .get ("InterPro","")or "No"
            go_terms =", ".join (info .get ("GO",[])[:8 ])if info .get ("GO")else "No"
            gene_context .append (
            f"- {rap_id } ({symbol })\n"
            f"  Description: {desc }\n"
            f"  InterPro: {interpro }\n"
            f"  GO: {go_terms }"
            )
        gene_context_text ="\n".join (gene_context )if gene_context else "NoAnnotationInformation。"

        if not abstracts_text .strip ():
            prompt =f"""
            you are a rigorous riceSystembiological expert。
            nowHasa set in specificExperimentConditionsignificant changes occurring belowExpression'sGenemodule：【{cluster_name }】。
            this module contains the following coreGene：{gene_list_text }。

            currently notSearchto available riceLiteratureabstract，so you cannot fabricateLiteratureConclusion，can only rely onProjectInternalAnnotationperform conservativeInference。

            【InternalAnnotationInformation】：
            {gene_context_text }

            【task requirements】：
            Please strictly follow the followingStructureOutput：
            **1. module coreFunctionInference**
            only based on Description、GO and InterPro，SummaryThe biological process most likely involved by this module。

            **2. key driver factor（Hub Genes）**
            indicating whichGenefromAnnotationlooks more likely to be coreRegulationfactor or key execution factor，and explain basis。

            **3. Functionleave blank prediction**
            correctLiteratureandAnnotationare weakerGene，can only make conservative inference；if evidence is too weak，please write clearly“currentlyNoreliable methodInference，NeedExperimentValidation”。

            【Constraint】：
            - must clearly write out“current lack of direct riceLiteratureEvidence”。
            - no fabrication PMID、paper title orExperimentResult。
            """
        else :
            prompt =f"""
            you are a rigorous riceSystembiological expert。
            nowHasa set in specificExperimentConditionsignificant changes occurring belowExpression'sGenemodule：【{cluster_name }】。
            this module contains the following coreGene：{gene_list_text }。

            【LiteratureDatapool】（strictly limited to rice research）：
            {abstracts_text [:300000 ]}

            【InternalAnnotationInformation】：
            {gene_context_text }

            【task requirements】：
            please do not analyze in isolationAnalysiseachGene，Instead, look for their biological intersections，follow strictly the followingStructureOutput：
            **1. module coreFunctionInference**
            this set of differencesGeneas a whole，What physiological process is most likely to be co-involved in rice、metabolismPathwayor environmental response？please use 2-3 this sentence highly condensed。

            **2. key driver factor（Hub Genes）**
            According toLiteratureevidence and internalAnnotation，indicating this set ofGenewhich are inHascore knowledgeRegulationfactor、transcribeRegulationfactors or key enzymes。

            **3. Functionleave blank prediction**
            ForLiteraturerare or poorly knownGene，please combine its same moduleGene'sFunctionas well as its own InterPro Domainand GO，Infer its possible collaborative role。

            【Constraint】：
            - only riceLiteratureConclusion。
            - if a judgment lacks directLiteraturesupport，must clearly write“Literatureinsufficient evidence，onlyAnnotationdrive speculation”。
            - do not fabricate PMID or paper title。
            """

        try :
            return chat_completion_text (
            self .client ,
            model =self .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            return f"> ❌ DEG cluster LiteratureSummaryInterfaceReport Error: {e }"

    def summarize_ncrna (self ,mirna_id ,targets_info ,abstracts_text ,rna_type ="miRNA",extra_context =""):
        """ Specialized logic for inferring function of miRNA/phasiRNA/lncRNA treatments"""
        if not abstracts_text .strip ():
            return "> ⚠️ Failed to fetch Literature for the target gene."

        target_names =[d .get ("Primary_Symbol")or d ["RAP_ID"]for d in targets_info ]

        # Determine regulation mechanism based on RNA type
        if rna_type in ["miRNA","phasiRNA"]:
            mechanism ="【Degradation/Inhibition】These target genes (via sequence complementarity guiding mRNA degradation or translation inhibition)"
            inference_rule ="When performing inference, you must rigorously apply logical inversion: if a target gene promotes a certain function, then overexpression of the RNA will inhibit that function."
        else :# lncRNA
            mechanism ="This relationship involves cis/trans association, co-expression, or state correlation with these genes, but the regulatory direction may not be consistent across different pairs."
            inference_rule =(
            "You cannot assume that lncRNA always promotes the same direction as its associated gene during inference."
            "For positive correlation or activation evidence, cautious forward inference can be made;"
            "For negative correlation, bivalent states, or context-dependent evidence, it should be clearly stated as antagonistic, dynamic, or condition-dependent."
            )

        prompt =f"""
        You are a rigorous expert in rice epigenetics and RNA RegulationExpert。
        youCurrentlyEvaluation {rna_type } '{mirna_id }' biologicalFunction。

        【HasknownInformation】：
        '{mirna_id }' is a {rna_type }，it works through following mechanismRegulationtargetGene：{mechanism }

        【responseGenepool and multi-omics evidence】：
        Sequenceprediction or co-ExpressionAnalysisshow，following core proteinsGeneis '{mirna_id }' target ofGene：{', '.join (target_names )}。
        following are theseGene'sHasknownFunctionLiteratureabstract：

        {abstracts_text [:500000 ]}

        【supplement contextual evidence】：
        {extra_context or "No additional contextual evidence."}

        【⚠️ inference rules】：
        {inference_rule }

        Please strictly follow the followingStructureOutputreport：
        **1. targetGeneFunctioninventory**
        briefSummaryeach targetGene'sFunction。

        **2. {rna_type } {mirna_id } FunctionInference**
        based on targetGeneFunctionand {rna_type } 'sRegulationmechanism，Inferencethis RNA biologicalFunction。
        if evidence suggests differentTissue、different partner Or inconsistent directions under different conditions，must clearly write out“ConditionDependency/TissueDependency/may be bidirectionalRegulation”，Do not force complex evidence into one-sided promotion。

        【⚠️ mandatory complete citation rules】：
        Must list complete references at the bottom of the reportLiterature，format as：* 《Literaturefull title》 - journal name
        """
        try :
            return chat_completion_text (
            self .client ,
            model =self .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            return f"> ❌ large model reasoningInterfacedisconnect orReport Error: {e }"

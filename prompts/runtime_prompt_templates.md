# RIGOR runtime prompt expressions

This file is a line-anchored export of the literal prompt expressions in the
released Python source. The executable source of truth is the cited code;
the formerly empty `src/prompts/system_prompts.py` is not imported by the
current runtime. Dynamic evidence fields are shown as Python f-string
expressions and are populated only during a run. This export contains no
patient, user, API-key, or dataset content.

Extracted prompt expressions: 20.

## `apps/protein_function_pipeline.py:97` — `prompt`

```python
f"""You are a rigorous rice molecular biologist. Based ONLY on the homology literature below, infer the biological function of gene '{target }'.

[Gene info]
- RAP ID: {target }
- Description: {note }
- InterPro: {interpro }

[Homology literature]
{abstracts [:100000 ]}

Rules:
- Prioritize the gene's own Description (e.g. "Similar to XXX") as the strongest evidence.
- If cross-species literature is cited, clearly state the species.
- Output 3-5 sentences. Cite PMIDs where applicable.
- Do NOT fabricate titles, journals, or PMIDs.
"""
```

## `apps/protein_function_pipeline.py:114` — `prompt`

```python
f"""You are a rigorous expert in rice molecular biology。please**just**based on the following homologsGeneLiterature，InferenceGene '{target }' biologicalFunction。

【GeneInformation】
- RAP ID: {target }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【homologousGeneLiterature】
{abstracts [:100000 ]}

【Inferencerule】
- highestPriority：Functionhomology describedInformation（if "Similar to XXX"）
- ifHascross-speciesLiterature，explicitly label species
- Output 3-5 sentence，reference PMID
- strictly prohibit fabricatingLiteraturetitle、journal name、PMID
"""
```

## `apps/protein_function_pipeline.py:144` — `prompt`

```python
f"""You are a rigorous rice molecular biologist. Based on the GO enrichment analysis and co-expression gene literature below, infer the function of gene '{target }'.

[Gene info]
- RAP ID: {target }
- Description: {note }
- InterPro: {interpro }

[GO enrichment (Top 10)]
{go_text or "No significant GO terms"}

[Co-expression gene literature]
{abstracts [:100000 ]if abstracts else "No co-expression literature"}

Rules:
- First summarize partner gene functions, then infer target gene function.
- Use Description and InterPro as primary evidence; co-expression as supporting.
- Output 3-5 sentences. Cite PMIDs where applicable.
- Do NOT fabricate titles, journals, or PMIDs.
"""
```

## `apps/protein_function_pipeline.py:164` — `prompt`

```python
f"""You are a rigorous expert in rice molecular biology。based on the following GO EnrichmentAnalysiswith commonExpressioncomradeGeneLiterature，InferenceGene '{target }' 'sFunction。

【GeneInformation】
- RAP ID: {target }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【GO EnrichmentAnalysis (Top 10)】
{go_text or "No significant GO terms found."}

【totalExpressionGeneLiterature】
{abstracts [:100000 ]if abstracts else "No co-expressed genes or literature."}

【Inferencerule】
- first briefly inventory comradesGene'sFunction，againInferencetargetGeneFunction
- Functiondescription and InterPro forPrimaryEvidence，totalExpressionforAuxiliaryValidation
- Output 3-5 sentence，reference PMID
- strictly prohibit fabricatingLiteraturetitle、journal name、PMID
"""
```

## `apps/protein_function_pipeline.py:193` — `prompt`

```python
f"""You are a rice epigenetics expert. Based on the promoter peak scan data below, assess the regulatory state of gene '{target }' across tissues.

[Promoter epigenomic data]
{epigenomic_summary }

[Interpretation guide]
- ATAC/FAIRE/MNase open peaks → chromatin accessibility
- H3K4me3, H3K27ac, H3K9ac → active promoter marks
- H3K27me3, H3K9me2 → repressive marks
- Co-existence of active and repressive marks → bivalent/dynamic state

Output: In which tissues is the promoter likely active, repressed, or bivalent? 2-3 sentences.
"""
```

## `apps/protein_function_pipeline.py:207` — `prompt`

```python
f"""You are an expert in rice epigenetics。please refer to the followingStartatom peak ScanData，JudgmentGene '{target }' differentTissueinRegulationstate。

【Startsub-epigenomeData】
{epigenomic_summary }

【explanation prompt】
- ATAC/FAIRE/MNase open peak → chromatin open
- H3K4me3, H3K27ac, H3K9ac → ActiveStartsub marker
- H3K27me3, H3K9me2 → Inhibitionmarker
- ActivitywithInhibitionmarkers coexist → bivalent/dynamic state

Output：thisGeneStartin which subTissuemoderately active、biasInhibitionor bivalent？2-3 sentence。
"""
```

## `apps/protein_function_pipeline.py:699` — `prompt`

```python
f"""You are a rigorous rice molecular biologist. Infer the biological function of rice gene '{target_gene }' based on five independent evidence sources, and assess potential promoter regulatory states across tissues.

[Gene info]
- RAP ID: {target_gene }
- Description: {note }
- InterPro: {interpro }

[Evidence 1: Homology literature, highest priority]
{part1_abstracts [:100000 ]if part1_abstracts else "No homology literature"}

[Evidence 2: GO enrichment]
{chr (10 ).join ([f"- {go } (count={count })"for go ,count in go_enrichment [:10 ]])}

[Evidence 3: Co-expression literature]
{part3_abstracts [:100000 ]if part3_abstracts else "No co-expression literature"}

[Evidence 4: Promoter epigenomic evidence]
Peak scans in promoter region (TSS upstream 2000 bp, downstream 500 bp):
{epigenomic_summary }

[Evidence 5: Multi-transcript expression]
{expression_summary }

[Citable references from system retrieval]
{references_for_prompt }

Output format:
1) Function inference from homology literature
2) Function inference from domain and GO enrichment
3) Function inference from co-expression network
4) Promoter regulatory-state judgment from epigenomic evidence
5) Judgment from multi-transcript expression patterns
Final conclusion in 2-4 sentences.

Do not create an extra references section inside the generated text. Do not fabricate title, journal, or PMID.
"""
```

## `apps/protein_function_pipeline.py:736` — `prompt`

```python
f"""You are a rigorous expert in rice molecular biology。Please base your analysis on the following five independent evidence sources，comprehensiveInferencericeGene '{target_gene }' biologicalFunction，and judge itsStartsub in differentTissuepotential inRegulationstate。

【GeneBasicInformation】
- RAP ID: {target_gene }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【Evidence1：homologousGeneLiterature（highestPriority）】
following is theGeneStudies of homologous proteins in different speciesLiterature：
{part1_abstracts [:100000 ]if part1_abstracts else "No homologous gene literature."}

【Evidence2：GOEnrichmentAnalysis】
based on targetGeneits10commonExpressionGene'sGOAnnotationstatistics，discover following high-frequencyGO terms：
{chr (10 ).join ([f"- {go } (appear{count }time)"for go ,count in go_enrichment [:10 ]])}

【Evidence3：totalExpressionGeneLiterature】
following are target-relatedGenetotalExpressioncomradeGene'sFunctionLiterature：
{part3_abstracts [:100000 ]if part3_abstracts else "No co-expressed genes or literature."}

【Evidence4：Startsub epigeneticRegulationEvidence】
following is based on targetGene TSS upstream 2000 bp、downstream 500 bp Startsubregion of peak ScanResult：
{epigenomic_summary }

【Evidence5：manyTranscriptExpressionpattern】
following are targetsGenein organ/period and differentTreatmentConditionunderExpressionmode summary：
{expression_summary }

【citableLiteraturelist（SystemSearch）】
{references_for_prompt }

【Inferencerule】
1. **highestPriority**：Evidence1（homologousGeneLiterature）- this is the most directFunctionEvidence
2. **timePriority**：Evidence2（Domain+GOEnrichment）- provide molecularFunctionstatistical support
 3. **AuxiliaryValidation**：Evidence3（totalExpressionGene）- used forValidationand supplement the aboveInference
 4. **epigeneticRegulationlayer**：Evidence4（promoter peak）- used to judge thisGenedifferentTissuepotential activation in、Inhibitionor dynamicRegulationstate
 5. **Expressionpattern**：Evidence5（organ/Period + TreatmentExpressionpattern）- used to judge differentTranscript/Expressionmodel in organ、period andTreatmentConditionunderExpressiondifference

【epigenome interpretation hints】
- ATAC/FAIRE/MNase Open peaks generally support open chromatin
- H3K4me3、H3K27ac、H3K9ac、H3K4ac usually support activityStartatom
- H3K27me3、H3K9me2、H3K9me1 commonly indicatesInhibitionstate
- ifActivitywithInhibitionmarkers coexist，should be consideredTissuespecificity、stage differences or bivalent/dynamic state
- promoter state score higher，more support thisTissuemiddle promoter preferential activation；if state for dynamic_or_bivalent，explain different markers coexist

please follow the followingStructureOutputreport：

**Result1：homology-basedGene'sFunctionInference**
evidence-based1'sLiterature，SummarythisGeneHomologous proteins in other speciesHasknownFunction。ifHascross-speciesLiterature，please clearly label species。

**Result2：basedDomainandGOenrichedFunctionInference**
CombinationInterProDomainandGOEnrichmentAnalysis，from molecularFunctionangleInferencethisGenePossible biological processes involved。

**Result3：based on co-ExpressionNetwork'sFunctionInference**
combining evidence3，from co-Expressioncomrade'sFunctionInferencetargetGenepotentialFunction。byResult1forPrimarybase comprehensive judgment on。

**Result4：basedStartevidence from offspring epigenomeRegulationstate judgment**
combining evidence4，this indicatesGene promoter in whichTissuemore like open activation、Inhibitionsilence or dynamic switch state，and explain supporting your peak type。

**Result5：based on multipleTranscriptExpressionpattern judgment**
SummarythisGenedifferentTranscript/Expressionmodel in organ、period andTreatmentConditionunderExpressionpattern，point out whether there is obviousTissuePreference、Differences in developmental stage preference or stress response。

**finalConclusion**
integrating fiveResult，gives thisGeneMost probable biologicalFunction，and supplement its possibleTissuespecificityRegulationcharacteristics andExpressionmode features（2-4sentence）。

【reference rule】
can mention in main text PMID as argument source，But do not create additional content at the end of the main text“referenceLiterature/Appendix”chapter。
strictly prohibit fabricatingLiteraturetitle、journal name、PMID。
"""
```

## `apps/protein_function_pipeline.py:830` — `prompt`

```python
f"""You are a rigorous rice epigenetics expert. The input is a genomic coordinate region, not a confirmed gene.

[Input region]
- Original coordinate: {coord }
- Strand: {strand or '+'}
- Estimated TSS: {promoter ['tss']}
- Promoter region: {promoter ['promoter_label']}

[Constraints]
1. No gene ID/domain/GO/co-expression/literature is provided.
2. Do not treat this region as a known gene for functional prediction.
3. Only infer potential regulatory state from promoter-region epigenomic peaks.
4. Do not fabricate gene names, literature, PMIDs, domains, or expression conclusions.

[Promoter epigenomic evidence]
{epigenomic_summary }

[Tissue score summary]
{chr (10 ).join ([f"- tissue={item ['tissue']}, score={item ['score']}, state={item ['state']}"for item in tissue_scores ])if tissue_scores else "- No tissue score"}

Output format:
1) Region-level epigenomic observations
2) Region-level regulatory-state judgment
3) What can and cannot be concluded
Final conclusion in 2-4 sentences. Subject must be "this region" or "this promoter region".
"""
```

## `apps/protein_function_pipeline.py:857` — `prompt`

```python
f"""You are a rigorous expert in rice epigenetics。now given to you is notHasknownGene，but rather aGenegroup coordinate region。

【Inputregion】
- original coordinates: {coord }
- chain direction: {strand or '+'}
- predicted TSS: {promoter ['tss']}
- promoter region: {promoter ['promoter_label']}

【important limitation】
1. here noHasprovide truthGene ID、GeneAnnotation、Domain、GO、totalExpressionorLiteratureInformation。
2. You cannot directly treat this coordinate region asHasknownGeneTo doFunctionpredict。
3. you can only base on promoter regional epigenome peak，judge the differences in this region acrossTissuepotential inRegulationstate。
4. strictly prohibit fabricatingGenename、Literature、PMID、DomainorExpressionConclusion。

【promoter epigenomic evidence】
{epigenomic_summary }

【Tissuescoring summary】
{chr (10 ).join ([f"- tissue={item ['tissue']}, score={item ['score']}, state={item ['state']}"for item in tissue_scores ])if tissue_scores else "- No tissue score"}

Please strictly follow the followingStructureOutput：

**Result1：regional level epigenome observation**
ObjectiveSummarythis promoter regions in differentTissueOpen chromatin and histone modification features。

**Result2：regional levelRegulationstate judgment**
judge the differences in this region acrossTissuemore like in active、poised_or_weakly_active、repressed Or dynamic_or_bivalent，and explain basis。

**Result3：can and cannot be drawnConclusion**
clarify which onesConclusioncan be derived from this regionDataget，whichConclusioncannot obtain。
for example：can determine promoter Activitytrend；cannot directly assert specificGeneFunction。

**finalConclusion**
use 2-4 sentenceSummaryThe most reasonable explanation for this coordinate region。subject must be“this region”or“this promoter region”，cannot be written as“thisGene”。
"""
```

## `src/agents/gwas_agent.py:299` — `prompt`

```python
f"""
            You are a rigorous rice quantitative genetics and breeding expert.
            Evaluate field-trait effects and breeding trade-offs for target gene '{target_gene }' and its partner network.

            [Background: macro-network phenotype context]
            {macro_evidence }

            [Core measured data: population GWAS quantitative effects]
            {gwas_desc }
            {data_markdown }
            Effect_Size(BETA) is retained as the raw regression coefficient from the GWAS model.
            Do not interpret the sign of BETA as positive or negative regulation, promotion, suppression,
            increase, or decrease. Use P_Value and -log10(P) to discuss association strength only.
            Trait_Normalized_Evidence is assigned within each trait's own genome-wide P-value distribution:
            Strong = top 0.1%, Moderate = top 1%, Weak = top 5%. Treat Strong/Moderate as the main GWAS evidence
            and Weak as suggestive background-level evidence.

            Requirements:
            Output concise, evidence-based Markdown with this structure:

            # {target_gene } GWAS Trait Association Assessment
            ## 1. Significant Trait Association Overview
            ## 2. Candidate Breeding Marker Value
            ## 3. Interpretation Limits and Validation Needs
            Use specific P_Value values and Trait_Normalized_Evidence levels. Do not infer trait direction.
            """
```

## `src/agents/gwas_agent.py:326` — `prompt`

```python
f"""
            You are a rigorous and logically sound expert in rice quantitative genetics and breeding。
            youCurrentlytargetedGene '{target_gene }' its comradesNetworkConduct GWAS trait association analysis。
            
            【backgroundInformation：totalExpressionNetworkphenotypic characteristics of communities】
            {macro_evidence }
            
            【core measuredData：natural population GWAS quantifying effect】
            {gwas_desc }
            {data_markdown }
            *(note：Effect_Size(BETA) is the raw GWAS regression coefficient. Do not interpret its sign as positive/negative regulation, promotion/inhibition, increase/decrease, or yield gain/loss. Use P_Value and -log10(P) to discuss association strength only. If BETA is mentioned, use only its absolute magnitude as an effect-size reference.)*
            *(Trait_Normalized_Evidence is computed within each trait's own genome-wide P-value distribution: Strong = top 0.1%, Moderate = top 1%, Weak = top 5%. Only Strong/Moderate should be treated as core GWAS evidence; Weak is only suggestive.)*
            
            【task requirements】：
            please discard filler，strictly based on the aboveData，Outputfollowing format is exquisite Markdown report：

            # {target_gene } GWAS trait association assessment
            ## 1. significant trait association overview
            ## 2. candidate breeding marker value
            ## 3. interpretation limits and validation needs
            """
```

## `src/agents/gwas_agent.py:365` — `prompt`

```python
f"""
        youCurrentlyfor non-Coding RNA '{mirna_id }' Evaluationits yield cost in the field。
        
        【microscopic GWAS empirical evidence】：
        weScandone '{mirna_id }' 's【targetGene】GWAS association signals in natural populations：
        {data_markdown }
        
        【important interpretation rule】：
        Effect_Size(BETA) is retained only as the raw regression coefficient from the GWAS model.
        Do not flip its sign and do not interpret it as positive/negative regulation, promotion/inhibition, increase/decrease, or yield gain/loss.
        Use P_Value and -log10(P) to discuss association strength only. If BETA is mentioned, use only the absolute magnitude.
        
        pleaseOutputreport：
        **1. target trait association summary**
        （briefly describe significant target-gene GWAS associations by P value）
        
        **2. {mirna_id } candidate marker value and validation needs**
        （do not infer phenotype direction; describe only associated traits and the need for validation）
        """
```

## `src/agents/literature_agent.py:258` — `prompt`

```python
f"""You are a bioinformatics expert. This is a functional description of a rice gene:

"{note }"

Identify:
1) core protein domains/modules (e.g., RING finger, C6HC, HSP70, G3BP)
2) major functional keywords (e.g., ubiquitin ligase, stress granule, chaperone)

Output only 2-3 core English search keywords, comma-separated, no explanation.
Example: RING finger, ubiquitin ligase
Example: G3BP, stress granule"""
```

## `src/agents/literature_agent.py:270` — `prompt`

```python
f"""you are biologicalInformationexpert in。This is a riceGene'sFunctionDescription：

"{note }"

please identify：
1. The core protein corresponding to this descriptionDomainorFunctionmodule（if RING finger、C6HC、HSP70、G3BP such as）
2. this protein'sPrimaryFunctionkeywords（if ubiquitin ligase、stress granule、chaperone such as）

onlyOutput 2-3 The top 5 core English search keywords，separated by commas，do not explain。
for example：RING finger, ubiquitin ligase
for example：G3BP, stress granule"""
```

## `src/agents/literature_agent.py:553` — `orphan_prompt`

```python
f"""
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
```

## `src/agents/literature_agent.py:637` — `prompt`

```python
f"""
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
```

## `src/agents/literature_agent.py:735` — `prompt`

```python
f"""
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
```

## `src/agents/literature_agent.py:761` — `prompt`

```python
f"""
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
```

## `src/agents/literature_agent.py:818` — `prompt`

```python
f"""
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
```

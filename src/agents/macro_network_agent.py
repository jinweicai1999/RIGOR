import pandas as pd 
import os 
from src .project_paths import EXPRESSION_DIR ,ANNOTATION_DIR ,NETWORK_DIR 

class MacroNetworkAgent :
    def __init__ (self ,
    edge_file =str (NETWORK_DIR /"coexpression_network_edges.csv"),
    module_file =str (EXPRESSION_DIR /"module_assignment.csv"),
    mapping_file =str (ANNOTATION_DIR /"RAP-MSU.txt")):
        self .edge_file =edge_file 
        self .module_file =module_file 
        self .mapping_file =mapping_file 

        self .rap2msu ={}
        self .msu2rap ={}
        self .gene2mod ={}
        self .edges_df =None 
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"

        self ._load_data ()

    def _load_data (self ):
        print ("Loading co-expression network base (8456 samples)..."if self .lang =="en"else "Currently loading 8,456-sample co-expression network base...")

        #1. Load absolute orthodox ID mapping dictionary
        with open (self .mapping_file ,'r')as f :
            for line in f :
                if line .startswith ("RAP")or line .startswith ("None"):continue 
                parts =line .strip ().split ('\t')
                if len (parts )>=2 and parts [0 ]!="None"and parts [1 ]!="None":
                    rap =parts [0 ].strip ()
                    #Take the first transcript from MSU
                    msu_raw =parts [1 ].split (',')[0 ].strip ()
                    #Only truncate transcript suffix (.1), absolutely preserve LOC_ prefix!
                    msu =msu_raw .split ('.')[0 ]

                    self .rap2msu [rap ]=msu 
                    self .msu2rap [msu ]=rap 

                    #2. Load macro module affiliation
        mod_df =pd .read_csv (self .module_file )
        #Dictionary keys are strictly in LOC_Os... format
        # self.gene2mod = dict(zip(mod_df['Gene'].astype(str).str.strip(), mod_df['Module']))
        #🌟 Upgrade to support one-to-many list dictionary
        self .gene2mod ={}
        for _ ,row in mod_df .iterrows ():
            gene =str (row ['Gene']).strip ()
            mod =str (row ['Module']).strip ()
            if gene not in self .gene2mod :
                self .gene2mod [gene ]=[]
            if mod not in self .gene2mod [gene ]:
                self .gene2mod [gene ].append (mod )
                #3. Load 610,000 truth edges
        self .edges_df =pd .read_csv (self .edge_file )

        #🌟 Replace function signature to pass user's original input and complete registry
    def get_top_hub_partners (self ,target_raw ,gene_info ,annotator ,top_n =5 ):
        """Core logic: Same lineage + high interaction weight, supports precise MSU ID localization"""

        #1. Ultimate precise mapping: Translate input to most accurate MSU ID
        target_msu =None 

        #Priority A: If user input starts with LOC_, absolutely respect the original!
        if target_raw .upper ().startswith ("LOC_"):
            target_msu =target_raw 
        else :
        #Priority B: Check its associated MSU_IDs—use the one present in the macro network
            for msu in gene_info .get ('MSU_IDs',[]):
                test_msu =f"LOC_{msu }"if not msu .startswith ("LOC_")else msu 
                if test_msu in self .gene2mod :
                    target_msu =test_msu 
                    break 
                    #Fallback C: Use old dictionary for forced conversion
            if not target_msu :
                target_msu =self .rap2msu .get (gene_info ['RAP_ID'])

        if not target_msu :
            print ("Cannot map the input gene to a standard LOC_ identifier."if self .lang =="en"else "No standard LOC_ format found for this Gene in the mapping table.")
            return None ,0 ,0 ,[]

            # 2. Check Sect
        target_mods =self .gene2mod .get (target_msu )# Currently, a list is obtained, such as ['M0035', 'M0371']
        if not target_mods :
            print (f"Target gene {target_msu } is not assigned to any macro module in the 8456-sample network."if self .lang =="en"else f"targetGene {target_msu } not belong to any 8456 Samplemacro module。")
            return None ,0 ,0 ,[]

        target_mod_str =" & ".join (target_mods )
        if self .lang =="en":
            print (f"Macro network mapping: target gene ({target_msu }) assigned to module(s) [{target_mod_str }].")
        else :
            print (f"macroNetworkLocalization：targetGene ({target_msu }) belonging to module [{target_mod_str }]。")

            # 3. Find Close Allies (find all Has edges containing it among 610K pairs)
        mask_a =self .edges_df ['Gene_A']==target_msu 
        mask_b =self .edges_df ['Gene_B']==target_msu 
        target_edges =self .edges_df [mask_a |mask_b ].copy ()

        if target_edges .empty :
            print ("No interaction records found for this gene in the network edge table."if self .lang =="en"else "No interaction record found for this Gene in the macro Network edge table.")
            #Note: 4 Parameters are returned here
            return target_mods [0 ],0 ,0 ,[]

            #Extract pure MSU ID of allies
        target_edges ['Partner_MSU']=target_edges .apply (
        lambda row :row ['Gene_B']if row ['Gene_A']==target_msu else row ['Gene_A'],axis =1 
        )

        # 4. Define [Same Sect] and [Other Sect]
        def check_same_mod (partner_msu ):
            partner_mods =self .gene2mod .get (partner_msu ,[])
            #As long as any intersection exists between ally's sect and target Gene's sect, they are considered same sect!
            return bool (set (target_mods )&set (partner_mods ))

        target_edges ['Is_Same_Mod']=target_edges ['Partner_MSU'].apply (check_same_mod )
        same_mod_edges =target_edges [target_edges ['Is_Same_Mod']].copy ()
        other_mod_edges =target_edges [~target_edges ['Is_Same_Mod']].copy ()

        total_same_mod =len (same_mod_edges )

        #Return the first module to the LLM as a representative to avoid changing too much subsequent code
        target_mod =target_mods [0 ]
        # 5. Sort by interaction weight from high to low respectively
        same_mod_edges =same_mod_edges .sort_values (by ='Partial_correlation_coefficient',ascending =False )
        other_mod_edges =other_mod_edges .sort_values (by ='Partial_correlation_coefficient',ascending =False )

        valid_partners =[]

        #Define an internal function to recruit elites to avoid duplicating code
        def recruit_elites (edges_df ,is_same_mod ):
            for _ ,row in edges_df .iterrows ():
                if len (valid_partners )>=top_n :
                    break #Stop immediately once enough people are gathered

                partner_msu =row ['Partner_MSU']
                weight =row ['Partial_correlation_coefficient']

                #Translate back to RAP ID
                partner_rap =self .msu2rap .get (partner_msu )
                if not partner_rap :continue 

                #Check household registration
                friend_info =annotator .query (partner_rap )
                if friend_info :
                    primary_sym =friend_info .get ('Primary_Symbol','')
                    syns =friend_info .get ('Synonyms',[])

                    #Golden filter: must have Symbol alias
                    if not primary_sym and not syns :
                        continue 

                    valid_partners .append ({
                    "RAP_ID":partner_rap ,
                    "MSU_ID":partner_msu ,
                    "Primary_Symbol":primary_sym ,
                    "Synonyms":syns ,
                    "PMIDs":friend_info .get ('PMIDs',[]),
                    "Note":friend_info .get ('Note','No description'),
                    "GO":friend_info .get ('GO',[]),
                    "Correlation":weight ,
                    "Is_Same_Mod":is_same_mod #🌟 Tag them to distinguish whether they are from the same sect or externalaid
                    })

                    #First tier: prioritize recruiting from the same sect
        recruit_elites (same_mod_edges ,is_same_mod =True )
        same_mod_found =len (valid_partners )

        #Second tier: if fewer than 5 elites in the same sect, immediately hire high-paid talent from other sects with direct interaction relationships!
        if len (valid_partners )<top_n :
            recruit_elites (other_mod_edges ,is_same_mod =False )

            #Note: now 4 values are returned! An additional same_mod_found is added to write a beautiful report
        return target_mod ,total_same_mod ,same_mod_found ,valid_partners 

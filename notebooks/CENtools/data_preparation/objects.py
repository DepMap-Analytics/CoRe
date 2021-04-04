"""
                                    ---------------------
                                     C E N  -  T O O L S

                                           Objects
                                    ---------------------
"""

########################################################################################################################
########################################################################################################################
# Importing packages and path info

from paths import raw_data_path, curated_data_path, object_path

from data_preparation import curation
from data_preparation.curation import essentiality, expression, mutation, drug, cnv, fusion, model_ID_conversion
from data_preparation.curation import BEG, BNEG, adam_core, adam_context, stem_core, known_core, only_stem
from data_preparation.curation import sanger_mapping, broad_mapping

import pandas, pickle

########################################################################################################################
########################################################################################################################
# PROJECT CLASSES #

# SANGER #

class Sanger(object):
    def __init__(self, model_name):
        self.model_name = model_name
        self.sanger = None
        self.EXPRESSION, self.MUTATION, self.CNV, self.FUSION, self.DRUG, self.ESSENTIALITY =\
            False,False,False,False,False,False
        self.mutations, self.oncogenic_mutations, self.hotspot_mutations = None, None, None
        self.mutated_genes, self.oncogenic_mutated_genes, self.hotspot_mutated_genes = None, None, None
        self.expression = None
        self.cnv = None
        self.fusion = None
        self.growth_property = None
        self.msi = None
        self.ploidy = None
        self.drug = None
        self.tissue, self.site = None, None
        self.cancer, self.cancer_subtype = None, None
        self.broad = None
        self.essentiality = None

    def model_info(self, model_df, conversion_table):

        if self.model_name in list(model_df.model_name):
            x = model_df[model_df.model_name == self.model_name]
            self.growth_property = x.growth_properties.values[0]
            self.ploidy = x.ploidy.values[0]
            self.site = x.sample_site.values[0]
            self.sanger = x.index.values[0]

        if self.sanger in list(conversion_table.sanger.values):
            print(self.sanger)
            y = conversion_table[(conversion_table.sanger == self.sanger) |
                                 (conversion_table.sanger_name == self.model_name)]
            self.msi = list(y.MSI_MSS.values)[0] if list(y.MSI_MSS.values)[0] != "None" else None
            self.broad = list(y.broad.values)[0] if list(y.broad.values)[0] != "None" else None
            self.tissue = list(y.tissue.values)[0]
            self.cancer = list(y.cancer.values)[0]
            self.cancer_subtype = list(y["cancer subtype"].values)[0] if list(y["cancer subtype"].values)[0] != "None" else None

    def expression_info(self, expression_df):

        if self.model_name in list(expression_df.columns):
            self.expression = expression_df[self.model_name]
            self.EXPRESSION = True

    def mutation_info(self, mutation_df):

        if self.model_name in list(mutation_df.Model_Name.values):
            x = mutation_df[mutation_df.Model_Name == self.model_name]
            self.mutations = list(x.Mutation)
            self.mutated_genes = [m.split(" p.")[0] for m in self.mutations]
            self.oncogenic_mutations = list(x[x.Oncogenic_Mutation].Mutation)
            self.oncogenic_mutated_genes = list(set([m.split(" p.")[0] for m in list(x[x.Oncogenic_Gene].Mutation)]))
            self.hotspot_mutations = list(x[x.HotSpot].Mutation)
            self.hotspot_mutated_genes = list(set([m.split(" p.")[0] for m in list(x[x.HotSpot].Mutation)]))
            self.MUTATION = True

    def cnv_info(self, cnv_df):

        if self.model_name in list(cnv_df.columns):
            self.cnv = cnv_df[self.model_name]
            self.CNV = True

    def fusion_info(self, fusion_df):

        if self.model_name in list(fusion_df.index):
            x = fusion_df.loc[self.model_name]
            if type(x) == pandas.DataFrame:
                self.fusion = list(set([row["5_Prime"] + "-" + row["3_Prime"] for ind, row in x.iterrows()]))
            else:
                self.fusion = list(x["5_Prime"] + "-" + x["3_Prime"])
            self.FUSION = True

    def drug_info(self, drug_df):

        if self.model_name in list(drug_df.CELL_LINE.values):
            self.drug = drug_df[drug_df.CELL_LINE == self.model_name][["DRUG_NAME", "IC50", "LN_IC50", "Z_SCORE"]]
            self.DRUG = True

    def essentiality_info(self, essentiality_df, normalisation):
        if self.model_name in list(essentiality_df.columns):
            self.essentiality = {normalisation: essentiality_df[self.model_name]}
            self.ESSENTIALITY = True

# BROAD #

class Broad(object):
    def __init__(self, broad):
        self.broad = broad
        self.model_name = None
        self.EXPRESSION, self.MUTATION, self.CNV, self.FUSION, self.DRUG, self.ESSENTIALITY =\
            False,False,False,False,False,False
        self.mutations, self.oncogenic_mutations, self.hotspot_mutations = None, None, None
        self.mutated_genes, self.oncogenic_mutated_genes, self.hotspot_mutated_genes = None, None, None
        self.expression = None
        self.cnv = None
        self.drug = None
        self.fusion = None
        self.site, self.tissue, self.growth_property = None, None, None
        self.cancer, self.cancer_subtype = None, None
        self.sanger = None
        self.essentiality = None

    def model_info(self, model_df, conversion_table):

        if self.broad in list(model_df.index):
            x = model_df.loc[self.broad]
            self.model_name = x.stripped_cell_line_name
            self.growth_property = x.culture_type
            self.site = x.sample_collection_site

        if self.broad in list(conversion_table.broad.values):
            y = conversion_table[conversion_table.broad == self.broad]
            self.sanger = list(y.sanger.values)[0] if list(y.sanger.values)[0] != "None" else None
            self.tissue = list(y.tissue.values)[0]
            self.cancer = list(y.cancer.values)[0]
            self.cancer_subtype = list(y["cancer subtype"].values)[0] if list(y["cancer subtype"].values)[0] != "None" else None

    def expression_info(self, expression_df):

        if self.broad in list(expression_df.columns):
            self.expression = expression_df[self.broad]
            self.EXPRESSION = True

    def mutation_info(self, mutation_df):

        if self.broad in list(mutation_df.index):
            x = mutation_df.loc[self.broad]
            self.mutations = list(x.Mutation)
            self.mutated_genes = [m.split(" p.")[0] for m in self.mutations]
            self.oncogenic_mutations = list(x[x.Oncogenic_Mutation].Mutation)
            self.oncogenic_mutated_genes = [m.split(" p.")[0] for m in list(x[x.Oncogenic_Gene].Mutation)]
            self.hotspot_mutations = list(x[x.HotSpot].Mutation)
            self.hotspot_mutated_genes = list(set([m.split(" p.")[0] for m in list(x[x.HotSpot].Mutation)]))
            self.MUTATION = True

    def cnv_info(self, cnv_df):

        if self.broad in list(cnv_df.columns):
            self.cnv = cnv_df[self.broad]
            self.CNV = True

    def fusion_info(self, fusion_df):

        if self.broad in list(fusion_df.index):
            x = fusion_df.loc[self.broad]
            if type(x) == pandas.DataFrame:
                self.fusion = list(set([row["5_Prime"] + "-" + row["3_Prime"] for ind, row in x.iterrows()]))
            else:
                self.fusion = list(x["5_Prime"] + "-" + x["3_Prime"])
            self.FUSION = True

    def drug_info(self, drug_df):

        if self.broad in list(drug_df.BROAD):
            self.drug = drug_df[drug_df.BROAD == self.broad][["DRUG_NAME", "IC50", "LN_IC50", "Z_SCORE"]]
            self.DRUG = True

    def essentiality_info(self, essentiality_df, normalisation):
        if self.broad in list(essentiality_df.columns):
            self.essentiality = {normalisation: essentiality_df[self.broad]}
            self.ESSENTIALITY = True

# INTEGRATED #

class Integrated(object):
    def __init__(self, model_name):
        self.model_name = model_name
        self.broad = None
        self.sampled_project = None
        self.EXPRESSION, self.MUTATION, self.CNV, self.FUSION, self.DRUG, self.ESSENTIALITY =\
            False,False,False,False,False,False
        self.mutations, self.oncogenic_mutations, self.hotspot_mutations = None, None, None
        self.mutated_genes, self.oncogenic_mutated_genes, self.hotspot_mutated_genes = None, None, None
        self.expression = None
        self.cnv = None
        self.fusion = None
        self.growth_property = None
        self.msi = None
        self.ploidy = None
        self.drug = None
        self.tissue, self.site = None, None
        self.cancer, self.cancer_subtype = None, None
        self.essentiality = None

    def model_info(self, model_df, conversion_table, sampled_project):

        self.sampled_project = sampled_project

        if self.sampled_project == "BROAD":
            self.broad = self.model_name
            if self.model_name in list(model_df.index):
                x = model_df.loc[self.model_name]
                self.growth_property = x.culture_type
                self.site = x.sample_collection_site

        elif self.sampled_project == "SANGER":
            self.broad = conversion_table[conversion_table.sanger_name == self.model_name]["broad"].values[0]
            if self.model_name in list(model_df.model_name):
                x = model_df[model_df.model_name == self.model_name]
                self.growth_property = x.growth_properties
                self.ploidy = x.ploidy
                self.site = x.sample_site

        if self.broad in list(conversion_table.broad.values):
            y = conversion_table[conversion_table.broad == self.broad]
            self.msi = list(y.MSI_MSS.values)[0]
            self.tissue = list(y.tissue.values)[0]
            self.cancer = list(y.cancer.values)[0]
            self.cancer_subtype = list(y["cancer subtype"].values)[0] if list(y["cancer subtype"].values)[0] != "None" else None

    def expression_info(self, expression_df):

        if self.model_name in list(expression_df.columns):
            self.expression = expression_df[self.model_name]
            self.EXPRESSION = True

    def mutation_info(self, mutation_df, sampled_project):
        if sampled_project == "SANGER":
            if self.model_name in list(mutation_df.Model_Name.values):
                x = mutation_df[mutation_df.Model_Name == self.model_name]
                self.mutations = list(x.Mutation)
                self.mutated_genes = [m.split(" p.")[0] for m in self.mutations]
                self.oncogenic_mutations = list(x[x.Oncogenic_Mutation].Mutation)
                self.oncogenic_mutated_genes = list(
                    set([m.split(" p.")[0] for m in list(x[x.Oncogenic_Gene].Mutation)]))
                self.hotspot_mutations = list(x[x.HotSpot].Mutation)
                self.hotspot_mutated_genes = list(set([m.split(" p.")[0] for m in list(x[x.HotSpot].Mutation)]))
                self.MUTATION = True
        elif sampled_project == "BROAD":
            if self.model_name in list(mutation_df.DepMap_ID.values):
                x = mutation_df[mutation_df.DepMap_ID == self.model_name]
                self.mutations = list(x.Mutation)
                self.mutated_genes = [m.split(" p.")[0] for m in self.mutations]
                self.oncogenic_mutations = list(x[x.Oncogenic_Mutation].Mutation)
                self.oncogenic_mutated_genes = list(
                    set([m.split(" p.")[0] for m in list(x[x.Oncogenic_Gene].Mutation)]))
                self.hotspot_mutations = list(x[x.HotSpot].Mutation)
                self.hotspot_mutated_genes = list(set([m.split(" p.")[0] for m in list(x[x.HotSpot].Mutation)]))
                self.MUTATION = True

    def cnv_info(self, cnv_df):

        if self.model_name in list(cnv_df.columns):
            self.cnv = cnv_df[self.model_name]
            self.CNV = True

    def fusion_info(self, fusion_df):

        if self.model_name in list(fusion_df.index):
            x = fusion_df.loc[self.model_name]
            if type(x) == pandas.DataFrame:
                self.fusion = list(set([row["5_Prime"] + "-" + row["3_Prime"] for ind, row in x.iterrows()]))
            else:
                self.fusion = list(x["5_Prime"] + "-" + x["3_Prime"])
            self.FUSION = True

    def drug_info(self, drug_df):

        if self.model_name in list(drug_df.CELL_LINE):
            self.drug = drug_df[drug_df.CELL_LINE == self.model_name][["DRUG_NAME", "IC50", "LN_IC50", "Z_SCORE"]]
            self.DRUG = True

    def essentiality_info(self, essentiality_df, normalisation):
        if self.model_name in list(essentiality_df.columns):
            self.essentiality = {normalisation: essentiality_df[self.model_name]}
            self.ESSENTIALITY = True


########################################################################################################################
########################################################################################################################
# PROJECT OBJECT DICTIONARIES & SERIALISATION #
# Expression, CNV, Fusion and Drug were not readed!

def main_project(project, normalisation):

    selected_class = Sanger if project == "SANGER" else (Broad if project == "BROAD" else Integrated)
    obj_dict = {}
    for model in essentiality(project, normalisation).columns:
        print(model)
        obj = selected_class(model)
        if project != "INTEGRATED":
            if project == "SANGER": obj.model_info(model_df=sanger_mapping(),
                                                   conversion_table = model_ID_conversion())
            elif project == "BROAD": obj.model_info(model_df=broad_mapping(),
                                                    conversion_table = model_ID_conversion())
        else:
            if model[:3] == "ACH":
                obj.model_info(model_df=broad_mapping(), conversion_table=model_ID_conversion(),
                               sampled_project = "BROAD")
            else:
                obj.model_info(model_df=sanger_mapping(), conversion_table=model_ID_conversion(),
                               sampled_project = "SANGER")
        print("mutation")
        #obj.expression_info(expression(project))
        if project != "INTEGRATED":
            obj.mutation_info(mutation(project))
        else:
            if model[:3] == "ACH":
                obj.mutation_info(mutation("INTEGRATED"), "BROAD")
            else:
                obj.mutation_info(mutation("INTEGRATED"), "SANGER")
        #obj.cnv_info(cnv(project))
        #obj.fusion_info(fusion(project))
        #obj.drug_info(drug(project))
        print("essentiality")
        #obj.essentiality_info(essentiality(project, normalisation), normalisation)
        print("dict")
        obj_dict[model] = obj
    print("pickle")
    pickle.dump(obj_dict, open(object_path + project + "_" + normalisation + "_OBJECT.pkl", "wb"))

    return obj_dict


########################################################################################################################
########################################################################################################################
# GENE CLASS #

# GENE #

class Gene(object):
    def __init__(self, gene):
        self.gene = gene
        self.BEG = False
        self.BNEG = False
        self.location = None
        self.ADAM_Core = False
        self.ADAM_Context = False
        self.pancancer_tractability = None
        self.pancancer_priority = None
        self.stem = False
        self.core_essential_process = False
        self.oncogenic, self.hotspot_carrier = False, False

    def add_bagel_info(self, essential_genes, non_essential_genes):

        self.BEG = True if self.gene in essential_genes else False
        self.BNEG = True if self.gene in non_essential_genes else False

    def add_location_info(self, location_df):

        self.location = location_df.loc[self.gene] if self.gene in location_df.index else None

    def add_adam_info(self, adam_core, adam_context, adam_df):

        self.ADAM_Core = True if self.gene in adam_core else False
        self.ADAM_Context = True if self.gene in adam_context else False
        self.pancancer_tractability = adam_df.loc[self.gene]["Bucket"]\
            if self.gene in adam_df.index else None
        self.pancancer_priority = adam_df.loc[self.gene]["Score"]\
            if self.gene in adam_df.index else None

    def add_stem_cell_info(self, only_stem_essential_genes):

        self.stem = True if self.gene in only_stem_essential_genes else False

    def add_core_process_genes(self, core_process_genes):

        self.core_essential_process = True if self.gene in core_process_genes else False


########################################################################################################################
########################################################################################################################
# OBJECT DICTIONARIES AND SERIALISATION #

def main_gene(gene_list, gene_list_name):
    gene_obj = {}
    for gene in gene_list:
        obj = Gene(gene)
        obj.add_bagel_info(essential_genes = BEG, non_essential_genes = BNEG)
        obj.add_location_info(location_df = annotation)
        obj.add_adam_info(adam_core = adam_core, adam_context = adam_context,
                          adam_df= adam_context_score)
        obj.add_stem_cell_info(only_stem_essential_genes = only_stem)
        obj.add_core_process_genes(core_process_genes=known_core)
        gene_obj[gene] = obj

    pickle.dump(gene_obj, open(object_path + "INITIAL_" + gene_list_name + "_GENE_OBJ.pkl", "wb"))

"""
main_gene(set(set(essentiality("SANGER", "MANUEL").index).intersection(
    set(essentiality("BROAD", "MANUEL").index))), "BS")
main_gene(set(essentiality("INTEGRATED", "MANUEL").index), "INTEGRATED")
"""

########################################################################################################################
########################################################################################################################
# DE-SERIALISATION OF THE OBJECTS #

def deserialisation_project(project):
    obj = pickle.load(open(object_path + project + '_MANUEL_OBJECT.pkl', "rb"))
    return obj

def deserialisation_gene(gene_list_gene):
    gene_obj = pickle.load(open(object_path + "INITIAL_" + gene_list_gene + "_GENE_OBJ.pkl", "rb"))
    return gene_obj


########################################################################################################################
########################################################################################################################




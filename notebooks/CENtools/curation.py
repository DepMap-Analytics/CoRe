"""
                                    ---------------------
                                     C E N  -  T O O L S

                                      Data  Preparation
                                    ---------------------
"""
########################################################################################################################
########################################################################################################################
# IMPORTING #

from paths import raw_data_path, curated_data_path, object_path
import necessary_packages

import pandas, pickle

########################################################################################################################
########################################################################################################################
# GENE ANNOTATIONS #

# BAGEL
BEG = [line.strip() for line in open(curated_data_path + "Curated_Bagel_Essential_Genes.txt").readlines()]
BNEG = [line.strip() for line in open(curated_data_path + "Curated_Bagel_Non_Essential_Genes.txt").readlines()]

# ADaM
adam_core = set(pandas.read_csv(curated_data_path + "ADAM_Core_Fitness_Genes.csv", index_col=0).index)
adam_context = set(pandas.read_csv(curated_data_path + "ADAM_Context_Genes.csv", header=None, index_col=1).index)
adam_context_score = pandas.read_csv(curated_data_path + "Adam_Context_Score_df.csv", index_col=0)
adam_pancancer_context = adam_context.intersection(set(adam_context_score[adam_context_score.Bucket <= 6].index))

# STEM
stem_core = set(pandas.read_csv(curated_data_path + "STEM_Cores.csv", index_col=0, header=None).index)
only_stem = set(pandas.read_csv(curated_data_path + "Cores_for_pluoripotency_only.csv", index_col=0, header=None).index)

# IHRY PAPER -- KNOWN CORE IN CORE ESSENTIAL PROCESSES
known_core = set(pandas.read_csv(curated_data_path + "Known_cores_for_all_Ihry_paper.csv", index_col=0).index)

# Dede et al.
dede_core = set(pandas.read_csv(curated_data_path + "Dede_etal_Core_genes.csv", index_col = 1).index)

# FAMILY
annotation = pandas.read_csv(curated_data_path + "Kinase_TF_Surface_SLCs.txt", sep="\t", index_col=0)

########################################################################################################################
########################################################################################################################
# ESSENTIALITY #

# Read curated essentiality tables

def essentiality(project, normalisation):
    if project == "SANGER":
        if normalisation == "BAYESIAN":
            sanger_essentiality = pandas.read_csv(curated_data_path + "sanger_essentiality.csv", index_col=0)
            return sanger_essentiality
        elif normalisation == "FC":
            sanger_cor_fc = pandas.read_csv(curated_data_path + "sanger_cor_FC_essentiality.csv", index_col=0)
            return sanger_cor_fc
        elif normalisation == "MANUEL":
            return pandas.read_csv(curated_data_path + "scaled_norm_sanger_cor_FC_essentiality.csv", index_col=0)

    elif project == "BROAD":
        if normalisation == "BAYESIAN":
            broad_essentiality = pandas.read_csv(curated_data_path + "broad_essentiality.csv", index_col=0)
            return broad_essentiality
        elif normalisation == "FC":
            broad_cor_fc = pandas.read_csv(curated_data_path + "broad_cor_FC_essentiality_19Q2.csv", index_col=0)
            return broad_cor_fc
        elif normalisation == "MANUEL":
            return pandas.read_csv(curated_data_path + "scaled_norm_broad_cor_FC_essentiality.csv", index_col=0)

    elif project == "INTEGRATED":
        if normalisation == "BAYESIAN":
            return "NO DATA"
        elif normalisation == "FC":
            integrated_cor_fc = pandas.read_csv(curated_data_path + "CERES_scaled_depFC.csv", index_col=0)
            return integrated_cor_fc
        elif normalisation == "MANUEL":
            return pandas.read_csv(curated_data_path + "scaled_norm_integrated_cor_FC_essentiality.csv", index_col=0)


# Coessentiality

def co_essentiality(project):
    if project == "SANGER":
        return sanger_co_essentiality
    elif project == "BROAD":
        return broad_co_essentiality
    elif project == "INTEGRATED":
        return integrated_co_essentiality


########################################################################################################################
########################################################################################################################
# CELL LINE MAPPING AND DISEASE INFO #

# Read curated conversion tables
def model_ID_conversion():
    ess_common_df = pandas.read_csv(curated_data_path + "mapping.csv")
    ess_common_df = ess_common_df.rename(columns={"Sanger_ID":"sanger", "BROAD": "broad",
                                                  "SANGER_NAME": "sanger_name",
                                                  "Broad_model_name": "broad_name",
                                                  "Tissue": "tissue", "Cancer_Type": "cancer",
                                                  "Cancer.Subtype": "cancer subtype"})
    return ess_common_df


def sanger_mapping():
    map = pandas.read_csv(raw_data_path + "sanger/sanger_model_list.csv", index_col=0)
    return map


def broad_mapping():
    map = pandas.read_csv(raw_data_path + "broad_19Q3/CCLE_sample_info_19Q3.csv", index_col=0)
    return map


########################################################################################################################
########################################################################################################################
# MUTATION INFO #

# Read curated mutation tables

def mutation(project):
    if project == "SANGER":
        mutation = pandas.read_csv(
            curated_data_path + "sanger_mutation.csv", index_col=0)
    elif project == "BROAD":
        mutation = pandas.read_csv(
            curated_data_path + "broad_19Q2_mutation.csv", index_col=0)
    elif project == "INTEGRATED":
        mutation = pandas.read_csv(
            curated_data_path + "integrated_mutation.csv", index_col=0)
    else:
        mutation = None
    return mutation


########################################################################################################################
########################################################################################################################
#  EXPRESSION - FPKM/TPM INFO #

# Read the curated expression files

def expression(project):
    if project == "SANGER":
        expression = pandas.read_csv(curated_data_path + "sanger_expression.csv", index_col=0)
    elif project == "BROAD":
        expression = pandas.read_csv(curated_data_path + "broad_19Q2_expression.csv", index_col=0)
    elif project == "INTEGRATED":
        expression = pandas.read_csv(curated_data_path + "integrated_expression.csv", index_col=0)
    return expression


########################################################################################################################
########################################################################################################################
# FUSION INFO #

# Read the curated data

def fusion(project):
    if project == "SANGER":
        fusion = pandas.read_csv(curated_data_path + "sanger_fusion.csv", index_col=0)
    elif project == "BROAD":
        fusion = pandas.read_csv(curated_data_path + "broad_19Q2_fusion.csv", index_col=0)
    elif project == "INTEGRATED":
        fusion = pandas.read_csv(curated_data_path + "integrated_fusion.csv", index_col=0)
    return fusion


########################################################################################################################
########################################################################################################################
# CNV INFO #

# Read the curated data

def cnv(project):
    if project == "SANGER":
        cnv = pandas.read_csv(curated_data_path + "sanger_cnv.csv", index_col=0)
    elif project == "BROAD":
        cnv = pandas.read_csv(curated_data_path + "broad_19Q2_cnv.csv", index_col=0)
    elif project == "INTEGRATED":
        cnv = pandas.read_csv(curated_data_path + "integrated_cnv.csv", index_col=0)
    return cnv


########################################################################################################################
########################################################################################################################
# DRUG INFO #

# Read the curated data

def drug(project):
    if project == "SANGER":
        drug = pandas.read_csv(curated_data_path + "sanger_drug.csv", index_col=0)
    elif project == "BROAD":
        drug = pandas.read_csv(curated_data_path + "broad_19Q2_drug.csv", index_col=0)
    elif project == "INTEGRATED":
        drug = pandas.read_csv(curated_data_path + "integrated_drug.csv", index_col=0)
    return drug

########################################################################################################################
########################################################################################################################

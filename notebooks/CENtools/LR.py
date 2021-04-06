"""
                                    ---------------------
                                     C E N  -  T O O L S

                                     Logistic Regression
                                    ---------------------
"""
########################################################################################################################
########################################################################################################################
# IMPORTING #

from paths import path, curated_data_path, object_path
from data_preparation.objects import Gene, Sanger, Broad, Integrated
from data_preparation.objects import deserialisation_project, deserialisation_gene
from data_preparation.curation import essentiality
from data_preparation.curation import BEG, BNEG, adam_core, adam_context, stem_core, known_core, only_stem

import os, pandas, networkx, pickle, matplotlib, numpy, sklearn, argparse
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.metrics import classification_report, confusion_matrix, \
    roc_curve, average_precision_score, precision_recall_curve, plot_precision_recall_curve
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import jsonpickle.ext.numpy as jsonpickle_numpy
jsonpickle_numpy.register_handlers()


numpy.random.seed(1234)

########################################################################################################################
########################################################################################################################
# FIGURE PROPERTIES #

matplotlib.rc('font', size=10)
medianprops = {"color": "red", "linewidth": 1.5}
boxprops = {"color": "black", "linestyle": "-", "linewidth": 1.5}
whiskerprops = {"color": "black", "linestyle": "--", "linewidth": 1.5}
flierprops = {"color": "black", "marker": "o", "markersize": 8}

########################################################################################################################
########################################################################################################################
# PATHS #

if "prediction_output" not in os.listdir(path): os.mkdir(path + "prediction_output")
prediction_path = path + "prediction_output/"

########################################################################################################################
########################################################################################################################
# CLASSES #

# CEN-tools Predictions

class CEN_Prediction(Gene):

    def __init__(self, gene):
        self.gene = gene
        self.LR_probabilities = dict()
        self.bin_vectors = dict()
        self.si_cluster = dict()
        self.classification = None

    def add_predicted_probabilities(self, mega_probability_dictionary, used_probability, init_gene_obj):

        if init_gene_obj[self.gene].BNEG != True and init_gene_obj[self.gene].BEG != True:
            for cell_line, probability_df in mega_probability_dictionary.items():
                self.LR_probabilities[cell_line] = probability_df.loc[self.gene, used_probability]

    def add_bin_vectors(self, bin_number):

        if self.LR_probabilities != {}:
            df = pandas.DataFrame([self.LR_probabilities.values()], index=["Probability"]).T
            histogram = numpy.histogram(df, bins=bin_number)
            self.bin_vectors[bin_number] = {i: histogram[0][i] for i in range(len(histogram[0]))}

    def add_classification(self, cluster_df):

        if pandas.notna(cluster_df.loc[self.gene, "Cluster_name"]):
            self.classification = cluster_df.loc[self.gene, "Cluster_name"]


    def add_cluster_info(self, bin_number, si_df, cluster_df):

        if pandas.notna(cluster_df.loc[self.gene, "Cluster_name"]) and pandas.notna(si_df.loc[self.gene, "sil_width"]):
            self.si_cluster[bin_number] = {cluster_df.loc[self.gene, "Cluster_name"]:
                                               si_df.loc[self.gene, "sil_width"]}


########################################################################################################################
########################################################################################################################
#  FUNCTIONS #

def construct_lr_model(df, project, gene_object, result_path):
    """
    Construct the Logistic Regression on Gold Standard Essential and Non-Essential Genes.
    :param df: The essentiality tables from the specified project.
    :param project: The project is from where the essentiality values are coming.
    :param gene_obj: The dictionary of gene annotation objects.
    :return: The best model (median AUC) and Probability table for each gene being essential or non-essential.
    """
    # TAKE ESSENTIALITY OF THE GOLD STANDARD GENES

    gold_genes = [gene for gene, obj in gene_object.items() if (obj.BEG) or (obj.BNEG)]

    gold_df = df.loc[gold_genes]

    # PERTURB THE DATA FRAME TO AVOID ANY BIAS

    gold_df = gold_df.reindex(numpy.random.permutation(gold_df.index))

    # LABEL THE GOLD DATA FRAME

    gold_df["Category"] = gold_df.apply(lambda x: "EG" if gene_object[x.name].BEG else (
        "NEG" if gene_object[x.name].BNEG else "NONE"), axis=1)

    # RESET THE INDEX

    # By using indeces we will reach the gene names in further steps.
    gold_df = gold_df.reset_index()

    # BINARISE THE DATA FRAME FOR LR

    # 0 for Non-Essential and 1 for Essential Genes (Success in lR)
    gold_df["Binarised"] = label_binarize(gold_df.loc[:, ["Category"]].values, ["NEG", "EG"])

    # MELT THE DATA FRAME

    gold_df = gold_df.drop(columns=["Category"])
    melt_df = gold_df.melt(id_vars=["index", "Binarised"])[["index", "Binarised", "value"]]
    melt_df.columns = ["Gene", "Binarised", "Essentiality_Score"]

    # GENE ESSENTIALITY DISTRIBUTION
    '''
    fig, ax = plt.subplots()
    ax.boxplot(x=[numpy.array(melt_df[melt_df.Binarised == 0].Essentiality_Score),
                  numpy.array(melt_df[melt_df.Binarised == 1].Essentiality_Score)],
               positions=[1, 2], showfliers=False, medianprops=medianprops, boxprops=boxprops,
               whiskerprops=whiskerprops, flierprops=flierprops)

    plt.xticks([1, 2], ["BAGEL Non Essential", "BAGEL Essential"])
    plt.xlabel("BAGEL GENES")
    plt.ylabel("ESSENTIALITY SCORES")
    plt.show()
    plt.savefig(result_path + project + "_BAGEL_ESS_DIST.pdf")
    plt.close()
    '''
    # PERTURB THE DATA FRAME TO AVOID ANY BIAS

    melt_df = melt_df.reindex(numpy.random.permutation(melt_df.index))

    # EXTRACT FEATURE AND VALUE

    X = melt_df[["Essentiality_Score"]].values
    y = melt_df["Binarised"].values

    # PREPARE CROSS VALIDATION

    cross_val = KFold(n_splits=5)

    seed = 0

    cols, rows = 3, 2
    ROC_lines, PR_lines, ROC_labels, PR_labels = [], [], [], []
    info = {}
    ROC_fig, ROC_axes = plt.subplots(rows, cols, constrained_layout=True)
    PR_fig, PR_axes = plt.subplots(rows, cols, constrained_layout=True)
    colors = ["lightslategray", "darkseagreen"]

    for train_index, validation_index in cross_val.split(X):

        # SEPARATE TRAIN AND VALIDATION SETS BOTH FROM FEATURE AND VALUE ARRAYS

        X_train, X_validation = X[train_index], X[validation_index]
        y_train, y_validation = y[train_index], y[validation_index]

        # LOGISTIC REGRESSION

        classifier = LogisticRegression(solver="lbfgs", random_state=seed)

        # FIT THE MODEL ON TO TRAIN SETS

        model = classifier.fit(X_train, y_train)
        classes = list(model.classes_)

        # FIGURE CORRECTION

        axes_ind = [0, seed] if seed < 3 else [1, seed - 3]

        # PREDICTION WITH THE TRAINED MODEL ON VALIDATION SET

        y_prediction, y_probability = model.predict(X_validation), model.predict_proba(X_validation)

        # CREATE PROBABILITY DF

        # Take the gene name index to create the probability data frame with gene annotations
        genes = melt_df.loc[validation_index]["Gene"]
        prob_DF = pandas.DataFrame(y_probability)
        prob_DF.index = genes

        # DRAW ROC CURVES ONE BY ONE

        fpr, tpr = {}, {}
        for element_ind in range(len(classes)):
            fpr[element_ind], tpr[element_ind], _ = sklearn.metrics.roc_curve(
                y_true=y_validation, y_score=y_probability[::, element_ind], pos_label=element_ind)

            #  ROC CURVE

            ROC_l, = ROC_axes[axes_ind[0], axes_ind[1]].plot(fpr[element_ind],
                                                             tpr[element_ind], color=colors[element_ind])
            ROC_lines.append(ROC_l)
            ROC_labels.append(classes[element_ind])

        info[seed + 1] = {"AUC": sklearn.metrics.auc(fpr[1], tpr[1]), "MODEL": model,
                          "MATRIX": confusion_matrix(y_validation, y_prediction),
                          "REPORT": classification_report(y_validation, y_prediction),
                          "AP_E": average_precision_score(y_validation, y_probability[:, 1]),
                          "AP_NE": average_precision_score(y_validation, y_probability[:, 0], pos_label=0),
                          "Probability_DF": prob_DF}

        ROC_axes[axes_ind[0], axes_ind[1]].plot([0, 1], [0, 1], color="silver", linestyle='dotted')
        ROC_axes[axes_ind[0], axes_ind[1]].set_xlabel("False Positive Rate")
        ROC_axes[axes_ind[0], axes_ind[1]].set_ylabel("True Positive Rate")
        ROC_axes[axes_ind[0], axes_ind[1]].set_title("KFold : %s - AUC : %.2f"
                                                     % (str(seed + 1), info[seed + 1]["AUC"]))

        precision, recall, _ = precision_recall_curve(y_validation, y_probability[:, 1].ravel())

        PR_axes[axes_ind[0], axes_ind[1]].set_xlabel("Recall")
        PR_axes[axes_ind[0], axes_ind[1]].set_ylabel("Precision")
        PR_axes[axes_ind[0], axes_ind[1]].set_ylim([0.1, 1.05])
        PR_l, = PR_axes[axes_ind[0], axes_ind[1]].plot(recall, precision, color=colors[1], lw=1)
        PR_lines.append(PR_l)
        PR_labels.append("EG")
        PR_axes[axes_ind[0], axes_ind[1]].set_title("KFold : %s - AP : %.2f" % (str(seed + 1),
                                                                                info[seed + 1]["AP_E"]))

        seed += 1

    # DELETE THE EMPTY AXES FROM THE FIGURES

    ROC_axes[1, 2].remove()
    PR_axes[1, 2].remove()

    eg_lines = [Line2D([0], [0], color=colors[0], lw=2), Line2D([0], [0], color=colors[1], lw=2)]

    # TAKE THE MEDIAN OF THE AUCS

    ROC_medians = [value["AUC"] for key, value in info.items()]
    '''
    ROC_fig.suptitle("SANGER - ROC Curve - Median AUC = %.3f" % (numpy.median(list(ROC_medians))))
    ROC_fig.legend(eg_lines, classes, loc="lower right")
    plt.show()
    ROC_fig.savefig(result_path + project + "_ROC_AUC_LR.pdf")
    plt.close()
    '''
    # PR CURVE

    PR_medians = [value["AP_E"] for key, value in info.items()]
    '''
    PR_fig.suptitle("SANGER - PR Curve - Median AP = %.3f" % (numpy.median(list(PR_medians))))
    PR_fig.legend(eg_lines, list(classes), loc="lower right")
    plt.show()
    PR_fig.savefig(result_path + project + "_PR_LR.pdf")
    plt.close()
    '''
    # BEST MODEL -- Median

    ROC_med = [key for key, value in info.items() if value["AUC"] == numpy.median(ROC_medians)]
    PR_med = [key for key, value in info.items() if value["AP_E"] == numpy.median(PR_medians)]

    if ROC_med == PR_med:
        print("ROC and Precision-Recall give the same median KFold.")

    med = PR_med[0]
    best = info[med]
    '''
    f = open(result_path + project + "_LR_Model_Scores.txt", "w")
    f.write(best["REPORT"] + "\n\n")
    f.write("Average Precision of the Median Model\n")
    f.write("EG:\t%.3f\n" % (info[med]["AP_E"]))
    f.write("NEG:\t%.3f\n" % (info[med]["AP_NE"]))

    f.close()
    '''
    prob_df = info[med]["Probability_DF"]

    return best["MODEL"], prob_df


def apply_models(model, df, gene_object):
    """
    Apply Pre-defined LR model on genes which have unknown essentiality annotation.
    :param model: The pre-defined model created with "construct_lr_model" function.
    :param df: The essentiality tables from the specified project.
    :param gene_object: The dictionary of gene annotation objects.
    :return: A dictionary having keys as cell lines and values as corresponding probability tables.
    """

    # ELIMINATE THE GOLD STANDARD GENES

    gold_genes = [gene for gene, obj in gene_object.items() if (obj.BEG) or (obj.BNEG)]
    df = df.drop(gold_genes)

    # PREDICT AND TAKE THE PROBABILITIES OF EACH GENE BEING ESSENTIAL OR NOT IN EACH CELL LINE

    results = {}
    for cell_line in df.columns:

        # Take the essentiality values of all genes in one cell line.
        inside_cl = df[[cell_line]].values

        # Take the probabilities of all cell in one cell line with the pre-defined model.
        inside_cl_probabilities = model.predict_proba(inside_cl)

        # Create Data Frame composed of probabilities.
        inside_cl_probabilities_df = pandas.DataFrame(inside_cl_probabilities,
                                                      index=df.index,columns=list(model.classes_))

        inside_cl_probabilities_df = inside_cl_probabilities_df.rename(columns={0: "NEG_LR", 1: "EG_LR"})
        results[cell_line] = inside_cl_probabilities_df

    return results


def main_lr(**args):

    # GENE OBJECT

    initial_gene_obj = deserialisation_gene(args["GENE_PROJECT"])

    # RESULT PATH

    if args["PROJECT"] not in os.listdir(prediction_path):
        os.mkdir(prediction_path + args["PROJECT"])

    result_path = prediction_path + args["PROJECT"] + "/"

    if args["APPLICATION"]:

        # LR will apply on to Corrected FCs

        essentiality_df = essentiality(args["PROJECT"], "FC")

        # CONSTRUCTION MODEL

        model, probability_df = construct_lr_model(
            df=essentiality_df, project=args["PROJECT"],
            gene_object=initial_gene_obj, result_path=result_path)

        # SAVE THE MODEL

        pickle.dump(model, open(result_path + args["PROJECT"] + "_MODEL.sav", "wb"))

        # APPLY MODEL TO PREDICT

        probabilities = apply_models(model=model, df=essentiality_df,
                                     gene_object=initial_gene_obj)

        # SAVE PROBABILITY TABLE AFTER PREDICTION

        pickle.dump(probabilities, open(result_path + args["PROJECT"] + "_Prob_DF.p", "wb"))


    else:

        probabilities = pickle.load(open(result_path + args["PROJECT"] + "_Prob_DF.p", "rb"))


    # After prediction probabilities or already predicted probabilities were read!

    # Add probability information from LR on to object

    cen_gene_obj = {}
    for gene in list(probabilities[list(probabilities.keys())[0]].index):
        obj = CEN_Prediction(gene)
        obj.add_predicted_probabilities(mega_probability_dictionary=probabilities,
                                        used_probability="EG_LR", init_gene_obj=initial_gene_obj)
        cen_gene_obj[gene] = obj

    pickle.dump(cen_gene_obj, open(object_path + args["PROJECT"] + "_PROB_CEN_OBJ.pkl", "wb"))

    prob_cen_obj = pickle.load(open(object_path + args["PROJECT"] + "_PROB_CEN_OBJ.pkl", "rb"))

    bin_vector_dicts = {}
    for gene in prob_cen_obj.keys():
        obj = prob_cen_obj[gene]
        obj.add_bin_vectors(bin_number=args["BIN_NUMBER"])
        prob_cen_obj[gene] = obj
        if gene in bin_vector_dicts.keys():bin_vector_dicts[gene][args["BIN_NUMBER"]] = obj.bin_vectors
        else: bin_vector_dicts[gene] = {bin_number: obj.bin_vectors}

    pickle.dump(prob_cen_obj, open(object_path + args["PROJECT"] + "_BIN_CEN_OBJ.pkl", "wb"))

    geneKeys = bin_vector_dicts.keys()
    bin_df = numpy.array([bin_vector_dicts[i][bin_number][bin_number][j] for i in geneKeys for j in range(0, bin_number)])
    bin_df = numpy.reshape(bin_df, (len(geneKeys), 20))
    bin_df = pandas.DataFrame(bin_df, index=geneKeys)
    bin_df.to_csv(result_path + args["PROJECT"] + "_Histogram_DF_" + str(args["BIN_NUMBER"]) +
                  "_BIN.txt", sep="\t")

    return probabilities, bin_df


def run(project, gene_project, bin_number, model_application = False):
    parsed_input = {"GENE_PROJECT": gene_project, "PROJECT": project,
                    "APPLICATION": model_application, "BIN_NUMBER" : bin_number}
    main_lr(**parsed_input)


########################################################################################################################
########################################################################################################################
# DE-SERIALISATION OF THE OBJECTS #

def deserialisation_cen_gene(project):
    global object_path
    cen_gene_obj = pickle.load(open(object_path + project + "_BIN_CEN_OBJ.pkl", "rb"))
    return cen_gene_obj

########################################################################################################################
########################################################################################################################

# Run from command line
gene_project = 'CUSTOM'
project = 'INTEGRATED'
bin_number = 20
model_application = True

essentiality_df = essentiality(project, "FC")

GenesDict = {}

for gene in essentiality_df.index:
    obj = Gene(gene)
    obj.add_bagel_info(essential_genes = BEG, non_essential_genes = BNEG)

    GenesDict[gene] = obj

pickle.dump(GenesDict, open(object_path + "INITIAL_CUSTOM_GENE_OBJ.pkl", "wb"))

run(project = "INTEGRATED", gene_project = "CUSTOM", bin_number = 20, model_application = True)

for f in os.listdir(object_path):
    os.remove(os.path.join(object_path, f))

print("pickle dictionary and output files have been deleted from " + object_path + " subfolder")

os.remove(curated_data_path + "CERES_scaled_depFC.csv")
print("CERES scaled dataset has been deleted from " + curated_data_path + " subfolder")

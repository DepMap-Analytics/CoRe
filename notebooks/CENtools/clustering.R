"
                                    ---------------------
                                     C E N  -  T O O L S

                                         Clustering
                                    ---------------------
                                 Created by: Sumana Sharma, PhD
"

# Install Packages
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(factoextra)
library(GGally)
library(plotly)
library(cluster)
library(dplyr)
library(enrichR)
library(CRISPRcleanR)

# User defined info

# core_marker <- ""
# context_marker <- ""
# rare_context_marker <- ""
# non_marker <- ""
# chosen_project <- ""
# path <- ""
# project_path <- paste0(path, "prediction_output/", chosen_project)


# Function for normalisation
normalize <- function(x) {
 return ((x - min(x)) / (max(x) - min(x)))
}


# To do the same with just CRISPR score
ClusterEssentiality<- function(Chosen_project, binPath, resultPath){
    set.seed(1234)

    Clusters_Project<- read.table(binPath, header=T, stringsAsFactors = F)
    Clusters_Project_matrix<- normalize(as.matrix(Clusters_Project))

    clustering_project <- kmeans(Clusters_Project_matrix,centers = 4, nstart = 20,  algorithm = "Hartigan-Wong")
    Combined.pca <- prcomp(t(Clusters_Project_matrix))

    Cluster_information<- as.data.frame(clustering_project$cluster)
    Data_to_plot<- data.frame(Combined.pca$rotation)
    Data_to_plot<- Data_to_plot[,(1:3)]
    Data_to_plot$Gene<- rownames(Data_to_plot)
    Data_to_plot$cluster<- Cluster_information$`clustering_project$cluster`[match(Data_to_plot$Gene, rownames(Cluster_information))]
    Data_to_plot$cluster<- as.factor(as.character(Data_to_plot$cluster))

    pdf(paste0(resultPath, Chosen_project, "_PCA_clusters.pdf"), width= 6, height =4)
    p <- ggplot(Data_to_plot, aes(PC1,PC2, color= cluster))
    p1 <- p + geom_point()+
      xlab("PC1")+
      ylab ('PC2')+
      theme_bw() +
      theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14 ),
            axis.text.x = element_text(angle = 0, size=14, hjust=0.5,vjust=0.5),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            plot.background = element_blank() ,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_color_manual(values=c('#003f5c', '#58508d', '#bc5090', '#ff6361', "#ffa600"))+
    theme(legend.title=element_blank())+
      theme(plot.title = element_text(size = 14))+
      theme(legend.position= 'right')
    print(p1)
    dev.off()

    sil <- silhouette(clustering_project$cluster, dist(Clusters_Project_matrix))
    sil_df_project<-as.data.frame(cbind(cluster=sil[,1],neighbor=sil[,2],sil_width=sil[,3]))
    sil_df_project$Gene<- rownames(Clusters_Project)
    write.csv(sil_df_project, file= paste0(resultPath, Chosen_project, "_clusters_sil_scores.csv"))
    out_sil <- fviz_silhouette(sil)

    Clusters_Project$cluster <- as.factor(clustering_project$cluster)

    # Marker_core<- core_marker
    # Marker_context_1<- context_marker
    # Marker_context_2<- rare_context_marker
    # Marker_non<- non_marker

    # Essential<- as.data.frame(as.numeric(Clusters_Project[rownames(Clusters_Project)==Marker_core, ncol(Clusters_Project)]), 'Essential_genes')
    # colnames(Essential)<- c('Cluster_name')
    # Context_1<- as.data.frame(as.numeric(Clusters_Project[rownames(Clusters_Project)==Marker_context_1, ncol(Clusters_Project)]), 'Context1')
    # colnames(Context_1)<- c('Cluster_name')
    # Context_2<- as.data.frame(as.numeric(Clusters_Project[rownames(Clusters_Project)==Marker_context_2, ncol(Clusters_Project)]), 'Context2')
    # colnames(Context_2)<- c('Cluster_name')
    # Non_essential<- as.data.frame(as.numeric(Clusters_Project[rownames(Clusters_Project)==Marker_non, ncol(Clusters_Project)]), 'Non_essential')
    # colnames(Non_essential)<- c('Cluster_name')
    #
    # Cluster_mapping<- rbind(Essential,Context_1, Context_2, Non_essential)
    #
    # Clusters_Project$Cluster_name<- rownames(Cluster_mapping)[match(Clusters_Project$cluster, Cluster_mapping$Cluster_name)]
    #
    # Clusters_project_essential<- Clusters_Project[Clusters_Project$Cluster_name=='Essential_genes',]
    # Essentials<- as.character(rownames(Clusters_project_essential))

    summary_sil <- out_sil$data %>% group_by(cluster) %>% summarise(avg = mean(sil_width))

    Essentials <- rownames(Clusters_Project)[which(Clusters_Project$cluster == which.max(summary_sil$avg))]

    return(Essentials)

    # write.csv(Clusters_Project, file= paste0(resultPath, Chosen_project, "_Clusters.csv"))
    # write.csv(Essentials, file= paste0(resultPath, Chosen_project, "_Essential_cluster.csv"))
}


# ClusterEssentiality(Chosen_project= 'INTEGRATED',
#                     binPath= paste0(project_path, "INTEGRATED_Histogram_DF_20_BIN.txt"),
#                     resultPath = project_path)

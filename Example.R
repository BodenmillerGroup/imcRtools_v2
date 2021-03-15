List_data = Load_IMC_data("Desktop/analysis/cpout/cell.csv",
                         Path_to_panel_file ="Desktop/analysis/cpout/panel.csv",
                         DNA_channel = c("Iridium191","Iridium193"),
                         Use_for_clustering = c("Myeloperoxidase.MPO","CD31","SMA","Cytokeratin.5",
                                                "Keratin.14","Vimentin","CD3","CD68", "Cytokeratin.8.18","CD45",
                                                "Carbonic.Anhydrase.IX","Fibronectin","Ki.67","CD20","CD44","pan.Cytokeratin"))

sce = Create_SCE(List_data)
sce = Poisson_normalization(sce)
sce = UMAP_embedding(sce)
sce = KNN_clustering(sce)
sce = Compute_neighborhood_graph(sce,graph_type = "KNN",graph_parameter = 5)
sce = Load_object_relationship(sce,Path_to_object_relationship = "Desktop/analysis/cpout/Object relationships.csv")


plot_UMAP(sce)
Plot_cluster_spatial(sce)
Plot_gene_expression_spatial(sce,Gene = "Iridium191")


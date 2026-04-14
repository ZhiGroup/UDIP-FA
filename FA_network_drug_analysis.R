required_packages <- c(
  "dplyr",
  "UpSetR",
  "data.table",
  "VennDiagram",
  "ggpubr",
  "tidyr",
  "ggplot2",
  "ggalluvial",
  "readr",
  "pheatmap",
  "stringr"
)

load_required_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

read_table_checked <- function(path, ..., data.table = FALSE) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path, call. = FALSE)
  }
  fread(path, ..., data.table = data.table)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

default_network_analysis_config <- function() {
  list(
    trn_network = "Three_network/TRN.csv",
    coexpression_network = "Three_network/gene_cor_expression.csv",
    ppi_network = "Three_network/PPI_network.csv",
    disorder_table = "Three_network/disease_gene.csv",
    genes_file = "../FA_5e-8/FUMA_job604167/genes.txt",
    lead_snp_file = "../FA_5e-8/FUMA_job604167/leadSNPs.txt",
    magma_file = "../FA_5e-8/FUMA_job604167/magma.genes.out",
    ppi_dir = "./PPI",
    coexpression_dir = "./Co_expression",
    trn_dir = "./TRN",
    drug_interactions = "drug/interactions.tsv",
    drug_classes = "drug_classes_pivoted.csv",
    selected_diseases = c(
      "Schizophrenia",
      "Parkinson Disease",
      "Attention Deficit and Disruptive Behavior Disorders",
      "Multiple Sclerosis",
      "Alzheimer Disease",
      "Depression",
      "Bipolar Disorder",
      "Autism Spectrum Disorder"
    )
  )
}

mappping_disorder_gene_PPI <- function(network_PPI, FA_gene, disorder_table, disorder_name, output_dir = ".") {
  disorder_gene <- disorder_table[, disorder_name]
  disorder_gene <- disorder_gene[disorder_gene != ""]
  rownames(network_PPI) <- seq_len(nrow(network_PPI))

  inter_gene <- intersect(FA_gene, disorder_gene)
  disorder_label <- apply(network_PPI, 1, function(x) {
    as.integer((x[1] %in% FA_gene & x[2] %in% disorder_gene) | (x[2] %in% FA_gene & x[1] %in% disorder_gene))
  })

  disorder_subnetwork <- cbind(network_PPI, disorder_label)
  disorder_subnetwork <- disorder_subnetwork[disorder_subnetwork$disorder_label == 1, ]

  write.table(
    disorder_subnetwork[, c(1, 2, 3)],
    file.path(output_dir, paste0(disorder_name, "net.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  all_gene_name <- unique(union(disorder_subnetwork[, 1], disorder_subnetwork[, 2]))
  FA_gene_ld <- setdiff(intersect(all_gene_name, FA_gene), inter_gene)
  disorder_ge_l <- setdiff(intersect(all_gene_name, disorder_gene), inter_gene)

  node_ano <- data.frame(
    gene = c(FA_gene_ld, disorder_ge_l, inter_gene),
    label_data = c(rep(1, length(FA_gene_ld)), rep(2, length(disorder_ge_l)), rep(3, length(inter_gene)))
  )
  colnames(node_ano) <- c("gene", paste0(disorder_name, "_anno"))

  write.table(
    node_ano,
    file.path(output_dir, paste0(disorder_name, "node.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  disorder_subnetwork
}

gene_data_tf_f <- function(disorder_net, disorder_gene, FA_gene, symbol_d) {
  colnames(disorder_net) <- c("node1", "node2", "weight")
  save_edge_d <- data.frame()

  if (symbol_d == "TRN") {
    inter_tfs <- intersect(disorder_net$node1, FA_gene)
    for (i in inter_tfs) {
      te_gene <- disorder_net[disorder_net$node1 == i, ]
      save_edge_d <- rbind(save_edge_d, te_gene)
    }
    return(unique(save_edge_d$node2))
  }

  bag_re_net <- disorder_net[disorder_net$node1 %in% FA_gene | disorder_net$node2 %in% FA_gene, ]
  union(bag_re_net$node1, bag_re_net$node2)
}

bag_upset <- function(FA_gene, root_path, symbol_d) {
  node_files <- list(
    AD = "ADnode.txt",
    PD = "PDnode.txt",
    ASD = "ASDnode.txt",
    BIP = "BIPnode.txt",
    SCZ = "SCZnode.txt",
    DEP = "DEPnode.txt",
    ADHD = "ADHDnode.txt",
    MS = "MSnode.txt"
  )
  edge_files <- list(
    AD = "ADnet.txt",
    PD = "PDnet.txt",
    ASD = "ASDnet.txt",
    BIP = "BIPnet.txt",
    SCZ = "SCZnet.txt",
    DEP = "DEPnet.txt",
    ADHD = "ADHDnet.txt",
    MS = "MSnet.txt"
  )

  read_named_table <- function(filename) {
    read.table(file.path(root_path, filename), header = TRUE)
  }

  AD_node <- read_named_table(node_files$AD)
  PD_node <- read_named_table(node_files$PD)
  ASD_node <- read_named_table(node_files$ASD)
  BIP_node <- read_named_table(node_files$BIP)
  SCZ_node <- read_named_table(node_files$SCZ)
  DEP_node <- read_named_table(node_files$DEP)
  ADHD_node <- read_named_table(node_files$ADHD)
  MS_node <- read_named_table(node_files$MS)

  AD_gene <- AD_node[AD_node$Alzheimer.s.Disease_anno %in% c(2, 3), ]$gene
  PD_gene <- PD_node[PD_node$Parkinson.s.Disease_anno == 2, ]$gene
  MS_gene <- MS_node[MS_node$MS_anno == 2, ]$gene
  ASD_gene <- ASD_node[ASD_node$Autism.Spectrum.Disorder_anno %in% c(2, 3), ]$gene
  BIP_gene <- BIP_node[BIP_node$Bipolar.Disorder_anno %in% c(2, 3), ]$gene
  SCZ_gene <- SCZ_node[SCZ_node$Schizophrenia_anno %in% c(2, 3), ]$gene
  DEP_gene <- DEP_node[DEP_node$Depression_anno %in% c(2, 3), ]$gene
  ADHD_gene <- ADHD_node[ADHD_node$Attention.deficit.hyperactivity.disorder_anno %in% c(2, 3), ]$gene

  AD_edge <- read_named_table(edge_files$AD)
  PD_edge <- read_named_table(edge_files$PD)
  ASD_edge <- read_named_table(edge_files$ASD)
  BIP_edge <- read_named_table(edge_files$BIP)
  SCZ_edge <- read_named_table(edge_files$SCZ)
  DEP_edge <- read_named_table(edge_files$DEP)
  ADHD_edge <- read_named_table(edge_files$ADHD)
  MS_edge <- read_named_table(edge_files$MS)

  upset_list <- list(
    AD = gene_data_tf_f(AD_edge, AD_gene, FA_gene, symbol_d),
    PD = gene_data_tf_f(PD_edge, PD_gene, FA_gene, symbol_d),
    ASD = gene_data_tf_f(ASD_edge, ASD_gene, FA_gene, symbol_d),
    BIP = gene_data_tf_f(BIP_edge, BIP_gene, FA_gene, symbol_d),
    SCZ = gene_data_tf_f(SCZ_edge, SCZ_gene, FA_gene, symbol_d),
    DEP = gene_data_tf_f(DEP_edge, DEP_gene, FA_gene, symbol_d),
    ADHD = gene_data_tf_f(ADHD_edge, ADHD_gene, FA_gene, symbol_d),
    MS = gene_data_tf_f(MS_edge, MS_gene, FA_gene, symbol_d)
  )

  whole_net <- rbind(AD_edge, PD_edge, ASD_edge, BIP_edge, SCZ_edge, DEP_edge, ADHD_edge, MS_edge)
  colnames(whole_net) <- c("node1", "node2", "weight")

  set_order <- c("MS", "ADHD", "DEP", "BIP", "PD", "AD", "SCZ", "ASD")
  set_colors <- c(
    "MS" = "#fdb462",
    "ADHD" = "#ffffb3",
    "DEP" = "#80b1d3",
    "BIP" = "#fb8072",
    "PD" = "#fccde5",
    "AD" = "#8dd3c6",
    "SCZ" = "#b3de69",
    "ASD" = "#bebada"
  )

  p <- upset(
    fromList(upset_list),
    sets = set_order,
    sets.bar.color = set_colors[set_order],
    nsets = 100,
    nintersects = 40,
    order.by = "freq",
    keep.order = TRUE,
    mb.ratio = c(0.6, 0.4),
    text.scale = 2
  )
  print(p)

  list(up_set = upset_list, whole_net = whole_net)
}

build_brain_drug_table <- function(FA_gene, drug_interactions_path, drug_classes_path, output_path = "brain_net_gene.txt") {
  drug_gene_all <- read_table_checked(drug_interactions_path, data.table = FALSE)
  FA_drug_gene <- drug_gene_all[drug_gene_all$gene_name %in% FA_gene, ]
  FA_drug_gene <- FA_drug_gene[
    FA_drug_gene$approved == "TRUE" &
      FA_drug_gene$immunotherapy == "FALSE" &
      FA_drug_gene$anti_neoplastic == "FALSE",
  ]
  FA_drug_gene$drug_concept_id <- vapply(strsplit(FA_drug_gene$drug_concept_id, ":"), `[`, character(1), 2)

  write.table(
    FA_drug_gene$drug_concept_id,
    "FA_interaction_drug.txt",
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE
  )

  drug_class_all <- read_table_checked(drug_classes_path, data.table = FALSE)
  brain_bone_df <- drug_class_all %>%
    filter(str_detect(
      tolower(DISEASE),
      "brain|alzheimer|parkinson|dementia|schizophrenia|epilepsy|cognitive|neurodegenerative|cerebral|stroke|multiple sclerosis|autism|bipolar|depression|anxiety|osteoporosis|bone density|fracture|osteopenia"
    ))
  write.table(brain_bone_df, "brain_drug.txt", row.names = FALSE, quote = FALSE)

  brain_bone_df$RxCUI <- as.character(brain_bone_df$RxCUI)
  result <- left_join(FA_drug_gene, brain_bone_df, by = c("drug_concept_id" = "RxCUI"))
  brain_drug_fa_sta <- result[result$drug_concept_id %in% brain_bone_df$RxCUI, ]
  brain_drug_fa_sta <- brain_drug_fa_sta[brain_drug_fa_sta$interaction_type != "NULL", ]

  df_expanded <- brain_drug_fa_sta %>%
    mutate(DISEASE = strsplit(as.character(DISEASE), ",")) %>%
    unnest(DISEASE) %>%
    mutate(DISEASE = trimws(DISEASE))

  brain_related_diseases <- c(
    "Schizophrenia", "Parkinson Disease", "Brain Ischemia", "Cerebral Hemorrhage", "Stroke",
    "Anxiety Disorders", "Seizures", "Unconsciousness", "Epilepsy", "Coma", "Psychotic Disorders",
    "Neuroleptic Malignant Syndrome", "Bipolar Disorder", "Autistic Disorder", "Alcohol Withdrawal Delirium",
    "Cerebral Palsy", "Huntington Disease", "Multiple Sclerosis", "Myelitis", "Spinal Cord Injuries",
    "Spinal Cord Ischemia", "Spinal Cord Neoplasms", "Cerebellar Ataxia", "Migraine Disorders",
    "Migraine with Aura", "Cerebral Infarction", "Alzheimer Disease", "Status Epilepticus",
    "Epilepsies, Partial", "Neuralgia, Postherpetic", "Phobic Disorders", "Depression, Postpartum",
    "Cluster Headache, Spontaneous", "Sleep Initiation and Maintenance Disorders", "Tonic-Clonic", "Tremor",
    "Dementia", "Restless Legs Syndrome", "Intracranial Hemorrhages", "Brain Neoplasms",
    "Lennox Gastaut Syndrome", "Spasms, Infantile", "Conduct Disorder", "Depressive Disorder",
    "Stress Disorders, Post-Traumatic", "Tourette Syndrome", "Cognitive Dysfunction", "Dysarthria",
    "Myoclonic", "Absence", "Panic Disorder", "Tic Disorders", "Catatonia", "Gait Disorders, Neurologic",
    "Myoclonic Epilepsy, Juvenile", "Headache", "Muscle Rigidity",
    "Attention Deficit and Disruptive Behavior Disorders", "Head Injuries, Closed",
    "Meningeal Neoplasms", "Obsessive-Compulsive Disorder", "Anorexia Nervosa",
    "Bulimia Nervosa", "Attention Deficit Disorder with Hyperactivity"
  )

  df_brain <- df_expanded %>% filter(DISEASE %in% brain_related_diseases)
  write.table(df_brain, output_path, row.names = FALSE, quote = FALSE, sep = "\t")
  df_brain
}

plot_brain_drug_alluvial <- function(df_brain, selected_diseases) {
  is_alluvia_form(df_brain, axes = 1:3, silent = TRUE)
  df_brain <- df_brain %>% filter(!is.na(interaction_score))
  df_brain_se <- df_brain[df_brain$DISEASE %in% selected_diseases, ]
  df_brain_se$freq <- 1

  ggplot(
    as.data.frame(df_brain_se),
    aes(axis1 = gene_name, axis2 = drug_name, axis3 = DISEASE, y = freq)
  ) +
    scale_x_discrete(limits = c("Gene", "Drug", "Disease"), expand = c(.1, .1)) +
    geom_alluvium(aes(fill = gene_name), width = 1 / 12, alpha = 0.8) +
    geom_stratum(width = 1 / 12, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    theme_minimal(base_size = 14) +
    theme(axis.title.x = element_blank()) +
    labs(title = "Gene-Drug-Brain Disease Alluvial Diagram", y = "Interaction Score")
}

run_fa_network_drug_analysis <- function(config = default_network_analysis_config()) {
  load_required_packages(required_packages)

  trn_network_data <- read_table_checked(config$trn_network, data.table = FALSE)
  co_expression_network <- read_table_checked(config$coexpression_network, data.table = FALSE)
  ppi_intersection_score_data <- read_table_checked(config$ppi_network, data.table = FALSE)
  genes_table <- read_table_checked(config$genes_file, data.table = FALSE)
  lead_snp <- read_table_checked(config$lead_snp_file, data.table = FALSE)
  magma <- read_table_checked(config$magma_file, data.table = FALSE)
  disorder_table <- read_table_checked(config$disorder_table, data.table = FALSE)

  ensure_dir(config$ppi_dir)
  ensure_dir(config$coexpression_dir)
  ensure_dir(config$trn_dir)

  FA_gene_pos <- genes_table[genes_table$posMapSNPs > 0, ]
  if ("symbol" %in% colnames(FA_gene_pos)) {
    write.table(FA_gene_pos$symbol, "FA_gene_anan.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if ("symbol" %in% colnames(FA_gene_pos) && "IndSigSNPs" %in% colnames(FA_gene_pos)) {
    SOX10 <- FA_gene_pos[FA_gene_pos$symbol == "SOX10", ]
    if (nrow(SOX10) > 0) {
      sox10_snps <- unlist(strsplit(SOX10$IndSigSNPs, ";"))
      message("SOX10 lead SNP overlap count: ", length(intersect(sox10_snps, lead_snp$rsID)))
    }
  }

  magma <- magma[magma$P < (0.05 / nrow(magma)), ]
  FA_gene <- union(unique(genes_table$symbol), magma$SYMBOL)
  disorder_names <- colnames(disorder_table)

  network_inputs <- list(
    TRN = list(data = trn_network_data, output_dir = config$trn_dir),
    Co_expression = list(data = co_expression_network, output_dir = config$coexpression_dir),
    PPI = list(data = ppi_intersection_score_data, output_dir = config$ppi_dir)
  )

  for (network_name in names(network_inputs)) {
    network_input <- network_inputs[[network_name]]
    for (disorder_name in disorder_names) {
      mappping_disorder_gene_PPI(
        network_PPI = network_input$data,
        FA_gene = FA_gene,
        disorder_table = disorder_table,
        disorder_name = disorder_name,
        output_dir = network_input$output_dir
      )
    }
  }

  FA_gene_PPI <- bag_upset(FA_gene, config$ppi_dir, "PPI")
  FA_gene_co_exp <- bag_upset(FA_gene, config$coexpression_dir, "EXP")
  FA_gene_TRN <- bag_upset(FA_gene, config$trn_dir, "TRN")

  df_brain <- build_brain_drug_table(
    FA_gene = FA_gene,
    drug_interactions_path = config$drug_interactions,
    drug_classes_path = config$drug_classes
  )

  write.csv(df_brain, "UFAG_gene.csv", row.names = FALSE, quote = FALSE)
  plot_brain_drug_alluvial(df_brain, config$selected_diseases)

  list(
    FA_gene = FA_gene,
    lead_snp = lead_snp,
    TRN = FA_gene_TRN,
    PPI = FA_gene_PPI,
    Co_expression = FA_gene_co_exp,
    brain_drug = df_brain
  )
}

if (sys.nframe() == 0) {
  run_fa_network_drug_analysis()
}

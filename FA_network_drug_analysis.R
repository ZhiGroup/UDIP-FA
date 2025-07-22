library(dplyr)
library(UpSetR)
library(data.table) # For UpSet plot (upset package, suitable for 2-7 sample sets)
library(VennDiagram)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(readr)
library(pheatmap)

## Network enrichment function ##
mappping_disorder_gene_PPI <- function(network_PPI, FA_gene, disorder_name) {
  # Extract disease gene set
  disorder_gene <- disorder_table[, disorder_name]
  rownames(network_PPI) <- c(1:nrow(network_PPI))
  disorder_gene <- disorder_gene[which(disorder_gene != "")]
  # Overlap genes
  inter_gene <- intersect(FA_gene, disorder_gene)
  # Identify disorder-gene and BAG-gene interactions
  disorder_label <- apply(network_PPI, 1, function(x) {
    if ((x[1] %in% FA_gene & x[2] %in% disorder_gene) | (x[2] %in% FA_gene & x[1] %in% disorder_gene)) {
      return(1)
    } else {
      return(0)
    }
  })
  disorder_subnetwork <- cbind(network_PPI, disorder_label)
  disorder_subnetwork <- disorder_subnetwork[which(disorder_subnetwork$disorder_label == 1), ]
  # Write subnetwork edge file
  write.table(disorder_subnetwork[, c(1, 2, 3)], paste(disorder_name, 'net.txt', sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)
  # Node annotation: 1 = BAG gene, 2 = disorder gene, 3 = shared gene
  all_gene_name <- union(disorder_subnetwork[, 1], disorder_subnetwork[, 2]) %>% unique()
  FA_gene_ld <- intersect(all_gene_name, FA_gene)
  disorder_ge_l <- intersect(all_gene_name, disorder_gene)
  FA_gene_ld <- setdiff(FA_gene_ld, inter_gene)
  disorder_ge_l <- setdiff(disorder_ge_l, inter_gene)
  label_data_gene <- c(rep(1, length(FA_gene_ld)), rep(2, length(disorder_ge_l)), rep(3, length(inter_gene)))
  gene_comb <- c(FA_gene_ld, disorder_ge_l, inter_gene)
  node_ano <- data.frame(gene_name = gene_comb, label_data = label_data_gene)
  colnames(node_ano) <- c('gene', paste(disorder_name, '_anno', sep = ""))
  write.table(node_ano, paste(disorder_name, 'node.txt', sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)
  return(disorder_subnetwork)
}

## Data loading ##
TRN_network_data <- fread('Three_network/TRN.csv', data.table = FALSE)
co_expression_network <- fread('Three_network/gene_cor_expression.csv', data.table = FALSE)
PPi_intersection_score_data <- fread('Three_network/PPI_network.csv', data.table = FALSE)
FA_gene <- fread('../FA_5e-8/FUMA_job604167/genes.txt', data.table = FALSE)
FA_gene_pos <- FA_gene[which(FA_gene$posMapSNPs > 0), ]
lead_SNP <- fread('../FA_5e-8/FUMA_job604167/leadSNPs.txt', data.table = FALSE)
SOX10 <- FA_gene_pos[which(FA_gene_pos$symbol == 'SOX10'), ]
str_split(SOX10$IndSigSNPs, ';') %>% unlist()
intersect(str_split(SOX10$IndSigSNPs, ';') %>% unlist(), lead_SNP$rsID)
write.table(FA_gene_pos$symbol, 'FA_gene_anan.txt', sep = '\t', row.names = FALSE, quote = FALSE)

TRN_FA <- FA_gene_TRN$whole_net
FA_lead_SNP <- fread('../FA_5e-8/FUMA_job604167/leadSNPs.txt', data.table = FALSE)

FA_magma <- fread('../FA_5e-8/FUMA_job604167/magma.genes.out')
FA_magma <- FA_magma[which(FA_magma$P < (0.05 / nrow(FA_magma))), ]
FA_gene <- FA_gene$symbol %>% unique()
FA_gene <- union(FA_gene, FA_magma$SYMBOL)
disorder_table <- fread('Three_network/disease_gene.csv', data.table = FALSE)
disorder_name <- colnames(disorder_table)

ADHD <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[1])
AD   <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[2])
MDD  <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[3])
SCZ  <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[4])
PD   <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[5])
ASD  <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[6])
BIP  <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[7])
MS   <- mappping_disorder_gene_PPI(TRN_network_data, FA_gene, disorder_name[8])

#### Co-expression network ###
ADHD <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[1])
AD   <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[2])
MDD  <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[3])
SCZ  <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[4])
PD   <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[5])
ASD  <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[6])
BIP  <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[7])
MS   <- mappping_disorder_gene_PPI(co_expression_network, FA_gene, disorder_name[8])

### PPI network ###
ADHD <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[1])
AD   <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[2])
MDD  <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[3])
SCZ  <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[4])
PD   <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[5])
ASD  <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[6])
BIP  <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[7])
MS   <- mappping_disorder_gene_PPI(PPi_intersection_score_data, FA_gene, disorder_name[8])

## Helper function for gene-TF enrichment ##
gene_data_tf_f <- function(disorder_net, disorder_gene, FA_gene, symbol_d) {
  colnames(disorder_net) <- c('node1', 'node2', 'weight')
  save_edge_d <- data.frame()
  if (symbol_d == 'TRN') {
    inter_tfs <- intersect(disorder_net$node1, FA_gene)
    for (i in inter_tfs) {
      te_gene <- disorder_net[which(disorder_net$node1 == i), ]
      save_edge_d <- rbind(save_edge_d, te_gene)
    }
    re_data <- unique(save_edge_d$node2)
  } else {
    bag_re_net <- disorder_net[which(disorder_net$node1 %in% FA_gene | disorder_net$node2 %in% FA_gene), ]
    re_data <- union(bag_re_net$node1, bag_re_net$node2)
  }
  return(re_data)
}

## UpSet and intersection for network ##
bag_upset <- function(FA_gene, root_path, symbol_d) {
  # Read node and edge files for all diseases
  AD_node   <- read.table(paste(root_path, 'ADnode.txt', sep = ""), header = TRUE)
  PD_node   <- read.table(paste(root_path, 'PDnode.txt', sep = ""), header = TRUE)
  ASD_node  <- read.table(paste(root_path, 'ASDnode.txt', sep = ""), header = TRUE)
  BIP_node  <- read.table(paste(root_path, 'BIPnode.txt', sep = ""), header = TRUE)
  scz_node  <- read.table(paste(root_path, 'SCZnode.txt', sep = ""), header = TRUE)
  dep_node  <- read.table(paste(root_path, 'DEPnode.txt', sep = ""), header = TRUE)
  ADHD_node <- read.table(paste(root_path, 'ADHDnode.txt', sep = ""), header = TRUE)
  MS_node   <- read.table(paste(root_path, 'MSnode.txt', sep = ""), header = TRUE)
  
  AD_gene   <- AD_node[which(AD_node$Alzheimer.s.Disease_anno == 2 | AD_node$Alzheimer.s.Disease_anno == 3), ]$gene
  pd_gene   <- PD_node[which(PD_node$Parkinson.s.Disease_anno == 2), ]$gene
  MS_gene   <- MS_node[which(MS_node$MS_anno == 2), ]$gene
  ASD_gene  <- ASD_node[which(ASD_node$Autism.Spectrum.Disorder_anno == 2 | ASD_node$Autism.Spectrum.Disorder_anno == 3), ]$gene
  BIP_gene  <- BIP_node[which(BIP_node$Bipolar.Disorder_anno == 2 | BIP_node$Bipolar.Disorder_anno == 3), ]$gene
  scz_gene  <- scz_node[which(scz_node$Schizophrenia_anno == 2 | scz_node$Schizophrenia_anno == 3), ]$gene
  dep_gene  <- dep_node[which(dep_node$Depression_anno == 2 | dep_node$Depression_anno == 3), ]$gene
  ADHD_gene <- ADHD_node[which(ADHD_node$Attention.deficit.hyperactivity.disorder_anno == 2 | ADHD_node$Attention.deficit.hyperactivity.disorder_anno == 3), ]$gene
  
  AD_edge   <- read.table(paste(root_path, 'ADnet.txt', sep = ""), header = TRUE)
  PD_edge   <- read.table(paste(root_path, 'PDnet.txt', sep = ""), header = TRUE)
  ASD_edge  <- read.table(paste(root_path, 'ASDnet.txt', sep = ""), header = TRUE)
  BIP_edge  <- read.table(paste(root_path, 'BIPnet.txt', sep = ""), header = TRUE)
  scz_edge  <- read.table(paste(root_path, 'SCZnet.txt', sep = ""), header = TRUE)
  dep_edge  <- read.table(paste(root_path, 'DEPnet.txt', sep = ""), header = TRUE)
  ADHD_edge <- read.table(paste(root_path, 'ADHDnet.txt', sep = ""), header = TRUE)
  MS_edge   <- read.table(paste(root_path, 'MSnet.txt', sep = ""), header = TRUE)
  
  bag_ad_re   <- gene_data_tf_f(AD_edge,   AD_gene,   FA_gene, symbol_d)
  bag_pd_re   <- gene_data_tf_f(PD_edge,   pd_gene,   FA_gene, symbol_d)
  bag_asd_re  <- gene_data_tf_f(ASD_edge,  ASD_gene,  FA_gene, symbol_d)
  bag_bip_re  <- gene_data_tf_f(BIP_edge,  BIP_gene,  FA_gene, symbol_d)
  bag_scz_re  <- gene_data_tf_f(scz_edge,  scz_gene,  FA_gene, symbol_d)
  bag_dep_re  <- gene_data_tf_f(dep_edge,  dep_gene,  FA_gene, symbol_d)
  bag_adhd_re <- gene_data_tf_f(ADHD_edge, ADHD_gene, FA_gene, symbol_d)
  bag_ms_re   <- gene_data_tf_f(MS_edge,   MS_gene,   FA_gene, symbol_d)
  
  whole_net <- rbind(AD_edge, PD_edge, ASD_edge, BIP_edge, scz_edge, dep_edge, ADHD_edge, MS_edge)
  colnames(whole_net) <- c('node1', 'node2', 'weight')
  upset_list <- list(bag_ad_re, bag_pd_re, bag_asd_re, bag_bip_re, bag_scz_re, bag_dep_re, bag_adhd_re, bag_ms_re)
  names(upset_list) <- c('AD', 'PD', 'ASD', 'BIP', 'SCZ', 'DEP', 'ADHD', 'MS')
  
  # Color vector for plotting
  disease_colors <- c(
    "AD" = "#8dd3c6",
    "ADHD" = "#ffffb3",
    "ASD" = "#bebada",
    "BIP" = "#fb8072",
    "DEP" = "#80b1d3",
    "MS" = "#fdb462",
    "SCZ" = "#b3de69",
    "PD" = "#fccde5"
  )
  set_order <- c("MS", "ADHD", "DEP", "BIP", "PD", "AD", "SCZ", "ASD")
  set_colors <- disease_colors[set_order]
  
  # Plot UpSet diagram
  p <- upset(
    fromList(upset_list),
    sets = set_order,
    sets.bar.color = set_colors,
    nsets = 100,
    nintersects = 40,
    order.by = "freq",
    keep.order = TRUE,
    mb.ratio = c(0.6, 0.4),
    text.scale = 2
  )
  print(p)
  return(list(up_set = upset_list, whole_net = whole_net))
}

FA_gene_PPI    <- bag_upset(FA_gene, './PPI/', 'PPI')
FA_gene_co_exp <- bag_upset(FA_gene, './Co_expression/', 'EXP')
FA_gene_TRN    <- bag_upset(FA_gene, './TRN/', 'TRN')

## Drug-gene interaction analysis ##
drug_gene_all <- fread('drug/interactions.tsv', data.table = FALSE)
FA_drug_gene <- drug_gene_all[which(drug_gene_all$gene_name %in% FA_gene), ]
FA_drug_gene <- FA_drug_gene[which((FA_drug_gene$approved == 'TRUE') & (FA_drug_gene$immunotherapy == 'FALSE') & (FA_drug_gene$anti_neoplastic == 'FALSE')), ]
FA_drug_gene$drug_concept_id <- lapply(FA_drug_gene$drug_concept_id, function(x) {
  x <- str_split(x, ':') %>% unlist()
  x <- x[2]
}) %>% unlist()

IRS1_drug <- FA_drug_gene[which(FA_drug_gene$gene_name == 'IRS1'), ]
target_count <- table(FA_drug_gene$gene_name) %>% as.data.frame()
target_count <- target_count[order(target_count$Freq, decreasing = TRUE), ]
write.table(FA_drug_gene$drug_concept_id, 'FA_interaction_drug.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)

## Drug class annotation for brain and bone diseases ##
drug_class_all <- fread('drug_classes_pivoted.csv', data.table = FALSE)
brain_bone_df <- drug_class_all %>%
  filter(str_detect(tolower(DISEASE),
                    "brain|alzheimer|parkinson|dementia|schizophrenia|epilepsy|cognitive|neurodegenerative|cerebral|stroke|multiple sclerosis|autism|bipolar|depression|anxiety|osteoporosis|bone density|fracture|osteopenia"
  ))
write.table(brain_bone_df, 'brain_drug.txt', row.names = FALSE, quote = FALSE)

brain_bone_df$RxCUI <- as.character(brain_bone_df$RxCUI)
result <- left_join(FA_drug_gene, brain_bone_df, by = c("drug_concept_id" = "RxCUI"))
brain_drug_fa_sta <- result[which(result$drug_concept_id %in% brain_bone_df$RxCUI), ]
brain_drug_fa_sta <- brain_drug_fa_sta[which(brain_drug_fa_sta$interaction_type != 'NULL'), ]

## Expand comma-separated DISEASE field to long format ##
df_expanded <- brain_drug_fa_sta %>%
  mutate(DISEASE = strsplit(as.character(DISEASE), ",")) %>%
  unnest(DISEASE) %>%
  mutate(DISEASE = trimws(DISEASE))

brain_related_diseases <- c(
  "Schizophrenia", "Parkinson Disease", "Brain Ischemia", "Cerebral Hemorrhage", "Stroke", "Anxiety Disorders", "Seizures", "Unconsciousness", "Epilepsy", "Coma", "Psychotic Disorders", "Neuroleptic Malignant Syndrome", "Bipolar Disorder", "Autistic Disorder", "Alcohol Withdrawal Delirium", "Cerebral Palsy", "Huntington Disease", "Multiple Sclerosis", "Myelitis", "Spinal Cord Injuries", "Spinal Cord Ischemia", "Spinal Cord Neoplasms", "Cerebellar Ataxia", "Migraine Disorders", "Migraine with Aura", "Cerebral Infarction", "Alzheimer Disease", "Status Epilepticus", "Epilepsies, Partial", "Neuralgia, Postherpetic", "Phobic Disorders", "Depression, Postpartum", "Cluster Headache, Spontaneous", "Sleep Initiation and Maintenance Disorders", "Tonic-Clonic", "Tremor", "Dementia", "Restless Legs Syndrome", "Intracranial Hemorrhages", "Brain Neoplasms", "Lennox Gastaut Syndrome", "Spasms, Infantile", "Conduct Disorder", "Depressive Disorder", "Stress Disorders, Post-Traumatic", "Tourette Syndrome", "Cognitive Dysfunction", "Dysarthria", "Myoclonic", "Absence", "Panic Disorder", "Tic Disorders", "Catatonia", "Gait Disorders, Neurologic", "Myoclonic Epilepsy, Juvenile", "Headache", "Muscle Rigidity", "Attention Deficit and Disruptive Behavior Disorders", "Head Injuries, Closed", "Meningeal Neoplasms", "Obsessive-Compulsive Disorder", "Anorexia Nervosa", "Bulimia Nervosa", "Attention Deficit Disorder with Hyperactivity"
)
df_brain <- df_expanded %>% filter(DISEASE %in% brain_related_diseases)
df_brain$gene_name %>% unique() %>% length()
df_brain$drug_concept_id %>% unique() %>% length()
df_brain$DISEASE %>% unique() %>% length()
write.table(df_brain, 'brain_net_gene.txt', row.names = FALSE, quote = FALSE, sep = '\t')

## Visualization - alluvial plot for selected diseases ##
df_brain <- fread("brain_net_gene.txt", data.table = FALSE, sep = "\t")
write.csv(df_brain, 'UFAG_gene.csv', row.names = FALSE, quote = FALSE)

is_alluvia_form(df_brain, axes = 1:3, silent = TRUE)
df_brain <- df_brain %>% filter(!is.na(interaction_score))
select_disease <- c("Schizophrenia", "Parkinson Disease", "Attention Deficit and Disruptive Behavior Disorders", "Multiple Sclerosis", "Alzheimer Disease", "Depression", "Bipolar Disorder", "Autism Spectrum Disorder")
df_brain_se <- df_brain[which(df_brain$DISEASE %in% select_disease), ]
df_brain_se$freq <- 1

ggplot(as.data.frame(df_brain_se),
       aes(axis1 = gene_name,
           axis2 = drug_name,
           axis3 = DISEASE,
           y = freq)) +
  scale_x_discrete(limits = c("Gene", "Drug", "Disease"), expand = c(.1, .1)) +
  geom_alluvium(aes(fill = gene_name), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank()) +
  labs(title = "Gene–Drug–Brain Disease Alluvial Diagram", y = "Interaction Score")

# (You can add more sections, such as pathway heatmap and other figures as needed.)
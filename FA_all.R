## 
library(dplyr)
library(data.table)
library(stringi)
library(stringr)
library(dplyr)
library(stringi)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(CMplot)
library(org.Hs.eg.db)
library(grid)
library(VennDiagram)
library(GenomicRanges)
library(Gviz)
#' Plot a genomic region on ideogram and axis using Gviz
#'
#' @param chr Chromosome number or name (e.g., "7")
#' @param start Start position (integer, 1-based)
#' @param end End position (integer)
#' @param region_name A label for the region (e.g., "7q31.31")
#'
#' @return A plot displaying the ideogram, genomic axis, and the annotated region
plot_genomic_region <- function(chr, start, end, region_name = "Region") {
  # Create a GRanges object for the target region
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = start, end = end),
                id = region_name)
  
  # Set up Gviz tracks
  ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  axisTrack <- GenomeAxisTrack()
  regionTrack <- AnnotationTrack(gr,
                                 name = region_name,
                                 genome = "hg19",
                                 chromosome = chr,
                                 fill = "red",
                                 col = NA)
  
  # Plot tracks with ±1Mb padding around region
  plotTracks(list(ideoTrack, axisTrack, regionTrack),
             from = start - 0,
             to = end + 0,
             chromosome = chr,
             genome = "hg19")
}
#
#’ Draw a Venn of two gene sets, run Fisher’s test, and annotate p‑value & OR
#’
#’ @param set1       Character vector of genes in set 1
#’ @param set2       Character vector of genes in set 2
#’ @param background Character vector of all background genes
#’ @param name1      Label for set 1
#’ @param name2      Label for set 2
#’ @param filename   Output file (e.g. "venn.png")
fisherVennPlot <- function(set1, set2, background,
                           name1 = "Set1", name2 = "Set2",
                           filename = "venn.png") {
  # 1. Counts
  inter_ab   <- length(intersect(set1, set2))
  only_a     <- length(set1) - inter_ab
  only_b     <- length(set2) - inter_ab
  neither_ab <- length(background) - length(union(set1, set2))
  
  # 2. Fisher table
  fisher_mat <- matrix(
    c(inter_ab, only_a,
      only_b,  neither_ab),
    nrow = 2, byrow = TRUE
  )
  
  # 3. Test
  ft   <- fisher.test(fisher_mat)
  pval <- ft$p.value
  or   <- as.numeric(ft$estimate)
  if (pval < (0.05/39)) {
    message("Significant set: ", name1)
    print(pval)
  }
  
  # # 4. Open device
  # png(filename, width = 900, height = 900, res = 300)
  # grid.newpage()
  # 
  # # 5. Draw (increase margin to make space at bottom)
  # venn.plot <- venn.diagram(
  #   x            = setNames(list(set1, set2), c(name1, name2)),
  #   filename     = NULL,
  #   fill         = c("#66C2A5", "#FC8D62"),
  #   alpha        = 0.5,
  #   cex          = 0.7,    # count label size
  #   cat.cex      = 0.7,    # category label size
  #   cat.fontface = "plain",
  #   cat.pos      = c(-20, 20),
  #   margin       = 0.2     # push diagram up a bit
  # )
  # grid.draw(venn.plot)
  # 
  # # 6. Annotate (move up to y = 0.12)
  # annotation <- sprintf(
  #   "Fisher's exact test:\np-value = %.3g\nOdds ratio = %.3g",
  #   pval, or
  # )
  # grid.text(
  #   annotation,
  #   x = unit(0.5, "npc"), y = unit(0.12, "npc"),
  #   gp = gpar(fontsize = 10)
  # )
  # 
  # # 7. Close
  # dev.off()
  message("Saved: ", filename)
}


plot_heart_manhattan<-function(heart_gwas,annotation_dir,save_name){
  #annotation_dir<-'./4ch_3Duin'
  #heart_gwas<-'HEART_4ich_inputation_3Duin.txt'
  heat_file<-fread(heart_gwas,data.table = FALSE)
  heat_file$BP<-heat_file$POS
  lead_snp<-fread(paste(annotation_dir,'/leadSNPs.txt',sep = ""),data.table = FALSE)
  lead_snp<-lead_snp[order(lead_snp$p),]
  lead_snp<-lead_snp[1:30,]
  ##
  gene_SNP<-fread(paste(annotation_dir,'/snps.txt',sep = ""),data.table = FALSE)
  gene_SNP<-gene_SNP[!duplicated(gene_SNP$rsID),]
  rownames(gene_SNP)<-gene_SNP$rsID
  gene_SNP<-gene_SNP[lead_snp$rsID,]
  ##
  #lead_snp<-lead_snp
  #heat_file_sig<-heat_file[which(heat_file$P<(5e-8)),]
  plot_heart<-heat_file[,c(2,1,9,8)]
  plot_heart<-plot_heart[which(plot_heart$P<0.5),]
  ## 1.482197e-323
  save_pdf<-paste(save_name,'.pdf',sep = "")
  CMplot(plot_heart,
         plot.type = "m",
         type = "p",
         LOG10 = TRUE,
         threshold = c(5e-8,5e-8/256),
         highlight.text=gene_SNP$nearestGene,
         highlight.text.cex =1.8,
         axis.cex = 1.5,
         lab.cex = 1.8,
         #  highlight.text.font=10,
         col = c("grey", "skyblue"),  # Alternating colors for different chromosomes
         highlight = lead_snp$rsID,  # Manually specify SNPs to highlight
         highlight.col = "red",       # Color for highlighted SNPs
         highlight.cex = 1.5,         # Point size for highlighted SNPs
         file = "pdf",
         file.name=save_pdf,
         dpi = 300)
  #####
 
}
##
plot_heart_manhattan('FA_JAGWAS_res/UDIP_FA.txt','./FA_JAGWAS_res/FUMA_res/Discovery/','JAGWAS_discover')
plot_heart_manhattan('FA_JAGWAS_res/UDIP_FA_rep.txt','./FA_JAGWAS_res/FUMA_res/Replication/','JAGWAS_discover_replicaiton')
plot_heart_manhattan('FA_JAGWAS_res/UDIP_FA_meta_JAGWAS_all.txt','./FA_meta_5e-8/','JAGWAS_meta')
## validation the result ##
lead_SNP_all<-fread('./FUMA_res/Discovery/leadSNPs.txt',data.table = FALSE)
all_rep_fa<-fread('./UDIP_FA_rep.txt',data.table = FALSE)
####
select_data<-all_rep_fa[which(all_rep_fa$SNP%in%lead_SNP_all$rsID),]
repliccated_snp<-select_data[which(select_data$P<=(0.05/128)),]%>%nrow()
ratio<-repliccated_snp/nrow(lead_SNP_all)  
## meta result analysis ##
lead_SNP_all_meta<-fread('Single_variant_FA_meta/leadSNPs.txt')
Genomic_loci_meta<-fread('Single_variant_FA_meta/GenomicRiskLoci.txt')
previous_zhao_snp<-fread('previous_science_zhao.txt')
previous_zhao_locus<-fread('previous_Locus.csv')
##
previous_zhao_locus <- previous_zhao_locus[grepl("FA$", ID)]
zhao_lead_snp<-previous_zhao_locus$LeadSNPs
flat_rsids <- unlist(strsplit(zhao_lead_snp, ";"))
##
intersect(flat_rsids,lead_SNP_all_meta$rsID)
## big40 GWAS ##. 1452-1526
Big_40_res<-fread('./big40_brain_study.csv')
Big_40_res$pheno <- as.numeric(gsub("V", "", Big_40_res$pheno))
FA_big40<-Big_40_res[which((Big_40_res$pheno>1452)&(Big_40_res$pheno<1526)),]
##
intersect(FA_big40$rsid,lead_SNP_all_meta$rsID)
#########
FA_big40$start<-FA_big40$pos-125000
FA_big40$end<-FA_big40$pos+125000

###
### found the overlap regions ##
buffer <- 125000  # ±125kb

# 
gr1 <- GRanges(seqnames = Genomic_loci_meta$chr,
               ranges = IRanges(start = Genomic_loci_meta$start - buffer,
                                end   = Genomic_loci_meta$end + buffer),
               locus_id = Genomic_loci_meta$GenomicLocus)

gr2 <- GRanges(seqnames = previous_zhao_locus$chr,
               ranges = IRanges(start = previous_zhao_locus$start - buffer,
                                  end   = previous_zhao_locus$end+ buffer),
               locus_id = previous_zhao_locus$GenomicLocus)
unique_ranges <- !duplicated(paste0(start(gr2), "-", end(gr2)))
gr2 <- gr2[unique_ranges]
# overlap regions 3
hits <- findOverlaps(gr1, gr2)

# 
overlap_df <- data.frame(
  fuma1_locus = mcols(gr1)$locus_id[queryHits(hits)],
  fuma1_chr   = as.character(seqnames(gr1))[queryHits(hits)],
  fuma1_start = start(gr1)[queryHits(hits)],
  fuma1_end   = end(gr1)[queryHits(hits)],
  fuma2_locus = mcols(gr2)$locus_id[subjectHits(hits)],
  fuma2_chr   = as.character(seqnames(gr2))[subjectHits(hits)],
  fuma2_start = start(gr2)[subjectHits(hits)],
  fuma2_end   = end(gr2)[subjectHits(hits)]
)

### our reult and big40 ##
# Function: Identify unique loci not overlapping with those from the BIG40 study
get_unique_loci <- function(big40, my_file, p_thresh = -log10(5e-8 / nrow(big40)), output_file = NULL) {
  # Load summary statistics
  #big40<-FA_big40
  buffer <- 125000
 # myres<-lead_SNP_all_meta
  colnames(big40)<-c("rsid","pheno","chr","pos", "a1","a2", "p","beta", "se","p-value(repro)",
                     "beta(repro)","se(repro)", "nominal" ,"start","end" )
  # Check that required columns are present
  required_cols <- c("chr", "pos", "p")
  if (!all(required_cols %in% names(big40))) stop("BIG40 file must contain: CHR, BP, P")
  if (!all(required_cols %in% names(myres))) stop("My summary file must contain: CHR, BP, P")
  
  # Filter genome-wide significant SNPs using the specified p-value threshold
  big40_sig <- big40[p<= p_thresh]
  my_sig <- myres
  ##
  print(intersect(big40_sig$rsid,my_sig$rsID))
  # Construct GRanges objects and expand each SNP region by ±125 kb
  big40_gr <- GRanges(seqnames = big40_sig$chr,
                      ranges = IRanges(start = big40_sig$pos - 125000,
                                       end = big40_sig$pos + 125000))
  my_gr <- GRanges(seqnames = Genomic_loci_meta$chr,
                   ranges = IRanges(start = Genomic_loci_meta$start - buffer,
                                    end   = Genomic_loci_meta$end + buffer),
                   locus_id = Genomic_loci_meta$GenomicLocus)
  
  # Merge overlapping/adjacent regions to avoid redundancy
  big40_gr <- reduce(big40_gr)
  my_gr <- reduce(my_gr)
  
  # Identify overlaps: any of 'my' loci overlapping with BIG40 loci
  hits <- findOverlaps(my_gr, big40_gr)
  overlapping_idx <- unique(queryHits(hits))
  
  # Extract non-overlapping loci — those at least 250 kb away from any BIG40 locus
  unique_my_gr <- my_gr[overlapping_idx]
  
  # Convert GRanges to a data frame and calculate the center of each interval
  unique_df <- as.data.frame(unique_my_gr)
  unique_df$CHR <- gsub("chr", "", unique_df$seqnames)
  unique_df$CENTER_BP <- round((unique_df$start + unique_df$end) / 2)
  
  # Optionally write the results to file
  if (!is.null(output_file)) {
    fwrite(unique_df[, c("CHR", "CENTER_BP", "start", "end")],
           file = output_file, sep = "\t")
  }
  # Return a data.frame of unique loci (non-overlapping with BIG40)
  return(unique_df[, c("CHR", "CENTER_BP", "start", "end")])
}
##
#### 4 GWAS overlap ##
FA_all_SNP<-fread('./FA_meta_5e-8/leadSNPs.txt',data.table = FALSE)
## compare with previous functions ##
GWAS_catelog_func<-fread('./FA_meta_5e-8/GWAS_catelog_lead_SNP.txt')
GWAS_catelog_genomic<-fread('./FA_meta_5e-8/gwascatalog.txt')
##
write.csv(GWAS_catelog_genomic,'./Supplementary/UDIP_FA_supplementary.csv',row.names = FALSE)

##
intersect(GWAS_catelog_func$SNPS,FA_all_SNP$rsID)%>%unique()
###   120634679   121222318 region visual
# two regions bone 120634679   121222318  7  7q31.31 ||| white 87026463  87544863 16  16q24.2
plot_genomic_region(chr = 7,
                    start = 120634679,
                    end = 121222318,
                    region_name = "7q31.31")
plot_genomic_region(chr = 16,
                    start = 87026463,
                    end = 87544863,
                    region_name = "16q24.2")
### classification of the GWAS catelag result #
brain_structure<-fread('FA_meta_5e-8/brain_structure_traits.txt')
brain_disease<-fread('FA_meta_5e-8/brain_disease_traits.txt')
others<-fread('FA_meta_5e-8/other_traits.txt')
###
intersect(brain_structure$SNPS,FA_all_SNP$rsID)%>%unique()%>%length()
intersect(brain_disease$SNPS,FA_all_SNP$rsID)%>%unique()%>%length()
intersect(others$SNPS,FA_all_SNP$rsID)%>%unique()%>%length()
bone_rows <- others[grepl("bone", others$`DISEASE/TRAIT`, ignore.case = TRUE), ]
##
bone_rows$SNPS%>%unique()
dd<-others$`DISEASE/TRAIT`%>%unique()
write.table(dd,'dd_other_sen.txt',row.names = FALSE,quote = FALSE)
##
plot_function_enrichment <- function(magma_file, output_file = "function_enrichment_plot.pdf") {
  
  magma_file<-'./FA_meta_5e-8/magma.gsa.out'
  magma_functions <- fread(magma_file, data.table = FALSE, skip = 4)
  magma_functions$FDR <- -log10(p.adjust(magma_functions$P, method = 'bonferroni'))
  magma_functions <- magma_functions[order(magma_functions$FDR, decreasing = TRUE), ]
  magma_functions <- magma_functions[magma_functions$FDR > 1.3, ]
  top_30 <- magma_functions[1:25,]
  print(nrow(magma_functions))
  top_30 <- top_30 %>%
    mutate(FULL_NAME = sapply(FULL_NAME, function(x) {
      parts <- stringi::stri_split_fixed(x, 'BP_', simplify = TRUE)
      if (ncol(parts) >= 2) parts[2] else NA
    }))
  top_30<-top_30[order(top_30$FDR),]
  p <- ggbarplot(top_30, x = 'FULL_NAME', y = 'FDR', fill = 'skyblue') +
    geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1) +
    xlab('GO Biological Process') + ylab('Enrichment (-log10 Adjusted P)') +
    coord_flip()
  
  ggsave(output_file, p, width = 10, height = 12)
  return(magma_functions)
}
## case highlight
library(locuscomparer)
##
FA_meta_GWAS<-fread('./FA_JAGWAS_res/UDIP_FA_meta_JAGWAS_all.txt',data.table = FALSE)
total_body_GWAS<-fread('./bone/Total_body_bone_mineral_format.txt',data.table = FALSE)
FA_meta_GWAS <- subset(FA_meta_GWAS, CHR == 7 & POS >= 120750000 & POS <= 121200000)
total_body_GWAS <- subset(total_body_GWAS, CHR == 7 & POS >= 120750000 & POS <= 121200000)
total_body_GWAS<-total_body_GWAS[,c('SNP','P')]
FA_meta_GWAS<-FA_meta_GWAS[,c('SNP','P')]
colnames(FA_meta_GWAS)<-c('rsid','pval')
colnames(total_body_GWAS)<-c('rsid','pval')
####
write.table(FA_meta_GWAS,'FA_gwas.tsv',sep = '\t',row.names = FALSE,quote = FALSE)
write.table(total_body_GWAS,'totalBMD_gwas.tsv',sep = '\t',row.names = FALSE,quote = FALSE)



###
locuscompare(in_fn1 = "./FA_gwas.tsv",  # 
             in_fn2 = "./totalBMD_gwas.tsv", 
             title1 = "UDIP-FA",
             title2 = "Total body BMD",
            # region = "chr7:120750000-121200000",
            # region = "chr9:120750000-121200000", 
             snp = "rs3779381")  # 主SNP
##
library(locuscomparer)
library(ggrepel)



highlight_snps <- c("rs3779381", "rs123456", "rs2295786")

# 绘图并标记
library(ggplot2)
ggplot(p$plot_data, aes(x = -log10(p1), y = -log10(p2))) +
  geom_point(aes(color = ld)) +
  geom_text_repel(data = subset(p$plot_data, SNP %in% highlight_snps),
                  aes(label = SNP), size = 3) +
  theme_bw() +
  labs(x = "-log10(P) in FA", y = "-log10(P) in BMD", color = "LD")


###

plot_disease_enrichment <- function(enrichment_file, output_file = "disease_enrichment_plot.pdf") {
  
  top_geneE <- fread(enrichment_file, sep = '\t', fill = TRUE, data.table = FALSE)
  disease <- top_geneE[top_geneE$Category == 'Disease', ]
  disease$FDR <- -log10(disease$`q-value Bonferroni`)
  disease <- head(disease[order(disease$FDR, decreasing = TRUE), ], 25)
  #disease=disease[1:25,]
  p <- ggbarplot(disease, x = 'Name', y = 'FDR', fill = 'skyblue') +
    geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1) +
    xlab('Phenotype') + ylab('Enrichment (-log10 Adjusted P)') +
    coord_flip()
  
  ggsave(output_file, p, width = 12, height = 6)
}
plot_tissue_enrichment <- function(tissue_file, output_file = "tissue_enrichment_plot.pdf", skip_lines = 5) {
  
  tissue <- fread(tissue_file, data.table = FALSE, fill = TRUE, skip = skip_lines)
  tissue$FDR <- -log10(p.adjust(tissue$P, method = 'bonferroni'))
  tissue <- tissue[order(tissue$FDR), ]
  
  p <- ggbarplot(tissue, x = 'FULL_NAME', y = 'FDR', fill = 'skyblue') +
    geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1) +
    xlab('Tissue') + ylab('Enrichment (-log10 Adjusted P)') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))+coord_flip()+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  print(p)
  ggsave(output_file, p, width = 6, height = 10)
}
plot_celltype_enrichment <- function(cell_file, output_file = "celltype_enrichment_plot.pdf", exclude_blood = TRUE) {
  
  #cell_file<-'./4ch_3D_GAN_heart_cell/magma_celltype_step1.txt'
  cell <- fread(cell_file, data.table = FALSE)
  
  if (exclude_blood) {
    cell <- cell[!(cell$Dataset %in% c(
      "GSE89232_Human_Blood", 
      "539_Travaglini_2020_Blood_level1", 
      "541_Xu_Human_2023_Blood_level1")), ]
  }
  
  enrichment_df <- tapply(cell$P.adj, cell$Cell_type, mean) %>% as.data.frame()
  enrichment_beta <- tapply(cell$BETA, cell$Cell_type, mean) %>% as.data.frame()
  enrichment_df$celltype <- rownames(enrichment_df)
  colnames(enrichment_df) <- c('adjust_p', 'celltype')
  enrichment_df$adjust_p <- -log10(enrichment_df$adjust_p)
  enrichment_df <- enrichment_df[order(enrichment_df$adjust_p), ]
  
  p <- ggbarplot(enrichment_df, x = 'celltype', y = 'adjust_p', fill = 'skyblue') +
    geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1) +
    xlab('Cell Types') + ylab('Enrichment (-log10 Adjusted P)') +
    coord_flip() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  print(p)
  ggsave(output_file, p, width = 12, height = 6)
}
# FA meta mvGWAS analysis #
view_4ch_function<-plot_function_enrichment('./FA_meta_5e-8/magma.gsa.sets.genes.out',output_file ='4ch_functional_enrichment.pdf')
plot_disease_enrichment('./heart_function.txt',output_file ='4ch_phenotype_enrichment.pdf' )
plot_tissue_enrichment('./FA_meta_5e-8/magma_exp_bs_age_avg_log2RPKM.gsa.out',output_file ='4ch_tissue_enrichment.pdf' )
plot_celltype_enrichment('./4ch_3D_GAN_heart_cell/magma_celltype_step1.txt',output_file ='4ch_cell_type_enrichment.pdf' )
###
FA_MAGAMA_gene<-fread('./FA_meta_5e-8/magma.genes.out')
FA_MAGAMA_gene_significant<-FA_MAGAMA_gene[which(FA_MAGAMA_gene$P<(0.05/19128)),]
preious_zhao_magma_gene_science<-fread('../compare_MAGMA_gene_zhao.csv')
preious_zhao_magma_gene_mp<-fread('../previous_zhao_MP.csv')
intersect(preious_zhao_magma_gene_science$GENE,FA_MAGAMA_gene_significant$GENE)%>%length()
intersect(preious_zhao_magma_gene_mp$GENE,FA_MAGAMA_gene_significant$GENE)%>%length()
Reduce(intersect,list(FA_MAGAMA_gene_significant$SYMBOL, preious_zhao_magma_gene_science$SYMBOL, preious_zhao_magma_gene_mp$SYMBOL))
###
write.csv(FA_MAGAMA_gene_significant,'./Supplementary/FA_MAGMA_gene_all.txt',row.names = FALSE,quote = FALSE)
####
FA_gene_phenotype_enrichment<-fread('Top_enrichment.txt')
FA_gene_phenotype_pheno<-FA_gene_phenotype_enrichment[which(FA_gene_phenotype_enrichment$Category=='Disease'),]
FA_gene_phenotype_pheno<-FA_gene_phenotype_pheno[which(FA_gene_phenotype_pheno$`q-value Bonferroni`<0.05),]
#######
FA_geneset_analysis<-fread('./FUMA_FA_meta/magma.gsa.out',skip=5,data.table = FALSE)
FA_geneset_analysis$Bonfer<-p.adjust(FA_geneset_analysis$V7, method = 'bonferroni')
FA_geneset_analysi_significant<-FA_geneset_analysis[which(FA_geneset_analysis$Bonfer<(0.05)),]
write.csv(FA_geneset_analysi_significant,'./Supplementary/MAGMA_gene_set_sig_enrichment.csv',row.names = FALSE,quote = FALSE)
##
FA_developmental_stage<-fread('./FA_meta_5e-8/magma_exp_bs_dev_avg_log2RPKM.gsa.out',skip = 5)
###
##
# Define a function to plot a Manhattan plot from MAGMA gene-level output
plot_magma_manhattan_with_labels <- function(magma_file, significance_level = 0.05, total_genes = 19128, output_name = "magma_manhattan") {
  # Load the MAGMA gene-level results
  #magma_file<-''
  magma_result <- fread(magma_file,data.table = FALSE)
  
  # Check necessary columns exist
  required_cols <- c("SYMBOL", "CHR", "START", "P")
  if (!all(required_cols %in% colnames(magma_result))) {
    stop("MAGMA output must contain the following columns: GENE, CHR, START, and P")
  }
  
  # Remove missing or invalid data
  magma_result <- magma_result %>%
    filter(!is.na(CHR), !is.na(START), !is.na(P))
  
  # Compute genome-wide significance threshold
  threshold <- significance_level / total_genes
  jk
  # Prepare data for CMplot
  cmplot_data <- magma_result %>%
    mutate(
      SNP = SYMBOL,
      Chromosome = CHR,
      Position = START,
      P.value = P
    ) %>%
    dplyr::select(SNP, Chromosome, Position, P.value)
  
  # Identify significant genes (make sure they actually exist in the data)
  significant_genes <- cmplot_data %>%
    filter(P.value < threshold) 
  
  significant_genes<-significant_genes[order(significant_genes$P.value),]
  # Confirm only genes that exist in cmplot_data are passed
  significant_genes<-significant_genes[1:30,1]
  significant_genes <- significant_genes[significant_genes %in% cmplot_data$SNP]
  
  # Generate Manhattan plot
  CMplot(
    cmplot_data,
    plot.type = "m",                # Manhattan plot
    LOG10 = TRUE,                   # -log10(P-value)
    threshold = threshold, # Significance threshold
    threshold.col = "red",           # Threshold line color
    threshold.lty = 1,               # Threshold line style (dashed)
    chr.den.col = NULL,              # No background density
    amplify = TRUE,                  # Enlarge significant points
    highlight = significant_genes,   # Highlight significant genes
    highlight.col = "blue",           # Highlight color
    highlight.cex = 1.5,             # Highlight point size
    highlight.text = significant_genes,            # Add gene labels
    highlight.text.cex = 0.8,         # Text size for labels
    signal.cex = 1.5,
    file.output = TRUE,              # Save plot to file
    file = "jpg",                    # Output format ("pdf" or "jpg")
    dpi = 300,                       # High resolution
    #outdir = "./",                   # Output directory (current folder)
    verbose = TRUE                   # Print progress
  )
  
  message("Manhattan plot saved as: ", output_name, ".jpg")
}
FA_maga_gene<-fread('./FA_meta_5e-8/magma.genes.out')
plot_magma_manhattan_with_labels('./FA_meta_5e-8/magma.genes.out')
view_4ch_magma<-view_4ch_magma[which(view_4ch_magma$P<(0.05/19291)),]
## 2ch view function ananlysis

view_2ch_magma<-view_2ch_magma[which(view_2ch_magma$P<(0.05/19291)),]
intersect(view_4ch_magma$SYMBOL,heart_gene$symbol)%>%length()
intersect(view_2ch_magma$SYMBOL,heart_gene1$symbol)%>%length()
###
venn_list<-list(MAGMA_4Ch=view_4ch_magma$SYMBOL,MAGMA_2Ch=view_2ch_magma$SYMBOL,eQTL_2ch=heart_gene1$symbol,eQTL_4ch=heart_gene$symbol)
save_name<-'heart_gene_overlap'
save_name<-paste(save_name,'.png',sep = "")
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("MAGMA_4Ch", "MAGMA_2Ch", "eQTL_2Ch", "eQTL_4Ch"),
  filename = save_name, 
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3500,
  resolution = 500,
  fill = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),  # 4 colors for 4 sets
  alpha = 0.6,
  cex = 1.5,            # Font size inside circles
  cat.cex = 1,        # Font size for category names
  cat.col = c("black", "black", "black", "black"),  # 4 colors for 4 labels
  lwd = 1.5,
  lty = "solid"
)


###
view_2Ch_gene_eqtl<-view_2Ch_gene[(!is.na(view_2Ch_gene$eqtlMapminP))|(view_2Ch_gene$ciMap=='Yes'),]
view_4Ch_gene_eqtl<-view_4Ch_view_gene[(!is.na(view_4Ch_view_gene$eqtlMapminP))|(view_4Ch_view_gene$ciMap=='Yes'),]
#
view_2Ch_eqtl_snp<-view_2Ch_gene_eqtl$IndSigSNPs%>%as.character()
view_4Ch_eqtl_snp<-view_4Ch_gene_eqtl$IndSigSNPs%>%as.character()
# clearn function
extract_rsids <- function(rs_column) {
  rs_column %>%
    gsub(" ", "", .) %>%                 
    gsub(":", ";", .) %>%                    
    paste(collapse = ";") %>%              
    strsplit(";") %>%                       
    unlist() %>%                            
    unique() %>%                             
    sort()                                   
}
# perform the filter
view_2Ch_eqtl_snp <- extract_rsids(view_2Ch_eqtl_snp)
view_4Ch_eqtl_snp<-extract_rsids(view_4Ch_eqtl_snp)
##
intersect(res_4ch_gan$rsID,view_4Ch_eqtl_snp)%>%length()
intersect(res_2ch_gan$rsID,view_2Ch_eqtl_snp)%>%length()
##
###
development_files<-fread('./magma_celltype_PsychENCODE_Developmental.gsa.out',skip = 5)
####
development_files<-development_files[order(development_files$P),]
development_files<-development_files[1:11,]
development_files<-development_files[-2,]
development_files$P<--log10(development_files$P)
ggbarplot(development_files,x='VARIABLE',y='P',x.text.angle=30,fill = 'skyblue')+
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)+xlab('cell types')+ylab('Enrichment(-log10(P))')
### PRS plot  ##
PRS2_FA_files<-fread('PRS_FA',sep = ',',data.table = FALSE)
PRS2_FA_files$FDR<--log10(p.adjust(PRS2_FA_files$p_value,method = 'bonferroni'))
PRS2_FA_files$adjust_p<-p.adjust(PRS2_FA_files$p_value,method = 'bonferroni')
PRS2_FA_files<-PRS2_FA_files[order(PRS2_FA_files$FDR),]
PRS2_FA_files$FDR_capped <- ifelse(PRS2_FA_files$FDR > 50, 50, PRS2_FA_files$FDR)
##

ggscatter(PRS2_FA_files, 
          x = "prs_file",                  
          y = "FDR_capped",  
          size = 4,
          color = "canonical_correlation",               
          gradient.color = c("red", "blue"),                
          alpha = 0.8,  
          #shape = 17
          repel = TRUE) +                  
  labs(title = "Bubble Plot of Canonical Correlation vs PRS",
       x = "PRS File",
       y = "-log10(Adjusted P)",
       color = "Canonical correlation",
       size = "-log10(Adjusted P)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  coord_flip() +
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)

# Plot using ggscatter from ggpubr package
ggscatter(PRS2_FA_files, 
          x = "prs_file",                  # x-axis: PRS file names
          y = "FDR",     # y-axis: canonical correlation
          #   size = "FDR",  # Bubble size based on canonical correlation
          color = "canonical_correlation",               # Bubble color based on p-value
          gradient.color = c("red", "blue"),                # Color palette (Red to Blue)
          alpha = 0.8,                     # Transparency of bubbles
          #     label = "prs_file",              # Label points with PRS file names
          repel = TRUE) +                  # Avoid overlapping labels
  labs(title = "Bubble Plot of Canonical Correlation vs PRS",
       x = "PRS File",
       y = "-log10(Adjusted P)",
       color = "Canonical correlation",
       size = "-log10(Adjusted P)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+coord_flip()+ geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)
##

### heart related gene analysis 
heart_disorder_gene<-fread('./Cardiac_Disease_Genes__Excluding_Traits_.csv',data.table = FALSE)
diorder_gene_sta<-table(heart_disorder_gene$`Disease/Trait ID`)%>%as.data.frame()
diorder_gene_choose<-diorder_gene_sta[which(diorder_gene_sta$Freq>200),]
###  ###
diorder_gene_selected<-heart_disorder_gene[which(heart_disorder_gene$`Disease/Trait ID`%in%diorder_gene_choose$Var1),]
### background gene ##
background_gene<-fread('NCBI37.3.gene.loc',data.table = FALSE)
##
# venn_list<-list(MAGMA_4Ch=view_4ch_magma$SYMBOL,MAGMA_2Ch=view_2ch_magma$SYMBOL,eQTL_2ch=heart_gene1$symbol,eQTL_4ch=heart_gene$symbol)
all_disorder_name<-diorder_gene_selected$`Disease/Trait`%>%unique()
gene_4Ch<-union(view_4ch_magma$SYMBOL,heart_gene$symbol)
gene_2Ch<-union(view_2ch_magma$SYMBOL,heart_gene1$symbol)
select_disease<-c('Dilated cardiomyopathy','Hypertrophic cardiomyopathy')
for(i in all_disorder_name){
  #i=all_disorder_name[1]
  heart_dis_temp<-diorder_gene_selected[which(diorder_gene_selected$`Disease/Trait`==i),]
  intersect()
  fisherVennPlot(heart_dis_temp$Gene,gene_2Ch,background_gene$V6,name1=i,name2='4Ch UDIP-heart',paste('./disoder_gene_overlap/4Ch/',i,'_2Ch_UDIP.png',sep = ""))
}

### single cell TRN ##
singlecell_heart_trn<-fread('single_cell_heart_TRN.csv')
count_degrer_tf<-table(singlecell_heart_trn$`TF name`)%>%as.data.frame()
UDIP_tf<-intersect(singlecell_heart_trn$`TF name`,gene_4Ch)

rownames(count_degrer_tf)<-count_degrer_tf$Var1
UDIP_tf_count<-count_degrer_tf[UDIP_tf,]
other_rd<-count_degrer_tf[setdiff(count_degrer_tf$Var1,UDIP_tf),]
###
group_df<-c(rep('UDIP-heart',nrow(UDIP_tf_count)),rep('Others',nrow(other_rd)))
combind_data<-rbind(UDIP_tf_count,other_rd)
combind_data$group<-group_df
t.test(UDIP_tf_count$Freq,other_rd$Freq)
ggboxplot(combind_data,x='group',y='Freq',xlab = 'Group',ylab = 'Degreee')
####
singlecell_heart_trn$`TF name`%>%unique()%>%length()
###
pos_4Ch<-heart_gene[which(heart_gene$posMapSNPs!=0),]
pos_2Ch<-heart_gene1[which(heart_gene1$posMapSNPs!=0),]
genomic_loc<-fread('2Ch_3D_GAN/GenomicRiskLoci.txt')
intersect(pos_4Ch$IndSigSNPs,res_4ch_gan$rsID)
####compare result #
Aung_res<-fread('./comparsion/Aung.csv',data.table = FALSE)
Bonazzola_res<-fread('./comparsion/Aung.csv',data.table = FALSE)
meyer_res<-fread('./comparsion/Meyer.csv',data.table = FALSE)
#
meyer_res<-meyer_res[which(meyer_res$Locus!='-'),]
Burns_res<-fread('./comparsion/PCA_left_MAGMA_gene.csv',data.table = FALSE)
Burns_res_lead_snp<-fread('./comparsion/PCA_left_lead_SNP.csv',data.table = FALSE)
PCA_fuma<-fread('./comparsion/PCA_FUMA.csv',data.table = FALSE)
Burns_res<-union(Burns_res$`Gene symbol`,PCA_fuma$Locus)
Burns_res<-Burns_res[which(Burns_res!="")]
Pirru_res<-fread('./comparsion/Pirruccello.csv',data.table = FALSE)
Pirru_res<-union(Pirru_res$Gene1,Pirru_res$Gene2)
Pirru_res<-Pirru_res[which(Pirru_res!="")]
####
Sooknah_log<-fread('./comparsion/Sooknah_et_al_log.csv')
Sooknah_short<-fread('./comparsion/Sooknah_short.csv')
##
Sooknah_lead_snp<-union(Sooknah_log$`Lead Variant`,Sooknah_short$`Lead Variant`)%>%unique()
Sooknah_gene<-union(Sooknah_log$`Gene Name`,Sooknah_short$`Gene Name`)%>%unique()


library(UpSetR)

plot_upset <- function(set_list, set_names = NULL, nsets = 5, nintersects = 20, 
                       order.by = "freq", mainbar.y.label = "Intersection Size",
                       sets.x.label = "Set Size", keep.order = FALSE) {
  # set_list: A named list where each element is a vector representing a set
  # set_names: Optional. A vector of names for the sets if set_list has no names
  # nsets: Number of individual sets to show on the x-axis
  # nintersects: Number of intersections to show in the main bar plot
  # order.by: How to order the intersections ("freq" or "degree")
  # mainbar.y.label: Label for the y-axis of the intersection size barplot
  # sets.x.label: Label for the x-axis (set sizes)
  # keep.order: Whether to keep the original input order of sets
  
  # Check if set names are provided
  if (is.null(names(set_list))) {
    if (is.null(set_names)) stop("Please provide either named set_list or set_names.")
    names(set_list) <- set_names
  }
  
  # Get all unique elements across all sets
  all_elements <- unique(unlist(set_list))
  
  # Create a binary presence/absence dataframe for UpSetR
  upset_df <- data.frame(element = all_elements)
  for (set_name in names(set_list)) {
    upset_df[[set_name]] <- as.integer(all_elements %in% set_list[[set_name]])
  }
  
  # Plot the UpSet diagram
  UpSetR::upset(upset_df, 
                nsets = nsets, 
                nintersects = nintersects, 
                order.by = order.by,
                mainbar.y.label = mainbar.y.label,
                sets.x.label = sets.x.label,
                keep.order = keep.order)
}
###
test_dd<-apply(upset_df[2:ncol(upset_df)], 1, sum)
###
names(test_dd)<-upset_df[,1]
test_dd<-as.data.frame(test_dd)
test_dd$gene<-rownames(test_dd)
##
test_dd<-test_dd[which(test_dd$test_dd==1),]
##
gene_one<-test_dd$gene
cc<-upset_df[which((upset_df$UIDP_2Ch==1)|(upset_df$UIDP_4Ch==1)),]
intersect(cc$element,gene_one)

###
plot_upset_advanced <- function(set_list, set_colors = NULL, top_n = 5, mainbar.y.label = "Intersection Size") {
  # set_list: named list where each element is a vector (one set)
  # set_colors: optional vector of colors (one color per set)
  # top_n: number of top intersection elements to label
  # mainbar.y.label: y-axis label for intersection size bar plot
  
  if (is.null(names(set_list))) {
    stop("set_list must be a named list with names for each set.")
  }
  
  # Prepare the data
  all_elements <- unique(unlist(set_list))
  upset_data <- data.frame(element = all_elements)
  
  for (set_name in names(set_list)) {
    upset_data[[set_name]] <- as.integer(all_elements %in% set_list[[set_name]])
  }
  
  # Find the intersection size
  upset_data_long <- upset_data %>%
    pivot_longer(cols = -element, names_to = "Set", values_to = "Present") %>%
    filter(Present == 1)
  
  intersection_counts <- upset_data_long %>%
    group_by(element) %>%
    summarise(sets_in = n()) %>%
    arrange(desc(sets_in))
  
  # Get top_n elements appearing in the most sets
  top_elements <- intersection_counts$element[1:min(top_n, nrow(intersection_counts))]
  
  # Check colors
  if (!is.null(set_colors)) {
    if (length(set_colors) != length(set_list)) {
      stop("Length of set_colors must match number of sets.")
    }
    names(set_colors) <- names(set_list)
  }
  
  # Create UpSet plot
  p <- upset(
    upset_data,
    intersect = names(set_list),
    name = "Genes",
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE
      ) +
        ylab(mainbar.y.label)
    ),
    set_sizes = (
      upset_set_size(aes(fill = after_stat(set))) +  # ⭐️ important: fill set size by set
        theme(axis.text.x = element_text(angle = 90))
    )
  )
  
  # Apply color if provided
  if (!is.null(set_colors)) {
    p <- p + scale_fill_manual(values = set_colors)
  }
  
  # Print plot
  print(p)
  
  # Print top elements
  cat("Top elements appearing in most sets:\n")
  print(top_elements)
}

nature_colors <- c(
  "#3C5488", # Blue
  "#E64B35", # Orange
  "#00A087", # Green
  "#B71C1C", # Red
  "#7E57C2", # Purple
  "#FFC107", # Yellow
  "#4DBBD5"  # Cyan
)


# Draw the UpSet plot
plot_upset(
  set_list = list(UIDP_2Ch = gene_2Ch, UIDP_4Ch = gene_4Ch,Sooknah_et_al=Sooknah_log$`Gene Name`,Aung_et_al = Aung_res$`Locus Name`, Bonazzola_et_al= Bonazzola_res$`Locus Name`,
                  Meyer_et_al=meyer_res$Locus,Burns_et_al=Burns_res, Pirru_et_al=Pirru_res),
  nsets = 8,
  nintersects = 200,
  order.by = "freq"
)
## lead SNP
plot_upset(
  set_list = list(UIDP_2Ch = res_2ch_gan$rsID, UIDP_4Ch = res_4ch_gan$rsID, Sooknah_et_al=Sooknah_log$`Lead Variant`,Aung_et_al = Aung_res$`Lead variant`,Bonazzola_et_al= Bonazzola_res$`Lead variant`,
                  Meyer_et_al=meyer_res$SNP, Burns_et_al=Burns_res_lead_snp$`Lead SNV`,Pirru_et_al=Pirru_res$`Lead SNP`),
  nsets = 8,
  nintersects = 200,
  order.by = "freq"
)
#
# loci overlap ##
T1_loci<-fread('./T1_T2_JAGWAS/T1_discovery_loci.csv')
T2_loci<-fread('./T1_T2_JAGWAS/T2_discovery_loci.csv')
FA_loci<-fread('./FA_meta/FA_meta_5e-8/GenomicRiskLoci.txt')
###
get_hit_res<-function(region_1,region_2){
#  region_1<-T1_loci
  #region_2<-T2_loci
  buffer <- 125000  # ±125kb
  # 
  gr1 <- GRanges(seqnames = region_1$chr,
                 ranges = IRanges(start = region_1$start - buffer,
                                  end   = region_1$end + buffer),
                 locus_id = region_1$GenomicLocus)
  
  gr2 <- GRanges(seqnames = region_2$chr,
                 ranges = IRanges(start = region_2$start - buffer,
                                  end   = region_2$end + buffer),
                 locus_id = region_2$GenomicLocus)
 # unique_ranges <- !duplicated(paste0(start(gr2), "-", end(gr2)))
 # gr2 <- gr2[unique_ranges]
  # overlap regions 3
  hits <- findOverlaps(gr1, gr2)
  overlap_df <- data.frame(
      fuma1_locus = mcols(gr1)$locus_id[queryHits(hits)],
      fuma1_chr   = as.character(seqnames(gr1))[queryHits(hits)],
      fuma1_start = start(gr1)[queryHits(hits)],
      fuma1_end   = end(gr1)[queryHits(hits)],
      fuma2_locus = mcols(gr2)$locus_id[subjectHits(hits)],
      fuma2_chr   = as.character(seqnames(gr2))[subjectHits(hits)],
      fuma2_start = start(gr2)[subjectHits(hits)],
      fuma2_end   = end(gr2)[subjectHits(hits)]
  )
  return(overlap_df)
}

T1_T2_region<-get_hit_res(T1_loci,T2_loci)
T1_FA_region<-get_hit_res(T1_loci,FA_loci)
T2_FA_region<-get_hit_res(FA_loci,T2_loci)
###

##
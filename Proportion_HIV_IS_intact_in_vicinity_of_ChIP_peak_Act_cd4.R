# 1. input dataframe
IS_intact_tous_H3K27ac <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K27ac_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K27ac_H3K4me1 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K27ac_H3K4me1_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K36me3 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K36me3_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K4me3 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K4me3_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K4me1 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K4me1_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K27me3 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K27me3_Act_peaks.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K9me3 <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/Einkauf_et_al_all_intact_IS_H3K9me3_Act_peaks.txt", header = F, stringsAsFactors = F)

# Calculation the proportion of IS in teh vicinity of different ChIP-seq condition
vicinity_to_ChIP <- function(df) {
  quant_untreat_loin <- dim(df %>% dplyr::filter(V5 == "untreat" & V8 == -1))[1]/dim(df %>% dplyr::filter(V5 == "untreat"))[1]
  quant_untreat_proche <- dim(df %>% dplyr::filter(V5 == "untreat" & V8 != -1))[1]/dim(df %>% dplyr::filter(V5 == "untreat"))[1]
  
  quant_treat_loin <- dim(df %>% dplyr::filter(V5 == "treat" & V8 == -1))[1]/dim(df %>% dplyr::filter(V5 == "treat"))[1]
  quant_treat_proche <- dim(df %>% dplyr::filter(V5 == "treat" & V8 != -1))[1]/dim(df %>% dplyr::filter(V5 == "treat"))[1]
  
  count <- c(quant_untreat_loin, quant_untreat_proche, quant_treat_loin, quant_treat_proche)
  condition <- c(rep("untreat", 2), rep("treat", 2))
  distance <- rep(c("Far", "Vicinity"), 2)
  
  df_fin <- data.frame(count, condition, distance)
  
  deposer_condition <- c("untreat", "treat")
  
  df_fin$condition <- factor(df_fin$condition, levels = deposer_condition)
  
  return(df_fin)
}

## all IS
Vicinity_H3K27ac <- vicinity_to_ChIP(IS_intact_tous_H3K27ac) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1 <- vicinity_to_ChIP(IS_intact_tous_H3K27ac_H3K4me1) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3 <- vicinity_to_ChIP(IS_intact_tous_H3K36me3) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1 <- vicinity_to_ChIP(IS_intact_tous_H3K4me1) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3 <- vicinity_to_ChIP(IS_intact_tous_H3K4me3) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3 <- vicinity_to_ChIP(IS_intact_tous_H3K27me3) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3 <- vicinity_to_ChIP(IS_intact_tous_H3K9me3) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP <- dplyr::bind_rows(Vicinity_H3K27ac, Vicinity_H3K27ac_H3K4me1, Vicinity_H3K36me3, Vicinity_H3K4me1, Vicinity_H3K4me3, Vicinity_H3K27me3, Vicinity_H3K9me3)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

##IS inactive-repressed in the vicinity of ChIP-seq signal
Vicinity_H3K27ac_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_IR <- dplyr::bind_rows(Vicinity_H3K27ac_IR, Vicinity_H3K27ac_H3K4me1_IR, Vicinity_H3K36me3_IR, Vicinity_H3K4me1_IR, Vicinity_H3K4me3_IR, Vicinity_H3K27me3_IR, Vicinity_H3K9me3_IR)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_IR.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_IR, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Reds")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

##IS elongated-repressed in the vicinity of ChIP-seq signal
Vicinity_H3K27ac_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_ER <- dplyr::bind_rows(Vicinity_H3K27ac_ER, Vicinity_H3K27ac_H3K4me1_ER, Vicinity_H3K36me3_ER, Vicinity_H3K4me1_ER, Vicinity_H3K4me3_ER, Vicinity_H3K27me3_ER, Vicinity_H3K9me3_ER)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_ER.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_ER, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Greens")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

## IS Fully-active in the vicinity of ChIP-seq signal
Vicinity_H3K27ac_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_FA <- dplyr::bind_rows(Vicinity_H3K27ac_FA, Vicinity_H3K27ac_H3K4me1_FA, Vicinity_H3K36me3_FA, Vicinity_H3K4me1_FA, Vicinity_H3K4me3_FA, Vicinity_H3K27me3_FA, Vicinity_H3K9me3_FA)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_FA.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_FA, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Purples")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

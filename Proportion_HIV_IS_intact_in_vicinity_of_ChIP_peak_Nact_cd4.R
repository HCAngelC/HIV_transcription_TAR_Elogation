# 1. Input bedtools intersect -wao dataframe
IS_intact_tous_H3K27ac_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K27ac_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K27ac_H3K4me1_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K27ac_H3K4me1_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K36me3_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K36me3_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K4me3_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K4me3_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K4me1_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K4me1_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K27me3_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K27me3_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)
IS_intact_tous_H3K9me3_nact <- read.table("/media/chen/DATA/HIV_transcription_initiation/overlap_ChIP/Pool_Replicate_peaks_hg38/cd4_nact/Einkauf_et_al_2022_H3K9me3_Nact_broadPeak_hg38.txt", header = F, stringsAsFactors = F)

# 2. Proportion of IS inactive-repressed in the vicinity of ChIP-seq signals
Vicinity_H3K27ac_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_nact_IR <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3_nact, V4 == "inactive_repressed")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_nact_IR <- dplyr::bind_rows(Vicinity_H3K27ac_nact_IR, Vicinity_H3K27ac_H3K4me1_nact_IR, Vicinity_H3K36me3_nact_IR, Vicinity_H3K4me1_nact_IR, Vicinity_H3K4me3_nact_IR, Vicinity_H3K27me3_nact_IR, Vicinity_H3K9me3_nact_IR)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_nact_IR.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_nact_IR, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Reds")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

# 3. Proportion of IS elongated-repressed in the vicinity of ChIP-seq signals
Vicinity_H3K27ac_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_nact_ER <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3_nact, V4 == "elongated_repressed")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_nact_ER <- dplyr::bind_rows(Vicinity_H3K27ac_nact_ER, Vicinity_H3K27ac_H3K4me1_nact_ER, Vicinity_H3K36me3_nact_ER, Vicinity_H3K4me1_nact_ER, Vicinity_H3K4me3_nact_ER, Vicinity_H3K27me3_nact_ER, Vicinity_H3K9me3_nact_ER)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_nact_ER.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_nact_ER, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Greens")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

# 3. Proportion of IS fully-active in the vicinity of ChIP-seq signals
Vicinity_H3K27ac_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27ac")
Vicinity_H3K27ac_H3K4me1_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27ac_H3K4me1_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27ac & H3K4me1")
Vicinity_H3K36me3_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K36me3_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K36me3")
Vicinity_H3K4me1_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me1_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K4me1")
Vicinity_H3K4me3_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K4me3_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K4me3")
Vicinity_H3K27me3_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K27me3_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K27me3")
Vicinity_H3K9me3_nact_FA <- vicinity_to_ChIP(dplyr::filter(IS_intact_tous_H3K9me3_nact, V4 == "fully_active")) %>% dplyr::mutate(ChIP = "H3K9me3")

Vicinity_tous_ChIP_nact_FA <- dplyr::bind_rows(Vicinity_H3K27ac_nact_FA, Vicinity_H3K27ac_H3K4me1_nact_FA, Vicinity_H3K36me3_nact_FA, Vicinity_H3K4me1_nact_FA, Vicinity_H3K4me3_nact_FA, Vicinity_H3K27me3_nact_FA, Vicinity_H3K9me3_nact_FA)

pdf("/media/chen/DATA/HIV_transcription_initiation/abb/Vicinity_tous_ChIP_nact_FA.pdf", height = 5, width = 6)
ggplot(Vicinity_tous_ChIP_nact_FA, aes(x = condition, y = count, fill = distance))+geom_bar(position = "stack", stat = "identity", color = "black")+facet_grid(. ~ ChIP)+theme_bw()+scale_fill_brewer(type = "seq", palette = "Purples")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 4))
dev.off()

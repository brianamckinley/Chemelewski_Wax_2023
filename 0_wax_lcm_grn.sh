########################################################
############# User defined parameters ##################
########################################################
Directory=/Users/brianmckinley/desktop/sorghum/144_wax_grn/
GRN_Directory=/Users/brianmckinley/desktop/bioinformatics/7_GRN/0_GRN_required_files/
Description=wax1_LCM_grn
filename=${Directory}/1_${Description}_summing_input.csv
sql=${Directory}/${Description}_select.sql
filename1=${Directory}/1.1_${Description}_summing_input.csv
##range in R for TPM threshold selection after summing
TPMcutoff=5
RmaxStartTPMdata=2
RmaxEndTPMdata=29
RmaxStartCountData=30
RmaxEndCountdata=57
group='c(1,1,1,1,2,2,2,2, 3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7)'
contrast='c(0.5, 0.5, -0.2, -0.2, -0.2, -0.2, -0.2)'
FC=5
FDR=0.05
enrichment=0
########################################################
############# Begin GRN calculation script #############
########################################################
cat > ${sql} << EOF
.mode tabs
.output ${filename}
.headers ON
.separator ","
SELECT DISTINCT
GeneIDV3 as Gene_ID,
"Wray_V-phase_Epidermis_R1_TPM",
"Wray_V-phase_Epidermis_R2_TPM",
"Wray_V-phase_Epidermis_R3_TPM",
"Wray_V-phase_Epidermis_R4_TPM",
"Wray_V-phase_Sub_Epidermis_R1_TPM",
"Wray_V-phase_Sub_Epidermis_R2_TPM",
"Wray_V-phase_Sub_Epidermis_R3_TPM",
"Wray_V-phase_Sub_Epidermis_R4_TPM",
"Wray_V-phase_Bundle_Sheath_R1_TPM",
"Wray_V-phase_Bundle_Sheath_R2_TPM",
"Wray_V-phase_Bundle_Sheath_R3_TPM",
"Wray_V-phase_Bundle_Sheath_R4_TPM",
"Wray_V-phase_Phloem_R1_TPM",
"Wray_V-phase_Phloem_R2_TPM",
"Wray_V-phase_Phloem_R3_TPM",
"Wray_V-phase_Phloem_R4_TPM",
"Wray_V-phase_Pith_R1_TPM",
"Wray_V-phase_Pith_R2_TPM",
"Wray_V-phase_Pith_R3_TPM",
"Wray_V-phase_Pith_R4_TPM",
"Wray_V-phase_Xylem_Parenchyma_R1_TPM",
"Wray_V-phase_Xylem_Parenchyma_R2_TPM",
"Wray_V-phase_Xylem_Parenchyma_R3_TPM",
"Wray_V-phase_Xylem_Parenchyma_R4_TPM",
"Wray_V-phase_Xylem_R1_TPM",
"Wray_V-phase_Xylem_R2_TPM",
"Wray_V-phase_Xylem_R3_TPM",
"Wray_V-phase_Xylem_R4_TPM",

"Wray_V-phase_Epidermis_R1_counts",
"Wray_V-phase_Epidermis_R2_counts",
"Wray_V-phase_Epidermis_R3_counts",
"Wray_V-phase_Epidermis_R4_counts",
"Wray_V-phase_Sub_Epidermis_R1_counts",
"Wray_V-phase_Sub_Epidermis_R2_counts",
"Wray_V-phase_Sub_Epidermis_R3_counts",
"Wray_V-phase_Sub_Epidermis_R4_counts",
"Wray_V-phase_Bundle_Sheath_R1_counts",
"Wray_V-phase_Bundle_Sheath_R2_counts",
"Wray_V-phase_Bundle_Sheath_R3_counts",
"Wray_V-phase_Bundle_Sheath_R4_counts",
"Wray_V-phase_Phloem_R1_counts",
"Wray_V-phase_Phloem_R2_counts",
"Wray_V-phase_Phloem_R3_counts",
"Wray_V-phase_Phloem_R4_counts",
"Wray_V-phase_Pith_R1_counts",
"Wray_V-phase_Pith_R2_counts",
"Wray_V-phase_Pith_R3_counts",
"Wray_V-phase_Pith_R4_counts",
"Wray_V-phase_Xylem_Parenchyma_R1_counts",
"Wray_V-phase_Xylem_Parenchyma_R2_counts",
"Wray_V-phase_Xylem_Parenchyma_R3_counts",
"Wray_V-phase_Xylem_Parenchyma_R4_counts",
"Wray_V-phase_Xylem_R1_counts",
"Wray_V-phase_Xylem_R2_counts",
"Wray_V-phase_Xylem_R3_counts",
"Wray_V-phase_Xylem_R4_counts"

FROM
(SELECT * FROM annotation
LEFT JOIN geneorder
ON geneorder.TranscriptIDV3 = annotation.TranscriptIDV3
--LEFT JOIN AnnotationV3
--I ON AnnotationV3.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN cytokinins
ON cytokinins.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN GA
ON GA.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN BR
ON BR.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN suberin
ON suberin.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN sugar_metabolism
ON sugar_metabolism.GeneIDV3 = annotation.GeneIDV3
LEFT JOIN internode_growth
ON internode_growth.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN buds
ON buds.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN della_internode
ON della_internode.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN della_tissue_specificity
ON della_tissue_specificity.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN SASV3_TPM
ON SASV3_TPM.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN cc100M
ON cc100M.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN ccSM100
ON ccSM100.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN ccBTx623
ON ccBTx623.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN DellaDiurnalCycling
ON DellaDiurnalCycling.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN DellaDevV3
ON DellaDevV3.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN KellerDiurnalCycling
ON KellerDiurnalCycling.TRANScriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN R07020_AER
ON R07020_AER.TRANScriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN Dw2dw2
ON Dw2dw2.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN nroots
ON nroots.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN RootDepth
ON RootDepth.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN mapman4
ON mapman4.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN TX08001_radial
ON TX08001_radial.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN temp
ON temp.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN SYM
ON SYM.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN dym
ON dym.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN ddym2
ON ddym2.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN pnnl
ON pnnl.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN dleaves
ON dleaves.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN wray_stem_dev
ON wray_stem_dev.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN tx_sas
ON tx_sas.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN x58M_BSI
ON x58M_BSI.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN x58M_diurnal
ON x58M_diurnal.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN wray_LCM
ON wray_LCM.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN xR07007_diurnal
ON xR07007_diurnal.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN x100M_diurnal
ON x100M_diurnal.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN x100M_BSI
ON x100M_BSI.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN TX08001_density
ON TX08001_density.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN BTx623_atlas
ON BTx623_atlas.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN seedling_MLG
ON seedling_MLG.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN LS_removal
ON LS_removal.TranscriptIDV3 = annotation.TranscriptIDV3
LEFT JOIN peptide_signaling
ON peptide_signaling.GeneIDV3 = annotation.GeneIDV3
LEFT JOIN nodal_buds
ON nodal_buds.transcriptIDV3 = annotation.transcriptIDV3
LEFT JOIN roots_field_water_stress
ON roots_field_water_stress.transcriptIDV3 = annotation.transcriptIDV3
LEFT JOIN roots_GH_water_stress
ON roots_GH_water_stress.transcriptIDV3 = annotation.transcriptIDV3
)
WHERE
NOT GeneIDV3 = "Sobic.K006100" AND
NOT GeneIDV3 = "Sobic.K000200" AND
NOT GeneIDV3 = "Sobic.K000400" AND
NOT GeneIDV3 = "Sobic.K000500" AND
NOT GeneIDV3 = "Sobic.K000700" AND
NOT GeneIDV3 = "Sobic.K000800" AND
NOT GeneIDV3 = "Sobic.K001000" AND
NOT GeneIDV3 = "Sobic.K001100" AND
NOT GeneIDV3 = "Sobic.K001200" AND
NOT GeneIDV3 = "Sobic.K001301" AND
NOT GeneIDV3 = "Sobic.K001401" AND
NOT GeneIDV3 = "Sobic.K001500" AND
NOT GeneIDV3 = "Sobic.K001650" AND
NOT GeneIDV3 = "Sobic.K001800" AND
NOT GeneIDV3 = "Sobic.K001900" AND
NOT GeneIDV3 = "Sobic.K002100" AND
NOT GeneIDV3 = "Sobic.K002200" AND
NOT GeneIDV3 = "Sobic.K002300" AND
NOT GeneIDV3 = "Sobic.K002400" AND
NOT GeneIDV3 = "Sobic.K002501" AND
NOT GeneIDV3 = "Sobic.K002600" AND
NOT GeneIDV3 = "Sobic.K002700" AND
NOT GeneIDV3 = "Sobic.K002800" AND
NOT GeneIDV3 = "Sobic.K002850" AND
NOT GeneIDV3 = "Sobic.K002900" AND
NOT GeneIDV3 = "Sobic.K003000" AND
NOT GeneIDV3 = "Sobic.K003100" AND
NOT GeneIDV3 = "Sobic.K003250" AND
NOT GeneIDV3 = "Sobic.K003400" AND
NOT GeneIDV3 = "Sobic.K003500" AND
NOT GeneIDV3 = "Sobic.K003600" AND
NOT GeneIDV3 = "Sobic.K003700" AND
NOT GeneIDV3 = "Sobic.K003700" AND
NOT GeneIDV3 = "Sobic.K003700" AND
NOT GeneIDV3 = "Sobic.K003800" AND
NOT GeneIDV3 = "Sobic.K004100" AND
NOT GeneIDV3 = "Sobic.K004201" AND
NOT GeneIDV3 = "Sobic.K004400" AND
NOT GeneIDV3 = "Sobic.K004500" AND
NOT GeneIDV3 = "Sobic.K004800" AND
NOT GeneIDV3 = "Sobic.K004900" AND
NOT GeneIDV3 = "Sobic.K005000" AND
NOT GeneIDV3 = "Sobic.K005101" AND
NOT GeneIDV3 = "Sobic.K005200" AND
NOT GeneIDV3 = "Sobic.K005300" AND
NOT GeneIDV3 = "Sobic.K005400" AND
NOT GeneIDV3 = "Sobic.K005500" AND
NOT GeneIDV3 = "Sobic.K005600" AND
NOT GeneIDV3 = "Sobic.K005700" AND
NOT GeneIDV3 = "Sobic.K005800" AND
NOT GeneIDV3 = "Sobic.K005900" AND
NOT GeneIDV3 = "Sobic.K006001" AND
NOT GeneIDV3 = "Sobic.K006100" AND
NOT GeneIDV3 = "Sobic.K006300" AND
NOT GeneIDV3 = "Sobic.K006400" AND
NOT GeneIDV3 = "Sobic.K006500" AND
NOT GeneIDV3 = "Sobic.K006600" AND
NOT GeneIDV3 = "Sobic.K006800" AND
NOT GeneIDV3 = "Sobic.K007000" AND
NOT GeneIDV3 = "Sobic.K007100" AND
NOT GeneIDV3 = "Sobic.K007150" AND
NOT GeneIDV3 = "Sobic.K007200" AND
NOT GeneIDV3 = "Sobic.K007500" AND
NOT GeneIDV3 = "Sobic.K007600" AND
NOT GeneIDV3 = "Sobic.K007900" AND
NOT GeneIDV3 = "Sobic.K008000" AND
NOT GeneIDV3 = "Sobic.K008300" AND
NOT GeneIDV3 = "Sobic.K008400" AND
NOT GeneIDV3 = "Sobic.K008450" AND
NOT GeneIDV3 = "Sobic.K008500" AND
NOT GeneIDV3 = "Sobic.K009100" AND
NOT GeneIDV3 = "Sobic.K009300" AND
NOT GeneIDV3 = "Sobic.K009500" AND
NOT GeneIDV3 = "Sobic.K009600" AND
NOT GeneIDV3 = "Sobic.K009800" AND
NOT GeneIDV3 = "Sobic.K009800" AND
NOT GeneIDV3 = "Sobic.K009900" AND
NOT GeneIDV3 = "Sobic.K010000" AND
NOT GeneIDV3 = "Sobic.K010100" AND
NOT GeneIDV3 = "Sobic.K010300" AND
NOT GeneIDV3 = "Sobic.K010500" AND
NOT GeneIDV3 = "Sobic.K012900" AND
NOT GeneIDV3 = "Sobic.K013000" AND
NOT GeneIDV3 = "Sobic.K013100" AND
NOT GeneIDV3 = "Sobic.K013200" AND
NOT GeneIDV3 = "Sobic.K014100" AND
NOT GeneIDV3 = "Sobic.K014200" AND
NOT GeneIDV3 = "Sobic.K015700" AND
NOT GeneIDV3 = "Sobic.K015800" AND
NOT GeneIDV3 = "Sobic.K016500" AND
NOT GeneIDV3 = "Sobic.K016700" AND
NOT GeneIDV3 = "Sobic.K016900" AND
NOT GeneIDV3 = "Sobic.K016900" AND
NOT GeneIDV3 = "Sobic.K017300" AND
NOT GeneIDV3 = "Sobic.K020400" AND
NOT GeneIDV3 = "Sobic.K020500" AND
NOT GeneIDV3 = "Sobic.K020800" AND
NOT GeneIDV3 = "Sobic.K021000" AND
NOT GeneIDV3 = "Sobic.K021301" AND
NOT GeneIDV3 = "Sobic.K021600" AND
NOT GeneIDV3 = "Sobic.K022000" AND
NOT GeneIDV3 = "Sobic.K022100" AND
NOT GeneIDV3 = "Sobic.K022200" AND
NOT GeneIDV3 = "Sobic.K022300" AND
NOT GeneIDV3 = "Sobic.K022500" AND
NOT GeneIDV3 = "Sobic.K022600" AND
NOT GeneIDV3 = "Sobic.K022800" AND
NOT GeneIDV3 = "Sobic.K022900" AND
NOT GeneIDV3 = "Sobic.K023150" AND
NOT GeneIDV3 = "Sobic.K023400" AND
NOT GeneIDV3 = "Sobic.K023450" AND
NOT GeneIDV3 = "Sobic.K023500" AND
NOT GeneIDV3 = "Sobic.K023600" AND
NOT GeneIDV3 = "Sobic.K023800" AND
NOT GeneIDV3 = "Sobic.K024000" AND
NOT GeneIDV3 = "Sobic.K024100" AND
NOT GeneIDV3 = "Sobic.K024300" AND
NOT GeneIDV3 = "Sobic.K024500" AND
NOT GeneIDV3 = "Sobic.K024600" AND
NOT GeneIDV3 = "Sobic.K024800" AND
NOT GeneIDV3 = "Sobic.K025800" AND
NOT GeneIDV3 = "Sobic.K025900" AND
NOT GeneIDV3 = "Sobic.K026200" AND
NOT GeneIDV3 = "Sobic.K026400" AND
NOT GeneIDV3 = "Sobic.K026800" AND
NOT GeneIDV3 = "Sobic.K026900" AND
NOT GeneIDV3 = "Sobic.K026900" AND
NOT GeneIDV3 = "Sobic.K027700" AND
NOT GeneIDV3 = "Sobic.K027900" AND
NOT GeneIDV3 = "Sobic.K028000" AND
NOT GeneIDV3 = "Sobic.K028000" AND
NOT GeneIDV3 = "Sobic.K028000" AND
NOT GeneIDV3 = "Sobic.K028200" AND
NOT GeneIDV3 = "Sobic.K028300" AND
NOT GeneIDV3 = "Sobic.K028600" AND
NOT GeneIDV3 = "Sobic.K028900" AND
NOT GeneIDV3 = "Sobic.K029400" AND
NOT GeneIDV3 = "Sobic.K029900" AND
NOT GeneIDV3 = "Sobic.K030000" AND
NOT GeneIDV3 = "Sobic.K030100" AND
NOT GeneIDV3 = "Sobic.K030200" AND
NOT GeneIDV3 = "Sobic.K030700" AND
NOT GeneIDV3 = "Sobic.K030700" AND
NOT GeneIDV3 = "Sobic.K031100" AND
NOT GeneIDV3 = "Sobic.K031200" AND
NOT GeneIDV3 = "Sobic.K031300" AND
NOT GeneIDV3 = "Sobic.K031400" AND
NOT GeneIDV3 = "Sobic.K031500" AND
NOT GeneIDV3 = "Sobic.K031500" AND
NOT GeneIDV3 = "Sobic.K031800" AND
NOT GeneIDV3 = "Sobic.K031800" AND
NOT GeneIDV3 = "Sobic.K031800" AND
NOT GeneIDV3 = "Sobic.K031800" AND
NOT GeneIDV3 = "Sobic.K031800" AND
NOT GeneIDV3 = "Sobic.K031900" AND
NOT GeneIDV3 = "Sobic.K032200" AND
NOT GeneIDV3 = "Sobic.K032300" AND
NOT GeneIDV3 = "Sobic.K032350" AND
NOT GeneIDV3 = "Sobic.K032400" AND
NOT GeneIDV3 = "Sobic.K032500" AND
NOT GeneIDV3 = "Sobic.K032600" AND
NOT GeneIDV3 = "Sobic.K032900" AND
NOT GeneIDV3 = "Sobic.K033000" AND
NOT GeneIDV3 = "Sobic.K040200" AND
NOT GeneIDV3 = "Sobic.K040800" AND
NOT GeneIDV3 = "Sobic.K041900" AND
NOT GeneIDV3 = "Sobic.K042000" AND
NOT GeneIDV3 = "Sobic.K042600" AND
NOT GeneIDV3 = "Sobic.K042700" AND
NOT GeneIDV3 = "Sobic.K043400" AND
NOT GeneIDV3 = "Sobic.K043500" AND
NOT GeneIDV3 = "Sobic.K043600" AND
NOT GeneIDV3 = "Sobic.K043700" AND
NOT GeneIDV3 = "Sobic.K043800" AND
NOT GeneIDV3 = "Sobic.K044200" AND
NOT GeneIDV3 = "Sobic.K044300" AND
NOT GeneIDV3 = "Sobic.K044400" AND
NOT GeneIDV3 = "Sobic.K044401" AND
NOT GeneIDV3 = "Sobic.K044402" AND
NOT GeneIDV3 = "Sobic.K044403" AND
NOT GeneIDV3 = "Sobic.K044404" AND
NOT GeneIDV3 = "Sobic.K044405" AND
NOT GeneIDV3 = "Sobic.K044406" AND
NOT GeneIDV3 = "Sobic.K044406" AND
NOT GeneIDV3 = "Sobic.K044407" AND
NOT GeneIDV3 = "Sobic.K044407" AND
NOT GeneIDV3 = "Sobic.K044407" AND
NOT GeneIDV3 = "Sobic.K044408" AND
NOT GeneIDV3 = "Sobic.K044409" AND
NOT GeneIDV3 = "Sobic.K044410" AND
NOT GeneIDV3 = "Sobic.K044411" AND
NOT GeneIDV3 = "Sobic.K044412" AND
NOT GeneIDV3 = "Sobic.K044413" AND
NOT GeneIDV3 = "Sobic.K044414" AND
NOT GeneIDV3 = "Sobic.K044415" AND
NOT GeneIDV3 = "Sobic.K044416" AND
NOT GeneIDV3 = "Sobic.K044417" AND
NOT GeneIDV3 = "Sobic.K044418" AND
NOT GeneIDV3 = "Sobic.K044419" AND
NOT GeneIDV3 = "Sobic.K044420" AND
NOT GeneIDV3 = "Sobic.K044505"
;
.exit
EOF
chmod u+x ${sql};
/Users/brianmckinley/desktop/sorghum/_01_SQLite3/sqlite-autoconf-3250000/SQLite3 /Users/brianmckinley/desktop/sorghum/Sorghum.db < ${sql};
rm ${sql};

#gcut --complement -d',' -f1 ${filename} > ${filename1};

cat > sum.R << EOF
library("dplyr")
library("matrixStats")

############# Sum transcripts to genes
transcript_tpm <- read.csv("${filename}", header = T, as.is = T)
ncol(transcript_tpm)
gene_tpm <- group_by(transcript_tpm, Gene_ID) %>% summarise(across(everything(), sum))


############# Apply TPM threshold
max <- rowMaxs(as.matrix(gene_tpm[,${RmaxStartTPMdata}:${RmaxEndTPMdata}]))
gene_tpm_max <- cbind(gene_tpm,max)
gene_tpm5up <- subset(gene_tpm_max,subset = max > ${TPMcutoff})
Gene_ID <- gene_tpm5up[,'Gene_ID']

############# Create TPM and counts datasets
${Description}_TPM_tpm5up <- cbind(Gene_ID,gene_tpm5up[,${RmaxStartTPMdata}:${RmaxEndTPMdata}])
${Description}_counts_tpm5up <- cbind(Gene_ID,gene_tpm5up[,${RmaxStartCountData}:${RmaxEndCountdata}])

############# Calculate differential expression
library("edgeR")
library("gtools")
group <- factor(${group})
DGEobject1 <- DGEList(counts=${Description}_counts_tpm5up[2:ncol(${Description}_counts_tpm5up)], genes=${Description}_counts_tpm5up[, 1], group=group)
levels(DGEobject1\$samples\$group)
keep <- rowSums(cpm(DGEobject1) > 0) > 0
DGEobject2 <- DGEobject1[keep, ,keep.lib.sizes=FALSE]
DGEobject3 <- calcNormFactors(DGEobject2, norm.method = "tmm", test.method = "edger", na.rm = TRUE)
design <- model.matrix(~0+group)
DGEobject4 <- estimateDisp(DGEobject3, design, robust=TRUE, tagwise=TRUE)
DGEobject5 <- glmQLFit(DGEobject4, design, winsor.tail.p=c(0.05, 0.1), robust=TRUE)
DGEobject6 <- glmQLFTest(DGEobject5, contrast = ${contrast})
DGEobject6\$table\$FC <- logratio2foldchange(DGEobject6\$table\$logFC, base=2)
topTags <- topTags(DGEobject6, n=47206, sort.by="none")
############# Assemble GRN input dataset of DE genes FC > threshold specified & FDR < threshold specified
FC <- as.matrix(topTags\$table\$FC)
FDR <- as.matrix(topTags\$table\$FDR)
TPMs <- ${Description}_TPM_tpm5up
combinedDGEobject <- cbind(TPMs,FC,FDR)
names(combinedDGEobject)[1] <- "Gene_ID"
write.csv(combinedDGEobject, "combinedDGEobject.csv", row.names = F)
PCCinput <- subset(combinedDGEobject, subset = FC > ${FC} & FDR < ${FDR})
write.csv(PCCinput, "DEgenes.csv", row.names = F)
print("Genes in grn input dataset with FC > ${FC} & FDR < ${FDR}")

############# Remove last two columns that were used for subsetting by FC and FDR
counts_ncol <- ncol(PCCinput)-2
PCCinput <- PCCinput[,1:counts_ncol]
nrow(PCCinput)
write.csv(PCCinput, "PCCinput.csv", row.names = F)

############# Calculate PCC matrix
cor.test.mk.R <- function(expvalues, cormethod = "pearson", corcutoff = 0.1, evalcutoff = 0.05) {

  correlation <- matrix(0, nrow = 1, ncol = 4)

  iminus1 <- nrow(expvalues) - 1
  corvalues <- numeric()

  for (i in 1:iminus1) {
    print(i)
    inext <- i + 1
    if (inext == nrow(expvalues)) {
      tmp.corvalues <- mycortest(expvalues[inext,], expvalues[i,], cormethod = cormethod)
      tmp.corvalues <- cbind(rownames(expvalues)[inext],
                             c(rownames(expvalues)[i]), t(tmp.corvalues))
    } else {
      tmp.corvalues <- apply(expvalues[inext:nrow(expvalues),], 1, mycortest,
                             expvalues[i, ], cormethod)
      tmp.corvalues <- cbind(rownames(expvalues)[inext:nrow(expvalues)],
                             c(rownames(expvalues)[i]), t(tmp.corvalues))
    }

    significantcorvalues <- which(as.numeric(tmp.corvalues[, 4]) <= evalcutoff)
    if (length(significantcorvalues) >= 1) {
      tmp.corvalues <- tmp.corvalues[significantcorvalues,]
      if (class(tmp.corvalues) != "matrix") {
        tmp.corvalues <- matrix(tmp.corvalues, ncol = 4, byrow = T)
      }
      poscorvalues <- which(as.numeric(tmp.corvalues[, 3]) >= corcutoff)
      if (length(poscorvalues) >= 1) {
        tmp.corvaluesup <- tmp.corvalues[poscorvalues,]
        if (class(tmp.corvaluesup) != "matrix") {
          tmp.corvaluesup <- matrix(tmp.corvaluesup, ncol = 4, byrow = T)
        }
        correlation <- rbind(correlation, tmp.corvaluesup)
        #write.table(tmp.corvaluesup, file = fileout, quote = F, row.names = F, append = T, col.names = F, sep = ",")
      }

      negcorvalues <- which(as.numeric(tmp.corvalues[, 3]) <= (-1 * corcutoff))
      if (length(negcorvalues) >= 1) {
        tmp.corvaluesdown <- tmp.corvalues[negcorvalues,]
        if (class(tmp.corvaluesdown) != "matrix") {
          tmp.corvaluesdown <- matrix(tmp.corvaluesdown, ncol = 4, byrow = T)
        }
        correlation <- rbind(correlation, tmp.corvaluesdown)
        #write.table(tmp.corvaluesdown, file = fileout, quote = F, row.names = F, append = T, col.names = F, sep = ",")
      }
    }
  }
  correlation <- correlation[-1, ]
  rownames(correlation) <- c()

  return(correlation)
}

#############################

mycortest<-function(vector1, vector2, cormethod="pearson") {

  vector1=as.numeric(vector1)
  vector2=as.numeric(vector2)
  tmp.cor=corr.test(vector1, vector2, method=cormethod, adjust="fdr")
  tmp.result=c(tmp.cor\$r, tmp.cor\$p)
  as.numeric(tmp.result)

}

#########################################################################################
library(psych)

print("Pairwise PCC calculation started")
for (i in 1:2){
  print(i)
  # cpm values of all genes for allcondirions
  Data <- PCCinput
  row.names(Data) <- Data[, 1]
  MRinput <- cor.test.mk.R(Data)
}
write.table(MRinput,"PCC.cor",quote=F,sep=",", row.names = FALSE)

print("Pairwise PCC calculation finished")


############ Calculate MR ###############################################################
mutual_Rank <- function (cor_data) {

  library(data.table)
  if (class(cor_data) == "character") {
    cor_data <- fread(cor_data, header = FALSE)
  }
  cor_data <- setDT(cor_data)

  setkey(cor_data)
  cor_data <- unique(cor_data)
  cor_data <- cor_data[!duplicated(cor_data[, c("V1", "V2")])]

  Neg_cor <- cor_data[cor_data\$V3 < 0, ]
  setkey(Neg_cor)
  Neg_cor <- unique(Neg_cor)

  Neg_cor_R1 <- data.frame(Neg_cor\$V1)
  Neg_cor_R1 <- data.frame(Neg_cor_R1[!duplicated(Neg_cor_R1), ])
  colnames(Neg_cor_R1)[1] <- "Genes"

  Neg_cor_R2 <- data.frame(Neg_cor\$V2)
  Neg_cor_R2 <- data.frame(Neg_cor_R2[!duplicated(Neg_cor_R2), ])
  colnames(Neg_cor_R2)[1] <- "Genes"

  Neg_cor_rank <- data.frame()

  for (i in Neg_cor_R1\$Genes)
  {
    Samp1 <- Neg_cor[Neg_cor\$V1 == i]
    Samp1 <- Samp1[order(V3), ]
    Samp1\$PCC_rank_G2toG1 <- seq.int(nrow(Samp1))
    Neg_cor_rank <- rbind(Neg_cor_rank, Samp1)
  }

  Neg_cor_rank_new <- data.frame()
  for (i in Neg_cor_R2\$Genes)
  {
    Samp1 <- Neg_cor_rank[Neg_cor_rank\$V2 == i]
    Samp1 <- Samp1[order(V3), ]
    Samp1\$PCC_rank_G1toG2 <- seq.int(nrow(Samp1))
    Neg_cor_rank_new <- rbind(Neg_cor_rank_new, Samp1)
  }

  Pos_cor <- cor_data[cor_data\$V3 >= 0, ]
  setkey(Pos_cor)
  Pos_cor <- unique(Pos_cor)

  Pos_cor_R1 <- data.frame(Pos_cor\$V1)
  Pos_cor_R1 <- data.frame(Pos_cor_R1[!duplicated(Pos_cor_R1), ])
  colnames(Pos_cor_R1)[1] <- "Genes"

  Pos_cor_R2 <- data.frame(Pos_cor\$V2)
  Pos_cor_R2 <- data.frame(Pos_cor_R2[!duplicated(Pos_cor_R2), ])
  colnames(Pos_cor_R2)[1] <- "Genes"

  Pos_cor_rank <- data.frame()

  for (i in Pos_cor_R1\$Genes)
  {
    Samp1 <- Pos_cor[Pos_cor\$V1 == i]
    Samp1 <- Samp1[order(-V3), ]
    Samp1\$PCC_rank_G2toG1 <- seq.int(nrow(Samp1))
    Pos_cor_rank <- rbind(Pos_cor_rank, Samp1)
  }

  Pos_cor_rank_new <- data.frame()
  for (i in Pos_cor_R2\$Genes)
  {
    Samp1 <- Pos_cor_rank[Pos_cor_rank\$V2 == i]
    Samp1 <- Samp1[order(-V3), ]
    Samp1\$PCC_rank_G1toG2 <- seq.int(nrow(Samp1))
    Pos_cor_rank_new <- rbind(Pos_cor_rank_new, Samp1)
  }

  CorNet_ranks <- rbind(Neg_cor_rank_new, Pos_cor_rank_new)
  CorNet_ranks\$Mutual_rank <- (CorNet_ranks\$PCC_rank_G2toG1 *
                                 CorNet_ranks\$PCC_rank_G1toG2) ^ (1 / 2)

  CorNet_ranks <- data.frame(CorNet_ranks)
  CorNet_ranks <- CorNet_ranks[, c("V1", "V2", "Mutual_rank", "V3", "V4",
                                   "PCC_rank_G2toG1", "PCC_rank_G1toG2")]
  CorNet_ranks\$New <- 0
  for (i in 1:nrow(CorNet_ranks)){
    CorNet_ranks[i,8] <- (1-(CorNet_ranks[i,3]/(max(CorNet_ranks[,3])-min(CorNet_ranks[,3]))))
  }

  return(CorNet_ranks)
}
############################################################

require(data.table)
print("Mutual Rank calculation starting")
for (x in 1:2){
  input <- as.data.table(MRinput)
  MR_Network <- mutual_Rank(input)
  write.table(MR_Network,"MRRank.cor",quote=F,sep=",", row.names = FALSE)
}
print("Mutual Rank calculation finished")
############################################################
############################################################

enrich_TFBS_Family <- function(SelValsPos, Gene, AveragesFile, FoldChange){

  Enrich <- data.frame(matrix(ncol = 3))
  colnames(Enrich) <- c("Gene_ID", "TFBS_ID", "FC")

  SelValsPos_Gene <- subset(SelValsPos, SelValsPos\$Sequence.ID == Gene)
  TF <- count(SelValsPos_Gene, SelValsPos_Gene\$TFBS.ID)
  colnames(TF) <- c("TFBS_ID", "Count")

  for (y in 1:nrow(TF)){
    FC <- (TF[y,2])/(AveragesFile[AveragesFile\$TFBS_ID == TF[y,1],][,2])
    Temp_n <- cbind(as.character(Gene),as.character(TF[y,1]),as.numeric(FC))
    colnames(Temp_n) <- c("Gene_ID", "TFBS_ID", "FC")
    Enrich <- rbind(Enrich, Temp_n)
  }

  t <-subset(Enrich, Enrich\$FC > FoldChange)
  return(t)

}

unloadNamespace("dplyr");
library("dplyr");
library("reshape2");
library("pheatmap");
library("RColorBrewer");
library("qvalue");


# Provide Parameters
SelValsPos <- read.table("${GRN_Directory}6.1_PlantPAN1KBGenesTFIDs_bothStrand_PlantTFDBFamilyUpdated.txt",header= TRUE, sep="\t")
AveragesFile <- read.table("${GRN_Directory}6.2_AvgFrequencyTFOccurencePlantPan1KbUpstream_bothStrand.txt", header= TRUE, sep="\t")
TFBS_Family <- read.csv("${GRN_Directory}6.3_TFBS_Family.csv", header = T, as.is = T)
TFBS_Family <- TFBS_Family[!duplicated(TFBS_Family\$Matrix_ID),]
AllDEG <- read.csv("${Directory}PCCinput.csv", header = T, as.is = T)
GeneList <- unique(AllDEG\$Gene_ID)

print("Motif enrichment starting; FC = ${enrichment}")
# Enriched CREs for all genes in the genome based on FC > 1.1
enrich_CRE <- data.frame(matrix(data = NA, ncol = 3))
colnames(enrich_CRE) <- c("Gene_ID", "TFBS_ID", "FC")
for (i in 1:length(GeneList)){
  print(i)
  t <- enrich_TFBS_Family(SelValsPos, GeneList[i], AveragesFile, ${enrichment})
  enrich_CRE <- rbind(enrich_CRE, t)
}

# CRE TF Family
enrich_CRE_Family <- merge(enrich_CRE, TFBS_Family, by.x = "TFBS_ID", by.y = "Matrix_ID")
# Write output
write.csv(enrich_CRE_Family, "${Directory}Enrichment.csv", row.names = F)
print("Motif enrichment finished")
############################################################
############################################################
############################################################
############################################################

print("TF-promoter binding calculation started")
# TFBS of all genes
input1 <- "${GRN_Directory}6.1_PlantPAN1KBGenesTFIDs_bothStrand_PlantTFDBFamilyUpdated.txt"
input2 <- "${Directory}Enrichment.csv"
input3 <- "${GRN_Directory}7.2_Sbi_TF_all.txt"
input4 <- "${Directory}PCCinput.csv"
output <- "${Directory}8.2_All_DEG_reg_intr_CRE_Enrichment_FC0.txt"

Motif <- read.table(input1, header=T, sep="\t")
Enriched_CRE <- read.csv(input2, header = T, as.is = T)
TF_Family <- read.table(input3, header = T, sep= "\t", as.is = T)
DEG <- read.csv(input4, header = T, as.is = T)

uniq <- unique(DEG\$Gene_ID)
network <- data.frame()

for (i in 1:length(uniq)){
  print(i)
  Enriched_CRE_sub <- subset(Enriched_CRE, Enriched_CRE\$Gene_ID == uniq[i])
  Motif_sub <- subset(Motif, Motif\$Sequence.ID == uniq[i] & Motif\$Family %in% Enriched_CRE_sub\$Family)
  TF <- unique(Motif_sub[Motif_sub\$Sequence.ID == uniq[i],]\$Family)
  tt <- TF_Family[TF_Family\$Family %in% TF,]
  df <- cbind(rep(uniq[i], times= nrow(tt)), tt)
  network <- rbind(network,df)
}

final_reg <- as.data.frame(cbind(as.character(network[,2]), rep("interacts", times= nrow(network)), as.character(network[,1])))
write.table(final_reg, output, sep="\t", quote=F, row.names = F)
print("TF-promoter binding finished")

############################################################
############################################################
############################################################
############################################################


# Correlation files
input1 <- data.frame()
input1 <- "${Directory}MRRank.cor"

# regulatory interaction files
input2 <- "${Directory}8.2_All_DEG_reg_intr_CRE_Enrichment_FC0.txt"

# Output files
output1 <- data.frame()
output1 <- gsub("Coexpression", "GRN", input1, fixed=TRUE)
output1 <- gsub(".cor", ".txt", output1, fixed=TRUE)
output1 <- gsub("5.2_Corr_MRRank", "${Description}", output1, fixed=TRUE)

final_reg <- read.csv(file = input2, header = T, sep= "\t", as.is = T)

##################################
print("GRN assembly starting")

for (i in 1:length(input1)){

  Cor1 <- read.table(input1[i], header = T, sep= ",", as.is = T)
  colnames(Cor1) <- c("Gene1", "Gene2", "MR", "PCC", "PValue", "PCC1", "PCC2", "MR_Scale")
  Cor_TS_All <- subset(Cor1, Cor1\$PCC >= 0.0 & Cor1\$PValue <= 0.05 & Cor1\$MR_Scale >= 0.0)

  colnames(final_reg) <- c("Gene1", "rg", "Gene2")
  temp1 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene1" = "Gene1", "Gene2" = "Gene2"))
  temp2 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene2" = "Gene1", "Gene1" = "Gene2"))
  Cor_TS_Final <- rbind(temp1,temp2)
  Cor_TS_Final <- Cor_TS_Final[!(duplicated(Cor_TS_Final)| duplicated(Cor_TS_Final, fromLast=TRUE)),]
  Cor_TS_Final <- Cor_TS_Final[which(Cor_TS_Final\$Gene1 != Cor_TS_Final\$Gene2),]
  Cor_TS_Final <- Cor_TS_Final[complete.cases(Cor_TS_Final), ]

  Cor_TS_Final\$rg <- "rg"
  Cor_TS_Final\$rg <- paste(Cor_TS_Final\$rg, round(as.numeric(Cor_TS_Final\$Corr),digits = 1), sep = "")
  Cor_TS_Final <- Cor_TS_Final[,c(1,2,3,5,6,9)]

  write.table(Cor_TS_Final, output1[i], sep="\t", quote=F, row.names = F)
}
write.table(Cor_TS_Final, "1_${Description}_${FC}_${TPMcutoff}.txt", sep="\t", quote=F, row.names = F)
print("GRN assembly finished")

temp3 <- t(unique(Cor_TS_Final\$Gene1))
temp4 <- t(unique(Cor_TS_Final\$Gene2))
GeneList <- t(cbind(temp3,temp4))
GeneList <- unique(GeneList)
write.csv(GeneList, "103_Genelist_for_selection.csv",  row.names=FALSE)

EOF

Rscript sum.R;
#rm sum.R;
#rm ${filename}
#rm ${filename1}
#rm PCCinput.csv
#rm MRRank.cor
#rm MRRank.txt
#rm Enrichment.csv
#rm 8.2_All_DEG_reg_intr_CRE_Enrichment_FC0.txt

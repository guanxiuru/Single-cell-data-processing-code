library(SeuratObject)
library(sp)
library(Seurat)
library(ggraph)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(Matrix)
library(usethis)
library(devtools)
library(Rcpp)
library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)






rm(list = ls())#清除所有环境变量

setwd('D:/D盘/实验数据/R/downloaded_packages/新单细胞')
# MTX文件读取 ---
h5_file_scc1 <- 'GSM5342787_scA_filtered_peak_bc_matrix.h5'
h5_file_scc2 <- 'GSM5342788_scB_filtered_peak_bc_matrix.h5'
h5_file_scc3 <- 'GSM5342789_scC_filtered_peak_bc_matrix.h5'
h5_file_scc4 <- 'GSM5342790_scD_filtered_peak_bc_matrix.h5'
h5_file_scc5 <- 'GSM5342791_scE_filtered_peak_bc_matrix.h5'
h5_file_scc6 <- 'GSM5342792_scF_filtered_peak_bc_matrix.h5'
h5_file_scc7 <- 'GSM5342793_scG_filtered_peak_bc_matrix.h5'
h5_file_scc8 <- 'GSM5342794_scH_filtered_peak_bc_matrix.h5'
h5_file_scc9 <- 'GSM5342795_scI_filtered_peak_bc_matrix.h5'
h5_file_scc10 <- 'GSM5342796_scJ_filtered_peak_bc_matrix.h5'
h5_file_scc11 <- 'GSM5342797_scK_filtered_peak_bc_matrix.h5'
h5_file_scc12 <- 'GSM5342798_scL_filtered_peak_bc_matrix.h5'
h5_file_scc13 <- 'GSM5342799_scM_filtered_peak_bc_matrix.h5'
h5_file_scc14 <- 'GSM5342800_scN_filtered_peak_bc_matrix.h5'
h5_file_scc15 <- 'GSM5342801_scO_filtered_peak_bc_matrix.h5'
h5_file_scc16 <- 'GSM5342802_scP_filtered_peak_bc_matrix.h5'
h5_file_scc17 <- 'GSM5342803_scQ_filtered_peak_bc_matrix.h5'
h5_file_scc18 <- 'GSM5342804_scR_filtered_peak_bc_matrix.h5'
h5_file_scc19 <- 'GSM5342805_scS_filtered_peak_bc_matrix.h5'
h5_file_scc20 <- 'GSM5342806_scT_filtered_peak_bc_matrix.h5'
h5_file_scc21 <- 'GSM5342807_scU_filtered_peak_bc_matrix.h5'
h5_file_scc22 <- 'GSM5342808_scV_filtered_peak_bc_matrix.h5'
h5_file_scc23 <- 'GSM5342809_scW_filtered_peak_bc_matrix.h5'
h5_file_scc24 <- 'GSM5342810_scX_filtered_peak_bc_matrix.h5'
h5_file_scc25 <- 'GSM5342811_scY_filtered_peak_bc_matrix.h5'
h5_file_scc26 <- 'GSM5342812_scZ_filtered_peak_bc_matrix.h5'
h5_file_scc27 <- 'GSM5342813_scAA_filtered_peak_bc_matrix.h5'
h5_file_scc28 <- 'GSM5342814_scAB_filtered_peak_bc_matrix.h5'
h5_file_scc29 <- 'GSM5342815_scAC_filtered_peak_bc_matrix.h5'
h5_file_scc30 <- 'GSM5342816_scAD_filtered_peak_bc_matrix.h5'
h5_file_scc31 <- 'GSM5342817_scAE_filtered_peak_bc_matrix.h5'
h5_file_scc32 <- 'GSM5342818_scAF_filtered_peak_bc_matrix.h5'
h5_file_scc33 <- 'GSM5342819_scAG_filtered_peak_bc_matrix.h5'
h5_file_scc34 <- 'GSM5342820_scAH_filtered_peak_bc_matrix.h5'
h5_file_scc35 <- 'GSM5342821_scAI_filtered_peak_bc_matrix.h5'
h5_file_scc36 <- 'GSM5342822_scAJ_filtered_peak_bc_matrix.h5'
h5_file_scc37 <- 'GSM5342823_scAK_filtered_peak_bc_matrix.h5'
h5_file_scc38 <- 'GSM5342824_scAL_filtered_peak_bc_matrix.h5'
h5_file_scc39 <- 'GSM5342825_scAM_filtered_peak_bc_matrix.h5'
h5_file_scc40 <- 'GSM5342826_scAN_filtered_peak_bc_matrix.h5'
h5_file_scc41 <- 'GSM5342827_scAO_filtered_peak_bc_matrix.h5'
h5_file_scc42 <- 'GSM5342828_scAP_filtered_peak_bc_matrix.h5'
h5_file_scc43 <- 'GSM5342829_scAQ_filtered_peak_bc_matrix.h5'
h5_file_scc44 <- 'GSM5342830_scAR_filtered_peak_bc_matrix.h5'

data_scc1 <- Read10X_h5(file=h5_file_scc1)
data_scc2 <- Read10X_h5(file=h5_file_scc2)
data_scc3 <- Read10X_h5(file=h5_file_scc3)
data_scc4 <- Read10X_h5(file=h5_file_scc4)
data_scc5 <- Read10X_h5(file=h5_file_scc5)
data_scc6 <- Read10X_h5(file=h5_file_scc6)
data_scc7 <- Read10X_h5(file=h5_file_scc7)
data_scc8 <- Read10X_h5(file=h5_file_scc8)
data_scc9 <- Read10X_h5(file=h5_file_scc9)
data_scc10 <- Read10X_h5(file=h5_file_scc10)
data_scc11 <- Read10X_h5(file=h5_file_scc11)
data_scc12 <- Read10X_h5(file=h5_file_scc12)
data_scc13 <- Read10X_h5(file=h5_file_scc13)
data_scc14 <- Read10X_h5(file=h5_file_scc14)
data_scc15 <- Read10X_h5(file=h5_file_scc15)
data_scc16 <- Read10X_h5(file=h5_file_scc16)
data_scc17 <- Read10X_h5(file=h5_file_scc17)
data_scc18 <- Read10X_h5(file=h5_file_scc18)
data_scc19 <- Read10X_h5(file=h5_file_scc19)
data_scc20 <- Read10X_h5(file=h5_file_scc20)
data_scc21 <- Read10X_h5(file=h5_file_scc21)
data_scc22 <- Read10X_h5(file=h5_file_scc22)
data_scc23 <- Read10X_h5(file=h5_file_scc23)
data_scc24 <- Read10X_h5(file=h5_file_scc24)
data_scc25 <- Read10X_h5(file=h5_file_scc25)
data_scc26 <- Read10X_h5(file=h5_file_scc26)
data_scc27 <- Read10X_h5(file=h5_file_scc27)
data_scc28 <- Read10X_h5(file=h5_file_scc28)
data_scc29 <- Read10X_h5(file=h5_file_scc29)
data_scc30 <- Read10X_h5(file=h5_file_scc30)
data_scc31 <- Read10X_h5(file=h5_file_scc31)
data_scc32 <- Read10X_h5(file=h5_file_scc32)
data_scc33 <- Read10X_h5(file=h5_file_scc33)
data_scc34 <- Read10X_h5(file=h5_file_scc34)
data_scc35 <- Read10X_h5(file=h5_file_scc35)
data_scc36 <- Read10X_h5(file=h5_file_scc36)
data_scc37 <- Read10X_h5(file=h5_file_scc37)
data_scc38 <- Read10X_h5(file=h5_file_scc38)
data_scc39 <- Read10X_h5(file=h5_file_scc39)
data_scc40 <- Read10X_h5(file=h5_file_scc40)
data_scc41 <- Read10X_h5(file=h5_file_scc41)
data_scc42 <- Read10X_h5(file=h5_file_scc42)
data_scc43 <- Read10X_h5(file=h5_file_scc43)
data_scc44 <- Read10X_h5(file=h5_file_scc44)

for(i in u){
  info <- sprintf("已完成 %d%%",round(i*100/length(u)))
  setTkProgressBar(pb, i*100/length(u), sprintf("进度 (%s)", info), info)
}




pbmc_1<- CreateSeuratObject(counts =data_scc1,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_2<- CreateSeuratObject(counts =data_scc2,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_3<- CreateSeuratObject(counts =data_scc3,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_4<- CreateSeuratObject(counts =data_scc4,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_5<- CreateSeuratObject(counts =data_scc5,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_6<- CreateSeuratObject(counts =data_scc6,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_7<- CreateSeuratObject(counts =data_scc7,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_8<- CreateSeuratObject(counts =data_scc8,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_9<- CreateSeuratObject(counts =data_scc9,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_10<- CreateSeuratObject(counts =data_scc10,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_11<- CreateSeuratObject(counts =data_scc11,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_12<- CreateSeuratObject(counts =data_scc12,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_13<- CreateSeuratObject(counts =data_scc13,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_14<- CreateSeuratObject(counts =data_scc14,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_15<- CreateSeuratObject(counts =data_scc15,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_16<- CreateSeuratObject(counts =data_scc16,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_17<- CreateSeuratObject(counts =data_scc17,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_18<- CreateSeuratObject(counts =data_scc18,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_19<- CreateSeuratObject(counts =data_scc19,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_20<- CreateSeuratObject(counts =data_scc20,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_21<- CreateSeuratObject(counts =data_scc21,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_22<- CreateSeuratObject(counts =data_scc22,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_23<- CreateSeuratObject(counts =data_scc23,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_24<- CreateSeuratObject(counts =data_scc24,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_25<- CreateSeuratObject(counts =data_scc25,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_26<- CreateSeuratObject(counts =data_scc26,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_27<- CreateSeuratObject(counts =data_scc27,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_28<- CreateSeuratObject(counts =data_scc28,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_29<- CreateSeuratObject(counts =data_scc29,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_30<- CreateSeuratObject(counts =data_scc30,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_31<- CreateSeuratObject(counts =data_scc31,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_32<- CreateSeuratObject(counts =data_scc32,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_33<- CreateSeuratObject(counts =data_scc33,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_34<- CreateSeuratObject(counts =data_scc34,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_35<- CreateSeuratObject(counts =data_scc35,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_36<- CreateSeuratObject(counts =data_scc36,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_37<- CreateSeuratObject(counts =data_scc37,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_38<- CreateSeuratObject(counts =data_scc38,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_39<- CreateSeuratObject(counts =data_scc39,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_40<- CreateSeuratObject(counts =data_scc40,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_41<- CreateSeuratObject(counts =data_scc41,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_42<- CreateSeuratObject(counts =data_scc42,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_43<- CreateSeuratObject(counts =data_scc43,project = 'pbmc_dm',min.cells=3,min.features=200)
pbmc_44<- CreateSeuratObject(counts =data_scc44,project = 'pbmc_dm',min.cells=3,min.features=200)

pbmc_1[['percent.mt']] <- PercentageFeatureSet(pbmc_1,pattern = '^MT-')
pbmc_2[['percent.mt']] <- PercentageFeatureSet(pbmc_2,pattern = '^MT-')
pbmc_3[['percent.mt']] <- PercentageFeatureSet(pbmc_3,pattern = '^MT-')
pbmc_4[['percent.mt']] <- PercentageFeatureSet(pbmc_4,pattern = '^MT-')
pbmc_5[['percent.mt']] <- PercentageFeatureSet(pbmc_5,pattern = '^MT-')
pbmc_6[['percent.mt']] <- PercentageFeatureSet(pbmc_6,pattern = '^MT-')
pbmc_7[['percent.mt']] <- PercentageFeatureSet(pbmc_7,pattern = '^MT-')
pbmc_8[['percent.mt']] <- PercentageFeatureSet(pbmc_8,pattern = '^MT-')
pbmc_9[['percent.mt']] <- PercentageFeatureSet(pbmc_9,pattern = '^MT-')
pbmc_10[['percent.mt']] <- PercentageFeatureSet(pbmc_10,pattern = '^MT-')
pbmc_11[['percent.mt']] <- PercentageFeatureSet(pbmc_11,pattern = '^MT-')
pbmc_12[['percent.mt']] <- PercentageFeatureSet(pbmc_12,pattern = '^MT-')
pbmc_13[['percent.mt']] <- PercentageFeatureSet(pbmc_13,pattern = '^MT-')
pbmc_14[['percent.mt']] <- PercentageFeatureSet(pbmc_14,pattern = '^MT-')
pbmc_15[['percent.mt']] <- PercentageFeatureSet(pbmc_15,pattern = '^MT-')
pbmc_16[['percent.mt']] <- PercentageFeatureSet(pbmc_16,pattern = '^MT-')
pbmc_17[['percent.mt']] <- PercentageFeatureSet(pbmc_17,pattern = '^MT-')
pbmc_18[['percent.mt']] <- PercentageFeatureSet(pbmc_18,pattern = '^MT-')
pbmc_19[['percent.mt']] <- PercentageFeatureSet(pbmc_19,pattern = '^MT-')
pbmc_20[['percent.mt']] <- PercentageFeatureSet(pbmc_20,pattern = '^MT-')
pbmc_21[['percent.mt']] <- PercentageFeatureSet(pbmc_21,pattern = '^MT-')
pbmc_22[['percent.mt']] <- PercentageFeatureSet(pbmc_22,pattern = '^MT-')
pbmc_23[['percent.mt']] <- PercentageFeatureSet(pbmc_23,pattern = '^MT-')
pbmc_24[['percent.mt']] <- PercentageFeatureSet(pbmc_24,pattern = '^MT-')
pbmc_25[['percent.mt']] <- PercentageFeatureSet(pbmc_25,pattern = '^MT-')
pbmc_26[['percent.mt']] <- PercentageFeatureSet(pbmc_26,pattern = '^MT-')
pbmc_27[['percent.mt']] <- PercentageFeatureSet(pbmc_27,pattern = '^MT-')
pbmc_28[['percent.mt']] <- PercentageFeatureSet(pbmc_28,pattern = '^MT-')
pbmc_29[['percent.mt']] <- PercentageFeatureSet(pbmc_29,pattern = '^MT-')
pbmc_30[['percent.mt']] <- PercentageFeatureSet(pbmc_30,pattern = '^MT-')
pbmc_31[['percent.mt']] <- PercentageFeatureSet(pbmc_31,pattern = '^MT-')
pbmc_32[['percent.mt']] <- PercentageFeatureSet(pbmc_32,pattern = '^MT-')
pbmc_33[['percent.mt']] <- PercentageFeatureSet(pbmc_33,pattern = '^MT-')
pbmc_34[['percent.mt']] <- PercentageFeatureSet(pbmc_34,pattern = '^MT-')
pbmc_35[['percent.mt']] <- PercentageFeatureSet(pbmc_35,pattern = '^MT-')
pbmc_36[['percent.mt']] <- PercentageFeatureSet(pbmc_36,pattern = '^MT-')
pbmc_37[['percent.mt']] <- PercentageFeatureSet(pbmc_37,pattern = '^MT-')
pbmc_38[['percent.mt']] <- PercentageFeatureSet(pbmc_38,pattern = '^MT-')
pbmc_39[['percent.mt']] <- PercentageFeatureSet(pbmc_39,pattern = '^MT-')
pbmc_40[['percent.mt']] <- PercentageFeatureSet(pbmc_40,pattern = '^MT-')
pbmc_41[['percent.mt']] <- PercentageFeatureSet(pbmc_41,pattern = '^MT-')
pbmc_42[['percent.mt']] <- PercentageFeatureSet(pbmc_42,pattern = '^MT-')
pbmc_43[['percent.mt']] <- PercentageFeatureSet(pbmc_43,pattern = '^MT-')
pbmc_44[['percent.mt']] <- PercentageFeatureSet(pbmc_44,pattern = '^MT-')

VlnPlot(pbmc_1,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_2,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_3,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_4,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_5,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_6,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_7,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_8,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_9,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_10,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_11,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_12,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_13,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_14,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_15,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_16,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_17,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_18,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_19,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_20,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_21,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_22,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_23,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_24,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_25,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_26,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_27,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_28,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_29,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_30,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_31,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_32,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_33,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_34,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_35,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_36,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_37,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_38,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_39,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_40,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_41,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_42,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_43,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(pbmc_44,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)


pbmc_1 <- subset(pbmc_1,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_2 <- subset(pbmc_2,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_3 <- subset(pbmc_3,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_4 <- subset(pbmc_4,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_5 <- subset(pbmc_5,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_6 <- subset(pbmc_6,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_7 <- subset(pbmc_7,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_8 <- subset(pbmc_8,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_9 <- subset(pbmc_9,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_10 <- subset(pbmc_10,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_11 <- subset(pbmc_11,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_12 <- subset(pbmc_12,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_13 <- subset(pbmc_13,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_14 <- subset(pbmc_14,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_15 <- subset(pbmc_15,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_16 <- subset(pbmc_16,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_17 <- subset(pbmc_17,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_18 <- subset(pbmc_18,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_19 <- subset(pbmc_19,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_20 <- subset(pbmc_20,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_21 <- subset(pbmc_21,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_22 <- subset(pbmc_22,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_23 <- subset(pbmc_23,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_24 <- subset(pbmc_24,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_25 <- subset(pbmc_25,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_26 <- subset(pbmc_26,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_27 <- subset(pbmc_27,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_28 <- subset(pbmc_28,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_29 <- subset(pbmc_29,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_30 <- subset(pbmc_30,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_31 <- subset(pbmc_31,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_32 <- subset(pbmc_32,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_33 <- subset(pbmc_33,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_34 <- subset(pbmc_34,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_35 <- subset(pbmc_35,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_36 <- subset(pbmc_36,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_37 <- subset(pbmc_37,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_38 <- subset(pbmc_38,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_39 <- subset(pbmc_39,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_40 <- subset(pbmc_40,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_41 <- subset(pbmc_41,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_42 <- subset(pbmc_42,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_43 <- subset(pbmc_43,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
pbmc_44 <- subset(pbmc_44,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)

table(pbmc_1@meta.data$orig.ident)
table(pbmc_2@meta.data$orig.ident)
table(pbmc_3@meta.data$orig.ident)
table(pbmc_4@meta.data$orig.ident)
table(pbmc_5@meta.data$orig.ident)
table(pbmc_6@meta.data$orig.ident)
table(pbmc_7@meta.data$orig.ident)
table(pbmc_8@meta.data$orig.ident)
table(pbmc_9@meta.data$orig.ident)
table(pbmc_10@meta.data$orig.ident)
table(pbmc_11@meta.data$orig.ident)
table(pbmc_12@meta.data$orig.ident)
table(pbmc_13@meta.data$orig.ident)
table(pbmc_14@meta.data$orig.ident)
table(pbmc_15@meta.data$orig.ident)
table(pbmc_16@meta.data$orig.ident)
table(pbmc_17@meta.data$orig.ident)
table(pbmc_18@meta.data$orig.ident)
table(pbmc_19@meta.data$orig.ident)
table(pbmc_20@meta.data$orig.ident)
table(pbmc_21@meta.data$orig.ident)
table(pbmc_22@meta.data$orig.ident)
table(pbmc_23@meta.data$orig.ident)
table(pbmc_24@meta.data$orig.ident)
table(pbmc_35@meta.data$orig.ident)
table(pbmc_36@meta.data$orig.ident)
table(pbmc_37@meta.data$orig.ident)
table(pbmc_38@meta.data$orig.ident)
table(pbmc_39@meta.data$orig.ident)
table(pbmc_40@meta.data$orig.ident)
table(pbmc_41@meta.data$orig.ident)
table(pbmc_42@meta.data$orig.ident)
table(pbmc_43@meta.data$orig.ident)
table(pbmc_44@meta.data$orig.ident)


























C1 <- Read10X(data.dir = 'c1')
C2<- Read10X(data.dir = 'c2')
C3<- Read10X(data.dir = 'c3')
C4 <- Read10X(data.dir = 'c4')
C5<- Read10X(data.dir = 'c5')
C6 <- Read10X(data.dir = 'c6')
D1 <- Read10X(data.dir = 'd1')
D2 <- Read10X(data.dir = 'd2')
D3 <- Read10X(data.dir = 'd3')
D4 <- Read10X(data.dir = 'd4')
D5 <- Read10X(data.dir = 'd5')
D6 <- Read10X(data.dir = 'd6')



C1 <- CreateSeuratObject(counts = C1,project = 'pbmc_con',min.cells=3,min.features=200)
C2 <- CreateSeuratObject(counts = C2,project = 'pbmc_con',min.cells=3,min.features=200)
C3 <- CreateSeuratObject(counts = C3,project = 'pbmc_con',min.cells=3,min.features=200)
C4 <- CreateSeuratObject(counts = C4,project = 'pbmc_con',min.cells=3,min.features=200)
C5 <- CreateSeuratObject(counts = C5,project = 'pbmc_con',min.cells=3,min.features=200)
C6 <- CreateSeuratObject(counts = C6,project = 'pbmc_con',min.cells=3,min.features=200)
D1<- CreateSeuratObject(counts = D1,project = 'pbmc_dm',min.cells=3,min.features=200)
D2<- CreateSeuratObject(counts = D2,project = 'pbmc_dm',min.cells=3,min.features=200)
D3<- CreateSeuratObject(counts = D3,project = 'pbmc_dm',min.cells=3,min.features=200)
D4<- CreateSeuratObject(counts = D4,project = 'pbmc_dm',min.cells=3,min.features=200)
D5<- CreateSeuratObject(counts = D5,project = 'pbmc_dm',min.cells=3,min.features=200)
D6<- CreateSeuratObject(counts = D6,project = 'pbmc_dm',min.cells=3,min.features=200)

C1[['percent.mt']] <- PercentageFeatureSet(C1,pattern = '^MT-')
C2[['percent.mt']] <- PercentageFeatureSet(C2,pattern = '^MT-')
C3[['percent.mt']] <- PercentageFeatureSet(C3,pattern = '^MT-')
C4[['percent.mt']] <- PercentageFeatureSet(C4,pattern = '^MT-')
C5[['percent.mt']] <- PercentageFeatureSet(C5,pattern = '^MT-')
C6[['percent.mt']] <- PercentageFeatureSet(C6,pattern = '^MT-')

D1[['percent.mt']] <- PercentageFeatureSet(D1,pattern = '^MT-')
D2[['percent.mt']] <- PercentageFeatureSet(D2,pattern = '^MT-')
D3[['percent.mt']] <- PercentageFeatureSet(D3,pattern = '^MT-')
D4[['percent.mt']] <- PercentageFeatureSet(D4,pattern = '^MT-')
D5[['percent.mt']] <- PercentageFeatureSet(D5,pattern = '^MT-')
D6[['percent.mt']] <- PercentageFeatureSet(D6,pattern = '^MT-')

VlnPlot(C1,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(C2,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(C3,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(C4,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(C5,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(C6,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)

VlnPlot(D1,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(D2,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(D3,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(D4,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(D5,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
VlnPlot(D6,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)








#数据质控
C1 <- subset(C1,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
C2 <- subset(C2,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
C3 <- subset(C3,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
C4 <- subset(C4,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
C5 <- subset(C5,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
C6 <- subset(C6,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)

D1<- subset(D1,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
D2 <- subset(D2,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
D3 <- subset(D3,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
D4 <- subset(D4,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
D5 <- subset(D5,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)
D6 <- subset(D6,subset = nFeature_RNA>200&nFeature_RNA<3000&percent.mt<15)

table(C1@meta.data$orig.ident)
table(C2@meta.data$orig.ident)
table(C3@meta.data$orig.ident)
table(C4@meta.data$orig.ident)
table(C5@meta.data$orig.ident)
table(C6@meta.data$orig.ident)
table(D1@meta.data$orig.ident)
table(D2@meta.data$orig.ident)
table(D3@meta.data$orig.ident)
table(D4@meta.data$orig.ident)
table(D5@meta.data$orig.ident)
table(D6@meta.data$orig.ident)

sample_list<- list(C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,D5,D6)
sample_list<- list(pbmc_1,pbmc_2,pbmc_3,pbmc_4,pbmc_5,pbmc_6,pbmc_7,pbmc_8,pbmc_9,pbmc_10,pbmc_11,pbmc_12,pbmc_13,pbmc_14,pbmc_15,pbmc_16,pbmc_17,pbmc_18)

sample_list<- list(pbmc_1,pbmc_2,pbmc_3,pbmc_4,pbmc_5,pbmc_6,pbmc_7,pbmc_8,pbmc_9,pbmc_10,pbmc_11,pbmc_12,pbmc_13,pbmc_14,pbmc_15,pbmc_16,pbmc_17,pbmc_18,pbmc_19,pbmc_20,pbmc_21,pbmc_22,pbmc_23,pbmc_24,pbmc_25,pbmc_26,pbmc_27,pbmc_28,
                   pbmc_29,pbmc_30,pbmc_31,pbmc_32,pbmc_33,pbmc_34,pbmc_35,pbmc_36,pbmc_37,pbmc_38,pbmc_39,pbmc_40,pbmc_41,pbmc_42,pbmc_43,pbmc_44)


pbmc_harmony <- merge(x=sample_list[[1]],y=sample_list[-1])
pbmc_harmony <- JoinLayers(pbmc_harmony)
all.genes <-rownames(pbmc_harmony )


all.genes <-rownames(pbmc_harmony)

print(memory.profile())






#标准化
pbmc_harmony <- NormalizeData(pbmc_harmony)%>%FindVariableFeatures()%>%ScaleData(features = all.genes)%>%RunPCA(verbose = FALSE)

library(harmony)
pbmc_harmony <- RunHarmony(pbmc_harmony,group.by.vars='orig.ident')#去除批次效应

mito_genes=rownames(pbmc_harmony)[grep("^MT-",rownames(pbmc_harmony))]
mito_genes
pbmc_harmony =PercentageFeatureSet(pbmc_harmony,"^MT-",col.name = "percent_mito")
fivenum(pbmc_harmony@meta.data$percent_mito)

ribo_genes=rownames(pbmc_harmony)[grep("^Rp[sl]",rownames(pbmc_harmony),ignore.case = T)]
ribo_genes
pbmc_harmony=PercentageFeatureSet(pbmc_harmony,"^RP[SL]",col.name = "percent_ribo")
fivenum(pbmc_harmony@meta.data$percent_ribo)

rownames(pbmc_harmony)[grep("^Hb[^(p)]",rownames(pbmc_harmony),ignore.case = T)]
pbmc_harmony=PercentageFeatureSet(pbmc_harmony,"^HB[^(P)]",col.name = "percent_hb")
fivenum(pbmc_harmony@meta.data$percent_hb)

feats<-c("nFeature_RNA","nCount_RNA")
p1=VlnPlot(pbmc_harmony,group.by = "orig.ident",features = feats,pt.size = 0.01,ncol = 2)+
  NoLegend()

p1
ggsave(filename = "vlnplot1.pdf",width = 250,height = 300,units = "mm",plot = p1)

feats<-c("percent_mito","percent_ribo","percent_hb")
p2=VlnPlot(pbmc_harmony,group.by = "orig.ident",features = feats,pt.size = 0.01,ncol = 3,same.y.lims = T)+
  scale_y_continuous(breaks = seq(0,100,5))+
  NoLegend()
p2
ggsave(filename = "vlnplot2.pdf",width = 250,height = 300,units = "mm",plot = p2)
library(tidyverse)
#过滤低质量细胞
View(pbmc_harmony)
selected_c <-WhichCells(pbmc_harmony,expression = nFeature_RNA >300)
selected_F <- rownames(pbmc_harmony)[Matrix::rowSums(pbmc_harmony@assays$RNA@counts>0) >3]
#提取部分基因
Sce.all.filt <-subset(pbmc_harmony,features=selected_F,cells=selected_c)
dim(pbmc_harmony)
dim(Sce.all.filt )
#过滤核糖体等
selected_mito<-WhichCells(Sce.all.filt,expression = percent_mito <15)
selected_ribo<-WhichCells(Sce.all.filt,expression = percent_ribo >3)
selected_hb<-WhichCells(Sce.all.filt,expression = percent_hb<0.1)
length(selected_mito)
length(selected_ribo)
length(selected_hb)
Sce.all.filt <-subset(Sce.all.filt,cells=selected_mito)
Sce.all.filt <-subset(Sce.all.filt,cells=selected_ribo)
Sce.all.filt <-subset(Sce.all.filt,cells=selected_hb)
dim(Sce.all.filt )
plot1 <-FeatureScatter(Sce.all.filt ,feature1 ="nCount_RNA",feature2 ="percent_mito"  )
plot2 <-FeatureScatter(Sce.all.filt ,feature1 ="nCount_RNA",feature2 ="nFeature_RNA"  )
plot1+plot2

table(Sce.all.filt $orig.ident)

feats<-c("nCount_RNA","nFeature_RNA")
p1_filtered=VlnPlot(Sce.all.filt,group.by = "orig.ident",features = feats,pt.size = 0.1,ncol = 2)+
  NoLegend()
p1
ggsave(filename = "vlnplot1_filtered.pdf",width = 250,height = 300,units = "mm",plot = p1)
feats<-c("nCount_RNA","nFeature_RNA")
feats<-c("percent_mito","percent_ribo","percent_hb")
p2=VlnPlot(pbmc_harmony,group.by = "orig.ident",features = feats,pt.size = 0.1,ncol = 3,same.y.lims = T)+
  scale_y_continuous(breaks = seq(0,100,5))+
  NoLegend()

p2
ggsave(filename = "vlnplot2_filtered.pdf",width = 250,height = 300,units = "mm",plot = p2)

ElbowPlot(Sce.all.filt,ndims = 50,reduction="harmony")#预分为50群
pbmc_harmony <- FindNeighbors(Sce.all.filt,reduction='harmony',dims = 1:9)%>%FindClusters(resolution =0.5)

table(pbmc_harmony@meta.data$seurat_clusters)
pbmc_harmony<-RunUMAP(pbmc_harmony,reduction ="harmony",dims = 1:9 )

plot1=DimPlot(pbmc_harmony,reduction = "umap",label = T)
plot2=DimPlot(pbmc_harmony,reduction = "umap",group.by = "orig.ident")

plotc<-plot1+plot2
plotc
ggsave(filename = "vlnplot2_16dim.pdf",width = 15,height = 7,plot = plotc)
saveRDS(pbmc_harmony,"pbmc_harmony.rds")

#tsne
plot1=DimPlot(pbmc_harmony,reduction ="tsne",label="T")
plot2=DimPlot(pbmc_harmony,reduction ="tsne",group.by = "orig.ident")
markers<-FindAllMarkers(object=pbmc_harmony,test.use="wilcox",
                        only.pos = TRUE,logfc.threshold = 0.25)


all.markers=markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)

top10<-all.markers%>%group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
write.csv(top10,"/12top10.csv",row.names = T)
#可视化
DoHeatmap(pbmc_harmony,features = top10$gene)+NoLegend()
VlnPlot(pbmc_harmony,features = top10$gene[1:10])
vlnplot(pbmc_harmony,Features=c("ifit1"))

table(pbmc_harmony@meta.data$seurat_clusters)

p<-DotPlot(pbmc_harmony,features = top10$gene[90:100])
genes_to_check=c('PTPRC','CD3D','CD3E','CD4','CD8A',#Tcell
                 'CSPG4',#PERCICYTE
                'ENPPS','FCER1A',
                'AQP3','MUC5AC','MUC5B','PIGR','SCGB1A1',
               'THBD','CD1C','CLEC4C','CD83','HLA.DOA2','HLA.DOA1',
                'KRT19','ITGAE','VCAM1','IL6R','ANPEP','CD24','CDH1',
                 'MUC1','ABCA3','LPCAT1','NAP5A','SFTPB','SFCTC','SLC3A2',
                 'ACTA2','PDGFRA','PDGFRB','THY1',
                 'CD19','CD79A','MS4A1',
                 'TAGLN2','CD5',
                 'CD27','CD38','LY9','LAIR1','ICAM1','KIT',
                 'NKG7','GNLY',
                 'CD8B','ZNF683','FCGR3A','FCGR3B','NCAM1','KLRB1',
                 'CXCR5','CCR7',
                 'CD6','IL7R','IL2RA','IKZF2',
                 'GP1BA','SELL','IFNG','CXCR3','IL17A','IL4','GATA3',
                 'CD33','ENTPD1',
                 'S100A8','S100A9','S100A12',
                 'CD68','CD163','MRC1','MSR1','CXCL10','CCL18',
                 'PECAM1','VWF','MCAM','CD34','ESM1',
                 'ALDH1A1','KRT18','PROM1'
)
genes_to_check=c('PTPRC','CD3D','CD3E','CD4','CD8A')
p<-DotPlot(pbmc_harmony,features = top10$gene[90:100])

DefaultAssay(sce) = "RNA"
DotPlot(sce,features = unique(genes_to_check),group.by = 'RNA_snn_res.1.5') + coord_flip()ggsave("har_gene_show.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show.png",width = 14,height = 12,dpi = 500)
celltype=data.frame(clusterID=0:16,
                    celltype="unkown")
celltype[celltype$clusterID%in%c(0,1),2]="T cell"
celltype[celltype$clusterID%in%c(2,13),2]="Smooth muscle cell"
celltype[celltype$clusterID%in%c(3,4),2]="Macrophage"
celltype[celltype$clusterID%in%c(5),2]="Fibroblast"
celltype[celltype$clusterID%in%c(6,15),2]="Dendritic cell"
celltype[celltype$clusterID%in%c(7,8,9,10),2]="Endothelial cell"
celltype[celltype$clusterID%in%c(11,14),2]="B cell"
celltype[celltype$clusterID%in%c(12),2]="Mast cell"
celltype[celltype$clusterID%in%c(16),2]="Cycling cell"


celltype[celltype$clusterID%in%c(1),2]="natural regulatory T cell"
celltype[celltype$clusterID%in%c(3),2]="Smooth muscle cell"


celltype[celltype$clusterID%in%c(6),2]="Megakaryocyte"

celltype[celltype$clusterID%in%c(8),2]="Myeloid cell"
celltype[celltype$clusterID%in%c(10),2]="Pro-tumor type-2 pericyte"




celltype[celltype$clusterID%in%c(16),2]="B cell"

celltype[celltype$clusterID%in%c(3),2]="Circulation fetal cell"


celltype[celltype$clusterID%in%c(6),2]="Natural killer cell"



celltype[celltype$clusterID%in%c(11),2]="Stem cell"

celltype[celltype$clusterID%in%c(16),2]="Monocyte"

celltype[celltype$clusterID%in%c(19),2]="Plasmacytoid dendritic cells"

celltype[celltype$clusterID%in%c(5),2]="activated effector cell"

celltype[celltype$clusterID%in%c(6,7),2]="megakaryocyte"
celltype[celltype$clusterID%in%c(11),2]="endothelial cell"

celltype[celltype$clusterID%in%c(13),2]="pro-tumor type-4 pericyte"
celltype[celltype$clusterID%in%c(14),2]="dendritic.cell"
celltype[celltype$clusterID%in%c(17),2]="plasmacytoid dendritic cells"
celltype[celltype$clusterID%in%c(16),2]="cycling cell"

celltype
table(celltype$celltype)

sce.in=pbmc_harmony
sce.in@meta.data$celltype="NA"

for (i in 1:nrow(celltype)) {
sce.in@meta.data[which(sce.in@meta.data$seurat_clusters==celltype$clusterID[i]),'celltype']<-celltype$celltype[i]  
}
table(sce.in@meta.data$celltype)

table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)
sce=sce.in

p.dim.cell=DimPlot(sce,reduction = "umap",group.by = "celltype",label = T,pt.size = 1)
p.dim.cell
ggsave(plot=p.dim.cell,filename = "diplot.pdf",width = 9,height = 7)
VlnPlot(sce,features = "IFIT1",group.by = "celltype",pt.size = 1)
DefaultAssay(sce)<-"RNA"

save(sce,file="jiyin")












save(sce,file="jiyin")

pbmc_harmony=sce
pbmc_harmony$celltype=Idents(pbmc_harmony)


library(msigdbr)
library(org.Hs.eg.db)
library(gplots)























































#差异分析+富集分析
pbmc_harmony=sce
names(pbmc_harmony@meta.data)
unique(pbmc_harmony$orig.ident)
DimPlot(pbmc_harmony,split.by = 'orig.ident')
Idents(pbmc_harmony)="celltype"
pbmc_harmony$celltype.group<-paste(pbmc_harmony$celltype,pbmc_harmony$orig.ident,sep = "_")
pbmc_harmony$celltype<-Idents(pbmc_harmony)
Idents(pbmc_harmony)<-"celltype.group"

cell.marker <- FindAllMarkers(object =pbmc_harmony,
                              only.pos = FALSE,
                              test.use = "wilcox",
                              slot = "data",
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
write.csv(cell.marker,"cell.csv")
table(pbmc_harmony$orig.ident)
prop.table(table(Idents(pbmc_harmony)))
table(Idents(pbmc_harmony),pbmc_harmony$orig.ident)
Cellratio<-prop.table(table(Idents(pbmc_harmony),pbmc_harmony$orig.ident),margin = 2)
Cellratio<-as.data.frame(Cellratio)
colnames(Cellratio)<-c("celltype","sample","ratio")
colourCount=length(unique(Cellratio$celltype))
ggplot(Cellratio)+
  geom_bar(aes(x=sample,y=ratio,fill = celltype),stat = "identity",width = 0.7,size=0.5,color='#222222')+
  theme_classic()+
  labs(x='sample',y='ratio')+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5,linetype = "solid"))






# 整个组织特定组别的差异分析 ---------------------------------------------------------------
#filter(p_val_adj<0.05)%>%可以在需要矫正p时加上。另外，group-by的值是分组那一列
degs =FindMarkers(pbmc_harmony,logfc.threshold=0.25,only.pos=FALSE,ident.1 = 'pbmc_dm',ident.2 = 'pbmc_con',group.by = 'orig.ident')%>%filter(p_val_adj<0.05)%>%mutate(gene=rownames(.))
#如果是样本之间的差异分析，则把ident1和2改成特定样本名，将group-by改到样本列名


# 两细胞类型间的差异分析（idnet1较ident2基础的差异分析基因） -------------------------------------------------------------
double_degs=FindMarkers(pbmc_harmony,logfc.threshold=0.25,only.pos=FALSE,
                        ident.1 = 'Tcell_',#分子
                        ident.2 = 'Macrophage_'#分母
)


# 组织中指定类型细胞内部不同分组的差异分析 -------------------------------------------------------------
#可以加上矫正p值(ident1与ident2的比值)
cells_degs =FindMarkers(subset(pbmc_harmony,celltype=='natural regulatory T cell'),
                        logfc.threshold=0.5,
                        ident.1 = 'pbmc_dm',#分子
                        ident.2 = 'pbmc_con',#分母
                        group.by = 'orig.ident')%>%mutate(gene=rownames(.))
# 保存数据
write.csv(cells_degs, file = "cardiomyocyte-cyfx.csv")

#差异分析的显著基因
top100gene <- degs%>%top_n(n=100,wt=avg_log2FC)%>%rownames()
bottom100gene <- degs%>%top_n(n=-100,wt=avg_log2FC)%>%rownames()
mygenest <- cbind(top100gene,bottom100gene)
mygenest <- as.data.frame(mygenest)
#保存差异分析基因数据
write.csv(mygenest, file = "cardiomyocyte-top-bottom.csv")





# 所有细胞进行笼统差异分析 ------------------------------------------------------------
#其含义是通过循环将每个类群与其他所有类群的基因进行差异分析，即一个类群中某基因（无论con与dm）与另外所有类群细胞（无论con与dm）的差异对比
markers_updown <- FindAllMarkers(object = pbmc_harmony,test.use = 'wilcox',only.pos = FALSE,logfc.threshold = 0.25)


# 分组细化后精确到细胞类型特定分组的差异分析 ---------------------------------------------------
#将每个细胞细化得到其特定细胞类型+分组处理的标签
pbmc_harmony$celltype.group <- paste(pbmc_harmony$celltype,pbmc_harmony$orig.ident,sep = '-')#注意这里的orig.ident代表的是分组，可能是group，根据情况定
pbmc_harmony$celltype.group
Idents(pbmc_harmony) <- 'celltype.group'
#报错原因：某组所含细胞数少
table(Idents(pbmc_harmony))#查看每组细胞数
#任何细胞类型+分组处理细胞之间差异分析
anydeg_cardiomyocyte <- FindMarkers(pbmc_harmony,
                                    ident.1 = 'natural regulatory T cell_pbmc_dm',#对比的细胞类型1（分子）
                                    ident.2 = 'natural regulatory T cell_pbmc_con',#对比的细胞类型2（分母）
                                    verbose=TRUE,#分析的进度条
                                    test.use='wilcox',#进行统计学分析的检验方法，wilcox为秩和检验
                                    min.pct=0.1)#基因在两组细胞间最低的表达比例，即纳入统计范围的最低比例


# 当细胞数量均较大：》3时，可通过循环分析 ----------------------------------------------------


#通过循环将dcm组与con组中所有各类细胞群单独差异分析
cellfordeg <- levels(pbmc_harmony$celltype)
cellfordeg[2]
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc_harmony,ident.1 = paste0(cellfordeg[i],'-pbmc_dm'),ident.2 = paste0(cellfordeg[i],'-pbmc_con'),verbose=TRUE)
  write.csv(CELLDEG,paste0(cellfordeg[i],'.csv'))
}


# 组织整体差异分析 ----------------------------------------------------------------

group_degs =FindMarkers(pbmc_harmony,
                        logfc.threshold=0.25,
                        only.pos=FALSE,
                        ident.1 = 'pbmc_dm',
                        ident.2 = 'pbmc_con',
                        group.by = 'orig.ident')%>%
  #filter(p_val_adj<0.05)%>%
  mutate(gene=rownames(.))





# 富集分析 --------------------------------------------------------------------
#获取各分组各细胞群的平均表达量的提取
library(GSVA)
genes <- unique(rownames(pbmc_harmony@assays$RNA$counts))





#表达矩阵提取
gene.expr <- as.matrix(pbmc_harmony[["RNA"]]$data)
dim(gene.expr)
library(GSVA)
gsva.result <- GSVA::gsva(gene.expr,mygenest,kcdf='Gaussian')

# 火山图 ---------------------------------------------------------------------
#参数
#参数设置
##,aesCol=c('',''))设置上调、下调点的颜色
#col.type='adjustP'设置矫正方式
#flip=T将火山图横置(特别慢)
#polar=T环形火山图



###火山图（组织差异分析）(该图含义是某个细胞类群相比于其他所有细胞类群加在一起，它所更明显表达 变化的基因类型，按照logfc排序的前五位)
library(scRNAtoolVis)
cell_type_cols <- c("#FEF295","#C9D79F","#6D93A9","#5C56A6","#B2258F","#FEF295","#C9D79F","#6D93A9","#5C56A6","#B2258F","#E90D71","#D23740","#B95B1C")#注意个数应与细胞群数对等
#该火山图是按照orig.ident种类数量进行分群的
#其含义是通过循环将每个类群与其他所有类群的基因进行差异分析，即一个类群中某基因（无论con与dm）与另外所有类群细胞（无论con与dm）的差异对比
jjVolcano(diffData = markers_updown,log2FC.cutoff = 0.25,size=3.5,fontface='italic',tile.col = cell_type_cols[1:9],topGeneN = 5)
ggsave(filename = '组织整合差异分析火山图.pdf',width = 12,height = 9)
#参数设置
##,aesCol=c('',''))设置上调、下调点的颜色
#col.type='adjustP'设置矫正方式

###火山图（特定组与组之间差异分析）
library(scRNAtoolVis)
cell_type_cols <- c("#FEF295","#C9D79F","#6D93A9","#5C56A6","#B2258F")#注意个数应与细胞群数对等
markers_fil <- markers_updown%>%filter(cluster%in%c('cardiomyocyte','Macrophage','fibroblast'))
jjVolcano(diffData = markers_fil,log2FC.cutoff = 0.25,size=3.5,fontface='italic',topGeneN = 5)
ggsave(filename = '特定组间差异分析火山图.pdf',width = 12,height = 9)

#火山图（展示特定基因在各细胞类群的表达情况）
cell_type_cols <- c("#FEF295","#C9D79F","#6D93A9","#5C56A6","#B2258F")#注意个数应与细胞群数对等
mygenes <- c('Tom20','Vmf','Mark4','Nkg7')
jjVolcano(diffData =markers_updown,tile.col = cell_type_cols[1:5],myMarkers = mygenes )
ggsave(filename = '特定基因在各细胞群中表达情况火山图.pdf',width = 12,height = 9)

###火山图（单独细胞群）
#火1
library(ggVolcano)
#注意log2FC_name和fdr_name参数需要对应准确
pbmc_harmony_cardiomyocyte_supper <- add_regulate(cells_degs,log2FC_name ="avg_log2FC",fdr_name = "p_val",log2FC=1,fdr=0.05)
write.csv(cells_degs,file = "afterdcm.csv")
#注意x="log2FoldChange",y="padj"固定不变
ggvolcano(pbmc_harmony_cardiomyocyte_supper,x="log2FoldChange",y="padj",fills=c("#e94234","#b4b4d8","#269846"),colors=c("#e94234","#b4b4d8","#269846"),label="gene",label_number = 10,output = FALSE)
ggsave(filename = '单个细胞群差异分析火山图1.pdf',width = 12,height = 9)

#火2（渐变）
library(ggVolcano)
library(ggVolcano)
library(RColorBrewer)
#注意log2FC_name和fdr_name参数需要对应准确
pbmc_harmony_cardiomyocyte_supper <- add_regulate(cells_degs,log2FC_name ="avg_log2FC",fdr_name = "p_val",log2FC=1,fdr=0.05)
#注意x="log2FoldChange",y="padj"固定不变
gradual_volcano(pbmc_harmony_cardiomyocyte_supper,x="log2FoldChange",y="padj",fills=brewer.pal(5,"RdYlBu"),colors=brewer.pal(8,"RdYlBu"),label="gene",label_number = 10,output = FALSE)
ggsave(filename = '单个细胞群差异分析火山图2.pdf',width = 12,height = 9)


# 热图展示差异基因 ----------------------------------------------------------------
library(scRNAtoolVis)
####多组间作图
#pbmc_harmony <- JoinLayers(pbmc_harmony)#将多组变量名整合
deg <- FindAllMarkers(pbmc_harmony,min.pct=0.25,logfc.threshold = 0.25)
#获取每种细胞类型的top10上调基因
deg_top5 <- deg%>%dplyr::group_by(cluster)%>%dplyr::top_n(n=5,wt=avg_log2FC)
###作图
p1 <- averageHeatmap(object = pbmc_harmony,
                     markerGene = deg_top5$gene,
                     group.by = 'celltype',
                     gene.order = deg_top5$gene)
p1
#调细胞类群分类颜色
colors <- c("#0000FF","#1C00E2","#3800C6","#5500AA","#71008D","#8D0071","#AA0055","#C60038",
            "#E2001C","#FF0000","#7FC97F","#9CBCA6","#B9AFCE","#D7B5B4","#F4BD90","#FDD58C",
            "#FEF295","#C9D79F","#6D93A9","#5C56A6","#B2258F","#E90D71","#D23740","#B95B1C",
            "#8f6041","#666666")
mycol1 <- colors[1:5]#有几个细胞类群就是多少
p2 <- averageHeatmap(object = pbmc_harmony,
                     markerGene = deg_top5$gene,
                     group.by = 'celltype',
                     gene.order = deg_top5$gene,
                     annoCol = T,
                     myanCol = mycol1)
p2
#仅展示部分感兴趣的基因
annogene <- c('Esam','Robo4','Shank3','Olfml2b')
p3 <-averageHeatmap(object = pbmc_harmony,
                    markerGene = deg_top5$gene,
                    group.by = 'celltype',
                    gene.order = deg_top5$gene,
                    showRowNames = F,#隐藏所有基因
                    markGenes = annogene)
p3


###单细胞群热图分析
#可以加上矫正p值
cells_degs =FindMarkers(subset(pbmc_harmony,celltype=='cardiomyocyte'),
                        logfc.threshold=0.5,
                        ident.1 = 'pbmc_dm',
                        ident.2 = 'pbmc_con',
                        group.by = 'orig.ident')%>%mutate(gene=rownames(.))

#观察cell_deg文件，只留下以基因名为行名，第一列为表达量


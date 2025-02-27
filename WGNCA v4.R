###############################################################################
#                  WGCNA 分析脚本（含候选关键基因筛选，基于 GS）
#                  可用于 RNA-Seq 或 Microarray 数据
#                  2025-xx-xx  by ChatGPT
###############################################################################

# ------------------- 0. 内存清理  -------------------
rm(list = ls())
gc()

# ------------------------------------------------------------
#      1. 参数设置 & 环境准备
# ------------------------------------------------------------
# 加载必要的库
library(WGCNA)
library(limma)           # Microarray 预处理 & voom 等
library(matrixStats)     # 用于 avereps
library(edgeR)           # 若使用 voom，需要先用 edgeR 的过滤函数
# 若使用 DESeq2，需要加载
# 若没有安装，请先运行 install.packages("BiocManager"); BiocManager::install("DESeq2")
library(DESeq2)

# 允许多线程（如不需要可注释）
# enableWGCNAThreads()

# ------------------- 用户可自由修改的参数 -------------------
# 工作目录（如果为 NULL 则自动使用脚本所在目录）
user_wd <- NULL  

# 指定数据类型（"RNA-Seq" 或 "Microarray"）
data_type <- "Microarray"   # "RNA-Seq" 或 "Microarray"

# 若 data_type == "RNA-Seq"，可指定变异稳定化处理方式
# 可选： "DESeq2" 或 "voom"
rna_seq_method <- "DESeq2"  

# 输入文件：表达矩阵（假设是count或microarray表达量） 和 分组信息
expFile   <- "merge.txt"
notesFile <- "notes.txt"

# 选择保留基因的变异百分比（即取 top X% 标准差最大的基因），默认 25%
topVarPercent <- 25

# 选择用于 pickSoftThreshold 的幂指数范围
powersToTry <- 1:20

# 样本层次聚类时判断离群样本的高度
sampleClustCutHeight <- 20000

# 动态剪切时最小模块基因数（例如 30）
minModuleSize <- 30

# 模块合并阈值（1 - 相关系数的剪切值）
MEDissThres <- 0.25

# ------------------- 关于 GS 和候选关键基因 -------------------
# 指定感兴趣的性状列名称（trait_data 中的列名），默认取第1列；如有需要可修改
trait_of_interest <- "MDD"  # 若为 NULL 则默认取 trait_data 的第1列

# 筛选候选基因的条件：GS > 0.5 且 p-value < 0.05
GS_threshold <- 0.5
GS_p_threshold <- 0.05

# ------------------------------------------------------------
# 存储数据处理过程的日志，用于输出到文件
steps_history <- character()

# ------------------------------------------------------------
# 函数1：路径处理函数（兼容 Windows）
processPath <- function(path) {
  if (.Platform$OS.type == "windows") {
    path <- gsub("\\\\", "/", path)
  }
  return(path)
}

# ------------------------------------------------------------
# 函数2：获取脚本所在目录
getScriptDir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, args)
  if (length(match) > 0) {
    script_path <- sub(file_arg, "", args[match])
    return(dirname(script_path))
  } else {
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      return(dirname(rstudioapi::getActiveDocumentContext()$path))
    } else {
      warning("无法确定脚本所在目录，使用当前工作目录。")
      return(getwd())
    }
  }
}

# ------------------------------------------------------------
# 函数3：设置工作目录
set_working_directory <- function(provided_wd = NULL) {
  if (is.null(provided_wd) || provided_wd == "") {
    script_dir <- getScriptDir()
    message("未提供工作目录，使用脚本所在目录：", script_dir)
    processed_wd <- processPath(script_dir)
  } else {
    processed_wd <- processPath(provided_wd)
    message("设置工作目录为：", processed_wd)
  }
  tryCatch({
    setwd(processed_wd)
    message("工作目录已设置为：", getwd())
  }, error = function(e) {
    stop("无法设置工作目录为：", processed_wd, "\n错误信息：", e$message)
  })
}

# ------------------------------------------------------------
# 执行：设置工作目录
set_working_directory(user_wd)

# ------------------------------------------------------------
#      2. 读取原始数据
# ------------------------------------------------------------
# 检查 data_type 参数
if (!data_type %in% c("RNA-Seq", "Microarray")) {
  stop('data_type 参数错误，请在脚本顶部将 data_type 设置为 "RNA-Seq" 或 "Microarray"')
}
steps_history <- c(steps_history, paste("数据类型选择:", data_type))

# 读取表达矩阵文件
expFile <- processPath(expFile)
if (!file.exists(expFile)) {
  stop("表达矩阵文件不存在，请检查路径: ", expFile)
}

expression_data <- read.table(expFile, sep="\t", header=TRUE, check.names=FALSE)
if (nrow(expression_data) < 1) {
  stop("表达矩阵文件中没有有效的行（基因），请检查文件格式。")
}
if (ncol(expression_data) < 2) {
  stop("表达矩阵文件中列数不足（至少包含基因名称和表达量），请检查文件格式。")
}

expression_data_matrix <- as.matrix(expression_data)
rownames(expression_data_matrix) <- expression_data_matrix[, 1]
expression_data <- expression_data_matrix[, -1, drop=FALSE]

# 转为数值矩阵
data_matrix <- matrix(as.numeric(as.matrix(expression_data)),
                      nrow=nrow(expression_data),
                      dimnames=dimnames(expression_data))

# 读取分组信息（notes.txt）
notesFile <- processPath(notesFile)
if (!file.exists(notesFile)) {
  stop("分组信息文件不存在，请检查路径: ", notesFile)
}
sample_info <- read.table(notesFile, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("sampleName", "group") %in% colnames(sample_info))) {
  stop("notes.txt 缺少必要的列：sampleName 或 group，请检查。")
}

# 只保留在分组信息中的样本
common_samples <- intersect(sample_info$sampleName, colnames(data_matrix))
if (length(common_samples) < 2) {
  stop("notes.txt 与表达矩阵的公共样本数 < 2，无法进行聚类。请检查文件匹配情况。")
}

# 仅保留有分组信息的样本的表达数据
data_matrix <- data_matrix[, common_samples, drop=FALSE]
sample_info_filtered <- sample_info[sample_info$sampleName %in% common_samples, ]

steps_history <- c(steps_history, paste("根据分组信息筛选样本：保留", length(common_samples), "个样本"))

# 绘制原始数据箱线图（用于手动判断是否需要进一步处理）
par(mfrow = c(1,1))
boxplot(data_matrix, main = "原始数据分布 (请检查是否需要进一步处理)",
        las=2, cex.axis=0.7)
cat("请在 RStudio 窗口查看箱线图: 原始数据分布。\n")

# ------------------------------------------------------------
#      3. 数据预处理：根据 data_type + rna_seq_method 选择不同流程
# ------------------------------------------------------------
if (data_type == "RNA-Seq") {
  
  # 用户可在此处选择 DESeq2 或 voom 两种方式
  if (rna_seq_method == "DESeq2") {
    steps_history <- c(steps_history, "执行RNA-Seq count数据标准化流程（DESeq2 + VST）。")
    
    # 构建分组因子
    group_factor <- factor(sample_info_filtered$group)
    # 构建 DESeq2 对象
    dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                  colData = data.frame(group = group_factor,
                                                       row.names = colnames(data_matrix)),
                                  design = ~ group)
    
    # 过滤低表达基因（简单过滤或由 DESeq2 自动处理）
    keep <- rowSums(counts(dds)) >= 1
    dds <- dds[keep, ]
    steps_history <- c(steps_history, paste("初步过滤低表达基因后剩余基因数:", nrow(dds)))
    
    # 运行 DESeq (估计离散度等)
    dds <- DESeq(dds)
    
    # 使用 VST 做变异稳定化变换
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    data_matrix <- assay(vsd)
    steps_history <- c(steps_history, "使用 DESeq2 的 VST 对count数据进行变异稳定化处理。")
    
  } else if (rna_seq_method == "voom") {
    steps_history <- c(steps_history, "执行RNA-Seq count数据标准化流程（edgeR + voom）。")
    
    # 先构建 DGEList
    dge <- DGEList(counts = data_matrix)
    
    # 过滤低表达基因
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    steps_history <- c(steps_history, paste("过滤低表达基因后剩余基因数:", nrow(dge$counts)))
    
    # 计算标准化因子 (TMM)
    dge <- calcNormFactors(dge)
    
    # 构建设计矩阵（假设单因素 group）
    group_factor <- factor(sample_info_filtered$group)
    design <- model.matrix(~group_factor)
    
    # 用 voom 进行变异稳定化处理
    voom_res <- voom(dge, design = design, plot = FALSE)
    data_matrix <- voom_res$E
    steps_history <- c(steps_history, "使用 limma::voom 对count数据进行变异稳定化处理。")
    
  } else {
    stop("rna_seq_method 参数错误，请设置为 'DESeq2' 或 'voom'")
  }
  
} else { 
  # --------------- Microarray 情况 ---------------
  steps_history <- c(steps_history, "执行Microarray数据处理流程。")
  
  # 如果最大值过大（>100），说明可能未做 log2 转换，则进行转换
  max_val <- max(data_matrix, na.rm=TRUE)
  if (max_val > 100) {
    data_matrix <- log2(data_matrix + 1)
    steps_history <- c(steps_history, "检测到表达值较大，执行 log2 转换 (data+1)。")
  } else {
    steps_history <- c(steps_history, "表达值较小，未进行 log2 转换。")
  }
  # 分位数归一化
  data_matrix <- normalizeBetweenArrays(data_matrix)
  steps_history <- c(steps_history, "使用 normalizeBetweenArrays 进行分位数归一化。")
}

# 绘制预处理后数据箱线图，观察标准化效果
par(mfrow = c(1,1))
boxplot(data_matrix, main = "预处理后数据分布 (检查标准化效果)",
        las=2, cex.axis=0.7)
cat("请在 RStudio 窗口查看箱线图: 预处理后数据分布。\n")

# ------------------------------------------------------------
#      4. 后续WGCNA流程（网络构建）
# ------------------------------------------------------------
# (1) 去重合并（若同一基因存在多个探针，则取平均表达）
data_matrix <- avereps(data_matrix)

# (2) 过滤低方差基因——这里替换为选择排名前 topVarPercent 的基因
all_sd <- apply(data_matrix, 1, sd)
num_genes <- length(all_sd)
num_keep <- floor(num_genes * (topVarPercent/100))
# 若保留数小于 1，则至少保留 1 个
num_keep <- max(num_keep, 1)
sorted_sd <- sort(all_sd, decreasing = TRUE)
cutoff_value <- sorted_sd[num_keep]
data_matrix <- data_matrix[all_sd >= cutoff_value, , drop=FALSE]
steps_history <- c(steps_history, paste("保留表达量标准差排名前", topVarPercent, "%的基因，剩余:", nrow(data_matrix)))
if (nrow(data_matrix) < 2) {
  stop("过滤后基因数 < 2，无法进行WGCNA，请调整 topVarPercent 或检查数据质量。")
}

# 转置表达矩阵：WGCNA要求行=样本，列=基因
dat_expr <- t(data_matrix)

# 检查缺失值
gsg <- goodSamplesGenes(dat_expr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("移除基因:", paste(colnames(dat_expr)[!gsg$goodGenes], collapse = ", "), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("移除样本:", paste(rownames(dat_expr)[!gsg$goodSamples], collapse = ", "), "\n")
  }
  dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes]
}
if (nrow(dat_expr) < 2) {
  stop("数据清理后有效样本数 < 2，无法进行聚类分析。请检查前面步骤。")
}

# ------------------------------------------------------------
# 样本聚类：检测离群样本
# ------------------------------------------------------------
sample_tree <- hclust(dist(dat_expr), method = "average")
pdf("1.样品聚类.pdf", width = 15, height = 10)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sample_tree, main = "样品聚类（检测离群样品）", xlab = "", sub = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = sampleClustCutHeight, col = "red")
dev.off()

clust <- cutreeStatic(sample_tree, cutHeight = sampleClustCutHeight, minSize = 2)
major_label <- which.max(table(clust))
keep_samples <- (clust == major_label)
if (sum(keep_samples) < 2) {
  stop("离群样品过滤后仅剩 <2 个样品，无法继续WGCNA，请调整 sampleClustCutHeight 或检查数据质量。")
}
dat_expr <- dat_expr[keep_samples, ]
sample_info_filtered <- sample_info_filtered[match(rownames(dat_expr), sample_info_filtered$sampleName), ]
steps_history <- c(steps_history, paste0("样品聚类（cutHeight=", sampleClustCutHeight,
                                        "）后保留样品数=", nrow(dat_expr)))

# ------------------------------------------------------------
# 构造性状矩阵（Trait Data）
# ------------------------------------------------------------
group_factor <- factor(sample_info_filtered$group)
trait_data <- model.matrix(~0 + group_factor)  # 哑变量矩阵
colnames(trait_data) <- levels(group_factor)
rownames(trait_data) <- sample_info_filtered$sampleName

# 绘制带性状信息的样品聚类热图
sample_tree2 <- hclust(dist(dat_expr), method = "average")
pdf("2.样品聚类热图.pdf", width = 15, height = 10)
plotDendroAndColors(sample_tree2,
                    colors = numbers2colors(trait_data, signed = FALSE),
                    groupLabels = colnames(trait_data),
                    main = "样品聚类与性状热图")
dev.off()

# ------------------------------------------------------------
# pickSoftThreshold：确定软阈值
# ------------------------------------------------------------
sft <- pickSoftThreshold(dat_expr, powerVector = powersToTry, verbose = 5)
soft_power <- sft$powerEstimate

pdf("3.power值散点图.pdf", width = 15, height = 10)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main="Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powersToTry, cex=cex1, col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main="Mean connectivity")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powersToTry, cex=cex1, col="red")
dev.off()

if (is.na(soft_power)) {
  warning("pickSoftThreshold 未能选择合适的 powerEstimate，默认使用6。")
  soft_power <- 6
}
steps_history <- c(steps_history, paste("pickSoftThreshold自动选择的 softPower =", soft_power))

# ------------------------------------------------------------
# 构建邻接矩阵与 TOM
# ------------------------------------------------------------
adjacency <- adjacency(dat_expr, power = soft_power, type = "signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
gene_tree <- hclust(as.dist(dissTOM), method = "average")

pdf("4.基因聚类.pdf", width = 15, height = 10)
plot(gene_tree, xlab="", sub="", main="基于 TOM 不相似性下的基因聚类",
     labels=FALSE, hang=0.04)
dev.off()

# ------------------------------------------------------------
# 动态剪切树：识别模块
# ------------------------------------------------------------
dynamicMods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamic_colors <- labels2colors(dynamicMods)

pdf("5.基因模块对应关系.pdf", width = 15, height = 10)
plotDendroAndColors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang=0.03,
                    addGuide = TRUE, guideHang=0.05,
                    main = "基因聚类及模块分配")
dev.off()
steps_history <- c(steps_history, paste("动态剪切树模块识别：minModuleSize =", minModuleSize))

# ------------------------------------------------------------
# 模块合并：合并相似模块
# ------------------------------------------------------------
MEList <- moduleEigengenes(dat_expr, colors = dynamic_colors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

pdf("6.模块聚类.pdf", width = 7, height = 6)
plot(METree, main="模块特征基因聚类", xlab="", sub="")
abline(h = MEDissThres, col="red")
dev.off()

merge <- mergeCloseModules(dat_expr, dynamic_colors, cutHeight = MEDissThres, verbose = 3)
merged_colors <- merge$colors
merged_MEs <- merge$newMEs

pdf("7.合并后的模块.pdf", width = 15, height = 10)
plotDendroAndColors(gene_tree, merged_colors, "Merged dynamic",
                    dendroLabels = FALSE, hang=0.03,
                    addGuide = TRUE, guideHang=0.05,
                    main = "基因聚类及合并后模块")
dev.off()

module_colors <- merged_colors
MEs <- merged_MEs
steps_history <- c(steps_history, paste("相似模块合并：cutHeight =", MEDissThres))

# ------------------------------------------------------------
# 模块-性状相关性分析
# ------------------------------------------------------------
module_trait_cor <- cor(MEs, trait_data, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(dat_expr))

pdf("8.模块与性状热图.pdf", width=10, height=10)
text_matrix <- paste0(signif(module_trait_cor, 2), "\n(",
                      signif(module_trait_pvalue, 1), ")")
dim(text_matrix) <- dim(module_trait_cor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(trait_data),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "模块与性状相关性热图")
dev.off()
steps_history <- c(steps_history, "计算并绘制模块-性状相关性热图。")

# ------------------------------------------------------------
# 基因显著性（GS）计算及候选关键基因筛选
# ------------------------------------------------------------
# 1. 选择感兴趣的性状向量
if (is.null(trait_of_interest)) {
  trait_of_interest <- colnames(trait_data)[1]
  steps_history <- c(steps_history, paste("未指定感兴趣性状，默认选择第一列:", trait_of_interest))
} else {
  if (!trait_of_interest %in% colnames(trait_data)) {
    stop("指定的感兴趣性状列名称在 trait_data 中不存在，请检查 trait_of_interest 参数。")
  }
  steps_history <- c(steps_history, paste("感兴趣性状列选择为:", trait_of_interest))
}

trait_vector <- trait_data[, trait_of_interest]

# 2. 对所有基因计算 GS（与感兴趣性状的相关系数的绝对值）及 p 值
geneTraitCor <- cor(dat_expr, trait_vector, use = "p")
geneGS <- abs(geneTraitCor)
geneGSp <- corPvalueStudent(as.matrix(geneTraitCor), nrow(dat_expr))

# 3. 选择关键模块：选取与感兴趣性状相关性绝对值最高的模块
absModuleCor <- abs(module_trait_cor[, trait_of_interest])
key_module_name <- rownames(module_trait_cor)[which.max(absModuleCor)]
key_module_color <- substring(key_module_name, 3)
steps_history <- c(steps_history,
  paste("关键模块选定为:", key_module_color,
        "（与性状", trait_of_interest, "的相关性=",
        round(module_trait_cor[key_module_name, trait_of_interest], 2), "）"))

# 4. 对该模块内的基因提取 GS 和 p 值
key_module_genes <- colnames(dat_expr)[module_colors == key_module_color]
GS_key <- geneGS[module_colors == key_module_color]
GSp_key <- geneGSp[module_colors == key_module_color]

# 5. 筛选候选基因：GS > GS_threshold 且 p-value < GS_p_threshold
candidate_idx <- which((GS_key > GS_threshold) & (GSp_key < GS_p_threshold))
candidate_genes <- key_module_genes[candidate_idx]

# 输出候选基因信息到文件
if (length(candidate_genes) > 0) {
  candidate_table <- data.frame(
    Gene = candidate_genes,
    GS = round(GS_key[candidate_idx], 3),
    p_value = round(GSp_key[candidate_idx], 3),
    stringsAsFactors = FALSE
  )
  write.table(candidate_table, file = "Candidate_Genes.txt", sep="\t",
              row.names = FALSE, quote = FALSE)
  steps_history <- c(steps_history,
    paste("候选关键基因筛选：在关键模块", key_module_color, 
          "中筛选出", nrow(candidate_table), "个候选基因。"))
  cat("候选关键基因已写入文件 'Candidate_Genes.txt'\n")
} else {
  steps_history <- c(steps_history,
    paste("候选关键基因筛选：在关键模块", key_module_color,
          "中未筛选到满足条件的候选基因。"))
  cat("在关键模块中未筛选到满足 GS > ", GS_threshold,
      " 且 p-value < ", GS_p_threshold, " 的候选基因。\n", sep="")
}

# ------------------------------------------------------------
# 输出各模块基因
# ------------------------------------------------------------
all_modules <- unique(module_colors)
for (mod in all_modules) {
  mod_genes <- colnames(dat_expr)[which(module_colors == mod)]
  out_file <- paste0("Module_", mod, ".txt")
  write.table(mod_genes, file = out_file, sep="\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
cat("WGCNA 分析完成，共识别到", length(all_modules), "个模块。\n")
steps_history <- c(steps_history,
                   paste0("WGCNA 分析完成，共识别到", length(all_modules), "个模块。"))

# ------------------------------------------------------------
# 将数据处理流程记录写入中文日志文件
# ------------------------------------------------------------
steps_history <- c(
  steps_history,
  "注意：若数据类型或预处理与实际不符，请在脚本开头修改 data_type & rna_seq_method。",
  "若候选关键基因筛选结果不理想，请调整 GS 与 p-value 的筛选阈值。",
  paste("保留基因采用的是表达量标准差排名前", topVarPercent, "%的筛选方法。")
)
writeLines(steps_history, con = "数据处理流程.txt")
cat("已将本次数据处理流程写入 '数据处理流程.txt'。\n")

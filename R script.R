install.packages(c("dplyr", "ggplot2", "readr", "BiocManager"))
BiocManager::install(c( 
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "matrixStats",
  "clusterProfiler",    
  "enrichplot",        
  "AnnotationDbi",      
  "DOSE",
  "biomaRt"
))

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(readr)
library(pheatmap)
library(matrixStats)
library(DOSE)
library(biomaRt)

countData <- read.table("counts_2.txt", header = TRUE, skip = 1, row.names = 1) # в зависимости от нахождения файла путь может отличаться 
# skip=1 пропускает первую строку с метаданными, header = TRUE - первый ряд - имена в таблице, row.names = 1 - первый столбец - названия генов 
countData <- countData[, 6:ncol(countData)] # обрезает ненужные для анализа данные
counts_data <- as.matrix(countData) # преобразовывает данные в матрицу
sampleNames <- colnames(counts_data) # создает отдельную переменную-вектор - имена колонок (а именно реплики и их порядковые номера в эксперименте и контроле)
condition <- factor(c("stress", "stress", "stress", "control", "control", "control")) # создаем вектор с информацией о состоянии 6 образцов - подвергшиеся стрессу и контроль
col_data <- data.frame(row.names = sampleNames, condition = condition) # связываем вектор имен колонок с состоянием, чтобы было проще отследить, где находится контроль, где нет
filtered_counts <- counts_data[rowSums(counts_data) >= 10, ] # фильтруем гены, у которых очень малые значения экспрессии, тем самым избавляясь от шума в статистических данных

dds <- DESeqDataSetFromMatrix( # создаем объект для анализа дифф. экспрессии с помощью DESeq2
  countData = filtered_counts, # выбираем отфильтрованную матрицу
  colData = col_data, # информация об образцах
  design = ~ condition # сравниваем образцы по переменной condition (стресс vs контроль)
)

dds <- DESeq(dds) # запуск анализа
res <- results(dds, contrast = c("condition", "stress", "control")) # Извлечение результатов: сравниваем "stress" против "control"

res_df <- as.data.frame(res) # преобразоание в data.frame
res_df$gene <- rownames(res_df) # добавление столбца с именами генов (имена строк)
deg <- subset(res_df, # фильтрация дифференциально экспрессированных генов
              !is.na(padj) & # избавление от NA
              padj < 0.05 & # p value должно быть меньше 0.05 для статистической значимости
              abs(log2FoldChange) > 1 # отбор биологически значимых изменений
              ) 

plotMA(res, ylim = c(-4, 4)) # построение графика зависимости логарифмического изменения экспрессии (M) от средней экспрессии (A)

vsd <- vst(dds, blind = FALSE) # стабилизирует дисперсию
topvar <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30) # выборка 30 генов с наибольшей вариабельностью
pheatmap( # построение heatmap
  assay(vsd)[topvar, ],
  scale = "row", # нормализация по строкам
  annotation_col = as.data.frame(colData(vsd))["condition", drop = FALSE], # добавляет полоску stress vs control, drop = FALSE — сохраняет формат датафрейм
  main = "Heatmap of top 30 variable genes", # название графика
  fontsize = 8 #размер шрифта
)

write.csv(res_df, "deseq2.csv", row.names = FALSE) # сохранение результатов в виде csv файлов, row.names = FALSE — не записывать имена строк как отдельную колонку
write.csv(deg, "deseq2_deg.csv", row.names = FALSE)
saveRDS(dds, "dds.rds") # сохранение результатов в бинарном виде

connect_ensembl <- function() { # при работе с biomart возникают проблемы с подключением, 
                # поэтому была создана отдельная функция которая пытается установить подключение с сервисом несколько раз
  
  mirrors <- c("www", "useast", "asia") # зеркала ensembl
  
  repeat {
    
    for (m in mirrors) {
      
      message(sprintf("[%s] Trying Ensembl mirror: %s",
                      Sys.time(), m))
      
      hamster <- tryCatch(
        useEnsembl(
          biomart = "genes",
          dataset = "cgpicr_gene_ensembl",
          mirror = m
        ),
        error = function(e) {
          message(sprintf("Failed (%s): %s", m, e$message))
          NULL
        }
      )
      
      if (!is.null(hamster)) {
        message(sprintf("Connected successfully via %s", m))
        return(hamster)
      }
    }
    Sys.sleep(3)
  }
}

hamster <- connect_ensembl()

mapping <- getBM( # мэппинг генов
  attributes = c('external_gene_name', 'entrezgene_id', 'ensembl_gene_id'),
  filters = 'external_gene_name',
  values = res_df$gene, 
  mart = hamster
)

res_df$ENTREZID <- mapping$entrezgene_id[match(res_df$gene, mapping$external_gene_name)] # связывание полученных ID генов с их названием в исходной переменной

# аналогичные операции проводим с deg переменной

mapping <- getBM(
  attributes = c('external_gene_name', 'entrezgene_id', 'ensembl_gene_id'),
  filters = 'external_gene_name',
  values = deg$gene, 
  mart = hamster
)

deg$ENTREZID <- mapping$entrezgene_id[match(deg$gene, mapping$external_gene_name)]

res_df <- res_df %>% filter(!is.na(ENTREZID)) # удаление всех пустых значений NA 
deg <- deg %>% filter(!is.na(ENTREZID))


# Формируем фоновый набор (все протестированные гены)
gene_universe <- unique(res_df$ENTREZID) # оставляет уникальные значения, убирая дубликаты
gene_deg <- unique(deg$ENTREZID)

gene_up <- deg %>% filter(log2FoldChange > 0) %>% pull(ENTREZID) %>% unique() # извлечение генов с повышенной экспрессии
gene_down <- deg %>% filter(log2FoldChange < 0) %>% pull(ENTREZID) %>% unique()

ekegg <- enrichKEGG( # KEGG enrichment
  gene = gene_deg,
  universe = gene_universe,
  organism = "cge", # выбор организма
  keyType = "kegg",  # идентификатор входных генов
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
)

write.csv(as.data.frame(ekegg), "KEGG_results.csv", row.names = FALSE) # таблица с результатами KEGG enrichment

saveRDS(ekegg, "ekegg.rds") # бинарный файл с результатами

######создание графиков#####

if(nrow(ekegg) > 0) {
  png("/home/bumblebee/rstudio/lemonhead/KEGG_cnetplot.png", width = 1800, height = 1400, res = 200)
  print(cnetplot(ekegg, showCategory = 10))
  dev.off()
}

if(nrow(ekegg) > 0) {
  png("/home/bumblebee/rstudio/lemonhead/KEGG_dotplot.png", width = 1200, height = 800, res = 150)
  print(dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment"))
  dev.off()
}

if(nrow(ekegg) > 0) {
  png("/home/bumblebee/rstudio/lemonhead/KEGG_barplot.png", width = 1200, height = 800, res = 150)
  print(barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment"))
  dev.off()
}






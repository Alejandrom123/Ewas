# Anotación de genes con BLAST (base de datos Drosophilidae + Culicidae) y resumen de calidad.
#
# 1) Comando de BLAST (ejecutado en el servidor, una vez por población):
#    blastp -num_threads 30 -db drosophila_culicidae_II_db \
#           -query gene_universe_completo_FASTA_HLvSP.fasta \
#           -max_target_seqs 1 -evalue 1e-6 \
#           -outfmt "10 qseqid sseqid sacc pident length mismatch qstart qend sstart send evalue bitscore salltitles" \
#           > resultadoblast_HLvSP.csv
#    BLAST+ v2.15.0. La base de datos se construyó con makeblastdb a partir de las proteínas de
#    Drosophilidae (txid 7214) y Culicidae (txid 7157).
#
# 2) Este script: calcula identidad, cobertura y e-value por proteína, completa la descripción de los
#    genes que VectorBase anota como "unspecified product", y resume las categorías asignadas manualmente.

library(openxlsx)

# ---- Archivos de entrada (modificar rutas) ----
blast_csv   <- "resultadoblast_HLvSP.csv"                 # salida de blastp (outfmt 10, sin encabezado)
query_fasta <- "gene_universe_completo_FASTA_HLvSP.fasta" # proteínas consultadas en Vectorbase
vectorbase  <- "LOC_genes_HLvSP.xlsx"                     # VectorBase: Gene ID, Product Description
categorias  <- "ACompleto_biomart_final_HLvSP_simple2026.xlsx" # columnas LOC y description_IR (manual)
salida      <- "blast_annotation_HL.csv"

# ---- 1. Resultados de BLAST ----
cols  <- c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore", "salltitles")
blast <- read.csv(blast_csv, header = FALSE, fill = TRUE, col.names = cols, sep = ",")[, 1:13]
blast <- blast[!is.na(suppressWarnings(as.numeric(blast$pident))), ]   # quita avisos de BLAST dentro del archivo
blast[c("pident", "length", "qstart", "qend", "evalue", "bitscore")] <-
  lapply(blast[c("pident", "length", "qstart", "qend", "evalue", "bitscore")], as.numeric)

# ---- 2. Longitud de cada proteína consultada y cobertura ----
fa  <- readLines(query_fasta)
id  <- sub("^>(\\S+).*", "\\1", fa[startsWith(fa, ">")])
len <- tapply(nchar(fa[!startsWith(fa, ">")]), cumsum(startsWith(fa, ">"))[!startsWith(fa, ">")], sum)
qlen <- setNames(as.numeric(len), id)

blast$qlength  <- qlen[blast$qseqid]
blast$coverage <- round((blast$qend - blast$qstart + 1) / blast$qlength * 100, 1)   # % de la proteína consulta alineada

# Un registro por gen (primer hit del primer transcrito), quitando el sufijo de isoforma (-RA, -RB...)
blast$Gene.ID <- sub("-R.$", "", blast$qseqid)
blast$salltitles <- trimws(sub("\\s*\\[.*$", "", blast$salltitles))
blast <- blast[!duplicated(blast$Gene.ID), ]

# ---- 3. Descripción final: VectorBase y, si es "unspecified product", BLAST ----
vb <- read.xlsx(vectorbase)
names(vb) <- make.names(names(vb))
ann <- merge(vb[, c("Gene.ID", "Product.Description")], blast, by = "Gene.ID", all.x = TRUE)

sin_desc <- is.na(ann$Product.Description) | ann$Product.Description == "unspecified product"
ann$source <- ifelse(!sin_desc, "VectorBase", ifelse(!is.na(ann$salltitles), "BLAST", "none"))
ann$description_final <- ifelse(ann$source == "BLAST", ann$salltitles,
                         ifelse(ann$source == "VectorBase", sub("\\s*\\[.*$", "", ann$Product.Description), NA))

# ---- 4. Categoría de resistencia asignada manualmente (description_IR) ----
cat <- read.xlsx(categorias)
cat <- cat[!is.na(cat$description_IR), ]
cat <- cat[!duplicated(cat$LOC), c("LOC", "description_IR")]
ann <- merge(ann, cat, by.x = "Gene.ID", by.y = "LOC", all.x = TRUE)

# ---- 5. Resumen ----
con_hit <- !is.na(ann$pident)
cat("Genes:", nrow(ann), "| con hit BLAST:", sum(con_hit), "\n")
cat("Fuente de la descripción:\n"); print(table(ann$source))
cat("Identidad (%): rango", range(ann$pident[con_hit]), "| mediana", median(ann$pident[con_hit]), "\n")
cat("Cobertura (%): rango", range(ann$coverage[con_hit], na.rm = TRUE), "| mediana", median(ann$coverage[con_hit], na.rm = TRUE), "\n")
cat("E-value máximo:", max(ann$evalue[con_hit]), "\n")
cat("Genes por categoría:\n"); print(table(ann$description_IR))

write.csv(ann[, c("Gene.ID", "Product.Description", "qseqid", "sacc", "pident", "length", "qlength",
                  "coverage", "evalue", "bitscore", "salltitles", "source", "description_final", "description_IR")],
          salida, row.names = FALSE)

# build_kinase_families.R
# Builds Kinase_Families.tsv: Manning-style group / family for every Kinase Library
# (TKL) matrix name shipped with OSMSS, from openly licensed (CC BY 4.0) sources only:
#   [J23] Johnson JL et al. (2023) Nature 613:759-766, Supplementary Table 1
#         (doi:10.1038/s41586-022-05575-3), CC BY 4.0
#   [Y24] Yaron-Barir TM et al. (2024) Nature 629:1174-1181, Supplementary Table 1
#         (doi:10.1038/s41586-024-07407-y), CC BY 4.0
#   [UP]  UniProt Consortium: pkinfam.txt ("Human and mouse protein kinases:
#         classification and index") and UniProtKB "protein_families" field, CC BY 4.0
# Changes made to the sources (CC BY "indicate changes"): J23 groups Alpha, PDHK, PIKK,
# FAM20 are filed under Group "Atypical" (Manning 2002 convention) and kept as Family;
# Y24 canonical Tyr kinases (RTK / nRTK) get Group "TK"; non-canonical Tyr activities
# (ncTK, matrix names *_TYR) get the Ser/Thr group of the same protein from J23/UniProt.
#
# Usage (from the OSMSS app folder): Rscript data/Kinase_Families_build.R <work_dir> <osmss_app_dir>
# Maintainer tool only; the app reads the finished data/Kinase_Families.tsv.

suppressPackageStartupMessages({ library(readxl); library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
SC  <- if (length(args) >= 1) args[1] else file.path(tempdir(), "kinase_families_build")
APP <- if (length(args) >= 2) args[2] else "."

# ---- 0. Downloads (official sources only; skipped if already present) -------------
dl <- function(url, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(dest)) download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}
f_j23 <- dl("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM3_ESM.xlsx",
            file.path(SC, "download_johnson2023", "st1.xlsx"))
f_y24 <- dl("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM3_ESM.xlsx",
            file.path(SC, "download_yaronbarir2024", "st1.xlsx"))
f_pk  <- dl("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pkinfam.txt",
            file.path(SC, "download_uniprot", "pkinfam.txt"))

# ---- 1. Kinase names the app needs --------------------------------------------------
st_names <- read.delim(file.path(APP, "data/matrices/ser_thr/ser_thr_all_log2_matrices.tsv"),
                       check.names = FALSE)[, 1]
ty_names <- read.delim(file.path(APP, "data/matrices/tyrosine/tyrosine_all_log2_matrices.tsv"),
                       check.names = FALSE)[, 1]

# ---- 2. Sources -----------------------------------------------------------------------
j23 <- as.data.frame(read_excel(f_j23, sheet = "Table S1 Data"))
j23 <- data.frame(Name = toupper(j23$Matrix_name), Gene = j23$Gene, UniProt = j23$`Uniprot id`,
                  J_Group = j23$Family, stringsAsFactors = FALSE)

y24 <- as.data.frame(read_excel(f_y24, sheet = "Table_S1_Data"))
y24 <- y24[!is.na(y24$SUBTYPE), ]                      # 93 human kinases (drops worm + mutants)
y24 <- data.frame(Name = toupper(y24$MATRIX_NAME), Gene = y24$GENENAME, UniProt = y24$UNIPROT_ID,
                  Subtype = y24$SUBTYPE, Y_Family = y24$FAMILY, stringsAsFactors = FALSE)

# pkinfam.txt: section headers sit between two lines of '='
pk_lines <- readLines(f_pk, warn = FALSE)
is_rule  <- grepl("^=+$", pk_lines)
hdr_idx  <- which(c(FALSE, is_rule[-length(is_rule)]) & c(is_rule[-1], FALSE) & !is_rule)
pk <- do.call(rbind, lapply(seq_along(hdr_idx), function(i) {
  h <- pk_lines[hdr_idx[i]]
  end <- if (i < length(hdr_idx)) hdr_idx[i + 1] - 2 else length(pk_lines)
  body <- pk_lines[(hdr_idx[i] + 2):end]
  body <- body[grepl("_HUMAN", body)]
  if (!length(body)) return(NULL)
  grp <- if (grepl("^Atypical", h)) "Atypical" else if (h == "Other") "Other" else
         if (grepl("^Tyr ", h)) "TK" else word(h, 1)
  fam <- if (grepl("^Atypical: ", h)) str_match(h, "^Atypical: (.+?) (Ser/Thr |protein )?(kinase )?family")[, 2] else NA
  data.frame(Gene = toupper(word(body, 1)),
             UniProt = str_match(body, "_HUMAN\\s*\\(([A-Z0-9]+)")[, 2],
             UP_Group = grp, UP_AtypFamily = fam, stringsAsFactors = FALSE)
}))

# UniProtKB protein_families for every accession involved (cached to disk)
f_upfam <- file.path(SC, "download_uniprot", "uniprotkb_protein_families.tsv")
acc_all <- unique(na.omit(c(j23$UniProt, y24$UniProt, pk$UniProt)))
if (!file.exists(f_upfam)) {
  chunks <- split(acc_all, ceiling(seq_along(acc_all) / 100))
  res <- lapply(chunks, function(a) {
    u <- paste0("https://rest.uniprot.org/uniprotkb/accessions?format=tsv&fields=accession,gene_primary,protein_families&accessions=",
                paste(a, collapse = ","))
    read.delim(url(u), stringsAsFactors = FALSE, check.names = FALSE)
  })
  write.table(do.call(rbind, res), f_upfam, sep = "\t", quote = FALSE, row.names = FALSE)
}
upfam <- read.delim(f_upfam, stringsAsFactors = FALSE, check.names = FALSE)
names(upfam) <- c("UniProt", "UP_Gene", "UP_Families")
up_sub <- function(acc) {
  s <- upfam$UP_Families[match(acc, upfam$UniProt)]
  sub <- str_match(s, "(?:^|, )([^,]+?) subfamily")[, 2]
  ifelse(is.na(sub), "", sub)
}

# ---- 3. Ser/Thr -------------------------------------------------------------------
ATYP <- c(Alpha = "Alpha", PDHK = "PDHK", PIKK = "PIKK", FAM20 = "FAM20")
# TKL names absent from J23 (added to the Kinase Library later): TKL name -> HGNC symbol
st_alias <- c(ALK1 = "ACVRL1", ALK7 = "ACVR1C")

st <- data.frame(Kinase = toupper(st_names), Type = "Ser/Thr", stringsAsFactors = FALSE)
m  <- match(st$Kinase, j23$Name)
st$Gene <- j23$Gene[m]; st$UniProt <- j23$UniProt[m]
st$Group <- j23$J_Group[m]; st$Source <- ifelse(is.na(m), NA, "Johnson2023_ST1")
miss <- is.na(m)
g  <- ifelse(st$Kinase[miss] %in% names(st_alias), st_alias[st$Kinase[miss]], st$Kinase[miss])
mp <- match(g, pk$Gene)
st$Gene[miss] <- g; st$UniProt[miss] <- pk$UniProt[mp]; st$Group[miss] <- pk$UP_Group[mp]
st$Source[miss] <- ifelse(is.na(mp), NA, "UniProt_pkinfam")
st$Family <- ""
at <- st$Group %in% names(ATYP)
st$Family[at] <- ATYP[st$Group[at]]; st$Group[at] <- "Atypical"
st$Group[st$Group == "NEK"] <- "Other"     # Manning 2002 files NEK under Other

# ---- 4. Tyr -----------------------------------------------------------------------
ty_alias <- c(MST1R = "RON")
ty <- data.frame(Kinase = toupper(ty_names), Type = "Tyr", stringsAsFactors = FALSE)
key <- sub("_TYR$", "", ty$Kinase); key <- ifelse(key %in% names(ty_alias), ty_alias[key], key)
m <- match(key, y24$Name)
ty$Gene <- y24$Gene[m]; ty$UniProt <- y24$UniProt[m]
canon <- y24$Subtype[m] %in% c("RTK", "nRTK")
ty$Type[!canon] <- "Tyr (non-canonical)"
ty$Group  <- ifelse(canon, "TK", NA)
ty$Family <- ifelse(canon, y24$Y_Family[m], "")
ty$Source <- ifelse(is.na(m), NA, "YaronBarir2024_ST1")
# non-canonical: take the Ser/Thr group of the same protein (J23 by accession, else pkinfam)
nc <- which(!canon & !is.na(m))
# match by accession, then by gene symbol (Y24 PINK1 is D6WMX4, a non-human ortholog)
jg <- j23$J_Group[match(ty$UniProt[nc], j23$UniProt)]
jg <- ifelse(is.na(jg), j23$J_Group[match(toupper(ty$Gene[nc]), toupper(j23$Gene))], jg)
pg <- pk$UP_Group[match(ty$UniProt[nc], pk$UniProt)]
pg <- ifelse(is.na(pg), pk$UP_Group[match(toupper(ty$Gene[nc]), pk$Gene)], pg)
grp <- ifelse(!is.na(jg), jg, pg)
ty$Family[nc] <- ifelse(grp %in% names(ATYP), ATYP[grp], "")
ty$Group[nc]  <- ifelse(grp %in% names(ATYP), "Atypical", ifelse(grp == "NEK", "Other", grp))
ty$Source[nc] <- paste0("YaronBarir2024_ST1+", ifelse(!is.na(jg), "Johnson2023_ST1", "UniProt_pkinfam"))

# ---- 5. Family / subfamily detail from UniProtKB (where not already set) ---------------
out <- rbind(st, ty)
sub_up <- up_sub(out$UniProt)
need <- out$Family == "" & sub_up != ""
out$Family[need] <- sub_up[need]
out$Source[need] <- paste0(out$Source[need], "+UniProtKB_family")
out$Subfamily <- ""   # no openly licensed 3rd level found; kept for schema stability
out <- out[, c("Kinase", "Type", "Group", "Family", "Subfamily", "Gene", "UniProt", "Source")]
out[is.na(out)] <- ""

# ---- 5b. Corrections: J23 group outvoted by UniProt pkinfam AND Manning 2002 ------------
# (Manning groups as distributed in PhosR::KinaseFamily.) HUNK: J23 CMGC -> CAMK;
# DSTYK: J23 TKL -> Other. Rule: majority of the three references.
fix <- c(HUNK = "CAMK", DSTYK = "Other")
hit <- out$Kinase %in% names(fix)
out$Group[hit]  <- fix[out$Kinase[hit]]
out$Source[hit] <- paste0(out$Source[hit], "+corrected_UniProt_Manning2002")

write.table(out, file.path(SC, "Kinase_Families.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 6. Sanity: J23 group vs UniProt pkinfam group, by accession ------------------------
jc <- merge(j23, pk[, c("UniProt", "UP_Group")], by = "UniProt")
jc$J_norm <- ifelse(jc$J_Group %in% names(ATYP), "Atypical", jc$J_Group)
jc$P_norm <- ifelse(jc$UP_Group == "NEK", "Other", jc$UP_Group)
dis <- jc[jc$J_norm != jc$P_norm, c("Name", "Gene", "J_Group", "UP_Group")]
cat("J23 vs pkinfam group disagreements (", nrow(dis), "of", nrow(jc), "):\n"); print(dis, row.names = FALSE)
cat("\nUnassigned rows:\n"); print(out[out$Group == "", ], row.names = FALSE)
cat("\nWrote", nrow(out), "rows to", file.path(SC, "Kinase_Families.tsv"), "\n")

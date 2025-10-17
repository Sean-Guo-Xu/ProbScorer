
# Handling tree structure
parse_brite_to_df <- function(lines, max_level = 3) {
  currentA <- NA_character_
  currentB <- NA_character_
  rows <- list()

  for (ln in lines) {
    if (nchar(trimws(ln)) == 0) next
    t <- sub("^\\s+", "", ln)        # Remove the leading spaces in each line.
    if (nchar(t) == 0) next

    first <- substr(t, 1, 1)
    # Only process the lines that start with A/B/C
    if (!first %in% c("A", "B", "C")) next

    if (first == "A") {
      currentA <- trimws(sub("^A\\s*", "", t))
      currentB <- NA_character_
    } else if (first == "B") {
      currentB <- trimws(sub("^B\\s*", "", t))
    } else if (first == "C") {
      # Remove the initial "C" and the spaces, and attempt to extract "code + name"
      rest <- trimws(sub("^C\\s*", "", t))
      m <- regexec("^([0-9]+)\\s+(.+)$", rest)
      mm <- regmatches(rest, m)[[1]]
      if (length(mm) >= 3) {
        code <- mm[2]
        cname <- mm[3]
      } else {
        # does not follow the "code + name" format, entire "rest" regarded as the name.
        code <- NA_character_
        cname <- rest
      }

      rows[[length(rows) + 1]] <- data.frame(
        A = ifelse(is.na(currentA), NA_character_, currentA),
        B = ifelse(is.na(currentB), NA_character_, currentB),
        C_code = code,
        C_name = cname,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(data.frame(A = character(0), B = character(0),
                      C_code = character(0), C_name = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out <- unique(out)   # Remove duplicates
  rownames(out) <- NULL
  return(out)
}



# auxiliary function: Convert entrezid to symbol and standardize the weights
convert_entrez_to_symbol <- function(weight_vec) {
  entrez_ids <- sub("ENTREZID:", "", names(weight_vec))

  mapping <- bitr(entrez_ids, fromType = "ENTREZID",
                  toType = "SYMBOL", OrgDb = org.Hs.eg.db)

  df <- data.frame(ENTREZID = entrez_ids, weight = weight_vec,
                   stringsAsFactors = FALSE)

  df <- inner_join(mapping, df, by = "ENTREZID")

  # divided equally when an entrezid corresponds to multiple symbols
  df <- df %>%
    group_by(ENTREZID) %>%
    mutate(weight = weight / n()) %>%
    ungroup()

  # Aggregate and accumulate by symbol
  df <- df %>%
    group_by(SYMBOL) %>%
    summarise(weight = sum(weight), .groups = "drop")

  # scale: Normalize to sum to 1
  df$weight <- df$weight / sum(df$weight)

  # Convert to named vectors
  weight_vec_symbol <- setNames(df$weight, df$SYMBOL)
  return(weight_vec_symbol)
}

# main
get_pathway_weights <- function(PathwaySet = "Metabolism") {
  # 1. Obtain human KEGG pathways
  all_human_pathways <- pathways("hsapiens", "kegg")

  # 2. Obtain KEGG BRITE hierarchical information
  human_brite <- keggGet("br:br08901")
  lines <- unlist(strsplit(human_brite, "\n"))
  df <- parse_brite_to_df(lines)   # 你已有的函数

  # 3. Select PathwaySet
  ids <- df[df$A %in% PathwaySet,  ]
  ids$id <- paste("hsa:", ids[,3], sep = "")

  # 4. Filtering pathway
  all_ids <- vapply(all_human_pathways, function(pw) pw@id, character(1))
  pathways <- all_human_pathways[all_ids %in% ids$id]

  # 5. Calculate the PageRank score for each path
  weighs <- list()
  for (i in seq_along(pathways)) {
    pathway <- pathways[[i]]
    if (length(pathway) == 0) next
    pathway_graph <- pathwayGraph(pathway)
    pathway_graph <- igraph::igraph.from.graphNEL(pathway_graph)
    score <- page.rank(pathway_graph)$vector
    if (length(score) == 0) next
    suppressMessages(  weighs[[names(pathways)[[i]]]] <- convert_entrez_to_symbol(score) )
  }

  return(weighs)
}

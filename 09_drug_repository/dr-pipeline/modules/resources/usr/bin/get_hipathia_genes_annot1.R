#!/usr/bin/env Rscript

# Purpose: extract information relating circuits and genes in HiPathia.
# No inputs required.
# Output 1: rds file with the metaginfo object from pathways
# Output 2: tsv file with extracted pathway/circuit/gene information

suppressMessages(library(hipathia))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))

# load resources

message("Loading HiPathia resources...")

hip_version <- packageVersion("hipathia")
pathways <- load_pathways("hsa")


### circuit/node/gene associations and node effect over circuit effector ###
# Activating/inhibitory signaling of nodes in RELation to EFFector is reflected in rel_eff

message("Extracting circuit-node-gene associations...")

rel_res <- tibble()

for (pathway in names(pathways$pathigraphs)) {
  for (circuit in names(pathways$pathigraphs[[pathway]]$effector.subgraphs)) {
    g <- pathways$pathigraphs[[pathway]]$effector.subgraphs[[circuit]]
    ## get table of relations between nodes
    rel_tab <- as_edgelist(g) %>% as.data.frame()
    names(rel_tab) <- c("tail", "head")
    rel_tab$relation <- E(g)$relation
    ## add effector row with column for relation to head node
    end_n <- names(V(g)[degree(g, mode="out") == 0])
    rel_tab <- full_join(tibble(tail = end_n, head = end_n, relation = 1, 
                                rel_head = 1), 
                         rel_tab, by = c("tail", "head", "relation"))
    ## get node relations over effector node recursively
    # Table rel_node will store results of calculated relations over effectors.
    # A temporary table will store previous results to check whether there is 
    # change between new relations table and previous one at each step.
    rel_node <- filter(rel_tab, !is.na(rel_head)) %>% select(head, rel_eff = rel_head)
    rel_prev <- tibble(head = "N", rel_eff = 1)
    while (!isTRUE(all.equal(rel_node, rel_prev))) {
      rel_prev <- rel_node
      new_relations <- full_join(rel_tab, rel_node, by = "head") %>% 
        filter(!is.na(tail)) %>% 
        .[.$head %in% rel_node$head & is.na(.$rel_head),] %>% 
        mutate(rel_head = relation * rel_eff) %>% 
        select(head = tail, rel_eff = rel_head)
      rel_node <- rbind(rel_node, new_relations) %>% unique()
      # if there are duplicates with different values, set value to 0
      dupes <- unique(rel_node$head[duplicated(rel_node$head)])
      if (length(dupes) > 0) {
        rel_node <- filter(rel_node, !head %in% dupes) %>% 
          rbind(tibble(head = dupes, rel_eff = 0))
      }
    }
    ## get gene ids for nodes
    l <- lapply(rel_node$head, function(x) V(g)[[x]]$genesList)
    names(l) <- rel_node$head
    # create table of genes per node
    l <- lapply(l, tibble::as_tibble_col, column_name = "gene") %>% 
      dplyr::bind_rows(.id = "node") %>% 
      dplyr::mutate(gene = suppressWarnings(as.numeric(gene))) %>%
      dplyr::filter(!is.na(gene)) %>% distinct()
    ## join results and add circuit data
    rel_node$node_name <- unlist(lapply(rel_node$head, function(z) V(g)[[z]]$label))
    node_relations <- inner_join(l, rel_node, by = c("node" = "head")) %>% 
      tibble::add_column("pathway" = pathway, "circuit" = circuit, 
                         circuit_name = get_path_names(pathways, circuit), .before = 1) %>% 
      dplyr::relocate(node_name, rel_eff, .after = "node")
    rel_res <- rbind(rel_res, node_relations)
  }
}


### assign last nodes with genes as circuit effectors ###

message("Selecting protein nodes as circuit effectors...")

circ_eff <- tibble()

for (pathway in names(pathways[["pathigraphs"]])) {
  circuits <- names(pathways[["pathigraphs"]][[pathway]][["effector.subgraphs"]])
  for (c in circuits) {
    g <- pathways[["pathigraphs"]][[pathway]][["effector.subgraphs"]][[c]]
    # circuit effector node
    n <- V(g)[degree(g, mode="out") == 0]
    ## in case the final node only contains a metabolite, search for genes in previous nodes
    i <- 1
    while (all(is.na(n$genesList)) & i < length(V(g))){
      n <- V(g)[adjacent_vertices(g, n, mode = c("in"))[[1]]]
      i <- i + 1
    }
    # get effector node(s) per circuit
    circ_eff <- dplyr::bind_rows(circ_eff, tibble("circuit" = c, node = names(n)))
  }
}
circ_eff$eff_node <- TRUE

# merge circuit/gene information with effector nodes

circ_genes <- rel_res %>% left_join(circ_eff, by = c("circuit", "node"))


### write rds of metaginfo and the table of hipathia circuit information ###

message("Saving pathway information...")
saveRDS(pathways, paste0("pathways-hipathia", hip_version, ".rds"))
out_name <- paste0("circuitGenes_hipathia", hip_version, ".tsv")
write_tsv(circ_genes, out_name, na = "")

message("Extraction of information from HiPathia pathways done.")

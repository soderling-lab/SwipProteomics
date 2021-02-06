#!/usr/bin/env Rscript

# title: SwipProteomics
# description: generate cytoscape networks for each module
# author: Tyler W Bradshaw

## ---- inputs

## ---- functions

is_connected <- function(graph, threshold) {
  # a helper function that checks if graph is connnected at a 
  # given edge weight threshold
  to_drop <- which(E(graph)$weight <= threshold)
  thresh_graph <- igraph::delete.edges(graph, to_drop)
  return(igraph::is.connected(thresh_graph))
}


mapParam <- function(visual_params) {
  # a helper function for RCy3::mapVisualProperty
  mapped_param <- do.call(RCy3::mapVisualProperty, visual_params)
  return(mapped_param)
}


myfun <- function(graph, threshold, pclust=1.0, too_small=1) {
  # just iterate through cutoffs and calc n_compoents, set threshold for
  # percent clustered to be high -- 95% 
	to_drop <- which(E(graph)$weight <= threshold)
	thresh_graph <- igraph::delete.edges(graph, to_drop)
	graph_comp <- igraph::components(thresh_graph)
	partition <- graph_comp$membership
	ncomp <- length(unique(partition))
	m_sizes <- table(partition)
	to_rm <- as.numeric(names(which(m_sizes <= too_small)))
	if (length(to_rm != 0)) {
		# set membership to 0
		partition[partition %in% to_rm] <- 0
	}
	percent_clust <- 1-sum(partition==0)/length(partition)
	return(percent_clust >= pclust)
} #EOF


createCytoscapeGraph <- function(module, netw_g, ppi_g, nodes) {
  # create a cytoscape graph of each module 

	## define Cytoscape layout
	netw_layout='force-directed edgeAttribute=weight'

	## subset graph: keep defined nodes
	graph <- netw_g
	idx <- match(nodes, names(V(graph)))
	if (sum(is.na(idx))>0) { warning("missing nodes") }
	vids <- idx[!is.na(idx)]
	subg <- igraph::induced_subgraph(graph, vids)

	## set node size ~ hubbiness or importance in its subgraph
	adjm <- as.matrix(as_adjacency_matrix(subg, attr="weight"))
	namen <- names(V(subg))
	node_importance <- apply(adjm,2,sum)
	subg <- igraph::set_vertex_attr(subg, "size", value=node_importance[namen])

	## Seq from min(edge weight) to max() to generate cutoffs
	#n_edges <- length(E(subg))
	min_weight <- min(E(subg)$weight)
	max_weight <- max(E(subg)$weight)
	cutoffs <- seq(min_weight, max_weight, length.out = 5000)

	## find a suitable cutoff for thresholding the graph
	doParallel::registerDoParallel(parallel::detectCores() - 1)
	is_single_component <- foreach(threshold=cutoffs) %dopar% {
		myfun(subg, threshold, pclust = 0.95) 
	} %>% unlist()

	## define limit as max(cutoff) at which the graph is still connected
	if (all(is_single_component)) { stop("Error thesholding graph.") }
	weight_limit <- cutoffs[max(which(is_single_component == TRUE))]

	## prune edges
	# NOTE: This removes all edge types connecting two nodes
	g <- delete.edges(subg, which(E(subg)$weight <= weight_limit))

	## write graph to file
	# NOTE: This is faster than sending to cytoscape via
	myfile <- paste0(module,".gml")

	igraph::write_graph(g, myfile, format = "gml")

	stopifnot(file.exists(myfile))

	# send to Cytoscape
	# NOTE: underscores in attribute names are removed
	# NOTE: if you keep getting errors with keys -- make sure they are there!
        # FIXME: this is OS specific!
	winfile <- gsub("/mnt/d/|\\~/", "D:/", file.path(getwd(), myfile))
	cysnetw <- RCy3::importNetworkFromFile(winfile)

	Sys.sleep(2)
	unlink(myfile)

	stopifnot(!file.exists(myfile))

	# for large networks, you man need to manually create a view
        #RCy3::getNetworkViews()
	#RCy3::commandsPOST("view create")

	## visual style defaults
	style.name <- "RCy3.style"
	visual_defaults <- list(
		 NETWORK_TITLE = "RCy3.network",
		 NODE_FILL_COLOR = "#BEBEBE", # gray
		 NODE_TRANSPARENCY = 255,
		 NODE_SIZE = 35,
		 NODE_SHAPE = "ellipse",
		 NODE_LABEL_TRANSPARENCY = 255,
		 NODE_LABEL_FONT_SIZE = 12,
		 NODE_LABEL_COLOR = "#000000", # black
		 NODE_BORDER_TRANSPARENCY = 200,
		 NODE_BORDER_WIDTH = 4,
		 NODE_BORDER_PAINT = "#000000", # black
		 EDGE_STROKE_UNSELECTED_PAINT = "#000000", # black
		 EDGE_WIDTH = 2,
		 NETWORK_BACKGROUND_PAINT = "#FFFFFF" # white
		)

	## mapped visual style parameters

	# some ranges for mapped params
        weight_range <- c(min(E(g)$weight), max(E(g)$weight))
	size_range <- c(min(V(g)$size),max(V(g)$size))
	edge_colors <- c("#BEBEBE", "#8B0000") # gray, 'dark red'
	mapped_params <- list()

	# EDGE TRANSPARENCY
        params <- list("edge transparency", "weight", "c", 
		       weight_range, c(155, 255))
        mapped_params[["EDGE_TRANSPARENCY"]] <- mapParam(params)
	# NODE FILL
	params <- list("node fill color","Color","p")
        mapped_params[["NODE_FILL_COLOR"]] <- mapParam(params)
	# NODE LABEL
        params <- list("node label","Protein", "p")
        mapped_params[["NODE_LABEL"]] <- mapParam(params)
	# EDGE TRANSPARENCY
	params <- list("edge transparency", "weight", "c", 
		       weight_range, c(155, 255))
	mapped_params[["EDGE_TRANSPARENCY"]] <- mapParam(params)
	# EDGE STROKE COLOR
        params <- list("edge stroke unselected paint", "weight", "c", 
		       weight_range, edge_colors)
	mapped_params[["EDGE_STROKE_UNSELECTED_PAINT"]] <- mapParam(params)
	# NODE SIZE
	params <- list("node size", "size", "c", size_range, c(35, 100))
	mapped_params[["NODE_SIZE"]] <- mapParam(params)

	## create a visual style
	RCy3::createVisualStyle(style.name, defaults = visual_defaults, 
			  mappings = mapped_params)

	# apply to graph
	RCy3::setVisualStyle(style.name)
	Sys.sleep(3)

	## Collect PPI edges
	idx <- match(nodes, names(V(ppi_g)))
		subg <- igraph::induced_subgraph(ppi_g, 
					 vids = V(ppi_g)[idx])
	edge_list <- apply(igraph::as_edgelist(subg, names = TRUE), 1, as.list)

	## If edge list is only of length 1, unnest it to avoid problems
	if (length(edge_list) == 1) {
		edge_list <- unlist(edge_list, recursive = FALSE)
	}

	## we add edges like this bc Cytoscape really only supports 1 type of
	## edge in a graph. We manually add additional PPI edges to create a
	## network with both co-variation and ppi edges.
	## NOTE: this can take several minutes if there are a large number of PPIs
	if (length(edge_list) > 0) {
		ppi_edges <- RCy3::addCyEdges(edge_list)
		# add PPIs and set to black
		selected_edges <- RCy3::selectEdges(ppi_edges, by.col = "SUID")
		# set edge bend to help distinguish ppi edges
		namen <- "EDGE_STROKE_UNSELECTED_PAINT"
		RCy3::setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = "#000000", # black 
				      visual.property = namen,
				      bypass = TRUE
				      )
		RCy3::setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = TRUE,
				      visual.property = "EDGE_BEND",
				      bypass = TRUE
				      )
	} # EIS

	##  clean-up
	RCy3::clearSelection()
	Sys.sleep(2) 

	## apply layout
	RCy3::layoutNetwork(netw_layout)
	Sys.sleep(2)

	## fit to screen
	RCy3::fitContent()
	Sys.sleep(2)

	## Mask color of non-significant nodes
	sig <- names(V(g))[V(g)$sigprot == 1]
	ns <- names(V(g))[names(V(g)) %notin% sig]
	if (length(ns) > 0) {
		RCy3::setNodeColorBypass(ns,new.colors="#BEBEBE")
	}

	## set bold border of BioID proteins
	# FIXME: report error in setNodeBorderWidthBypass 
        # Error in .cyFinally(res) : object 'res' not found
	#wash_nodes <- names(V(g))[V(g)$isWASH==1]
	#if (length(wash_nodes) > 0) {
	#	RCy3::setNodeBorderWidthBypass(wash_nodes,new.sizes=10)
	#}

	## free up some memory
	RCy3::cytoscapeFreeMemory()

} #EOF


## ---- set-up the workspace

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

# global imports
suppressPackageStartupMessages({
  library(RCy3) 
  library(dplyr) 
  library(igraph) 
  library(data.table) 
  library(doParallel)
})


# project functions and data
devtools::load_all(root, quiet=TRUE)

# project directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Output directory for cytoscape networks
netwdir <- file.path(root,"figs","Networks")
if (!dir.exists(netwdir)) {
	dir.create(netwdir)
}


## ---- Load the data in root/data

# Load the data from root/data
data(swip_tmt)
data(wash_interactome)
wash_prots <- wash_interactome

data(swip_partition) # partition
data(swip_gene_map) # gene_map
data(swip_colors) # module_colors

# we label sig_prots from MSstatsTMT
data(msstats_results)
data(msstats_sig_prots)



## ---- Load networks in root/rdata

# adjm
myfile <- file.path(root,"rdata","adjm.rda")
load(myfile) 

# ne_adjm
myfile <- file.path(root,"rdata","ne_adjm.rda")
load(myfile) 

# ppi_adjm
myfile <- file.path(root,"rdata","ppi_adjm.rda")
load(myfile) 


## ---- Create igraph graph to be passed to Cytoscape

# create a list of all modules
modules <- split(names(partition),partition)[-1] # drop M0
names(modules) <- paste0("M",names(modules))

# insure that matrices are in matching order
check <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!check) { stop("input adjacency matrices should be of matching dimensions") }

# create igraph graph objects
# NOTE: graph edge weight is enhanced(pearson.cor)
netw_g <- graph_from_adjacency_matrix(ne_adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
ppi_g <- graph_from_adjacency_matrix(ppi_adjm,mode="undirected",diag=FALSE,
				     weighted=TRUE)


# annotate graphs with protein NAMES
proteins <- toupper(gene_map$symbol[match(names(V(netw_g)),gene_map$uniprot)])
netw_g <- set_vertex_attr(netw_g,"protein",value = proteins)
proteins <- toupper(gene_map$symbol[match(names(V(ppi_g)),gene_map$uniprot)])
ppi_g <- set_vertex_attr(ppi_g,"protein",value = proteins)


## ---- annotate graphs with additional metadata

# collect meta data from msstats_results as noa data.table
tmp_dt <- data.table(Protein = names(V(netw_g)),
		  Module = paste0("M",partition[names(V(netw_g))]))
noa <- left_join(tmp_dt, msstats_results, by = "Protein") %>% 
	filter(Contrast == "Mutant-Control")

# add module colors
noa$Color <- module_colors[noa$Module]

stopifnot(!any(is.na(noa$Color)))

# add WASH iBioID annotation
noa$isWASH <- as.numeric(noa$Protein %in% wash_prots)

# add sig prot annotations
noa$sigprot <- as.numeric(noa$Protein %in% sig_prots)

# NOTE: seems like all attributes should be strings to make it into Cytoscape
# convert each col to character
noa <- as.data.frame(apply(noa, 2, function(x) as.character(x)))

# Loop to add node attributes to igraph netw_graph
# NOTE: this can take a couple seconds
for (i in c(1:ncol(noa))) {
	namen <- colnames(noa)[i]
	col_data <- setNames(noa[[i]],nm=noa$Protein)
	netw_g <- set_vertex_attr(netw_g,namen,value=col_data[names(V(netw_g))])
}


fwrite(noa,"noa.csv")

quit()

## ---- Create Cytoscape graphs

# all modules
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# Check that we are connected to Cytoscape
cytoscapePing()

## Loop to create graphs:
message("\nCreating Cytoscape graphs!")
pbar <- txtProgressBar(max=length(modules),style=3)

for (module in names(modules)){ 

	nodes <- modules[[module]]

	createCytoscapeGraph(module, netw_g, ppi_g, nodes)

	setTxtProgressBar(pbar, value = match(module,names(modules)))

}
close(pbar)

## When done, save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.

myfile <- file.path(netwdir,"Modules.cys")
winfile <- gsub("/mnt/d/|\\~/","D:/",myfile) 

RCy3::saveSession(winfile)

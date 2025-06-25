wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(dplyr)
library(igraph)

set.seed(121212)

Ranges <- read.table("./SE_Global.txt", header = TRUE)

# Function to create the adjacency matrix
synon <- function(Ranges, p, q) {
  M <- matrix(0, nrow = nrow(Ranges), ncol = nrow(Ranges))
  
  for (i in 1:(nrow(Ranges)-1)) {
    for (j in (i+1):nrow(Ranges)) {
      genus_i <- strsplit(Ranges$species[i], "_")[[1]][1]
      genus_j <- strsplit(Ranges$species[j], "_")[[1]][1]
      
      if (genus_i != genus_j) {
        M[i, j] <- 0
        M[j, i] <- 0
        next
      }
      
      overlap_time <- min(Ranges$ts[i], Ranges$ts[j]) - max(Ranges$te[i], Ranges$te[j])
      
      check_overlap <- function(dist_i, dist_j) {
        dist_i_set <- unlist(strsplit(dist_i, ";"))
        dist_j_set <- unlist(strsplit(dist_j, ";"))
        return(length(intersect(dist_i_set, dist_j_set)) > 0)
      }
      
      if (overlap_time > 0) {
        if (check_overlap(Ranges$distribution[i], Ranges$distribution[j])) {
          M[i, j] <- rbinom(1, 1, p)
          M[j, i] <- M[i, j]
        } else {
          M[i, j] <- rbinom(1, 1, q)
          M[j, i] <- M[i, j]
        }
      }
    }
  }
  rownames(M) <- Ranges$species
  colnames(M) <- Ranges$species
  
  return(M)
}

# Clustering behavior
library(ggplot2)

p_vals <- seq(0, 0.05, by = 0.005)
q_vals <- seq(0, 0.05, by = 0.005)

size <- matrix(0, nrow = length(p_vals), ncol = length(q_vals))
clust <- matrix(0, nrow = length(p_vals), ncol = length(q_vals))

for (i in seq_along(p_vals)) {
  for (j in seq_along(q_vals)) {
    
    p <- p_vals[i]
    q <- q_vals[j]
    
    if (p > q) {
      temp1 <- numeric()
      temp2 <- numeric()
      for (k in 1:10) {
        M <- synon(Ranges, p, q)
        graph <- graph_from_adjacency_matrix(M, mode = "undirected", weighted = NULL, diag = FALSE)
        temp1[[k]] <- length(components(graph)$csize)
        temp2[[k]] <- max(components(graph)$csize)
        cat(sprintf("\rCalculating... i=%.3f, j=%.3f", p, q))
        flush.console()
      }
      size[i, j] <- mean(temp1)
      clust[i, j] <- mean(temp2)
    } else {
      size[i, j] <- NA
      clust[i, j] <- NA
    }
  }
}

save(size, file = "size_results.RData")
save(clust, file = "clust_results.RData")

load(file="./size_results.RData")
df1 <- expand.grid(p = p_vals, q = q_vals)
df1$max_group_size <- as.vector(size)

ggplot(df1, aes(x = p, y = q, fill = max_group_size)) +
  geom_tile() +
  geom_text(aes(label = round(max_group_size, 1)), color = "white", size = 12) +
  scale_fill_gradient(low = "pink", high = "black", guide = "none") +
  scale_x_continuous(breaks = seq(0, 0.05, by = 0.01)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01)) +
  theme_minimal(base_size = 36) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 36),
    axis.title = element_text(size = 64),
    plot.title = element_text(size = 72, face = "bold"),
  ) +
  labs(x = "p", y = "q", fill = "Size")

load(file="./clust_results.RData")
df2 <- expand.grid(p = p_vals, q = q_vals)
df2$max_group_size <- as.vector(clust)

ggplot(df2, aes(x = p, y = q, fill = max_group_size)) +
  geom_tile() +
  geom_text(aes(label = round(max_group_size, 1)), color = "white", size = 12) +
  scale_fill_gradient(low = "skyblue", high = "black", guide = "none") +
  scale_x_continuous(breaks = seq(0, 0.05, by = 0.01)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01)) +
  theme_minimal(base_size = 36) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 36),
    axis.title = element_text(size = 64),
    plot.title = element_text(size = 72, face = "bold"),
  ) +
  labs(x = "p", y = "q")



# Generate random synonymization scenarios for p=0.01, 0.02, 0.03 and q=p/10, p/5, p/2

p_values <- c(0.01, 0.02, 0.03)

PyRate_sc <- function(dataset) {
  for (p in p_values) {
    for (i in c(10, 5, 2)) {
      for (k in 1:5) {
        q <- p/i
        M <- synon(Ranges, p, q)
        g <- graph_from_adjacency_matrix(M, mode = "undirected", weighted = NULL, diag = FALSE)
        
        components <- components(g)
        group_labels <- paste("sp", 1:length(components$csize), sep = "")
        
        newdataset <- dataset
        
        for (group_index in 1:length(group_labels)) {
          group_species <- c()
          group_species <- names(components$membership)[components$membership == group_index]
          newdataset$Taxon_name[newdataset$Taxon_name %in% group_species] <- group_labels[group_index]
        }
        
        dataset_name <- deparse(substitute(dataset))
        output_filename <- sprintf("./%s_p%.2fq%.3f_%d.txt", dataset_name, p, q, k)
        write.table(newdataset, file = output_filename, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  }
}



## PyRate (spatially standardized)
### 1. Global
Global_std <- read.table("./Global_standardized.txt", header = TRUE)
PyRate_sc(Global_std)
### 2. NAfEu
NAfEu_std <- read.table("./NAfEu_standardized.txt", header = TRUE)
PyRate_sc(NAfEu_std)
### 3. WAfSAm
WAfSAm_std <- read.table("./WAfSAm_standardized.txt", header = TRUE)
PyRate_sc(WAfSAm_std)

# Plot
plotgraph <- function(syngraph) {
  
  syn <- graph_from_adjacency_matrix(syngraph, mode = "undirected", diag = FALSE)
  
  regions <- c("Eu", "NAf", "NAf_Som", "NAm", "WAf", "SAm", "SAm_Arg")
  region_colors <- c("blue", "skyblue", "purple", "yellow", "orange", "red", "pink")
  
  pie <- function(dist) {
    region <- unlist(strsplit(dist, ";"))
    pie_colors <- region_colors[match(names(table(region)), regions)]
    return(list(sizes = as.vector(table(region)), colors = pie_colors))
  }
  
  vertex_pie <- lapply(Ranges$distribution, pie)
  
  plot(syn, vertex.label = NA, vertex.size = 5, vertex.shape = "pie",
    vertex.pie = lapply(vertex_pie, function(x) x$sizes),
    vertex.pie.color = lapply(vertex_pie, function(x) x$colors),
    vertex.pie.lty = 1, edge.width = 1, edge.color = "black")
  
  legend("bottomleft", legend = regions, fill = region_colors,  cex = 1, text.font = 2, bty = "n")
}

# Example plots
ex <- synon(Ranges, 0.02, 0.01)
M <- ex
plotgraph(M)


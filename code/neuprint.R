# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("natverse/neuprintr")

rm(list = ls())
setwd("~/Documents/projects/DB-connectomics/code/")

library(neuprintr)
library(reshape2)
library(ggplot2)
library(dplyr)
theme_set(theme_bw())

# usethis::edit_r_environ()
neuprint_dataset = "hemibrain:v1.0"
datasets_neuprint <- neuprint_datasets()

re_load_dataset <- T
if(!re_load_dataset){

  datasets_neuprint$`hemibrain:v1.2.1`$ROIs
  datasets_neuprint$`hemibrain:v1.2.1`$superLevelROIs
  
  neurons_per_region <- lapply(unique(datasets_neuprint$`hemibrain:v1.2.1`$ROIs),
                               function(i) try(neuprintr::neuprint_find_neurons(input_ROIs = i)))
  saveRDS(neurons_per_region, "../data/neurons_per_region.RDS")
}else{
  cat('Reading existent RDS\n')
  neurons_per_region <- readRDS("../data/neurons_per_region.RDS")
}

ROIs <- unique(datasets_neuprint$`hemibrain:v1.2.1`$ROIs)[!sapply(neurons_per_region, function(i) is.null(dim(i)))]
# datasets_neuprint <- datasets_neuprint[!sapply(neurons_per_region, function(i) is.null(dim(i)))]
neurons_per_region <- neurons_per_region[!sapply(neurons_per_region, function(i) is.null(dim(i)))]

neurons_per_region_summary <- sapply(neurons_per_region, function(i) try(table(i$bodytype)))
names(neurons_per_region_summary) <- ROIs
# neuprint_dump(dir = "~/Desktop/", roi = datasets_neuprint$`hemibrain:v1.2.1`$ROIs[1])

table(sapply(neurons_per_region_summary, typeof))

neurons_per_region_summary <- neurons_per_region_summary[sapply(neurons_per_region_summary, typeof) != 'character']

first_col_to_rownames <- function(i){
  rownames(i) <- i[,1]
  i[,-1]
}

NA_are_zero <- function(i){
  i[is.na(i)] <- 0
  i
}

normalise_rw <- function(i){
  sweep(i, 1, rowSums(i), '/')
}

neurons_per_region_matrix <- NA_are_zero(first_col_to_rownames(dcast(melt(neurons_per_region_summary), Var1~L1, value.var='value')))
neurons_per_region_summary_ph <- pheatmap::pheatmap(normalise_rw(1+neurons_per_region_matrix))
## binary: if there are any of each type of neuron in the area
neurons_per_region_summary_ph_bin <- pheatmap::pheatmap(apply((neurons_per_region_matrix)>0, 1, as.numeric))
neurons_per_region_summary_ph
# neurons_per_region_summary_matrix_collapse <- 

neurons_per_region_matrix

# neurons_per_region_summary_ph_bin_collapse=

ggplot(melt(neurons_per_region_summary), aes(x=factor(Var1, levels=neurons_per_region_summary_ph$tree_row$labels[neurons_per_region_summary_ph$tree_row$order]),
                                             y=factor(L1, neurons_per_region_summary_ph$tree_col$labels[neurons_per_region_summary_ph$tree_col$order]), fill=value))+
  geom_tile()
ggsave("~/Documents/projects/connectomics/results/number-cell-types-per-region.pdf", height = 10, width = 10)

neurons_per_region_norm <- melt(as(normalise_rw(neurons_per_region_matrix), 'matrix')) %>% filter(value > 0)

ggplot(neurons_per_region_norm,
       aes(x=factor(Var1, levels=neurons_per_region_summary_ph$tree_row$labels[neurons_per_region_summary_ph$tree_row$order]),
           y=factor(Var2, neurons_per_region_summary_ph$tree_col$labels[neurons_per_region_summary_ph$tree_col$order]),
           fill=value))+
  geom_tile()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_fill_stepsn(colours = c("blue","black"), breaks=c(0.01, 0.05, 0.1, 0.7))
# scale_fill_stepsn(colours = c("#fbf8cc","#fde4cf","#ffcfd2","#f1c0e8"), breaks=c(0.05, 0.1, 0.2, 0.3, 0.7))
ggsave("~/Documents/projects/connectomics/results/number-cell-types-per-region_norm.pdf", height = 10, width = 10)


ggplot(neurons_per_region_norm %>% filter(Var1 %in% names(sort(table(neurons_per_region_norm$Var1), decreasing = T)[1:50])) %>%
         filter(Var2 %in% names(sort(table(neurons_per_region_norm$Var2), decreasing = T)[1:50])),
       aes(x=factor(Var1, levels=neurons_per_region_summary_ph$tree_row$labels[neurons_per_region_summary_ph$tree_row$order]),
           y=factor(Var2, neurons_per_region_summary_ph$tree_col$labels[neurons_per_region_summary_ph$tree_col$order]),
           fill=value))+
  geom_tile()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank(),
  #       axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank())+
  scale_fill_stepsn(colours = c("blue","black"), breaks=c(0.01, 0.05, 0.1, 0.7))


ggplot(neurons_per_region_norm %>% filter(grepl('MB', Var2)),
       aes(x=factor(Var1, levels=neurons_per_region_summary_ph$tree_row$labels[neurons_per_region_summary_ph$tree_row$order]),
           y=factor(Var2, neurons_per_region_summary_ph$tree_col$labels[neurons_per_region_summary_ph$tree_col$order]),
           fill=value))+
  geom_tile()+
  scale_fill_stepsn(colours = c("blue","black"), breaks=c(0.01, 0.05, 0.1, 0.7))


neurons_per_region

neurons_per_region_summary_melt <- melt(neurons_per_region_summary)
to_named_vector <- function(i){
  .x <- unlist(i[,2])
  names(.x) <- unlist(i[,1])
  .x
}

neurons_per_region_summary_melt_abundance <- to_named_vector(neurons_per_region_summary_melt %>% group_by(Var1) %>% dplyr::summarise(count_neuron_type=sum(value)))
neurons_per_region_summary_melt_abundance_norm <- to_named_vector(neurons_per_region_summary_melt %>% group_by(Var1) %>% dplyr::summarise(count_neuron_type=sum(value)))/sum(neurons_per_region_summary_melt$value)
neurons_per_region_summary_melt_abundance_norm <- sort(neurons_per_region_summary_melt_abundance_norm, decreasing = T)[1:100]
ggplot(cbind.data.frame(count=neurons_per_region_summary_melt_abundance_norm, name=names(neurons_per_region_summary_melt_abundance_norm)),
       aes(x=factor(name, levels=names(sort(neurons_per_region_summary_melt_abundance_norm))), y=count))+geom_line(aes(group=1))+
  scale_y_continuous(trans = 'log2')+  theme(axis.text.x = element_text(angle = 30, hjust=1))
          

##' check compatibility between clusters and annotated cell types from janelia. we want something scale-invariant - possibly clr transformation? although that is
##' is not subcomposition-invariant

## get a clustered dataset, or sets of clusters
clusters <- readRDS("~/Desktop/clusters_npc10_res0p2.RDS")


plot(cbind(clusters=as.vector(table(clusters)/length(clusters)),
     neuprint_cell_types=as.vector(neurons_per_region_summary_melt_abundance_norm[1:length(unique(clusters))])), type='l')
abline(coef=c(0,1))


vec1 <- as.vector(table(clusters)/length(clusters))
vec2 <- as.vector(neurons_per_region_summary_melt_abundance_norm[1:length(unique(clusters))])

n=1
m=length(vec1)
h=length(vec2)
max_h_m <- max(c(m,h))

adj_matrix <- matrix(0, nrow=m, ncol=h)
diag(adj_matrix) <- 1 
adj_matrix

# nits <- 1000
# dists <- rep(NA, nits)
# for(it in 1:nits){
#   adj_matrix <- matrix(0, nrow=m, ncol=h)
#   suggested_matrix_passes_test <- F
#   count_does_not_pass <- 0
#   while(!suggested_matrix_passes_test){
#     cat('Does not pass ', count_does_not_pass, '\n')
#     suggested_matrix <- adj_matrix
#     # suggested_matrix[sample(1:length(suggested_matrix), size = sample(1:length(suggested_matrix), size = 1), replace = F)] <- 1
#     suggested_matrix[sample(1:length(suggested_matrix), size = sample(1:30, size = 1), replace = F)] <- 1
#     suggested_matrix_passes_test <- (sum(colSums(suggested_matrix)) <= max_h_m)
#     ##additional constraints
#     if(sum(adj_matrix)==0) suggested_matrix_passes_test <- F
#     count_does_not_pass <- count_does_not_pass+1
#   }  
#   dists[it] <- as.vector(dist(rbind((vec1 %*% suggested_matrix), vec2)))
# }

plot(dists)
suggested_matrix

suggested_matrix
as.vector(dist(rbind((vec1 %*% adj_matrix), vec2)))
rbind((vec1 %*% adj_matrix))
vec2
plot(cbind((vec1 %*% adj_matrix), vec2))

aa <- matrix(0, ncol=3, nrow=3)
diag(aa) <- 1
aa[1,2] <- 1
aa

##---------------------------------------------------------------------------##
## optimal transport (it hasn't worked)
library(transport)
res_OT <- transport::wasserstein(c(1-sum(vec1), vec1), c(1-sum(vec2), vec2),
                       costm = matrix(1, nrow=m+1, ncol=h+1))
res_OT


# total_weight_of_matching <- function(arg_A, arg_T, adjacency_matrix){
#   abs((arg_A %*% adjacency_matrix) - arg_T)
# }

linear_programming_assignment_problem <- function(arg_A, arg_T){
  ## https://en.wikipedia.org/wiki/Assignment_problem
  ## edge (i,j) where i is an index for 1,...|arg_A| and j is an index for 1,...|arg_T|
  x <- matrix(NA, length(arg_A), length(arg_T))
  ## where xij is {0,1} and is the adjacency_matrix
  
  ## LP:
  ## maximise (arg_A %*% adjacency_matrix) - arg_T
  # sum_j xij = 1 for all i (not necessarily!)
  # sum_i xij = 1 for all j (not necessarily!)
  adjacency_matrix
  ## 

}

library(lpSolve)


f.obj <- c(5, 7)

# Set matrix corresponding to coefficients of constraints by rows
# Do not consider the non-negative constraint; it is automatically assumed
f.con <- matrix(c(1, 0,
                  2, 3,
                  1, 1), nrow = 3, byrow = TRUE)

# Set unequality signs
f.dir <- c("<=",
           "<=",
           "<=")

# Set right hand side coefficients
f.rhs <- c(16,
           19,
           8)

# Final value (z)
lp("max", f.obj, f.con, f.dir, f.rhs)

# Variables final values
lp("max", f.obj, f.con, f.dir, f.rhs)$solution


#------
# Set assignment costs matrix
# costs <- matrix(c(15, 10, 9,
#                   9, 15, 10,
#                   10, 12 ,8), nrow = 3, byrow = TRUE)
# 
# # Print assignment costs matrix
# costs
# 
# # Final value (z)
# lp.assign(costs)
# 
# # Variables final values
# lp.assign(costs)$solution


arg_A <- c(4,2,1,1)
arg_B <- c(2,1,2,3)
num_factors <- 1:max(c(length(arg_A), length(arg_B)))
latent_factors <- lapply(num_factors, function(num_fact){
  combinat::combn(x = 1:num_fact, m = num_fact)
})

expand.grid(1:4, 1:4)

combinat::combn(1:10, 4)

utils::combn(10, 5)

expand.grid(do.call('cbind', lapply(1:3, function(i) 1:num_fact)))

## start with all in separate grousps and then amalgamate?

arg_A
arg_B

##
initial_factors_A <- c(1,2,3,3)
initial_factors_B <- c(1,2,3,4)
# initial_factors_A <- c(1,2,3,2) ## true solution
# initial_factors_B <- c(1,3,1,2) ## true solution
initial_factors_A <- factor(initial_factors_A, levels=unique(c(unique(initial_factors_A),
                                                               unique(initial_factors_B))))
initial_factors_B <- factor(initial_factors_B, levels=unique(c(unique(initial_factors_A),
                                                               unique(initial_factors_B))))

compute_loss <- function(suggested_factors_list){
  FA <- factor(suggested_factors_list[[1]], levels=unique(unlist(suggested_factors_list)))
  FB <- factor(suggested_factors_list[[2]], levels=unique(unlist(suggested_factors_list)))
  
  sum(( sapply(split(suggested_factors_list[[1]], f = FA), function(i) sum(as.numeric(i))) -   
          sapply(split(suggested_factors_list[[2]], f = FB), function(i) sum(as.numeric(i))))**2)
  
}

compute_loss_vec <- function(suggested_factors_list){
  FA <- factor(suggested_factors_list[[1]], levels=unique(unlist(suggested_factors_list)))
  FB <- factor(suggested_factors_list[[2]], levels=unique(unlist(suggested_factors_list)))
  
  rbind( unlist(sapply(split(suggested_factors_list[[1]], f = FA), function(i) sum(as.numeric(i)))),
         unlist(sapply(split(suggested_factors_list[[2]], f = FB), function(i) sum(as.numeric(i)))))
}

suggest_new_factors <- function(list_of_vec_factors){
  sampled_list <- sample(x = 1:2, 1)
  idx_to_change <- sample(1:length(list_of_vec_factors[[sampled_list]]), 1)
  ## i below we might be sampling the same again!
  list_of_vec_factors[[sampled_list]][idx_to_change] <- sample(levels(list_of_vec_factors[[sampled_list]]), 1)
  return(list_of_vec_factors)
}

nits <- 100
loss <- rep(NA, nits)
loss[1] <- sum(( sapply(split(arg_A, f = initial_factors_A), sum) -   sapply(split(arg_B, f = initial_factors_B), sum))**2)
prev_factors <- replicate(nits, list())
pass <- F
additional_constraints <- F
for(it in 2:nits){
  prev_factors[[it]] <- list(A=initial_factors_A, B=initial_factors_B)
  ct <- 0
  while(!pass){
    cat('Count in trial: ', ct, '\n')
    suggested_factors <- prev_factors[[it]]
    cat(sapply(suggested_factors, paste0, collapse='-'), sep=' | ', '\n')
    suggested_factors <- suggest_new_factors(suggested_factors)
    suggested_loss <- compute_loss(suggested_factors)
    cat('Suggested loss=', suggested_loss, '\n')
    if(suggested_loss <= loss[it-1]){
      cat('Accepted\n')
      ## accept
      pass <- T
    }else{
      cat('Failed\n')
      pass <- F
    }
    ## make sure a) that if we are sampling the first vector, the elements never seen before are ordered [TO DO]
    if(additional_constraints){
      if(!all(as.numeric(unique(prev_factors[[it]][[1]])) == levels(prev_factors[[it]][[1]]))) pass <- F
    }
    # lapply(2:length(prev_factors[[it]][[1]]), function(i) lapply(sort(unique(prev_factors[[it]][[1]])), function(i){.x <- which(prev_factors[[it]][[1]] == i); sapply(2:length(.x), function(j) .x[j]>.x[j-1])}))
    
    ## b) that we haven't lost any elements in both lists, and if we have remove them [TO DO]
    ct <- ct+1
  }
  loss[it] <- suggested_loss
  pass <- F
}


loss


compute_loss(suggested_factors)
compute_loss_vec(suggested_factors)

cowplot::plot_grid(
ggplot(melt(do.call('rbind', sapply(prev_factors, unlist))),
       aes(x=Var1, y=Var2, fill=factor(value)))+labs(x='Iteration', y='Labels')+
  geom_tile()+jcolors::scale_fill_jcolors(palette = 'pal4')+theme(legend.position = 'top'),
ggplot(melt(matrix(loss)), aes(x=Var1, y=value))+geom_line()+labs(y='Loss', y='Iteration'),
rel_heights = c(3,1), ncol=1)

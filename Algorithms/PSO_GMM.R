###############################################################################-
############################     PSO--GMM   ####################################
###############################################################################-

# This code : Use GMM for clustering and PSO for learning the best DEG list
# library packages
library(Seurat)
library(cluster)
library(mclust)
library(dplyr)
### read in data
epicell <- readRDS("~/wyc_20231021/single_cell_SEQ/project/chemotherapy/epi_0403_final.rds")

DEG <- read.csv("~/wyc_20231021/Python/GMM_mixture/DEG_new.csv",row.names = 1)

### set the parameters
seurat_object <- epicell
size <- 30
poolrange <- 500
generation <- 30
toprange <- 50
dim <- toprange*2
vmax <- 2
c1 <- 2
c2 <- 2
r1 <- runif(1)
r2 <- runif(1)
DEG_list <- DEG$genename
DEG_range <- length(DEG$genename)
certainty <- 0.01


###############################################################################-
###############    creat the initial population and velocities   ###############

# create individual 
create_individual <- function(poolrange, DEG_range) {
  numbers1 <- 1:poolrange
  weights1 <- exp(-1e-5 * (numbers1 - 1)^2)
  weights1 <- weights1 / sum(weights1) # scale the weight
  sample1 <- sample(numbers1, size = 50, replace = FALSE, prob = weights1)
  sample1 <- sort(sample1)
  
  numbers2 <- DEG_range:(DEG_range - poolrange + 1)
  weights2 <- exp(-1e-5 * (DEG_range - numbers2)^2)
  weights2 <- weights2 / sum(weights2) # scale the weight
  sample2 <- sample(numbers2, size = 50, replace = FALSE, prob = weights1)
  sample2 <- sort(sample2)
  
  individual <- c(sample1, sample2)
  return(individual)
}

# creat the initial population 
initial_population <- function(poolrange, DEG_range, size) {
  population_matrix <- matrix(nrow = size, ncol = 100)
  for (i in 1:size) {
    individual <- create_individual(poolrange, DEG_range)
    population_matrix[i, ] <- individual
  }
  return(population_matrix)
}

initial_population_matrix <- initial_population(poolrange, DEG_range, size)

###############################################################################-
######################    creat the initial velocities #########################

# initial velocities
velocities <- matrix(runif(dim * size, -vmax, vmax), ncol = dim)


###############################################################################-
######################    get fitness function  ################################

#evaluate class by WGSS & BGSS
eva_ss <- function(data) {
  # every cluster's center
  centroids <-
    data.frame(sapply(unique(data$class), function(class) {
      colMeans(data[data$class == class, 1:2])
    }))
  colnames(centroids) <- c("centroid_x", "centroid_y")
  # WGSS
  WGSS <- sum(sapply(unique(data$class), function(class) {
    sum((data[data$class == class, 1] - centroids[class, "centroid_x"]) ^ 2 +
          (data[data$class == class, 2] - centroids[class, "centroid_y"]) ^2)
  }))
  # overall center
  overall_centroid <- colMeans(data[, 1:2])
  # BGSS
  BGSS <- sum(sapply(unique(data$class), function(class) {
    n_cluster <- sum(data$class == class)
    sum((centroids[class, "centroid_x"] - overall_centroid[1]) ^ 2 +
          (centroids[class, "centroid_y"] - overall_centroid[2]) ^ 2) * n_cluster
  }))
  fitness <-  BGSS-WGSS 
  return(fitness)
}

get_fitness <- function(x,DEG_list,seurat_object,toprange,dim,certainty){
  list1 <- list(DEG_list[x[1:toprange]])
  list2 <- list(DEG_list[x[toprange:dim]])
  seurat_object <- AddModuleScore(seurat_object,features = list1,name = "score1_")
  seurat_object <- AddModuleScore(seurat_object,features = list2,name = "score2_")
  score_matrix<- seurat_object@meta.data[, c("score1_1", "score2_1")]
  GMM <- Mclust(score_matrix,G = 2)
  score_matrix$uncertainty <- GMM$uncertainty
  score_matrix$class <- GMM$classification
  certain_SM <- subset(score_matrix,score_matrix$uncertainty<certainty)
  fitness <- eva_ss(data = certain_SM)
  return(fitness)
}

###############################################################################-
######################      MAIN FUNCTION       ################################
###############################################################################-


################ pBest and gBest (personal and group)  #########################

pBest <- initial_population_matrix

initial_fitness <- apply(initial_population_matrix, MARGIN = 1, get_fitness,
                         DEG_list = DEG_list, seurat_object = seurat_object, 
                         toprange = toprange, dim = dim, certainty = certainty)

gBest <- initial_population_matrix[which.max(initial_fitness), ]

pBest_fitness <- initial_fitness
gBest_fitness <- max(initial_fitness)

population_matrix <- initial_population_matrix


##############################   show result   #################################
result_matrix <- matrix(ncol = 2,nrow = generation)

# PSO MAIN FUNCTION
for (gen in 1:generation) {
  for (i in 1:size) {
    # update individual position
    population_matrix[i,] <- round(population_matrix[i,] + velocities[i,])
    
    # get fitness for every individual
    fitness <- get_fitness(population_matrix[i,], DEG_list, seurat_object,toprange,dim,certainty)
    
    # update individual best fitness
    if (fitness > pBest_fitness[i]) {
      pBest_fitness[i] <- fitness
      pBest[i,] <- population_matrix[i,]
    }
    
    # update globle best fitness
    if (fitness > gBest_fitness) {
      gBest_fitness <- fitness
      gBest <- population_matrix[i,]
    }
  }
  
  gBest_matrix <- rep(gBest,times=size)
  gBest_matrix <- matrix(gBest_matrix,nrow = size,byrow = T)
  
  velocities <- velocities + c1 * r1 * (pBest - population_matrix) + c2 * r2 * (gBest_matrix - population_matrix)
  velocities[velocities > vmax] <- vmax
  velocities[velocities < (-vmax)] <- (-vmax)

  cat("Generation:", gen, " gBestFitness:", gBest_fitness, "\n")
  result_matrix[gen,] <- c(gen,gBest_fitness)
  return(result_matrix)
}

cat("Global Best Solution:\n")
print(gBest)
cat("Global Best Fitness:\n")
print(gBest_fitness)

















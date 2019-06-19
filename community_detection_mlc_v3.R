#### Run as : Rscript community_detection_mlc.R <"LAML.filtered.positive.sorted"> <output_filename>

## Import packages

library(igraph)
library(dplyr)

## Set script arguments

args = commandArgs(trailingOnly = TRUE)
file_name = args[1]
out_file = args[2]
out_file2 = args[3]
## Data import and cleaning

# Import gene list with scores
parsed_genes = read.csv(file_name, sep='\t',header=TRUE,stringsAsFactors = FALSE)

## Creating igraph object

edgelist = as.matrix(parsed_genes[,1:2])
g = graph.edgelist(edgelist, directed=FALSE)

## Community detection : Multilevel algorithm
mlc <- multilevel.community(g)

## Writing the results to a file

mlc_community_list <- as.data.frame(as.matrix(membership(mlc)))
mlc_community_list$gene <- rownames(mlc_community_list)
mlc_community_member <- arrange(mlc_community_list, V1) %>%
  select(gene, V1)
colnames(mlc_community_member)[2] <- 'community'

community_info = mlc_community_member

# Lookup Table
vals <- community_info[,2]
keys <- community_info[,1]
lookup <- setNames(vals, keys)

index = which(lookup[parsed_genes[,1]] != lookup[parsed_genes[,2]])
filtered_edgelist <- parsed_genes[-index,]

write.table(mlc_community_member, row.names=FALSE, col.names=TRUE, file=out_file, sep='\t', quote=FALSE)
write.table(filtered_edgelist, row.names=FALSE, col.names=TRUE, file=paste(out_file2,sep=''), sep='\t', quote=FALSE)

print(paste('Number of Edges Before filtered: ', as.character(dim(parsed_genes)[1]), sep=""))
print(paste('Number of Eges After filtered: ', as.character(dim(filtered_edgelist)[1]), sep=""))

# Check to make sure the number of nodes did not change
check1 = c(parsed_genes[,1],parsed_genes[,2])
uniq1 = unique(check1)
a = length(uniq1)
print(paste('Number of Nodes Before filtered: ', as.character(a)))

check2 = c(filtered_edgelist[,1], filtered_edgelist[,2])
uniq2 = unique(check2)
b = length(uniq2)
print(paste('Number of Nodes After filtered: ', as.character(b)))

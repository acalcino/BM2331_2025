# Phylogenetics

Phylogenetics is the study of the relatedness of genes and species. The first branched phylogenetic tree was created by Darwin as a rough sketch in a notebook in 1837 which has since become known as the "I think" sketch.

![think](images/think.jpg)

Beside this he scribbled:

>Case must be that one generation should have as many living as now.
>To do this & to have many species in same genus (as is) requires extinction.
>Thus between A and B immense gap of relation.
>C and B the finest gradation, B and D rather greater distinction.
>Thus genera would be formed. — bearing relation to ancient types with several extinct forms.

This was the first articulation of the idea of `descent with modification` which underpins our understanding of evolution by natural selection to this day. What the figure and the caption illustrate is the idea that all species alive today share common ancestors and that over the immense time since the divergence of all extant species, an incalculable number of other species which *fill the gaps* between modern species have gone extinct. This process of extinction of less well suited organisms and the survival of more well suited ones is the engine that drives the evolutionary process. This was a revolutionary insight by Darwin is made even more remarkable by the fact that he had no idea about genetics or the process by which novel character emergence could come to pass.

While early trees were based entirely on phenotypic character trait analysis, today we also have the option of building relationship trees based on genetic traits ie. genetic sequences. When we do this for whole genomes to ascertain the relatedness of species, this is called `phylogenomics`.

![phylogenomics](images/phylogenomics.png)

But when we do this for single genes/proteins, this is called `phylogenetics`, although you do see people using this term for the study of species too.

![phylogeneitcs](images/phylogenetics.png)

Today we are going to build phylogenetic trees of our gene family of interest, based on the multiple sequence alignment we created in the previous tutorial. As with all things bioinformatics, there is more than one way to skin a cat. In the phylogenetics lectures you learnt about multiple methods to infer relatedness between aligned sequences including Neighbour Joining (NJ), Maximum Likelihood (ML), Bayesian Inference etc. Today we'll try to build two trees, one using a simple NJ algorithm and another using a more sophisticated ML algorthim to compare their performance in terms of accuracy and computational demand.


## Neighbour Joining tree

In RStudio, set your working directory to your `phylogenetic_project` directory, create a new folder inside called something like `trees` and then inside this, create `NJ`.

```R

setwd("~/working-directory/phylogenetic_project")

dir.create("trees/NJ", recursive = TRUE)

```

Now we are going to load all the appropriate libraries including the `phangorn` package for tree construction and we'll also load our MSA.

```R

library(phangorn)
library(Biostrings)

# Load the alignment file
alignment <- readAAStringSet("msa/msa_muscle_aligned_sequences.fa")

# Convert Convert to phangorn format
phyDat_alignment <- as.phyDat(alignment, type = "AA")

# Create a distance matrix and initial NJ tree for model testing. This calculates the maximum likelihood distances between all pairs of sequences in your alignment, which is the expected number of amino acid substitutions per site.
dm <- dist.ml(phyDat_alignment)
nj_tree <- NJ(dm)

# Visualise this initial tree
plot(nj_tree, main = "Neighbor-Joining Tree (midpoint rooted)", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

```

This creates an `unrooted` tree which you will be able to see in your RStudio Plot pane. Have a look at the tree you've created to see if it makes evolutionary sense ie. do related sequences share a common ancestor? Next up, we'll root the tree by defining an outgroup which we spoke about in our previous tutorial. 

```R

# Root the tree on your outgroup sequence
outgroup <- "Cin_Pax6_NP_001027641.1 homeobox transcription factor Pax6 [Ciona intestinalis]"
nj_tree_outgroup <- root(nj_tree, outgroup = outgroup, resolve.root = TRUE)

# Plot your rooted tree
plot(nj_tree_outgroup, main = "Neighbor-Joining Tree (outgroup rooted)", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

# Save your rooted NJ tree in Newick format
write.tree(nj_tree, file = "trees/NJ/NJ_tree.nwk"

```

We can also calculate the confidence level of each node in the tree by `bootstrapping` it. Bootstrapping is a statistical resampling method used to estimate the reliability of a result by repeatedly sampling from your original data with replacement. In this case, *replacement* means to randomly sample columns from your MSA, rather than using the original sequence order to see if these purturbations influence the final topology of the tree. The theory is that if your tree is **robust**, this process of resampling with replacement should produce a tree with the same topology more often than not. Topology referes to how each sample relates to one another. Here's an example:

Original aligment of 8 amino acid sequeces
>        Position: 1  2  3  4  5  6  7  8
>Sequence A:        M  K  R  D  L  M  R  E
>Sequence B:        M  K  H  D  V  M  R  Q
>Sequence C:        L  S  R  D  L  I  G  E
>Sequence D:        L  S  H  N  V  I  G  Q

This creates a tree with the topology:

>Original tree topology: ((A,B),(C,D))

Now we do bootstrap replicate 1:
Randomly selected positions: 1, 1, 3, 5, 2, 7, 8, 4
>Position:         1  1  3  5  2  7  8  4
>Sequence A:        M  M  R  L  K  R  E  D
>Sequence B:        M  M  H  V  K  R  Q  D
>Sequence C:        L  L  R  L  S  G  E  D
>Sequence D:        L  L  H  V  S  G  Q  N

Now we recalculate topology:

>Bootstrap tree 1: ((A,B),(C,D)) ✓ Same topology

Now we do bootstrap replicate 2:

Randomly selected positions: 2, 4, 4, 6, 8, 3, 1, 1
>Position:         2  4  4  6  8  3  1  1
>Sequence A:        K  D  D  M  E  R  M  M
>Sequence B:        K  D  D  M  Q  H  M  M
>Sequence C:        S  D  D  I  E  R  L  L
>Sequence D:        S  N  N  I  Q  H  L  L

Calculate topology:

>Bootstrap tree 2: (A,(B,(C,D))) ✗ Different topology!

After doing this 100 or 10000 times, we see how many trees we build from these resampled MSAs have the same or different topology.

>Clade (A,B):    appeared in 85 trees → Bootstrap support = 85%
>Clade (C,D):    appeared in 92 trees → Bootstrap support = 92%
>Clade ((A,B),(C,D)): appeared in 78 trees → Bootstrap support = 78%

And now we can add these bootstrap values to our original tree:

>                         78
>            ┌────────────┴────────────┐
>            │                         │
>         85 │                      92 │
>    ┌───────┴───────┐        ┌────────┴────────┐
>    │               │        │                 │
>Sequence A       Sequence B   Sequence C     Sequence D

Let's do this for our tree.

```R

# Convert phyDat to AAbin format
alignment_AAbin <- as.AAbin(phyDat_alignment)

# Bootstrap using manual resampling
set.seed(123)  # For reproducibility
boot_trees <- list()

# Run bootstrapping with 1000 replicates
for(i in 1:1000) {

  # Create a vector of random column numbers
  boot_indices <- sample(ncol(alignment_AAbin), replace = TRUE)
  
  # Use those column numbers to create new modified alignment
  boot_alignment <- alignment_AAbin[, boot_indices]
  
  # Convert from AAbin format to phyDat (phangorn) format
  boot_phyDat <- as.phyDat(boot_alignment)
  
  # Calculate distance matrix of evolutionary distances between every pair of sequences
  boot_dm <- dist.ml(boot_phyDat)
  
  # Build tree
  boot_tree <- NJ(boot_dm)
  
  # Root with the same outgroup
  boot_tree <- root(boot_tree, outgroup = outgroup, resolve.root = TRUE)
  boot_trees[[i]] <- boot_tree
  
  if(i %% 10 == 0) cat(i, " ")
}
cat("Done!\n")

# Calculate bootstrap support for all nodes
nj_bs <- prop.clades(nj_tree_outgroup, boot_trees)
nj_bs_percent <- round(nj_bs / length(boot_trees) * 100, 0)

# Plot NJ tree with bootstrap values
plot(nj_tree_outgroup, main = "NJ Tree with Bootstrap Support (outgroup rooted)", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
nodelabels(nj_bs_percent, cex = 0.6, frame = "none", adj = c(1.2, -0.5))
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

```

Neighbour Joining doesn't care about the evolutionary process - it just takes a distance matrix and clusters accordingly. NJ doesn't need to know how distances 



# Test evolutionary models
cat("Testing evolutionary models for protein sequences...\n")
model_test <- modelTest(phyDat_alignment, tree = tree_init, 
                        model = c("JTT", "WAG", "LG", "Dayhoff", "cpREV", "Blosum62"),
                        G = TRUE, I = TRUE)

# Display model test results
print(model_test)

# Get best model based on AIC
best_model <- model_test$Model[which.min(model_test$AIC)]
cat("\nBest model based on AIC:", best_model, "\n\n")
# Phylogenetics I

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

>Phylogenetic trees represent the evolutionary relationships among sequences. Closely related sequences are grouped together, and branch lengths indicate the amount of evolutionary change or genetic distance between them.

This means that to go from a multiple sequence alignment to a phylogenetic tree, we first need a way to compute the *pairwise* distances between each sequence in our alignment. There are many ways to create a distance matrix, and we will try a few, but the simplest one is to simply calculate a p-value describing the number of differences at each position in an alignment between two sequences. We will try this out in the following neighbour joining tree example.

## Neighbour Joining tree

Neighbour joining builds a tree that iteratively joins pairs of taxa that minimise the total branch length of the tree. To do this, NJ algorithms take in a distance matrix which contains pairwise evolutionary distances between taxa, calculated from a multiple sequence alignment. Let's have a go at creating one so we can see what this means along the way.

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

```

We are now at the point where we need to take our MSA and use this to calculate a distance matrix. We are going to use the dist.hamming() function from the `phangorn` package that calculates the number of differing positions in a pairwise alignment divided by the total number of positions - essentially a p-value.

```R

# Create a p-value based distance matrix
dm <- dist.hamming(phyDat_alignment)

# To inspect dm, convert to a matrix
dm_matrix <- as.matrix(dm)

```

Now we can build and inspect our first tree.

```R

# Build an NJ tree from our distance matrix
nj_tree <- NJ(dm)

# Visualise this initial tree
plot(nj_tree, main = "Neighbor-Joining Tree", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

```

This creates an `unrooted` tree which you will be able to see in your RStudio Plot pane. Have a look at the tree you've created to see if it makes evolutionary sense ie. do related sequences share a common ancestor? Next up, we'll root the tree by defining an outgroup which we spoke about in our previous tutorial. Remember that rooting a tree uses prior knowledge of your sequences to declare which is most divergent, prior to plotting. Rooting at different nodes preserves topology but changes the interpretation of evolutionary relationships.

![root](images/root.png)

```R

# Root the tree on your outgroup sequence
outgroup <- "Cin_Pax6_NP_001027641.1 homeobox transcription factor Pax6 [Ciona intestinalis]"
nj_tree_outgroup <- root(nj_tree, outgroup = outgroup, resolve.root = TRUE)

# Plot your rooted tree
plot(nj_tree_outgroup, main = "Neighbor-Joining Tree (outgroup rooted)", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

# Save your rooted NJ tree in Newick format
write.tree(nj_tree_outgroup, file = "trees/NJ/NJ_tree_outgroup.nwk"

```

We can also calculate the confidence level of each node in the tree by `bootstrapping` it. Bootstrapping is a statistical resampling method used to estimate the reliability of a result by repeatedly sampling from your original data with replacement. In this case, *replacement* means to randomly sample columns from your MSA, rather than using the original sequence order to see if these purturbations influence the final topology of the tree. The theory is that if your tree is **robust**, this process of resampling with replacement should produce a tree with the same topology more often than not. Topology referes to how each sample relates to one another. Here's an example:

Original aligment of 8 amino acid sequeces
```text

Position:          1  2  3  4  5  6  7  8
Sequence A:        M  K  R  D  L  M  R  E
Sequence B:        M  K  H  D  V  M  R  Q
Sequence C:        L  S  R  D  L  I  G  E
Sequence D:        L  S  H  N  V  I  G  Q

```

This creates a tree with the topology:

```text

Original tree topology: ((A,B),(C,D))

```
Now we do bootstrap replicate 1:
Randomly selected positions: 1, 1, 3, 5, 2, 7, 8, 4

```text

Position:          1  1  3  5  2  7  8  4
Sequence A:        M  M  R  L  K  R  E  D
Sequence B:        M  M  H  V  K  R  Q  D
Sequence C:        L  L  R  L  S  G  E  D
Sequence D:        L  L  H  V  S  G  Q  N

```

Now we recalculate topology:

```text

Bootstrap tree 1: ((A,B),(C,D)) ✓ Same topology

```

Now we do bootstrap replicate 2:

Randomly selected positions: 2, 4, 4, 6, 8, 3, 1, 1

```text

Position:         2  4  4  6  8  3  1  1
Sequence A:        K  D  D  M  E  R  M  M
Sequence B:        K  D  D  M  Q  H  M  M
Sequence C:        S  D  D  I  E  R  L  L
Sequence D:        S  N  N  I  Q  H  L  L

```

Calculate topology:

```text

>Bootstrap tree 2: (A,(B,(C,D))) ✗ Different topology!

```

After doing this 100 or 10000 times, we see how many trees we build from these resampled MSAs have the same or different topology.

```text

Clade (A,B):    appeared in 85 trees → Bootstrap support = 85%
Clade (C,D):    appeared in 92 trees → Bootstrap support = 92%
Clade ((A,B),(C,D)): appeared in 78 trees → Bootstrap support = 78%

```

And now we can add these bootstrap values to our original tree:

```text

                         78
            ┌────────────┴────────────┐
            │                         │
         85 │                      92 │
    ┌───────┴───────┐        ┌────────┴────────┐
    │               │        │                 │
Sequence A       Sequence B   Sequence C     Sequence D

```

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
  boot_dm <- dist.hamming(boot_phyDat)
  
  # Build tree
  boot_tree <- NJ(boot_dm)
  
  # Root with the same outgroup
  boot_tree <- root(boot_tree, outgroup = outgroup, resolve.root = TRUE)
  boot_trees[[i]] <- boot_tree
  
  # Monitor progress and report every time ten new bootstraps have been completed
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

### Neighbour Joining tree with Maximum Likelihood calculated distance matrix

This is all good and well, but we can do better than just using a p-value to calculate distances. While the p-distance is just the proportion of mismatched sites, we can use more sophisticated models to correct for multiple substitutions at the same site, unequal base frequencies, transition/transversion biases amongst other parameters. Prior to building our tree above, we created a distance matrix using `dist.hamming()`. Another option is to use `dist.ml()` which applies a *maximum likelihood* approach. This method asks: `Given a certain evolutionary distance, what's the probability of observing the differences we see in the sequences?` 

The reason that the answer to that question is not 100% is because not all substitutions are equally probable given an amount of evolutionary time. As we saw in our BLAST lecture and tutorial, the probability of one amino acid being substituted for another amino acid with similar chemical properties is is more likely than it being substituted for an amino acid with different chemical properties. These probabilities were described by the substitution matrix that we chose, the most common of which is `BLOSUM62`. 

Like the BLAST algorithm, phylogenetic analysis can take advantage of sophisticated substitution matrices to more accurately model evolutionary change between two sequences. When we run `dist.ml()`, at every position in the pairwise alignment, it calculates the probability that we see the two amino acids observed given an evolutionary distance of `d`, the units of which are `expected subsitutions per site`. It then sums up the log of these likelihood values for each position to calculate a `log likelihood` value for the whole pairwise alignment. The algorithm then recalculates this value multiple times using a range of vaulues of `d` so that in the end it can determine the value of `d` at which the log likelihood is maximised. It are these `maximum likelihood distance values (d)` that are used to populate the distance matrix.

Here's a visual representation:

```text

# Alignment
Position:    1  2  3  4  5
Sequence A:  M  K  R  D  L
Sequence B:  M  K  H  D  V
             ✓  ✓  ✗  ✓  ✗

# Key equations
P(i→j | distance d) = probability that amino acid i changes to j 
                      over evolutionary distance d

P(i=i | distance d) = probability amino acid i stays unchanged
                      over evolutionary distance d
                      
# Calculate likelihood for distance d = 0.1 (0.1 changes per site)
Position 1: M → M (no change)
P(M stays M | d=0.1) = 0.95  (high - short time, unlikely to change)

Position 2: K → K (no change)
P(K stays K | d=0.1) = 0.94

Position 3: R → H (changed)
P(R→H | d=0.1) = 0.008  (low - not much time for change)

Position 4: D → D (no change)
P(D stays D | d=0.1) = 0.95

Position 5: L → V (changed)
P(L→V | d=0.1) = 0.012  (low - similar amino acids but still unlikely in short time)

Total likelihood = 0.95 × 0.94 × 0.008 × 0.95 × 0.012 = 0.0000081
Log-likelihood = ln(0.0000081) = −9.42

# Calculate likelihood for distance d = 0.5
Position 1: M → M (no change)
P(M stays M | d=0.5) = 0.65  (still likely to stay same)

Position 2: K → K (no change)
P(K stays K | d=0.5) = 0.63

Position 3: R → H (changed)
P(R→H | d=0.5) = 0.042  (more reasonable - had time to change)

Position 4: D → D (no change)
P(D stays D | d=0.5) = 0.64

Position 5: L → V (changed)
P(L→V | d=0.5) = 0.055  (reasonable for similar amino acids)

Total likelihood = 0.65 × 0.63 × 0.042 × 0.64 × 0.055 = 0.000061
Log-likelihood = ln(0.000061) = -7.41

# Calculate likelihood for distance d = 1.0
Position 1: M → M (no change)
P(M stays M | d=1.0) = 0.35  (lower - had lots of time to change)

Position 2: K → K (no change)
P(K stays K | d=1.0) = 0.33

Position 3: R → H (changed)
P(R→H | d=1.0) = 0.068  (high - had time)

Position 4: D → D (no change)
P(D stays D | d=1.0) = 0.34

Position 5: L → V (changed)
P(L→V | d=1.0) = 0.082  (high - had time)

Total likelihood = 0.35 × 0.33 × 0.068 × 0.34 × 0.082 = 0.000022
Log-likelihood = ln(0.000022) = -8.42

# Plot the log-likelihood distance
Log-likelihood

     |
-7.5 |            ★ (maximum at d ≈ 0.5)
     |         /     \
-8.0 |        /       \
     |      /           \
-8.5 |     /             \
     |    /               \
-9.0 |   /                  \
     |  ●                    ●
-9.5 |/________________________\_
     0   0.2  0.4  0.6  0.8  1.0    d (substitutions/site)

     d=0.1        d=0.5        d=1.0
     ln(L)=-9.42  ln(L)=-7.41  ln(L)=-8.42
     
```

Ok, let's try create a maximum likelihood distance matrix and then use this to recreate our NJ tree. In this example, as we are using the [JTT](https://academic.oup.com/bioinformatics/article/8/3/275/193076?login=true) model which is an empirically derived model based on observed substitution rates.

```R

# Create a maximum likelihood distance matrix
dm <- dist.ml(phyDat_alignment, model = "JTT")

```

You can now go back to the `NJ(dm)` step above to use this new distance matrix to compare it to the old p-distance matrix we created first.

## Maximum Likelihood tree

Rather than using ML to calculate a distance matrix for a neighbour joining tree, we can use ML to produce the tree itself. This uses a completely different approach in which, rather than first calculating a distance matrix and then clustering branches to minimise branch lengths, we instead use a ML approach to evaluate tree topologies directly against the multiple sequence alignment. In the end, this selects the tree with the topology that maximises the likelihood of the alignment data. This is the workflow of producing a NJ tree using a ML distance matrix:

>Alignment  
>    ↓  
>Calculate ML distances (pairwise)  
>    ↓  
>Distance Matrix  
>    ↓  
>NJ clustering algorithm  
>    ↓  
>Tree  

and this is what the workflow of producing a ML tree:

>Alignment  
>    ↓  
>Propose tree topology  
>    ↓  
>Calculate likelihood of entire alignment given tree  
>    ↓  
>Try different topologies  
>    ↓  
>Find tree with maximum likelihood  

A logical question to ask is, what tree topology should we initially propose? There are several options including using a random topology or using a *star* topology in which every sequence radiates out of the centre, but most of the time, the best option is to first produce an NJ tree and use this to kick start your ML tree production. As we already have a NJ tree, we can just use that!

Once the initial tree is definied, the ML approach will iterate over small changes to this topology and each time, will assess the maximum likelihood of the entire tree. It will iterate over these potential topologies until it reaches **convergence**. That is to say, after reaching a particular point, more iterations does not lead to improved ML scores. Here's a visual representation of convergence:

```text

Log-likelihood over iterations:

ln L
     |
-1400|________________________________  ← Converged (flat)
     |                      ___/
-1420|                 ____/
     |            ____/
-1440|       ____/
     |   ___/
-1460|  /
     | /
-1480|/
     |
-1500|
     └────────────────────────────────> Iteration
      0    2    4    6    8   10   12

     Rapid          Slower        Flat → STOP
     improvement    gains         (converged)
     
```

Before we get to this though, we need to choose a model. Models differ because not every set of genes or proteins evolve in the same way and at the same rate. Let's run a test to see which model fits our dataset best.


```R

# Run the model test
model_test <- modelTest(phyDat_alignment, tree = nj_tree_outgroup, 
                        model = c("JTT", "WAG", "LG", "Dayhoff", "cpREV", "Blosum62"),
                        G = TRUE, I = TRUE)

# Select the model with the lowest AIC score as the best
best_model <- model_test$Model[which.min(model_test$AIC)]

```

This tests the following six models, however many, many more than this have been developed:

>JTT: General protein evolution model  
WAG: Emphasizes different amino acid frequencies  
LG: More recent, based on larger datasets  
Dayhoff: Older model, based on closely related proteins  
cpREV: Specifically for chloroplast proteins  
Blosum62: Based on conserved blocks in alignments  

It also tests if allowing rate variation (meaning it won't assume every site evolves at the same rate) will improve the result (G) or if allowing some sites to be invariant will improve the result (I). The most common metric to look at to determine the best model for your alignment is `AIC` which stands for *Akaike Information Criterion*. You don't really need to know the deep mathematics behind how it works, just know that lower is better!

Now that we know the best model for our data, we can calculate our ML tree.

```R

# Create pml object using best model
pml_obj <- pml(nj_tree_outgroup, phyDat_alignment, model = best_model)

# Optimise everything
ml_tree <- optim.pml(pml_obj, 
                     optNni = TRUE,      # Find best topology
                     optBf = TRUE,       # Optimize frequencies
                     optQ = TRUE,        # Optimize rates
                     optGamma = TRUE,    # Optimize gamma
                     optInv = TRUE,      # Optimize invariant
                     control = pml.control(trace = 1),
                     multicore = TRUE,   # Enable multicore
                     mc.cores = 4)       # Use 4 cores
                     
```

As with our NJ tree, we need to calculate bootstrap values to see how robust our tree is. After this, we can plot it.

```R

# Calculate bootstrap values
bs <- bootstrap.pml(ml_tree, 
                    bs = 1000,          # 100 bootstrap replicates
                    optNni = TRUE,      # Optimize topology for each
                    multicore = TRUE,   # Use multiple cores
                    mc.cores = 4)       # 4 cores

# Convert these to a percentage
ml_bs_percent <- round(prop.clades(ml_tree$tree, bs) / length(bs) * 100, 0)

# Plot tree with bootstrap values
plot(ml_tree$tree, main = "ML Tree with Bootstrap Support", 
     cex = 0.7, label.offset = 0.01, direction = "rightwards")
nodelabels(ml_bs_percent, cex = 0.6, frame = "none", adj = c(1.2, -0.5))
add.scale.bar(x = 0, y = 0.5, cex = 0.7, lwd = 2)

# Save tree with bootstrap values
ml_tree_bs <- ml_tree$tree
ml_tree_bs$node.label <- ml_bs_percent
write.tree(ml_tree_bs, file = "trees/ML/ML_tree_bootstrap.nwk")
cat("\nML tree with bootstrap values saved\n")

```
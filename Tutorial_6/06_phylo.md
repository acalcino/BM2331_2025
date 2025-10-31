# Phylogenetics II

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

```R
# Set working directory and load libraries
library(phangorn)
library(Biostrings)
library(ggtree)
library(ggplot2)

setwd("~/working-directory/phylogenetic_project")

# Load the NJ tree you created last week
load("trees/NJ/nj_ml_tree.RData")

# Also load your alignment and convert to phangorn format
alignment <- readAAStringSet("msa/msa_muscle_aligned_sequences.fa")
phyDat_alignment <- as.phyDat(alignment, type = "AA")

```

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
model_test <- modelTest(phyDat_alignment, tree = nj_ml_tree, 
                        model = c("JTT", "WAG", "LG", "Dayhoff", "cpREV", "Blosum62"),
                        G = TRUE, I = TRUE)

```

This tests the following six models, however many, many more than this have been developed:

>JTT: General protein evolution model  
WAG: Emphasizes different amino acid frequencies  
LG: More recent, based on larger datasets  
Dayhoff: Older model, based on closely related proteins  
cpREV: Specifically for chloroplast proteins  
Blosum62: Based on conserved blocks in alignments  

It also tests if allowing rate variation (meaning it won't assume every site evolves at the same rate) will improve the result (G) or if allowing some sites to be invariant will improve the result (I). The most common metric to look at to determine the best model for your alignment is `AIC` which stands for *Akaike Information Criterion*. You don't really need to know the deep mathematics behind how it works, just know that lower is better!

Now that we know the best model for our data, we can calculate our ML tree. I am also specifying `k = 4` and `inv = TRUE` here to reflect that my model test found that using `G` and `I` improved results.

```R

# If your tree is rooted, you first need to unroot it
unrooted_tree <- unroot(nj_ml_tree)

# Create pml object using best model
pml_obj <- pml(unrooted_tree, phyDat_alignment, model = "Dayhoff", k = 4, inv = 0.1)

# Optimise everything
ml_pml_obj <- optim.pml(pml_obj, 
                     optNni = TRUE,      # Find best topology
                     optBf = TRUE,       # Optimize frequencies
                     optGamma = TRUE,    # Optimize gamma
                     optInv = TRUE,      # Optimize invariant
                     control = pml.control(trace = 1),
                     multicore = TRUE,   # Enable multicore
                     mc.cores = 4)       # Use 4 cores
                     
```

As with our NJ tree, we need to calculate bootstrap values to see how robust our tree is. After this, we can plot it.

```R

# Calculate bootstrap values using pbmcapply for multithreaded support progress bar output

library(pbmcapply)

# Start of with 1000 bootstraps and 1 cores
bs_results <- pbmclapply(1:1000, function(i) {
              bootstrap.pml(ml_pml_obj, bs = 1, 
              optNni = TRUE)
              }, 
              mc.cores = 1)
              
# If this is going to take too long, kill the job and think about how you might make it more efficient.

# Combine the results and calculate percentage
bs_combined <- do.call(c, bs_results)
ml_bs_percent <- round(prop.clades(ml_pml_obj$tree, bs_combined) / length(bs_combined) * 100, 0)

# Extract tree from pml object and attache bootstrap labels
ml_tree <- ml_pml_obj$tree
ml_tree$node.label <- as.character(ml_bs_percent)

# Plot tree with bootstrap values
ml_tree_plot <- ggtree(ml_tree) +
  geom_tiplab(size = 3, offset = 0.01) +
  geom_nodelab(size = 2.5, hjust = -0.2, vjust = 0.5) +
  geom_treescale(x = 0, y = -0.5, width = 1, linesize = 2, fontsize = 3) +
  hexpand(.4) +
  ggtitle("ML Tree with Bootstrap Support")
    
ml_tree_plot

# Root the tree and plot again

outgroup <- "Cin_Pax6_NP_001027641.1 homeobox transcription factor Pax6 [Ciona intestinalis]"

# Root the tree using the outgroup
rooted_tree <- root(ml_tree, outgroup = outgroup, resolve.root = TRUE)

# Plot the rooted tree
ml_bootstrapped_rooted <- ggtree(rooted_tree) +
    geom_tiplab(size = 3, offset = 0.01) +
    geom_nodelab(aes(label = label),
                 size = 2.5,
                 hjust = -0.2,
                 vjust = 0.5) +
    geom_treescale(x = 0, y = -0.5, 
                   width = 1, 
                   linesize = 2,
                   fontsize = 3) +
    hexpand(.4) +
    ggtitle("ML Tree with Bootstrap Support (outgroup rooted)")

ml_bootstrapped_rooted

# If you find your tip labels too long, shorten them by removing everything after the first space
rooted_tree$tip.label <- gsub(" .*", "", rooted_tree$tip.label)

# Save as png
ggsave("trees/ML/ml_bootstrapped_rooted.png", plot = ml_bootstrapped_rooted, width = 12, height = 8, dpi = 300)

# Save tree with bootstrap values in Newick format
write.tree(rooted_tree, file = "trees/ML/ml_bootstrapped_rooted_tree.nwk")

```

ggtree is a powerful package that lets you modify your trees is many interesting ways. Feel free to explore everything ggtree offers [here](https://yulab-smu.top/treedata-book/). If you'd prefer a graphical interface, there are stand alone programs you can use One of these is [TreeViewer](https://treeviewer.org/) which is available for Linux, Mac and, ughh, Windows too. Once you've loaded your tree, you'll need to click the `Modules` tab, then `Plot actions`. Under `Labels`, change `Attribute` to `Support` and make sure `Attribute type` is set as `Number`.

![treeviewer](images/treeviewer.jpg)

TreeViewer also has a bajillion ways to customise your tree and you can read about it all [here](https://github.com/arklumpus/TreeViewer/wiki). Once you have finished making changes, click `File`, `Export` to save it as a pdf, png or svg file. You can then import this into your RMarkdown document if you wish.
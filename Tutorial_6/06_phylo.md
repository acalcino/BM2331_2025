# Phylogenetics II

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
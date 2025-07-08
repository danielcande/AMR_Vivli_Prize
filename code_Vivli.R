#loading packages
library(tidyverse)
library(poLCA)
library(nnet)
library(gt)
library(gtsummary)

set.seed(123)

#loading data
data <- read.csv("AMR_Atlas_data.csv")

data <- data %>%
  left_join(country_to_region_lookup, by = "Country")

e_faecium <- data %>%
  filter(Species == "Enterococcus faecium")

s_aureus <- data %>%
  filter(Species == "Staphylococcus aureus")

k_pneumoniae <- data %>%
  filter(Species == "Klebsiella pneumoniae")

a_baumannii <- data %>%
  filter(Species == "Acinetobacter baumannii")

a_baumannii_genetics <- a_baumannii %>% #Extracting the genetic variables
  dplyr::select(Isolate.Id, 112:135) %>%
  mutate(across(where(is.character), ~na_if(.,"")))

p_aeruginosa <- data %>%
  filter(Species == "Pseudomonas aeruginosa")

p_aeruginosa_genetics <- data %>%
  filter(Species == "Pseudomonas aeruginosa") %>%
  dplyr::select(Isolate.Id, 112:134)

e_spp <- data %>%
  filter(Species == "Enterobacter spp")

#e_coli <- data %>%
#  filter(Species == "Escherichia coli")

# Replace all empty strings "" with NA across all columns
a_baumannii <- a_baumannii %>%
  mutate(across(where(is.character), ~na_if(.,"")))

e_faecium <- e_faecium %>%
  mutate(across(where(is.character), ~na_if(., "")))

e_spp <- e_spp %>%
  mutate(across(where(is.character), ~na_if(., "")))

k_pneumoniae <- k_pneumoniae %>%
  mutate(across(where(is.character), ~na_if(., "")))

p_aeruginosa <- p_aeruginosa %>%
  mutate(across(where(is.character), ~na_if(., "")))

s_aureus <- s_aureus %>%
  mutate(across(where(is.character), ~na_if(., "")))

#e_coli <- e_coli %>%
#  mutate(across(where(is.character), ~na_if(., "")))

#Removing columns with >95% missingness

a_baumannii <- a_baumannii %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))

e_faecium <- e_faecium %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))

e_spp <- e_spp %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))

k_pneumoniae <- k_pneumoniae %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))


p_aeruginosa <- p_aeruginosa %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))

s_aureus <- s_aureus %>%
  dplyr::select(where(~ mean(is.na(.)) < 0.95))

#e_coli <- e_coli %>%
#  dplyr::select(where(~ mean(is.na(.)) < 0.95))

#Remove the antibiotic columns that don't show resistance
a_baumannii <- a_baumannii %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

e_faecium <- e_faecium %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

e_spp <- e_spp %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

k_pneumoniae <- k_pneumoniae %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

p_aeruginosa <- p_aeruginosa %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

s_aureus <- s_aureus %>%
  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

#e_coli <- e_coli %>%
#  dplyr::select( !where(~is.character(.) && any(str_detect(., "<|>"), na.rm = TRUE)))

#Trying running poLCA per bug - let's do 5 clusters and see what happens - I'll do model selection later

#STep 1: defining antibiotic names
a_baumannii_antibiotics <- a_baumannii %>%
  dplyr::select(14:23,) %>%
  colnames() %>%
  gsub("_I$","",.)

e_faecium_antibiotics <- e_faecium %>%
  dplyr::select(13:23,) %>%
  colnames() %>%
  gsub("_I$","",.)


e_spp_antibiotics <- e_spp %>%
  dplyr::select(13:35,) %>%
  colnames() %>%
  gsub("_I$","",.)

k_pneumoniae_antibiotics <- k_pneumoniae %>%
  dplyr::select(14:40,) %>%
  colnames() %>%
  gsub("_I$","",.)

p_aeruginosa_antibiotics <- p_aeruginosa %>%
  dplyr::select(14:26,) %>%
  colnames() %>%
  gsub("_I$","",.)

s_aureus_antibiotics <- s_aureus %>%
  dplyr::select(14:27,) %>%
  colnames() %>%
  gsub("_I$","",.)

#e_coli_antibiotics <- e_coli %>%
#  dplyr::select(14:40,) %>%
#  colnames() %>%
#  gsub("_I$","",.)

#Removing suffix from the main df for each bug
colnames(a_baumannii) <- gsub("_I$","",colnames(a_baumannii))
colnames(e_faecium) <- gsub("_I$","",colnames(e_faecium))
colnames(e_spp) <- gsub("_I$","",colnames(e_spp))
colnames(k_pneumoniae) <- gsub("_I$","",colnames(k_pneumoniae))
colnames(p_aeruginosa) <- gsub("_I$","",colnames(p_aeruginosa))
colnames(s_aureus) <- gsub("_I$","",colnames(s_aureus))
#colnames(e_coli) <- gsub("_I$","",colnames(e_coli))

#Changing all columns (besides ID) to factor
a_baumannii <- a_baumannii %>%
  mutate(across(where(is.character), ~as.factor(.)))

e_faecium <- e_faecium %>%
  mutate(across(where(is.character), ~as.factor(.)))

e_spp <- e_spp %>%
  mutate(across(where(is.character), ~as.factor(.)))

k_pneumoniae <- k_pneumoniae %>%
  mutate(across(where(is.character), ~as.factor(.)))

p_aeruginosa <- p_aeruginosa %>%
  mutate(across(where(is.character), ~as.factor(.)))

s_aureus <- s_aureus %>%
  mutate(across(where(is.character), ~as.factor(.)))

#e_coli <- e_coli %>%
#  mutate(across(where(is.character), ~as.factor(.)))

#Creating formula objects for poLCA
f_a_baumannii <- as.formula(paste("cbind(", paste(a_baumannii_antibiotics, collapse = ", "), ") ~ 1"))
f_e_faecium <- as.formula(paste("cbind(", paste(e_faecium_antibiotics, collapse = ", "), ") ~ 1"))
f_e_spp <- as.formula(paste("cbind(", paste(e_spp_antibiotics, collapse = ", "), ") ~ 1"))
f_k_pneumoniae <- as.formula(paste("cbind(", paste(k_pneumoniae_antibiotics, collapse = ", "), ") ~ 1"))
f_p_aeruginosa <- as.formula(paste("cbind(", paste(p_aeruginosa_antibiotics, collapse = ", "), ") ~ 1"))
f_s_aureus <- as.formula(paste("cbind(", paste(s_aureus_antibiotics, collapse = ", "), ") ~ 1"))
#f_e_coli <- as.formula(paste("cbind(",paste(e_coli_antibiotics, collapse = ", "), ") ~ 1"))
#Extracting the antibiograms
a_baumannii_abx <- a_baumannii %>%
  dplyr::select(13:27)

e_faecium_abx <- e_faecium %>%
  dplyr::select(13:23)

e_spp_abx <- e_spp %>%
  dplyr::select(13:35)

k_pneumoniae_abx <- k_pneumoniae %>%
  dplyr::select(14:40)

p_aeruginosa_abx <- p_aeruginosa %>%
  dplyr::select(14:26)

s_aureus_abx <- s_aureus %>%
  dplyr::select(14:27)

#e_coli_abx <- e_coli %>%
#  dplyr::select(14:40)


##THIS BLOCK WAS USED FOR TESTING THE SUMMER PROJECT MODEL
# s_aureus_abx <- s_aureus_abx %>%
#   rename(Clindamycine = Clindamycin_I, Erythromycine = Erythromycin_I, Linezolide = Linezolid_I,
#          TigÃ©cycline = Tigecycline_I, Vancomycine = Vancomycin_I, Daptomycine = Daptomycin_I, 
#          Gentamicine = Gentamicin_I, Oxacilline = Oxacillin_I, Teicoplanine = Teicoplanin_I, Trimethoprime = Trimethoprim.sulfa_I)
# 
# s_aureus_abx <- s_aureus_abx %>%
#   dplyr::select(-c(Levofloxacin_I, Minocycline_I, Ceftaroline_I, Moxifloxacin_I))
# 
# #keep complete cases
# s_aureus_abx <- s_aureus_abx %>%
#   filter(complete.cases(.))
# 
# s_aureus_core_abx_name <- colnames(s_aureus_abx)
# 
# convert_to_numeric_core  <- function(x) {
#   recode(x, "Resistant" = 5L, "Intermediate" = 4L, "Susceptible" = 3L,"NT"=2L,"Autre"=1L)
# }
# 
# s_aureus_abx[] <- lapply(s_aureus_abx, convert_to_numeric_core)
# 
# #Running new poLCA on original data only containing the core antibiotics that are common to both datasets
# f_core_staph <- as.formula(paste("cbind(", paste(s_aureus_core_abx_name, collapse = ", "), ") ~ 1"))
# 
# LCA_core_staph_3 <- poLCA(f_core_staph, drug_data_staph_clean, na.rm = FALSE, nclass = 3)
# 
# #run validation on Atlas data
# 
# posterior_probs <- poLCA.posterior(lc = LCA_core_staph_3,
#                                    y = s_aureus_abx)
# 
# # --- STEP 2: Assign hard classifications ---
# 
# # Find the most likely class for each isolate in the external data
# predicted_classes <- apply(posterior_probs, 1, which.max)
# 
# # Add this new classification to the external dataset
# s_aureus_abx <- s_aureus_abx %>%
#   mutate(Predicted_Cluster = cbind(predicted_classes))
# 
# s_aureus_abx <- s_aureus_abx %>%
#   filter(Predicted_Cluster == (1 | 2 | 3))
# 
# s_aureus_abx$Predicted_Cluster <- unlist(s_aureus_abx$Predicted_Cluster)
# 
# # --- STEP 3: EVALUATION ---
# 
# # 1. Check the distribution of the new classifications
# #    This shows you how many isolates from the external data fall into each of your original classes.
# print("Distribution of predicted classes in external data:")
# print(table(s_aureus_abx$Predicted_Cluster))
# 
# 
# # 2. Check the coherence of the class profiles
# #    Do the resistance patterns still match?
# 
# # Get the original model's resistance profiles (the "ground truth")
# original_profiles <- LCA_core_staph_3$probs
# 
# # Calculate the OBSERVED resistance profiles for the newly classified external data
# # This dplyr pipe calculates the % resistance for each antibiotic, for each new cluster
# validation_profiles <- s_aureus_abx %>%
#   # Group by the new cluster assignments
#   group_by(Predicted_Cluster) %>%
#   
#   # The fix is to place the function inside across() as the .fns argument
#   summarise(
#     across(
#       all_of(s_aureus_core_abx_name),
#       .fns = ~ mean(.x == 5, na.rm = TRUE)
#     )
#   )
# 
# # Now, view the result
# print(validation_profiles)
# 
# print("Original model profiles:")
# print(original_profiles)
# 
# print("Observed profiles in the externally validated data:")
# print(validation_profiles)
# 
# #----------SAME FOR ECOLI-------------
# 
# e_coli_abx <- e_coli_abx %>%
#            rename(
#              Amikacine = Amikacin,
#              `Amoxicilline.ac..clav.` = Amoxycillin.clavulanate, # Use backticks for special characters
#              Ampicilline = Ampicillin,
#              `CÃ©fÃ©pime` = Cefepime,
#              Ceftazidime = Ceftazidime,
#              `CÃ©ftriaxone` = Ceftriaxone,
#              `ImipÃ©nÃ¨me` = Imipenem,
#              `LÃ©vofloxacine` = Levofloxacin,
#              `PipÃ©racilline.tazob.` = `Piperacillin.tazobactam`,
#              `AztrÃ©onam` = Aztreonam,
#              Cefixime = Cefixime,
#              `Ceftazidime.avibactam` = `Ceftazidime.avibactam`,
#              Ciprofloxacine = Ciprofloxacin,
#              Colistine = Colistin,
#              Ertapeneme = Ertapenem,
#              Gentamicine = Gentamicin,
#              Cotrimoxazole = `Trimethoprim.sulfa`
#            )
# 
# #Remove columns with no match
# 
# e_coli_abx <- e_coli_abx %>%
#   dplyr::select(-c(Meropenem, Minocycline, Tigecycline, Ampicillin.sulbactam,
#             Ceftaroline, Doripenem, Ceftolozane.tazobactam, Meropenem.vaborbactam,
#             Cefpodoxime, Ceftibuten))
# 
# 
# #keep complete cases
# #e_coli_abx <- e_coli_abx %>%
# #  filter(complete.cases(.))
# 
# e_coli_abx_name <- colnames(e_coli_abx)
# 
# convert_to_numeric_core  <- function(x) {
#   recode(x, "Resistant" = 5L, "Intermediate" = 4L, "Susceptible" = 3L,"NT"=2L,"Autre"=1L)
# }
# 
# e_coli_abx[] <- lapply(e_coli_abx, convert_to_numeric_core)
# 
# #Running new poLCA on original data only containing the core antibiotics that are common to both datasets
# f_core_ecoli <- as.formula(paste("cbind(", paste(e_coli_abx_name, collapse = ", "), ") ~ 1"))
# 
# LCA_core_ecoli_5 <- poLCA(f_core_ecoli, drug_data_ecoli_clean, na.rm = FALSE, nclass = 5)
# 
# #run validation on Atlas data
# 
# posterior_probs_ecoli <- poLCA.posterior(lc = LCA_core_ecoli_5,
#                                    y = e_coli_abx)
# 
# # --- STEP 2: Assign hard classifications ---
# 
# # Find the most likely class for each isolate in the external data
# predicted_classes_ecoli <- apply(posterior_probs_ecoli, 1, which.max)
# 
# # Add this new classification to the external dataset
# e_coli_abx <- e_coli_abx %>%
#   mutate(Predicted_Cluster = cbind(predicted_classes_ecoli))
# 
# e_coli_abx <- e_coli_abx %>%
#   filter(Predicted_Cluster == (1 | 2 | 3 | 4 | 5))
# 
# e_coli_abx$Predicted_Cluster <- unlist(e_coli_abx$Predicted_Cluster)
# 
# # --- STEP 3: EVALUATION ---
# 
# # 1. Check the distribution of the new classifications
# #    This shows you how many isolates from the external data fall into each of your original classes.
# print("Distribution of predicted classes in external data:")
# print(table(e_coli_abx$Predicted_Cluster))
# 
# 
# # 2. Check the coherence of the class profiles
# #    Do the resistance patterns still match?
# 
# # Get the original model's resistance profiles (the "ground truth")
# original_profiles <- LCA_core_ecoli_5$probs
# 
# # Calculate the OBSERVED resistance profiles for the newly classified external data
# # This dplyr pipe calculates the % resistance for each antibiotic, for each new cluster
# validation_profiles <- e_coli_abx %>%
#   # Group by the new cluster assignments
#   group_by(Predicted_Cluster) %>%
#   
#   # The fix is to place the function inside across() as the .fns argument
#   summarise(
#     across(
#       all_of(e_coli_abx_name),
#       .fns = ~ mean(.x == 5, na.rm = TRUE)
#     )
#   )
# 
# # Now, view the result
# print(validation_profiles)
# 
# print("Original model profiles:")
# print(original_profiles)
# 
# print("Observed profiles in the externally validated data:")
# print(validation_profiles)

#validation_clusters <- predict(LCA_core_staph_3, newdata = s_aureus_abx)

# Extract parameters from your fitted LCA model
# lca_params <- LCA_core_staph_3$probs  # Item-response probabilities
# class_probs <- LCA_core_staph_3$P     # Class probabilities
# 
# # Function to calculate posterior probabilities for new data
# # More robust function with error checking
# assign_clusters_robust <- function(new_data, lca_model) {
#   # Extract basic info with checks
#   n_classes <- 3
#   n_obs <- nrow(new_data)
#   
#   print(paste("Classes:", n_classes, "Observations:", n_obs))
#   
#   if (is.na(n_classes) || n_classes <= 0) {
#     stop("Invalid number of classes")
#   }
#   
#   # Get variable names from original model
#   manifest_vars <- names(lca_model$probs)
#   
#   # Check which variables are available
#   available_vars <- intersect(manifest_vars, names(new_data))
#   print(paste("Available variables:", length(available_vars)))
#   
#   # Initialize matrix
#   posterior_probs <- matrix(0, nrow = n_obs, ncol = n_classes)
#   
#   # Calculate posteriors
#   for (i in 1:n_obs) {
#     for (k in 1:n_classes) {
#       likelihood <- 1
#       
#       for (var in available_vars) {
#         obs_value <- new_data[i, var]
#         if (!is.na(obs_value)) {
#           # poLCA stores probs as [response_level, class]
#           # Assuming binary coding (0/1), so add 1 for indexing
#           response_idx <- as.numeric(obs_value) + 1
#           if (response_idx %in% 1:2) {  # Valid response
#             likelihood <- likelihood * lca_model$probs[[var]][response_idx, k]
#           }
#         }
#       }
#       
#       posterior_probs[i, k] <- likelihood * lca_model$P[k]
#     }
#   }
#   
#   # Normalize
#   row_sums <- rowSums(posterior_probs)
#   posterior_probs <- posterior_probs / row_sums
#   
#   # Assign classes
#   predicted_class <- apply(posterior_probs, 1, which.max)
#   
#   return(list(
#     predicted_class = predicted_class,
#     posterior_probs = posterior_probs,
#     certainty = apply(posterior_probs, 1, max)
#   ))
# }
# 
# # Try the robust version
# validation_results <- assign_clusters_robust(s_aureus_abx, LCA_core_staph_3)
# 

#Running LCA (standard - no covariates)
#Starting with 5 clusters (arbitrary) WILL NEED BIC/AIC check
a_baumannii_LCA <- poLCA(f_a_baumannii, a_baumannii_abx, na.rm = FALSE, nclass = 4)
e_faecium_LCA <- poLCA(f_e_faecium, e_faecium_abx, nclass = 5, na.rm = FALSE)
e_spp_LCA <- poLCA(f_e_spp, e_spp_abx, nclass = 5, na.rm = FALSE)
k_pneumoniae_LCA <- poLCA(f_k_pneumoniae, k_pneumoniae_abx, nclass = 5, na.rm = FALSE)
p_aeruginosa_LCA <- poLCA(f_p_aeruginosa, p_aeruginosa_abx, na.rm = FALSE, nclass = 5)
s_aureus_LCA <- poLCA(f_s_aureus, s_aureus_abx, na.rm=FALSE, nclass = 4)


#Properly checking which model is a best fit - function initialisation

# Function to calculate VLMR test
calculate_vlmr <- function(model_k_minus_1, model_k) {
  # Log-likelihoods
  ll_k_minus_1 <- model_k_minus_1$llik
  ll_k <- model_k$llik
  
  # Number of parameters
  npar_k_minus_1 <- model_k_minus_1$npar
  npar_k <- model_k$npar
  
  # Sample size
  n <- model_k$N
  
  # VLMR test statistic
  vlmr_stat <- -2 * (ll_k_minus_1 - ll_k)
  
  # Degrees of freedom (difference in parameters)
  df <- npar_k - npar_k_minus_1
  
  # Approximate p-value (chi-square distribution)
  p_value <- 1 - pchisq(vlmr_stat, df)
  
  return(list(
    statistic = vlmr_stat,
    df = df,
    p_value = p_value,
    ll_k_minus_1 = ll_k_minus_1,
    ll_k = ll_k
  ))
}


# Function to run complete model comparison
compare_lca_models <- function(data, formula, max_classes = 6) {
  
  results <- data.frame(
    classes = 1:max_classes,
    loglik = NA,
    aic = NA,
    bic = NA,
    vlmr_stat = NA,
    vlmr_p = NA,
    entropy = NA
  )
  
  models <- list()
  
  # Fit models
  for (k in 1:max_classes) {
    cat("Fitting", k, "class model...\n")
    models[[k]] <- poLCA(formula, data, nclass = k, verbose = FALSE, na.rm= FALSE)
    
    results$loglik[k] <- models[[k]]$llik
    results$aic[k] <- models[[k]]$aic
    results$bic[k] <- models[[k]]$bic
    
    # Calculate entropy (classification quality)
    if (k > 1) {
      probs <- models[[k]]$posterior
      entropy <- -sum(probs * log(probs), na.rm = TRUE) / nrow(probs)
      results$entropy[k] <- 1 - (entropy / log(k))
    }
  }
  
  # Calculate VLMR tests
  for (k in 2:max_classes) {
    vlmr <- calculate_vlmr(models[[k-1]], models[[k]])
    results$vlmr_stat[k] <- vlmr$statistic
    results$vlmr_p[k] <- vlmr$p_value
  }
  
  return(list(results = results, models = models))
}

# Run for E. coli and S. aureus
a_baumannii_comparison <- compare_lca_models(a_baumannii_abx, f_a_baumannii, max_classes = 4) #4 classes
#e_coli_comparison <- compare_lca_models(e_coli_abx, f_e_coli, max_classes = 20)
e_faecium_comparison <- compare_lca_models(e_faecium_abx, f_e_faecium, max_classes = 5) #5 clusters
e_spp_comparison <- compare_lca_models(e_spp_abx, f_e_spp, max_classes = 5) #5 clusters
k_pneumoniae_comparison <- compare_lca_models(k_pneumoniae_abx, f_k_pneumoniae, max_classes = 7) #5 clusters
p_aeruginosa_comparison <- compare_lca_models(p_aeruginosa_abx, f_p_aeruginosa, max_classes = 7)
s_aureus_comparison <- compare_lca_models(s_aureus_abx, f_s_aureus, max_classes = 7) #5 clusters


# Visualise results

# Plot information criteria
a_baumannii_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - A. baumannii",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

e_faecium_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - E. faecium",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

e_spp_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - E. spp",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

k_pneumoniae_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - K. pneumoniae",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

p_aeruginosa_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - P. aeruginosa",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

s_aureus_comparison$results %>%
  dplyr::select(classes, aic, bic) %>%
  pivot_longer(cols = c(aic, bic), names_to = "criterion", values_to = "value") %>%
  ggplot(aes(x = classes, y = value, color = criterion)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Selection Criteria - S. aureus",
       x = "Number of Classes", y = "Information Criterion") +
  theme_minimal()

# Plot VLMR p-values
a_baumannii_comparison$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - E. coli",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

e_faecium_comparsion$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - S. aureus",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

e_spp_comparison$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - S. aureus",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

k_pneumoniae_comparison$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - S. aureus",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

p_aeruginosa_comparsion$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - S. aureus",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

s_aureus_comparsion$results %>%
  filter(classes > 1) %>%
  ggplot(aes(x = classes, y = vlmr_p)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "VLMR Test p-values - S. aureus",
       x = "Number of Classes", y = "VLMR p-value") +
  theme_minimal()

a_baumannii_model <- a_baumannii_comparison$models[[4]]
e_faecium_model <- e_faecium_comparison$models[[5]]
e_spp_model <- e_spp_comparison$models[[5]]
k_pneumoniae_model <- k_pneumoniae_comparison$models[[5]]
p_aeruginosa_model <- p_aeruginosa_comparison$models[[5]]
s_aureus_model <- s_aureus_comparison$models[[4]]

#Adding class to bug dfs
a_baumannii <- a_baumannii %>%
  cbind(as.factor(a_baumannii_model$predclass)) %>%
  rename(Cluster= `as.factor(a_baumannii_model$predclass)`)

e_faecium <- e_faecium %>%
  cbind(as.factor(e_faecium_model$predclass)) %>%
  rename(Cluster= `as.factor(e_faecium_model$predclass)`)

e_spp <- e_spp %>%
  cbind(as.factor(e_spp_model$predclass)) %>%
  rename(Cluster= `as.factor(e_spp_model$predclass)`)

k_pneumoniae <- k_pneumoniae %>%
  cbind(as.factor(k_pneumoniae_model$predclass)) %>%
  rename(Cluster=`as.factor(k_pneumoniae_model$predclass)`)

p_aeruginosa <- p_aeruginosa %>%
  cbind(as.factor(p_aeruginosa_model$predclass)) %>%
  rename(Cluster=`as.factor(p_aeruginosa_model$predclass)`)

s_aureus <- s_aureus %>%
  cbind(as.factor(s_aureus_model$predclass)) %>%
  rename(Cluster=`as.factor(s_aureus_model$predclass)`)


#Changing clusters to factors with descriptive names
a_baumannii$Cluster <- factor(a_baumannii$Cluster, levels = c(1,2,3,4),
                              labels = c("Pan-susceptible","Non-CRAB MDRAB","CRAB MDRAB","Low-Level Resistance"))


e_faecium$Cluster <- factor(e_faecium$Cluster, levels = c(1,2,3,4,5), 
                            labels = c("VS-MDR (Vancomycin-susceptible)","VR-MR (Vancomycin-resistant",
                                       "Amp-S, Van-S, E-R (Ampicillin-Susceptible, Erythromycin-Resistant Cluster",
                                       "Pan-susceptible","High-Level Penicillin-Resistant VSE"))

e_spp$Cluster <- factor(e_spp$Cluster, levels = c(1,2,3,4,5), 
                        labels= c("XDR-CRE","ESBL, Carbapenem-S","MDR, Carbapenem-S, (AmpC-like)",
                                  "Amp-R, ESBL-Negative","Pan-Susceptible"))

k_pneumoniae$Cluster <- factor(k_pneumoniae$Cluster, levels=c(1,2,3,4,5),
                               labels=c("CRE XDR","ESBL-like MDR","AmpC-like MDR","Susceptible(Intrinsic Amp-R)",
                                        "non-CRE XDR"))  #Critical cluster is CRE XDR, high risk is non-CRE XDR, ref is Intrinsic Amp-R 

p_aeruginosa$Cluster <- factor(p_aeruginosa$Cluster, levels=c(1,2,3,4,5),
                               labels=c("Pan-S","Fluoroq-R","MDR","CRP","XDR-CRP"))

s_aureus$Cluster <- factor(s_aureus$Cluster, levels=c(1,2,3,4),
                           labels=c("MRSA(Gentamicin-R","MRSA(Fluoroq-R)","MSSA","MRSA(Gentamicin-S)"))

#Adding genetic data
a_baumannii <- a_baumannii %>%
  left_join(a_baumannii_genetics, by="Isolate.Id") %>%
  janitor::remove_empty(which = "cols") %>% #removing columns with only NAs
  mutate(across(c(SHV, TEM, VEB, PER, GES, KPC, OXA, NDM, IMP, VIM, SPM, GIM),as.factor))

#Running regression with class as output

#Set reference class
a_baumannii$Cluster <- relevel(a_baumannii$Cluster, ref = "Pan-susceptible")
e_faecium$Cluster <- relevel(e_faecium$Cluster, ref = "Pan-susceptible")
e_spp$Cluster <- relevel(e_spp$Cluster, ref = "Pan-Susceptible")
p_aeruginosa$Cluster <- relevel(p_aeruginosa$Cluster, ref = "Pan-S")
k_pneumoniae$Cluster <- relevel(k_pneumoniae$Cluster, ref = "Susceptible(Intrinsic Amp-R)")
s_aureus$Cluster <- relevel(s_aureus$Cluster, ref = "MSSA")

a_baumannii_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                                   Speciality + Source + In...Out.Patient + Year, data = a_baumannii) #Removed genetic vars with <2 lvl

e_faecium_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                                 Speciality + Source + In...Out.Patient + Year, data = e_faecium)

e_spp_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                             Speciality + Source + Year, data = e_spp) #removed in/outpatient as they're ll "Non-Given"

s_aureus_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                             Speciality + Source + In...Out.Patient + Year, data = s_aureus)

k_pneumoniae_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                                Speciality + Source + In...Out.Patient + Year, data = k_pneumoniae)

p_aeruginosa_multinom <- multinom(Cluster ~ Super_Region + Gender + Age.Group + 
                                    Speciality + Source + In...Out.Patient + Year, data = p_aeruginosa)

# 
# #Extracting odds ratios and p-values
# 
# calc_odds_ratios <- function(model) {
#   coeffs <- model$coefficients
#   std_errs <- model$standard.errors
#   
#   # Calculate z-scores
#   z_scores <- coeffs / std_errs
#   
#   # Calculate two-tailed p-values
#   p_values <- (1 - pnorm(abs(z_scores))) * 2
#   
#   
#   # Calculate the RRRs by exponentiating the coefficients
#   OR <- exp(coef(model))
# return(OR)
#   }
# 
# a_baumannii_OR <- calc_odds_ratios(a_baumannii_multinom)
# e_faecium_OR <- calc_odds_ratios(e_faecium_multinom)
# e_spp_OR <- calc_odds_ratios(e_spp_multinom)
# p_aeruginosa_OR <- calc_odds_ratios(p_aeruginosa_multinom)
# k_pneumoniae_OR <- calc_odds_ratios(k_pneumoniae_multinom)
# s_aureus_OR <- calc_odds_ratios(s_aureus_multinom)


get_multinom_results_fast <- function(model) {
  
  # --- Step 1: Extract coefficients and standard errors ---
  # This is the core information from the model summary
  summary_model <- summary(model)
  coeffs <- summary_model$coefficients
  std_errs <- summary_model$standard.errors
  
  # --- Step 2: Manually calculate p-values and confidence intervals ---
  
  # Z-scores
  z_scores <- coeffs / std_errs
  
  # Two-tailed p-values
  p_values <- (1 - pnorm(abs(z_scores))) * 2
  
  # Confidence intervals on the log-odds scale
  conf_low_log <- coeffs - 1.96 * std_errs
  conf_high_log <- coeffs + 1.96 * std_errs
  
  # --- Step 3: Convert matrices to tidy data frames and combine ---
  
  # Helper function to convert a matrix to a long-format tibble
  reshape_matrix <- function(mat, value_name) {
    as.data.frame(mat) %>%
      tibble::rownames_to_column(var = "y.level") %>%
      pivot_longer(
        cols = -y.level,
        names_to = "term",
        values_to = value_name
      )
  }
  
  # Reshape all our calculated matrices
  tidy_coeffs <- reshape_matrix(coeffs, "estimate_log")
  tidy_pvals <- reshape_matrix(p_values, "p.value")
  tidy_conf_low <- reshape_matrix(conf_low_log, "conf.low_log")
  tidy_conf_high <- reshape_matrix(conf_high_log, "conf.high_log")
  
  # Join them all together into one final data frame
  results_df <- tidy_coeffs %>%
    left_join(tidy_pvals, by = c("y.level", "term")) %>%
    left_join(tidy_conf_low, by = c("y.level", "term")) %>%
    left_join(tidy_conf_high, by = c("y.level", "term"))
  
  # --- Step 4: Exponentiate and finalize ---
  
  final_results <- results_df %>%
    # Exponentiate the log-odds and CIs to get Odds Ratios
    mutate(
      estimate = exp(estimate_log),
      conf.low = exp(conf.low_log),
      conf.high = exp(conf.high_log),
      is_significant = p.value < 0.05
    ) %>%
    # Filter out the intercept term for cleaner plotting
    filter(term != "(Intercept)") %>%
    # Select and reorder columns to match the broom output
    select(y.level, term, estimate, conf.low, conf.high, p.value, is_significant)
  
  return(final_results)
}

# Use the new, faster function on your model
k_pneumoniae_results_fast <- get_multinom_results_fast(k_pneumoniae_multinom)

ggplot(
  data = k_pneumoniae_results_fast, 
  aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = is_significant)
) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  facet_wrap(~ y.level, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Significant (p < 0.05)") +
  labs(
    title = "Forest Plot of Multinomial Regression for K. pneumoniae",
    subtitle = "Odds Ratios and 95% Confidence Intervals (Generated with Fast Function)",
    x = "Odds Ratio (log scale)",
    y = "Predictor Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
# View the first few rows of the tidy results
# print(head(k_pneumoniae_results))

# Create the forest plot
ggplot(
  data = k_pneumoniae_results, 
  aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = is_significant)
) +
  
  # Add the points for the odds ratio estimates
  geom_point(size = 3) +
  
  # Add the horizontal error bars for the confidence intervals
  geom_errorbarh(height = 0.2, linewidth = 1) +
  
  # Add a vertical line at 1.0, which is the line of "no effect"
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  
  # Use a log scale for the x-axis, which is standard for odds ratios
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")
  ) +
  
  # Separate the plot into panels for each outcome cluster
  # This is crucial for interpreting a multinomial model
  facet_wrap(~ y.level, ncol = 2, scales = "free_y") +
  
  # Manually set the colors to make significance clear (e.g., Red for significant)
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     name = "Significant (p < 0.05)",
                     labels = c("Yes", "No")) +
  
  # Add labels and a clean theme
  labs(
    title = "Forest Plot of Multinomial Regression for K. pneumoniae",
    subtitle = "Odds Ratios and 95% Confidence Intervals",
    x = "Odds Ratio (log scale)",
    y = "Predictor Variable"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank(), # Clean up grid lines
    panel.border = element_rect(colour = "grey80", fill=NA), # Add border to facets
    strip.text = element_text(face = "bold") # Make facet titles bold
  )

a_baumannii_1 <- a_baumannii %>%
  filter(Cluster== "Pan-susceptible")

a_baumannii_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

a_baumannii_2 <- a_baumannii %>%
  filter(Cluster=="Non-CRAB MDRAB")

a_baumannii_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

a_baumannii_3 <- a_baumannii %>%
  filter(Cluster=="CRAB MDRAB")

a_baumannii_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

a_baumannii_4 <- a_baumannii %>%
  filter(Cluster=="Low-Level Resistance")

a_baumannii_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()



#E. faecium

e_faecium_1 <- e_faecium %>%
  filter(Cluster=="VS-MDR (Vancomycin-susceptible)")

e_faecium_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_faecium_2 <- e_faecium %>%
  filter(Cluster=="VR-MR (Vancomycin-resistant")

e_faecium_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_faecium_3 <- e_faecium %>%
  filter(Cluster=="Amp-S, Van-S, E-R (Ampicillin-Susceptible, Erythromycin-Resistant Cluster")

e_faecium_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_faecium_4 <- e_faecium %>%
  filter(Cluster=="Pan-susceptible")

e_faecium_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_faecium_5 <- e_faecium %>%
  filter(Cluster=="High-Level Penicillin-Resistant VSE")

e_faecium_5 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

#E. spp
e_spp_1 <- e_spp %>%
  filter(Cluster=="XDR-CRE")

e_spp_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_spp_2 <- e_spp %>%
  filter(Cluster == "ESBL, Carbapenem-S")

e_spp_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_spp_3 <- e_spp %>%
  filter(Cluster=="MDR, Carbapenem-S, (AmpC-like)")

e_spp_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_spp_4 <- e_spp %>%
  filter(Cluster=="Pan-Susceptible")

e_spp_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

e_spp_5 <- e_spp %>%
  filter(Cluster=="Amp-R, ESBL-Negative")

e_spp_5 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

#S. aureus

s_aureus_1 <- s_aureus %>%
  filter(Cluster == "MRSA(Gentamicin-R")

s_aureus_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

s_aureus_2 <- s_aureus %>%
  filter(Cluster == "MRSA(Fluoroq-R)")

s_aureus_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

s_aureus_3 <- s_aureus %>%
  filter(Cluster == "MSSA")

s_aureus_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

s_aureus_4 <- s_aureus %>%
  filter(Cluster == "MRSA(Gentamicin-S)")

s_aureus_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()



k_pneumoniae_1 <- k_pneumoniae %>%
  filter(Cluster == "CRE XDR")

k_pneumoniae_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

k_pneumoniae_2 <- k_pneumoniae %>%
  filter(Cluster == "ESBL-like MDR")

k_pneumoniae_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


k_pneumoniae_3 <- k_pneumoniae %>%
  filter(Cluster == "AmpC-like MDR")

k_pneumoniae_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


k_pneumoniae_4 <- k_pneumoniae %>%
  filter(Cluster == "Susceptible(Intrinsic Amp-R)")

k_pneumoniae_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


k_pneumoniae_5 <- k_pneumoniae %>%
  filter(Cluster == "non-CRE XDR")

k_pneumoniae_5 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


p_aeruginosa_1 <- p_aeruginosa %>%
  filter(Cluster == "Pan-S")

p_aeruginosa_1 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()

p_aeruginosa_2 <- p_aeruginosa %>%
  filter(Cluster == "Fluoroq-R")

p_aeruginosa_2 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


p_aeruginosa_3 <- p_aeruginosa %>%
  filter(Cluster == "MDR")

p_aeruginosa_3 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


p_aeruginosa_4 <- p_aeruginosa %>%
  filter(Cluster == "CRP")

p_aeruginosa_4 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


p_aeruginosa_5 <-  p_aeruginosa %>%
  filter(Cluster == "XDR-CRP")

p_aeruginosa_5 %>%
  dplyr::select(Super_Region, Gender, Age.Group, Speciality,
                Source, In...Out.Patient, Year) %>%
  tbl_summary()


#Plotting variables in relationship to cluster

a_baumannii_age <- ggplot(data = a_baumannii, mapping = aes(x = Age.Group, 
                                                            fill = as.factor(a_baumannii_model$predclass))) + geom_bar(position = "fill")

a_baumanni_region <- ggplot(data = a_baumannii, mapping = aes(x = Super_Region,
                                                              fill = as.factor(a_baumannii_model$predclass))) +
  geom_bar(position = "fill", orientation = .75) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

e_faecium_age <- ggplot(data = e_faecium, mapping = aes(x = Age.Group, 
                                                        fill = as.factor(e_faecium_model$predclass))) + geom_bar(position = "fill")

e_faecium_region <- ggplot(data = e_faecium, mapping = aes(x = Super_Region,
                                                           fill = as.factor(e_faecium_model$predclass))) +
  geom_bar(position = "fill", orientation = .75) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



e_spp_age <- ggplot(data = e_spp, mapping = aes(x = Age.Group, 
                                                fill = as.factor(e_spp_model$predclass))) + geom_bar(position = "fill")

e_spp_region <- ggplot(data = e_spp, mapping = aes(x = Super_Region,
                                                   fill = as.factor(e_spp_model$predclass))) +
  geom_bar(position = "fill", orientation = .75) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



s_aureus_age <- ggplot(data = s_aureus, mapping = aes(x = Age.Group, 
                                                fill = as.factor(Cluster))) + geom_bar(position = "fill")

s_aureus_region <- ggplot(data = s_aureus, mapping = aes(x = Super_Region,
                                                   fill = as.factor(Cluster))) +
  geom_bar(position = "fill", orientation = .75) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


k_pneumoniae_age <- ggplot(data = k_pneumoniae, mapping = aes(x = Age.Group, 
                                                      fill = as.factor(Cluster))) + geom_bar(position = "fill")

k_pneumoniae_region <- ggplot(data = k_pneumoniae, mapping = aes(x = Super_Region,
                                                         fill = as.factor(Cluster))) +
  geom_bar(position = "fill", orientation = .75) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Creating resistance probability object
a_baumannii_resistance_prob <- a_baumannii_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  mutate(Class = levels(a_baumannii$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")

#Removing medication with 0 resistance
e_faecium_model$probs$Tigecycline <- NULL

e_faecium_resistance_prob <- e_faecium_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  # Add class names
  mutate(Class = levels(e_faecium$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")

e_spp_resistance_prob <- e_spp_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  # Add class names
  mutate(Class = levels(e_spp$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")

#Remove two attributes without any resistance registered
s_aureus_model$probs$Daptomycin <- NULL
s_aureus_model$probs$Tigecycline <- NULL

s_aureus_resistance_prob <- s_aureus_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  # Add class names
  mutate(Class = levels(s_aureus$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")


k_pneumoniae_resistance_prob <- k_pneumoniae_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  # Add class names
  mutate(Class = levels(k_pneumoniae$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")


p_aeruginosa_resistance_prob <- p_aeruginosa_model$probs %>%
  # Keep only the probability for "Resistant"
  lapply(., function(x) x[, "Resistant"]) %>% 
  # Convert the list to a dataframe
  as.data.frame() %>%
  # Add class names
  mutate(Class = levels(p_aeruginosa$Cluster)) %>%
  # Reshape data for ggplot
  pivot_longer(cols = -Class, names_to = "Antibiotic", values_to = "Probability")

#################################



#Making class lookup objects
a_baumannii_antibiotic_class_lookup <- c(
  # --- Beta-lactams ---
  # Cephalosporins (often grouped by generation or just as 'Cephalosporin')
  "Cefepime" = "Cephalosporin",      # Typically 4th generation
  "Ceftazidime" = "Cephalosporin",   # Typically 3rd generation
  "Ceftriaxone" = "Cephalosporin",   # Typically 3rd generation
  
  # Carbapenems
  "Imipenem" = "Carbapenem",
  "Meropenem" = "Carbapenem",
  
  # Penicillin + Beta-Lactamase Inhibitor Combinations
  "Piperacillin.tazobactam" = "Penicillin + Inhibitor", # Often called Piperacillin/Tazobactam
  "Ampicillin.sulbactam" = "Penicillin + Inhibitor",   # Often called Ampicillin/Sulbactam
  
  # --- Fluoroquinolones ---
  "Levofloxacin" = "Fluoroquinolone",
  "Ciprofloxacin" = "Fluoroquinolone",
  
  # --- Tetracyclines ---
  "Minocycline" = "Tetracycline"
)

e_faecium_class_lookup <- c(
  "Ampicillin" = "Penicillin",
  "Erythromycin" = "Macrolide",
  "Levofloxacin" = "Fluoroquinolone",
  "Linezolid" = "Oxazolidinone",
  "Minocycline" = "Tetracycline",
  "Penicillin" = "Penicillin", # Generic Penicillin
  "Vancomycin" = "Glycopeptide",
  "Daptomycin" = "Lipopeptide",
  "Quinupristin.dalfopristin" = "Streptogramin", # Often called Quinupristin/Dalfopristin
  "Teicoplanin" = "Glycopeptide"
)

e_spp_class_lookup <- c(
  "Amikacin" = "Aminoglycoside",
  "Amoxycillin.clavulanate" = "Penicillin + Beta-lactamase inhibitor",
  "Ampicillin" = "Penicillin",
  "Cefepime" = "Cephalosporin (4th gen)",
  "Ceftazidime" = "Cephalosporin (3rd gen)",
  "Imipenem" = "Carbapenem",
  "Levofloxacin" = "Fluoroquinolone",
  "Meropenem" = "Carbapenem",
  "Piperacillin.tazobactam" = "Penicillin + Beta-lactamase inhibitor",
  "Tigecycline" = "Glycylcycline",
  "Ampicillin.sulbactam" = "Penicillin + Beta-lactamase inhibitor",
  "Aztreonam" = "Monobactam",
  "Cefixime" = "Cephalosporin (3rd gen)",
  "Ceftaroline" = "Cephalosporin (5th gen)",
  "Ceftazidime.avibactam" = "Cephalosporin + Beta-lactamase inhibitor",
  "Ciprofloxacin" = "Fluoroquinolone",
  "Colistin" = "Polymyxin",
  "Gentamicin" = "Aminoglycoside",
  "Trimethoprim.sulfa" = "Sulfonamide combination",
  "Ceftolozane.tazobactam" = "Cephalosporin + Beta-lactamase inhibitor",
  "Meropenem.vaborbactam" = "Carbapenem + Beta-lactamase inhibitor",
  "Cefpodoxime" = "Cephalosporin (3rd gen)",
  "Ceftibuten" = "Cephalosporin (3rd gen)"
)

s_aureus_class_lookup <- c(
  "Clindamycin" = "Lincosamide",
  "Erythromycin" = "Macrolide",
  "Levofloxacin" = "Fluoroquinolone",
  "Linezolid" = "Oxazolidinone",
  "Minocycline" = "Tetracycline",
  "Vancomycin" = "Glycopeptide",
  "Ceftaroline" = "Cephalosporin (5th gen)", 
  "Gentamicin" = "Aminoglycoside",
  "Moxifloxacin" = "Fluoroquinolone",
  "Oxacillin" = "Penicillin", 
  "Teicoplanin" = "Glycopeptide",
  "Trimethoprim.sulfa" = "Sulfonamide combination"
)

k_pneumoniae_class_lookup  <- tibble(
  Antibiotic = c(
    "Amikacin", "Amoxycillin.clavulanate", "Ampicillin", "Cefepime", 
    "Ceftazidime", "Ceftriaxone", "Imipenem", "Levofloxacin", 
    "Meropenem", "Minocycline", "Piperacillin.tazobactam", "Tigecycline", 
    "Ampicillin.sulbactam", "Aztreonam", "Cefixime", "Ceftaroline", 
    "Ceftazidime.avibactam", "Ciprofloxacin", "Colistin", "Doripenem", 
    "Ertapenem", "Gentamicin", "Trimethoprim.sulfa", "Ceftolozane.tazobactam", 
    "Meropenem.vaborbactam", "Cefpodoxime", "Ceftibuten"
  ),
  Antibiotic_Class = c(
    "Aminoglycoside", "Penicillin + Beta-lactamase inhibitor", "Penicillin", "Cephalosporin (4th gen)",
    "Cephalosporin (3rd gen)", "Cephalosporin (3rd gen)", "Carbapenem", "Fluoroquinolone",
    "Carbapenem", "Tetracycline", "Penicillin + Beta-lactamase inhibitor", "Glycylcycline",
    "Penicillin + Beta-lactamase inhibitor", "Monobactam", "Cephalosporin (3rd gen)", "Cephalosporin (5th gen)",
    "Cephalosporin + Beta-lactamase inhibitor", "Fluoroquinolone", "Polymyxin", "Carbapenem",
    "Carbapenem", "Aminoglycoside", "Sulfonamide combination", "Cephalosporin + Beta-lactamase inhibitor",
    "Carbapenem + Beta-lactamase inhibitor", "Cephalosporin (3rd gen)", "Cephalosporin (3rd gen)"
  )
)

p_aeruginosa_class_lookup <- tibble(
  Antibiotic = c(
    "Amikacin", "Cefepime", "Ceftazidime", "Imipenem", "Levofloxacin",
    "Meropenem", "Piperacillin.tazobactam", "Aztreonam", "Ceftazidime.avibactam",
    "Ciprofloxacin", "Colistin", "Doripenem", "Ceftolozane.tazobactam"
  ),
  Antibiotic_Class = c(
    "Aminoglycoside", "Cephalosporin (4th gen)", "Cephalosporin (3rd gen)", "Carbapenem", "Fluoroquinolone",
    "Carbapenem", "Penicillin + Beta-lactamase inhibitor", "Monobactam", "Cephalosporin + Beta-lactamase inhibitor",
    "Fluoroquinolone", "Polymyxin", "Carbapenem", "Cephalosporin + Beta-lactamase inhibitor"
  ))


a_baumannii_antibiotic_class_lookup <- a_baumannii_antibiotic_class_lookup %>%
  tibble(
    Antibiotic = names(a_baumannii_antibiotic_class_lookup),
    Antibiotic_Class = unname(a_baumannii_antibiotic_class_lookup)
  )

e_faecium_class_lookup <- e_faecium_class_lookup %>%
  tibble(
    Antibiotic = names(e_faecium_class_lookup),
    Antibiotic_Class = unname(e_faecium_class_lookup)
  )

e_spp_class_lookup <- e_spp_class_lookup %>%
  tibble(
    Antibiotic = names(e_spp_class_lookup),
    Antibiotic_Class=unname(e_spp_class_lookup)
  )

s_aureus_class_lookup <- s_aureus_class_lookup %>%
  tibble(
    Antibiotic = names(s_aureus_class_lookup),
    Antibiotic_Class=unname(s_aureus_class_lookup)
  )

a_baumannii_resistance_prob <- a_baumannii_resistance_prob %>%
  left_join(a_baumannii_antibiotic_class_lookup, by = "Antibiotic")

e_faecium_resistance_prob <- e_faecium_resistance_prob %>%
  left_join(e_faecium_class_lookup, by="Antibiotic")

e_spp_resistance_prob <- e_spp_resistance_prob %>%
  left_join(e_spp_class_lookup, by="Antibiotic")

s_aureus_resistance_prob <- s_aureus_resistance_prob %>%
  left_join(s_aureus_class_lookup, by="Antibiotic")

k_pneumoniae_resistance_prob <- k_pneumoniae_resistance_prob %>%
  left_join(k_pneumoniae_class_lookup, by="Antibiotic")

p_aeruginosa_resistance_prob <- p_aeruginosa_resistance_prob %>%
  left_join(p_aeruginosa_class_lookup, by="Antibiotic")

ggplot(a_baumannii_resistance_prob, aes(x = Antibiotic_Class, y = Probability, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() + #
  labs(
    title = "Average Resistance Probability by Antibiotic Class",
    subtitle = "Comparing resistance profiles for each of the 4 latent classes",
    x = "Antibiotic Class",
    y = "Average Probability of Resistance",
    fill = "Latent Class"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"))

ggplot(e_faecium_resistance_prob, aes(x = Antibiotic_Class, y = Probability, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() + #
  labs(
    title = "Average Resistance Probability by Antibiotic Class",
    subtitle = "Comparing resistance profiles for each of the 5 latent classes",
    x = "Antibiotic Class",
    y = "Average Probability of Resistance",
    fill = "Latent Class"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"))

ggplot(e_spp_resistance_prob, aes(x = Antibiotic_Class, y = Probability, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() + #
  labs(
    title = "Average Resistance Probability by Antibiotic Class",
    subtitle = "Comparing resistance profiles for each of the 5 latent classes",
    x = "Antibiotic Class",
    y = "Average Probability of Resistance",
    fill = "Latent Class"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"))

ggplot(s_aureus_resistance_prob, aes(x = Antibiotic_Class, y = Probability, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() + #
  labs(
    title = "Average Resistance Probability by Antibiotic Class",
    subtitle = "Comparing resistance profiles for each of the 5 latent classes",
    x = "Antibiotic Class",
    y = "Average Probability of Resistance",
    fill = "Latent Class"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"))
# 
# 
# ###Checking ESKAPE Pathogens and finding new ways of grouping data####
# eskape_pathogens <- data %>%
#   filter(Species == c("Enterococcus faecium","Staphylococcus aureus",
#                       "Klebsiella pneumoniae","Acinetobacter baumannii",
#                       "Pseudomonas aeruginosa","Enterobacter spp"))
# 
# 
# #Checking stats per country and bug
# samples_per_country <- eskape_pathogens %>%
#   group_by(Species, Country) %>%
#   summarise(sample_nr = n())
# 
# samples_per_country_over_100 <- samples_per_country %>%
#   filter(sample_nr > 100)

country_to_region_lookup <- tibble::tribble(
  ~Country, ~Super_Region,
  # East Asia and Pacific
  "Australia", "East Asia and Pacific",
  "China", "East Asia and Pacific",
  "Hong Kong", "East Asia and Pacific",
  "Indonesia", "East Asia and Pacific",
  "Japan", "East Asia and Pacific",
  "Malaysia", "East Asia and Pacific",
  "New Zealand", "East Asia and Pacific",
  "Philippines", "East Asia and Pacific",
  "Singapore", "East Asia and Pacific",
  "Korea, South", "East Asia and Pacific", # Republic of Korea in UNICEF's list
  "Taiwan", "East Asia and Pacific", # Not explicitly listed by UNICEF, but geographically fits here.
  "Thailand", "East Asia and Pacific",
  "Vietnam", "East Asia and Pacific",
  
  # Europe and Central Asia - Eastern Europe and Central Asia (EECA)
  "Bulgaria", "Eastern Europe and Central Asia",
  "Croatia", "Eastern Europe and Central Asia",
  "Czech Republic", "Eastern Europe and Central Asia", # Czechia in UNICEF's list
  "Hungary", "Eastern Europe and Central Asia",
  "Poland", "Eastern Europe and Central Asia",
  "Romania", "Eastern Europe and Central Asia",
  "Russia", "Eastern Europe and Central Asia", # Russian Federation in UNICEF's list
  "Serbia", "Eastern Europe and Central Asia",
  "Slovak Republic", "Eastern Europe and Central Asia", # Slovakia in UNICEF's list
  "Slovenia", "Eastern Europe and Central Asia",
  "Turkey", "Eastern Europe and Central Asia", # Türkiye in UNICEF's list. Also appears in MENA context. Prioritizing direct country listing under EECA.
  "Ukraine", "Eastern Europe and Central Asia",
  
  # Europe and Central Asia - Western Europe (WE)
  "Austria", "Western Europe",
  "Belgium", "Western Europe",
  "Denmark", "Western Europe",
  "Estonia", "Western Europe",
  "Finland", "Western Europe",
  "France", "Western Europe",
  "Germany", "Western Europe",
  "Greece", "Western Europe",
  "Ireland", "Western Europe",
  "Italy", "Western Europe",
  "Latvia", "Western Europe",
  "Lithuania", "Western Europe",
  "Netherlands", "Western Europe", # Kingdom of the Netherlands in UNICEF's list
  "Norway", "Western Europe",
  "Portugal", "Western Europe",
  "Spain", "Western Europe",
  "Sweden", "Western Europe",
  "Switzerland", "Western Europe",
  "United Kingdom", "Western Europe",
  
  # Latin America and Caribbean
  "Argentina", "Latin America and Caribbean",
  "Brazil", "Latin America and Caribbean",
  "Chile", "Latin America and Caribbean",
  "Colombia", "Latin America and Caribbean",
  "Costa Rica", "Latin America and Caribbean",
  "Dominican Republic", "Latin America and Caribbean",
  "El Salvador", "Latin America and Caribbean",
  "Guatemala", "Latin America and Caribbean",
  "Honduras", "Latin America and Caribbean",
  "Jamaica", "Latin America and Caribbean",
  "Mexico", "Latin America and Caribbean",
  "Nicaragua", "Latin America and Caribbean",
  "Panama", "Latin America and Caribbean",
  "Puerto Rico", "Latin America and Caribbean", # Not explicitly listed by UNICEF, but geographically fits here.
  "Venezuela", "Latin America and Caribbean", # Bolivarian Republic of Venezuela in UNICEF's list
  
  # Middle East and North Africa
  "Egypt", "Middle East and North Africa",
  "Israel", "Middle East and North Africa",
  "Jordan", "Middle East and North Africa",
  "Kuwait", "Middle East and North Africa",
  "Lebanon", "Middle East and North Africa",
  "Morocco", "Middle East and North Africa",
  "Oman", "Middle East and North Africa",
  "Qatar", "Middle East and North Africa",
  "Saudi Arabia", "Middle East and North Africa",
  "Tunisia", "Middle East and North Africa",
  
  # North America
  "Canada", "North America",
  "United States", "North America",
  
  # South Asia
  "India", "South Asia",
  "Pakistan", "South Asia",
  
  # Sub-Saharan Africa - Eastern and Southern Africa
  "Kenya", "Eastern and Southern Africa",
  "Malawi", "Eastern and Southern Africa",
  "Mauritius", "Eastern and Southern Africa",
  "Namibia", "Eastern and Southern Africa",
  "South Africa", "Eastern and Southern Africa",
  "Uganda", "Eastern and Southern Africa",
  
  # Sub-Saharan Africa - West and Central Africa
  "Cameroon", "West and Central Africa",
  "Ghana", "West and Central Africa",
  "Ivory Coast", "West and Central Africa", 
  "Nigeria", "West and Central Africa",
  "Senegal", "West and Central Africa" 
)
# 
# eskape_pathogens <- eskape_pathogens %>%
#   left_join(country_to_region_lookup, by = "Country")
# 
# samples_per_region <- eskape_pathogens %>%
#   group_by(Super_Region, Species) %>%
#   summarise(sample_nr = n())
# 
# eskape_pathogens_country_NA <- eskape_pathogens %>%
#   filter(is.na(Super_Region))
# 
# ###Preparing data for LCA
# 
# #Separating by region and bug
# list_of_dfs_by_bug_region <- eskape_pathogens %>%
#   group_by(Species, Super_Region) %>% 
#   group_nest()
# 
# #saveRDS(Vivli + data validation, "amr_processed_data.rds")
# 
# ####TESTING OUT PREDICTIVE MODELLING OF CLUSTER MEMBERSHIP
# a_baumannii_prediction_clean <- a_baumannii %>%
#   count(Super_Region, Year, Cluster) %>% 
#   group_by(Super_Region, Year) %>%                             
#   mutate(proportion = n / sum(n)) %>%                          
#   ungroup() 
# 
# # --- Step-by-Step Guide to Building a Cluster Prediction Model ---
# # This script implements a multi-model XGBoost approach to forecast future
# # cluster proportions, including bootstrapped confidence intervals.
# 
# library(xgboost)
# 
# # We'll focus on a single region to build our first model.
# model_data <- a_baumannii_prediction_clean %>%
#   filter(Super_Region == "Western Europe")
# 
# # --- 1. Feature Engineering ---
# # Pivot the data so each cluster's proportion gets its own column.
# feature_data <- model_data %>%
#   dplyr::select(Year, Cluster, proportion) %>%
#   pivot_wider(names_from = Cluster,
#               values_from = proportion,
#               names_prefix = "Cluster_",
#               values_fill = 0) %>%
#   arrange(Year)
# 
# # Create "lagged" features (predictors)
# feature_data_lagged <- feature_data %>%
#   mutate(across(starts_with("Cluster_"), ~lag(.x, 1), .names = "{.col}_lag1"))
# 
# # --- 2. Multi-Model Training with Bootstrapping ---
# # We will train an ensemble of models for EACH cluster.
# 
# predictor_variables <- colnames(feature_data_lagged)[grepl("_lag1$", colnames(feature_data_lagged))]
# target_variables <- colnames(feature_data)[grepl("^Cluster_", colnames(feature_data))]
# 
# # Create a complete dataset for modeling (remove NAs from lagging)
# modeling_df <- feature_data_lagged %>% na.omit()
# 
# # This list will hold the bootstrap model ensembles for each target cluster
# all_bootstrap_models <- list()
# 
# # Loop through each target cluster to create its own set of models
# for (target_var in target_variables) {
#   cat(paste("\n--- Training models for:", target_var, "---\n"))
#   
#   # Train an initial model to get residuals
#   initial_model <- xgboost(
#     data = as.matrix(modeling_df[, predictor_variables]),
#     label = modeling_df[[target_var]],
#     nrounds = 100, objective = "reg:squarederror", verbose = 0
#   )
#   residuals <- modeling_df[[target_var]] - predict(initial_model, as.matrix(modeling_df[, predictor_variables]))
#   
#   # Bootstrap loop to create an ensemble for this specific cluster
#   n_bootstrap <- 100
#   bootstrap_models_for_target <- list()
#   for (i in 1:n_bootstrap) {
#     bootstrap_target <- predict(initial_model, as.matrix(modeling_df[, predictor_variables])) + sample(residuals, size = nrow(modeling_df), replace = TRUE)
#     
#     bootstrap_model <- xgboost(
#       data = as.matrix(modeling_df[, predictor_variables]),
#       label = bootstrap_target,
#       nrounds = 100, objective = "reg:squarederror", verbose = 0
#     )
#     bootstrap_models_for_target[[i]] <- bootstrap_model
#   }
#   all_bootstrap_models[[target_var]] <- bootstrap_models_for_target
# }
# 
# 
# # --- 3. Iterative Forecasting with a Multi-Model Approach ---
# n_forecast_years <- 5
# last_known_year <- max(feature_data$Year)
# last_known_features <- feature_data %>%
#   filter(Year == last_known_year) %>%
#   dplyr::select(all_of(target_variables)) %>%
#   as.numeric()
# 
# # This array will store all predictions for all clusters from all bootstrap models
# # Dims: bootstrap_run, forecast_year, cluster_id
# future_predictions_array <- array(NA, dim = c(n_bootstrap, n_forecast_years, length(target_variables)))
# 
# # Loop through each bootstrap simulation
# for (i in 1:n_bootstrap) {
#   # Start with the last known real data
#   current_features <- last_known_features
#   
#   # Iteratively predict future years
#   for (j in 1:n_forecast_years) {
#     # This vector will hold predictions for all clusters for a single year
#     yearly_predictions <- numeric(length(target_variables))
#     
#     # Prepare input matrix (it's the same for all cluster models in a given year)
#     input_matrix <- matrix(current_features, nrow = 1)
#     colnames(input_matrix) <- predictor_variables
#     
#     # Predict each cluster's proportion for the year
#     for (k in 1:length(target_variables)) {
#       target_var <- target_variables[k]
#       current_model <- all_bootstrap_models[[target_var]][[i]]
#       yearly_predictions[k] <- predict(current_model, input_matrix)
#     }
#     
#     # --- Normalize predictions to sum to 1 ---
#     # Clip at 0 to avoid negative proportions
#     yearly_predictions[yearly_predictions < 0] <- 0
#     normalized_predictions <- yearly_predictions / sum(yearly_predictions)
#     
#     # Store the normalized predictions in our results array
#     future_predictions_array[i, j, ] <- normalized_predictions
#     
#     # The new features for the next iteration are the predictions we just made
#     current_features <- normalized_predictions
#   }
# }
# 
# # --- 4. Summarize Forecasts and Plot a Specific Cluster ---
# # Let's choose which cluster we want to visualize
# cluster_to_plot <- "Cluster_3"
# cluster_index <- which(target_variables == cluster_to_plot)
# 
# # Extract the predictions for our chosen cluster
# predictions_for_one_cluster <- future_predictions_array[, , cluster_index]
# 
# # Calculate the mean, lower, and upper bounds
# forecast_summary <- data.frame(
#   Year = (last_known_year + 1):(last_known_year + n_forecast_years),
#   Point_Forecast = apply(predictions_for_one_cluster, 2, mean),
#   Lower_CI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.025),
#   Upper_CI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.975)
# )
# 
# print("Forecast Summary with Confidence Intervals:")
# print(forecast_summary)
# 
# # --- Prepare data for a seamless plot ---
# plot_data_historical <- feature_data %>%
#   dplyr::select(Year, Actual_Proportion = !!sym(cluster_to_plot))
# 
# last_actual_point <- plot_data_historical %>%
#   filter(Year == last_known_year)
# 
# forecast_line_data <- bind_rows(
#   data.frame(Year = last_actual_point$Year, Point_Forecast = last_actual_point$Actual_Proportion),
#   forecast_summary %>% 
#     dplyr::select(Year, Point_Forecast)
# )
# 
# # Plot the results
# ggplot(plot_data_historical, aes(x = Year, y = Actual_Proportion)) +
#   geom_line(aes(color = "Actual"), linewidth = 1.2) +
#   geom_ribbon(data = forecast_summary, aes(x = Year, ymin = Lower_CI, ymax = Upper_CI),
#               fill = "skyblue", alpha = 0.5, inherit.aes = FALSE) +
#   geom_line(data = forecast_line_data, aes(x = Year, y = Point_Forecast, color = "Forecast"),
#             linewidth = 1.2, linetype = "dashed") +
#   labs(
#     title = "XGBoost Forecast with 95% Confidence Interval",
#     subtitle = paste("Forecasting", cluster_to_plot, "Proportion in Western Europe"),
#     y = "Proportion of Isolates",
#     color = "Legend"
#   ) +
#   scale_color_manual(values = c("Actual" = "black", "Forecast" = "red")) +
#   theme_minimal(base_size = 14)
# 
# 
# 
# 
# e_faecium$`e_faecium_model$predclass` <- as.factor(e_faecium$`e_faecium_model$predclass`)
# 
# e_spp$`e_spp_model$predclass` <- as.factor(e_spp$`e_spp_model$predclass`)
# 
# a_baumannii$Country <- factor(a_baumannii$Country, levels=unique(a_baumannii$Country), labels = unique(a_baumannii$Country))
# 

#Making HEATMAPS for resistance to antibiotic classes

#Plotting resistance grouper by antibiotic class and cluster

generate_resistance_heatmap <- function(resistance_prob){
  aggregated_probs <- resistance_prob %>%
    group_by(Class, Antibiotic_Class) %>%
    summarise(
      mean_probs=mean(Probability, na.rm=TRUE)
    ) %>%
    ungroup()
  
  class_probs_table <- aggregated_probs %>%
    pivot_wider(
      names_from = Antibiotic_Class,
      values_from = mean_probs
    )
  
  
  # Create the table
  heatmap_table <- class_probs_table %>%
    gt(rowname_col = "Cluster") %>% 
    
    # Add a main title and subtitle
    tab_header(
      title = md("**Antibiotic Resistance Profiles by Cluster**"),
      subtitle = "Mean probability of resistance for each antibiotic class"
    ) %>%
    
    # Add a spanner header over the antibiotic columns
    tab_spanner(
      label = md("**Antibiotic Class**"),
      columns = everything() # Apply to all columns except the row labels
    ) %>%
    
    # Format the numbers as percentages with one decimal place
    fmt_percent(
      columns = everything(),
      decimals = 1
    ) %>%
    data_color(
      columns = where(is.numeric),
      colors = scales::col_numeric(
        palette = c("#63BE7B", "#FFEB84", "#F8696B"), # Green, Yellow, Red
        domain = c(0, 1) 
      )
    )
return(heatmap_table)
  }

a_baumannii_heatmap <- generate_resistance_heatmap(a_baumannii_resistance_prob)
e_faecium_heatmap <-  generate_resistance_heatmap(e_faecium_resistance_prob)
e_spp_heatmap <- generate_resistance_heatmap(e_spp_resistance_prob)
k_pneumoniae_heatmap <- generate_resistance_heatmap(k_pneumoniae_resistance_prob)
p_aeruginosa_heatmap <- generate_resistance_heatmap(p_aeruginosa_resistance_prob)
s_aureus_heatmap <- generate_resistance_heatmap(s_aureus_resistance_prob)

# 
# #Create dendograms to see cluster distance
# library(ggdendro)
# library(patchwork)
# 
# # Step 1: Perform Hierarchical Clustering to get the order
# # First, we need the data in a 'wide' format for clustering
# cluster_profiles_wide <- a_baumannii_resistance_prob %>%
#   pivot_wider(names_from = Antibiotic, values_from = Probability) %>%
#   # Use Class as row names for the matrix
#   tibble::column_to_rownames("Class")
# 
# # Now, perform the clustering
# dist_matrix <- dist(cluster_profiles_wide, method = "euclidean")
# hc <- hclust(dist_matrix, method = "ward.D2")
# 
# # This is the crucial part: get the order of clusters from the dendrogram
# cluster_order <- hc$labels[hc$order]
# 
# 
# # Step 2: Reorder your original (long) dataframe to match the dendrogram
# k_pneumoniae_prob_reordered <- k_pneumoniae_resistance_prob %>%
#   mutate(Class = factor(Class, levels = cluster_order))
# 
# 
# # Step 3: Create the two separate plots
# 
# # Plot 1: The Dendrogram (rotated)
# dendro_plot <- ggdendrogram(hc, rotate = TRUE) + 
#   theme(axis.text.y = element_blank()) # Hides the now-vertical labels

# Load all libraries
library(Hmisc)       # for read.xport
library(foreign)     # legacy SAS import
library(brms)
library(cmdstanr)
library(dplyr)       # data wrangling
library(tidyr)       # for drop_na, fill, etc.
# library(mirt)      # for IRT modeling
library(psych)       # for PCA / FA
library(Matrix)      # for SVD
library(caret)       # for preprocessing
library(corrr)       # for correlation filtering
library(readr)       # for read_csv
library(mgcv)
library(glmnet)
library(latex2exp)
library(pbapply)
library(xgboost)
library(mclust)
library(SHAPforxgboost)

# Data wrangling and response modeling with IRT BRMS
# Load and Merge NHANES Files
files <- list(
  BIOPRO = "P_BIOPRO.XPT",
  BMX = "P_BMX.XPT",
  BPQ = "P_BPQ.XPT",
  BPXO = "P_BPXO.XPT",
  CBC = "P_CBC.XPT",
  CDQ = "P_CDQ.XPT",
  COT = "P_COT.XPT",
  DEMO = "P_DEMO.XPT",
  DPQ = "P_DPQ.XPT",
  FASTQX = "P_FASTQX.XPT",
  GLU = "P_GLU.XPT",
  FETIB = "P_FETIB.XPT",
  GHB = "P_GHB.XPT",
  HDL = "P_HDL.XPT",
  HEPA = "P_HEPA.XPT",
  HEPB_S = "P_HEPB_S.XPT",
  HSCRP = "P_HSCRP.XPT",
  IHGEM = "P_IHGEM.XPT",
  INQ = "P_INQ.XPT",
  INS = "P_INS.XPT",
  MCQ = "P_MCQ.XPT",
  PAQ = "P_PAQ.XPT",
  PBCD = "P_PBCD.XPT",
  TST = "P_TST.XPT",
  TRIGLY = "P_TRIGLY.XPT",
  UVOC = "P_UVOC.XPT",
  UVOC2 = "P_UVOC2.XPT"
)
dfs <- lapply(files, function(f) as.data.frame(read.xport(f)))
df <- Reduce(function(x, y) full_join(x, y, by = "SEQN"), dfs)

# Clean and filter
df <- df[rowMeans(is.na(df)) <= 0.5, ]
df <- df[, colMeans(is.na(df)) <= 0.3]
df <- df %>% mutate(across(where(is.numeric), ~ ifelse(abs(.) < 1e-50, 0, .)))

skip_codes <- c(7, 9, 77, 99, 777, 999)

suspect_vars <- df %>%
  select(where(is.numeric)) %>%
  select(where(~ all(na.omit(.) == floor(na.omit(.))))) %>%  # only true integers
  select(where(~ n_distinct(na.omit(.)) <= 20)) %>%           # low-cardinality
  summarise(across(everything(), ~ any(. %in% skip_codes))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "has_skip_code") %>%
  filter(has_skip_code) %>%
  pull(variable)

# Recode skip codes to NA
df <- df %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 7))) %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 9))) %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 77))) %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 99))) %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 777))) %>%
  mutate(across(all_of(suspect_vars), ~ na_if(., 999)))

# Bayesian Graded Response Model
## Step 1: Extract and clean PHQ-9 item responses
dep_questions <- c('DPQ010','DPQ020','DPQ030','DPQ040','DPQ050','DPQ060','DPQ070','DPQ080','DPQ090')

phq_data <- df[dep_questions]
phq_data <- lapply(phq_data, function(col) {
  col <- as.numeric(col)
  col[col %in% c(7, 9)] <- NA  # Handle PHQ skip codes
  return(col)
}) %>% as.data.frame()

# Only keep rows with any non-missing responses
phq_data_valid <- phq_data[rowSums(!is.na(phq_data)) > 0, ]
phq_data_valid$id <- as.integer(rownames(phq_data_valid))  # preserve linkage to df

# Step 2: Convert to long format for brms
phq_long <- phq_data_valid %>%
  pivot_longer(cols = -id, names_to = "item", values_to = "response") %>%
  filter(!is.na(response)) %>%
  mutate(
    response = ordered(response)  # required for ordinal model
  )

# Fit the 1D IRT model

set.seed(42)

# Fit the Bayesian Graded Response Model

brm_mod <- brm(
  response ~ 1 + item + (1 | id),
  data = phq_long,
  family = cumulative("logit"),
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options = list("O1")),
  refresh = 10
)

saveRDS(brm_mod, file="part_hier_13jun.rds")

# load brm_model
brm_model <- readRDS("part_hier_13jun.rds")
# started at 2320

# Extract posterior samples of person-level latent traits
theta_draws <- ranef(brm_model)$id[, , "Intercept"]  # person-by-draw matrix

# Summarize
theta_summary <- data.frame(
  id = as.integer(rownames(theta_draws)),
  theta_mean = rowMeans(theta_draws),
  theta_sd = apply(theta_draws, 1, sd),
  theta_q025 = apply(theta_draws, 1, quantile, 0.025),
  theta_q975 = apply(theta_draws, 1, quantile, 0.975)
)

df$theta_dep <- NA
df$theta_dep[theta_summary$id] <- theta_summary$theta_mean
df <- df %>% filter(!is.na(theta_dep))

ggplot(df, aes(theta_dep)) + geom_density()

ggplot(df, aes(theta_dep)) + geom_histogram(bins=100)

write.csv(theta_summary, file="theta_summary.csv")
write.csv(df, file="brm_response_full.csv")

# ------------ START HERE, ABOVE IS RESPONSE MODELING (ALREADY DONE)------------
df <- read_csv('brm_response_full.csv') %>% select(-1)

df %>% colnames
dep_questions <- c('DPQ010','DPQ020','DPQ030','DPQ040','DPQ050','DPQ060','DPQ070','DPQ080','DPQ090')

y <- df$theta_dep
X <- df %>%
  select(-all_of(c("SEQN", "MCQ520", "MCQ366A", "MCQ366B", "MCQ366C", "MCQ366D", "MCQ080", "RIAGENDR", "RIDRETH3", "DMDMARTZ", "DMDBORN4", dep_questions, "theta_dep"))) %>%
  select(where(is.numeric))

X_imputed <- X %>% mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
nzv <- nearZeroVar(X_imputed, saveMetrics = TRUE)
X_reduced <- X_imputed[, !nzv$nzv]

cor_matrix <- cor(X_reduced, use = "pairwise.complete.obs")
high_corr <- findCorrelation(cor_matrix, cutoff = 0.98)
X_deduped <- X_reduced[, -high_corr]

X_centered <- scale(X_deduped, center = TRUE, scale = TRUE)
svd_out <- svd(X_centered)
rank <- sum(svd_out$d > 1e-10)

X_fullrank <- svd_out$u[, 1:rank] %*% diag(svd_out$d[1:rank]) %*% t(svd_out$v[, 1:rank]) %>% as.matrix
X_fullrank <- data.frame(X_fullrank + matrix(rnorm(prod(dim(X_fullrank)), 0, 1e-10), nrow = nrow(X_fullrank)))
colnames(X_fullrank) <- colnames(X_centered)

# ------------------------------------------------

# GMM clusters to determine slices and breaks
gmm <- Mclust(y)
summary(gmm)

# H <- 10

breaks <- unique(quantile(y, probs = seq(0, 1, length.out = H + 1), na.rm = TRUE))
# y_sliced <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE) using quantile
y_sliced <- gmm$classification
H <- length(unique(y_sliced))

table(y_sliced)

p_h <- as.numeric(table(y_sliced)) / length(y_sliced)
x_bar <- colMeans(X_fullrank)

M <- matrix(0, ncol = ncol(X_fullrank), nrow = ncol(X_fullrank))
for (h in 1:H) {
  idx <- which(y_sliced == h)
  if (length(idx) > 1) {
    xh <- colMeans(X_fullrank[idx, , drop = FALSE])
    diff <- xh - x_bar
    M <- M + p_h[h] * (diff %*% t(diff))
  }
}

Sigma_x <- cov(X_fullrank)
Sigma_x_reg <- Sigma_x + diag(1e-8, ncol(Sigma_x))
eig <- eigen((solve(Sigma_x_reg) %*% M + t(solve(Sigma_x_reg) %*% M)) / 2)

# Scree plot
plot(eig$values[1:20], type = "b", pch = 19,
     xlab = "SIR Component", ylab = "Eigenvalue",
     main = "Scree Plot of SIR Directions")
abline(h = 1e-3, col = "red", lty = 2)

d <- 10
sir_directions <- eig$vectors[, 1:d]
X_fullrank <- X_fullrank %>% as.matrix
X_sir <- X_fullrank %*% sir_directions

sir_loadings <- data.frame(Variable = colnames(X_centered),
                           SIR1 = sir_directions[, 1]) %>%
                          arrange(desc(abs(SIR1)))

# Compute loadings for both directions
sir_loadings <- data.frame(
  Variable = colnames(X_centered),
  SIR1 = sir_directions[, 1],
  SIR2 = sir_directions[, 2],
  SIR3 = sir_directions[, 3]
)

# Top 10 absolute loadings for SIR1
top_sir1 <- sir_loadings %>% select(Variable, SIR1) %>%
  arrange(desc(abs(SIR1))) %>%
  dplyr::slice(1:10)

# Top 10 absolute loadings for SIR2
top_sir2 <- sir_loadings %>% select(Variable, SIR2) %>%
  arrange(desc(abs(SIR2))) %>%
  dplyr::slice(1:10)

# Top 10 absolute loadings for SIR2
top_sir3 <- sir_loadings %>%
  arrange(desc(abs(SIR3))) %>%
  dplyr::slice(1:10)

# Optional: print for checking
top_sir1 %>% select(Variable, SIR1)
top_sir2 %>% select(Variable, SIR2)
top_sir3 %>% select(Variable, SIR3)

# (1) Make sure labels are the correct length
labels <- c("Very Low", "Low", "Mild", "Moderate", "Moderate-High", "High", "Very High", "Extreme")[1:H]

# (2) Assign as factor, making sure to include NA as its own level if you want
phq_cat <- factor(y_sliced, levels = 1:H, labels = labels)

# (3) Create sir_df in a way that preserves row order and length
sir_df <- as.data.frame(X_sir)
colnames(sir_df) <- c("SIR1", "SIR2", "SIR3")
sir_df$Severity <- phq_cat

# Plot using ggplot2
sir_df %>%
  ggplot(aes(x = SIR1, y = SIR2, color = Severity)) +
  geom_point(alpha = 0.4, size = 3) +
  scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  labs(
    title = "SIR Projection by PHQ-9 Depression Severity",
    x = "SIR1",
    y = "SIR2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + xlim(min(sir_df$SIR1), max(sir_df$SIR1)) +
  ylim(min(sir_df$SIR2), max(sir_df$SIR2))

colnames(sir_df) <- make.names(colnames(sir_df), unique=TRUE)
sir_df_plot <- sir_df %>%
  filter(!is.na(as.character(Severity)))

# Now plot safely!
sir_scatter <- ggplot(sir_df_plot, aes(x = SIR1, y = SIR2, color = Severity)) +
  geom_point(alpha = 0.4, size = 3) +
  scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  labs(
    title = "SIR Projection by PHQ-9 Depression Severity",
    x = "SIR1",
    y = "SIR2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("sir_scatter.png", plot=sir_scatter, width = 8, height = 6, dpi=300)

summary(aov(SIR1 ~ Severity, data = sir_df))
summary(aov(SIR2 ~ Severity, data = sir_df))
summary(aov(SIR3 ~ Severity, data = sir_df))

# Pivot long, THEN filter out NA
sir_long <- tidyr::pivot_longer(
  sir_df, 
  cols = starts_with("SIR"), 
  names_to = "SIR_Component", 
  values_to = "SIR_Value"
)

# Filter for SIR1 and SIR2 only, then filter out NA Severity
sir_long_plot <- sir_long %>%
  filter(SIR_Component %in% c("SIR1", "SIR2")) %>%
  filter(!is.na(as.character(Severity)))

sir_long_plot$Severity_Num <- as.integer(droplevels(sir_long_plot$Severity)) # recode only non-NA

sir_boxplot <- ggplot(sir_long_plot, aes(x = factor(Severity_Num), y = SIR_Value, fill = factor(Severity_Num))) +
  geom_boxplot() +
  facet_grid(. ~ SIR_Component) +
  labs(
    title = "SIR1 and SIR2 vs Depression Severity",
    x = "Severity Cluster",
    y = "SIR Value"
  ) +
  scale_fill_viridis_d(option = "plasma", direction = 1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

sir_boxplot
ggsave("sir_boxplot.png", plot=sir_boxplot, width = 8, height = 6, dpi=300)

# Define complete subset
idx <- complete.cases(y, X_sir[,1])

# Re-plot
plot(y[idx], X_sir[idx,1], pch = 19, col = rgb(0,0,0,0.1),
     xlab = "IRT Theta", ylab = "SIR1", main = "SIR1 vs Latent Depression")

# Add smoothed line
lines(lowess(y[idx], X_sir[idx,1]), col = "red", lwd = 2)

# Re-plot
plot(y[idx], X_sir[idx,2], pch = 19, col = rgb(0,0,0,0.1), 
     xlab = "IRT Theta", ylab = "SIR2", main = "SIR2 vs Latent Depression")

# Add smoothed line
lines(lowess(y[idx], X_sir[idx,2]), col = "red", lwd = 2)

# Re-plot
plot(y[idx], X_sir[idx,3], pch = 19, col = rgb(0,0,0,0.1), 
     xlab = "IRT Theta", ylab = "SIR2", main = "SIR3 vs Latent Depression")

# Add smoothed line
lines(lowess(y[idx], X_sir[idx,3]), col = "red", lwd = 2)

# Use complete cases
idx <- complete.cases(y, X_sir[,1:3])

# Fit GAM with smooth terms
gam_sir <- gam(y[idx] ~ s(X_sir[idx,3]) + s(X_sir[idx,2]) + s(X_sir[idx,3]) + s(X_sir[idx,4]), method = "REML")

# Summary
summary(gam_sir)

plot(gam_model, pages = 1, shade = TRUE)

y_gam <- predict(gam_model)

# Plot predictions vs actual
plot(y[idx], y_gam, pch = 20, col = rgb(0, 0, 0, 0.3),
     xlab = "True IRT Theta", ylab = "GAM Prediction",
     main = "GAM Fit using SIR1 + SIR2")

# final SIR workflow output --------------------
eig_values <- eig$values[1:3]
sir_loadings$SCombined <- sqrt(
  eig_values[1] * sir_loadings$SIR1^2 +
    eig_values[2] * sir_loadings$SIR2^2 +
    eig_values[3] * sir_loadings$SIR3^2
)
# Confirm the new column exists:
head(sir_loadings)
# Now arrange:
top_combined <- sir_loadings %>% arrange(desc(SCombined)) %>% dplyr::slice(1:10)

top_vars <- top_combined$Variable  # character vector of your top 10 variable names
gam_df <- df %>%
  select(theta_dep, all_of(top_vars)) %>%
  drop_na()
gam_formula <- as.formula(
  paste("theta_dep ~", paste(sprintf("s(%s)", top_vars), collapse = " + "))
)
gam_sir <- gam(gam_formula, data = gam_df, method = "REML")

summary(gam_sir)
plot(gam_sir, pages = 1, shade = TRUE)

# ------------------------------------------


# PCA
pca_out <- prcomp(X_fullrank, center = FALSE, scale. = FALSE)
X_pca <- pca_out$x[, 1:4]  # First two PCs

# Create a data frame for plotting
pca_df <- as.data.frame(X_pca, phq_cat)
colnames(pca_df) <- c("PCA1", "PCA2")
pca_df$Severity <- phq_cat

# Filter out NA severity values
pca_df <- pca_df %>% filter(!is.na(Severity))

# Plot using ggplot2
ggplot(pca_df, aes(x = PCA1, y = PCA2, color = Severity)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  labs(
    title = "PCA Projection by PHQ-9 Depression Severity",
    x = "PCA1",
    y = "PCA2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Scree plot
plot(pca_out$sdev^2, type = "b", pch = 19,
     xlab = "PCA Component", ylab = "Eigenvalue",
     main = "Scree Plot of PCA Directions")
abline(h = 1e-3, col = "red", lty = 2)

# Use complete cases
idx <- complete.cases(y, X_pca[,1:4])

# Fit GAM with smooth terms
gam_pca <- gam(y[idx] ~ s(X_pca[idx,1]) + s(X_pca[idx,2]) + s(X_pca[idx,3]) + s(X_pca[idx,4]), method = "REML")

# Summary
summary(gam_pca)

plot(gam_pca, pages = 1, shade = TRUE)

# Get loadings
loadings <- pca_out$rotation

# Top 10 variables for PC1
top_pc1 <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:10]
top_pc1_vars <- names(top_pc1)

# Top 10 variables for PC2
top_pc2 <- sort(abs(loadings[, 2]), decreasing = TRUE)[1:10]
top_pc2_vars <- names(top_pc2)

# Top 10 variables for PC3
top_pc3 <- sort(abs(loadings[, 3]), decreasing = TRUE)[1:10]
top_pc3_vars <- names(top_pc3)

pc1_vars <- data.frame(
  PC1_Loading = loadings[names(top_pc1), 1]
)

pc2_vars <- data.frame(
  PC2_Loading = loadings[names(top_pc2), 2]
)

data.frame(
  PC3_Loading = loadings[names(top_pc3), 3]
)

plot(y[idx], X_pca[idx,1], pch = 19, col = rgb(0,0,0,0.1),
     xlab = "IRT Theta", ylab = "PC1", main = "PC1 vs Latent Depression")
lines(lowess(y[idx], X_pca[idx,1]), col = "red", lwd = 2)
plot(y[idx], X_pca[idx,2], pch = 19, col = rgb(0,0,0,0.1),
     xlab = "IRT Theta", ylab = "PC1", main = "PC2 vs Latent Depression")
lines(lowess(y[idx], X_pca[idx,2]), col = "red", lwd = 2)

# LASSO
# Fit LASSO with cross-validation
cv_lasso <- cv.glmnet(X_fullrank, y, alpha = 1, standardize = FALSE)  # alpha = 1 for LASSO
best_lambda <- cv_lasso$lambda.min

# Extract non-zero coefficients
lasso_coef <- coef(cv_lasso, s = best_lambda)
active <- lasso_coef[lasso_coef[, 1] != 0, , drop = FALSE]
active <- as.data.frame(as.matrix(active))
active$Variable <- rownames(active)
colnames(active)[1] <- "Coefficient"
active <- active[order(-abs(active$Coefficient)), ]

# View top selected variables (excluding intercept)
active_vars <- active[active$Variable != "(Intercept)", ]
lasso_variables <- head(active_vars, 10)

lasso_vars <- active_vars$Variable[1:10]

# Create the modeling dataframe with theta_dep and selected predictors
lasso_df <- df %>%
  select(all_of(c("theta_dep", lasso_vars))) %>%
  drop_na()

smooth_terms <- sapply(lasso_vars, function(var) {
  n_unique <- length(unique(lasso_df[[var]]))
  if (is.numeric(lasso_df[[var]]) && n_unique > 10) {
    paste0("s(", var, ")")
  } else {
    var  # just use as linear/factor
  }
})
formula_str <- paste("theta_dep ~", paste(smooth_terms, collapse = " + "))
gam_lasso <- gam(
  formula = as.formula(formula_str),
  data = lasso_df,
  method = "REML"
)
summary(gam_lasso)
plot(gam_lasso, pages=1, shade=TRUE)

# xgboost
set.seed(1991)

dtrain <- xgb.DMatrix(data = as.matrix(X_fullrank), label = as.numeric(y))

# Fit XGBoost model
xgb_mod <- xgboost(
  data = dtrain,
  objective = "reg:squarederror",
  nrounds = 200,
  eta = 0.05,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8,
  verbose = FALSE
)

# Create explainer
shap_explainer <- shap.prep(xgb_model = xgb_mod, X_train = X_fullrank)

xg_vars <- shap_explainer %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value))) %>%
  arrange(desc(mean_abs_shap)) %>% head(10)

# Step 1: Safely pull the top 20 variable names
top_xg_vars <- xg_vars %>%
  slice_head(n = 10) %>%     # works even if fewer than 20
  pull(variable) %>%         # pull variable names as character
  as.character()             # ensure it's a character vector

# Step 2: Create the modeling dataframe
gam_df_xg <- df %>%
  select(all_of(c("theta_dep", top_xg_vars))) %>%
  drop_na()

# Step 3: Create formula
# Threshold for splines (change as needed)
spline_thresh <- 10

# Count unique values for each variable in top_vars
unique_counts <- sapply(gam_df_xg[top_vars], function(x) length(unique(x)))

# Split variables into spline and linear groups
vars_spline <- names(unique_counts[unique_counts >= spline_thresh])
vars_linear <- names(unique_counts[unique_counts < spline_thresh])
# Build terms
spline_terms <- paste0("s(", vars_spline, ")")
all_terms <- c(spline_terms, vars_linear)
formula_str <- paste("theta_dep ~", paste(all_terms, collapse = " + "))

# Step 4: Fit GAM
gam_xg <- gam(
  formula = as.formula(formula_str),
  data = gam_df_xg,
  method = "REML"
)

# Step 5: Review summary
summary(gam_xg)


plot(gam_xg, pages = 1, shade=TRUE)

# Data
model_df <- data.frame(
  Model = c("SIR GAM", "PCA GAM", "LASSO GAM", "XGBoost GAM"),
  R2_adj = c(0.109, 0.0421, 0.114, 0.116),
  Deviance_Explained = c(11.1, 4.3, 11.6, 12)
)

# Rescale R2_adj
max_deviance <- max(model_df$Deviance_Explained)
max_r2 <- max(model_df$R2_adj)
scale_factor <- max_deviance / max_r2
model_df$R2_adj_scaled <- model_df$R2_adj * scale_factor

# Pivot for plotting
model_long <- model_df %>%
  select(Model, Deviance_Explained, R2_adj_scaled) %>%
  pivot_longer(-Model, names_to = "Metric", values_to = "Value")

# Assign colors and metric labels
metric_labels <- c("Deviance_Explained" = "Deviance Explained (%)", 
                   "R2_adj_scaled" = "Adjusted R²")
metric_colors <- c("Deviance_Explained" = "#377eb8", 
                   "R2_adj_scaled" = "#e41a1c")

# Plot (position_dodge to offset bars)
sum_plot <- ggplot(model_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = metric_colors, labels = metric_labels, name = "Metric") +
  scale_y_continuous(
    name = "Deviance Explained (%)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Adjusted R²")
  ) +
  labs(title = "Model Comparison: Deviance Explained and Adjusted R²") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title.y.left = element_text(color = "#377eb8", face = "bold"),
    axis.title.y.right = element_text(color = "#e41a1c", face = "bold")
  )
sum_plot
ggsave("sum_plot.png", plot = sum_plot, width = 8, height = 6, dpi = 300)


# Calculate scaling
max_count <- max(hist(df$theta_dep, plot=FALSE)$counts)
max_density <- max(density(df$theta_dep)$y)
scale_factor <- max_count / max_density

dep_plot <- ggplot(df, aes(x = theta_dep)) +
  # Histogram
  geom_histogram(
    bins = 50, 
    fill = "gold", 
    color = "black", 
    alpha = 0.6, 
    linewidth = 0.4
  ) +
  # Density (scaled to match counts)
  geom_density(
    aes(y = ..density.. * scale_factor), 
    color = "#377eb8", 
    linewidth = 1.5, 
    alpha = 0.8
  ) +
  labs(
    title = expression("Distribution of Depression Severity Scores " ~ (hat(theta)[i])),
    x = "Depression Severity",
    y = "Count"
  ) +
  scale_y_continuous(
    name = "Count",
    sec.axis = sec_axis(~ . / scale_factor, name = "Density")
  ) +
  theme_bw(base_size = 17) +
  theme(
    plot.title = element_text(hjust = 0, size = 22, face = "bold", margin=margin(b=16)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 15),
    axis.title.y.right = element_text(color = "#377eb8", face = "bold"),
    axis.title.y.left = element_text(color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("dep_plot.png", plot = dep_plot, width = 8, height = 6, dpi = 300)


table_data <- data.frame(
  SIR1_Var = c("LBDSGBSI", "LBDSTPSI", "LBDSALSI", "BPXOSY3", "BPXOSY2", "RIDAGEYR", "BPXOSY1", "LBDBPBSI", "LBXNEPCT", "DMDEDUC2"),
  SIR1_Val = c(0.450, -0.431, 0.259, -0.184, -0.177, -0.169, -0.152, -0.134, -0.124, 0.119),
  
  SIR2_Var = c("LBXLYPCT", "LBXWBCSI", "LBXNEPCT", "LBDSALSI", "LBXHGB", "LBDLYMNO", "INDFMPIR", "MCQ160A", "LBXMOPCT", "PAQ650"),
  SIR2_Val = c(0.362, -0.347, 0.297, 0.205, 0.201, 0.174, 0.145, 0.139, 0.137, -0.137),
  
  PC1_Var = c("BMXWAIST", "BMXBMI", "BMXWT", "BMXARMC", "BMXHIP", "BPAOCSZ", "BPQ020", "BPXODI1", "BPXOSY2", "BPXODI2"),
  PC1_Val = c(-0.278, -0.256, -0.248, -0.237, -0.237, -0.206, 0.173, -0.169, -0.164, -0.164),
  
  PC2_Var = c("LBXHGB", "LBDUIBSI", "LBXHCT", "LBDPCT", "LBXMCHSI", "LBDIRNSI", "LBXMCVSI", "LBXMC", "LBXRDW", "LBDTIB"),
  PC2_Val = c(0.259, -0.243, 0.243, 0.236, 0.205, 0.204, 0.184, 0.174, -0.174, -0.161),
  
  LASSO_Var = c("RIDAGEYR", "MCQ160A", "MCQ160P", "BPQ020", "MCQ160L", "INDFMPIR", "PAQ620", "MCQ300A", "PAQ650", "MCQ160M"),
  LASSO_Val = c(-0.215, -0.125, -0.108, -0.097, -0.072, -0.069, -0.066, -0.062, 0.059, -0.056),
  
  SHAP_Var = c("MCQ160A", "RIDAGEYR", "LBXHCOT", "MCQ160P", "INDFMPIR", "LBDBPBSI", "PAQ620", "BPQ020", "BMXHIP", "LBXSCK"),
  SHAP_Val = c(0.0930, 0.0788, 0.0768, 0.0602, 0.0593, 0.0590, 0.0576, 0.0537, 0.0536, 0.0411)
)


table_data$Rank <- 1:10
 
# Move Rank to first column
table_data <- table_data %>% select(Rank, everything())

tab <- gt(table_data) %>%
  tab_header(
    title = md("**Top 10 Variables and Magnitudes for Each Method**")
  ) %>%
  fmt_number(columns = ends_with("_Val"), decimals = 3) %>%
  cols_label(
    SIR1_Var = "Variable", SIR1_Val = "Value",
    SIR2_Var = "Variable", SIR2_Val = "Value",
    PC1_Var = "Variable", PC1_Val = "Value",
    PC2_Var = "Variable", PC2_Val = "Value",
    LASSO_Var = "Variable", LASSO_Val = "Value",
    SHAP_Var = "Variable", SHAP_Val = "Value"
  ) %>%
  tab_spanner(label = "SIR1", columns = c(SIR1_Var, SIR1_Val)) %>%
  tab_spanner(label = "SIR2", columns = c(SIR2_Var, SIR2_Val)) %>%
  tab_spanner(label = "PC1", columns = c(PC1_Var, PC1_Val)) %>%
  tab_spanner(label = "PC2", columns = c(PC2_Var, PC2_Val)) %>%
  tab_spanner(label = "LASSO", columns = c(LASSO_Var, LASSO_Val)) %>%
  tab_spanner(label = "SHAP", columns = c(SHAP_Var, SHAP_Val)) %>%
  cols_align(align = "center", columns = everything()) %>%
  # Bold column labels
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  # Bold method spanner labels
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) %>%
  # Bold the table title
  tab_style(
    style = cell_text(weight = "bold", size = "large"),
    locations = cells_title(groups = "title")
  ) %>%
  # Bold all table body cells (data)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body()
  ) %>%
  # Thicker/darker borders
  tab_options(
    table.border.top.width = px(3),
    table.border.bottom.width = px(3),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.width = px(3),
    heading.border.bottom.color = "black",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    row_group.border.bottom.width = px(2),
    row_group.border.bottom.color = "black",
    data_row.padding = px(4)
  )
# Helper: which variable names are in each group
blood_cell_vars <- c("LBDSGBSI", "LBXWBCSI", "LBXLYPCT", "LBNEXPCT", "LBXHGB", "LBDUIBSI", "LBXHCT", "LBDPCT", "LBXMCHSI", "LBDIRNSI", "LBXMCVSI", "LBXMC", "LBXRDW", "LBDTIB", "LBXMOPCT", "LBDLYMNO", "LBXNEPCT")
blood_chem_vars <- c("LBDSALSI", "LBDBPBSI", "LBXHCOT", "LBXSCK", "LBDSTPSI")
anthro_vars <- c("BMXWAIST", "BMXBMI", "BMXWT", "BMXARMC", "BMXHIP", "BPAOCSZ")
demo_vars <- c("RIDAGEYR", "INDFMPIR", "DMDEDUC2", "PAQ620", "PAQ650", "MCQ160A", "MCQ160L", "MCQ300A", "MCQ160M")
survey_vars <- c("BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPQ020", "MCQ160P", "DMDMARTZ")

# Returns row indices in a column matching your group
rows_with_vars <- function(df, colname, vargroup) {
  which(df[[colname]] %in% vargroup)
}
tab <- tab %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = SIR1_Var, rows = rows_with_vars(table_data, "SIR1_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = SIR2_Var, rows = rows_with_vars(table_data, "SIR2_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = PC1_Var, rows = rows_with_vars(table_data, "PC1_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = PC2_Var, rows = rows_with_vars(table_data, "PC2_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = LASSO_Var, rows = rows_with_vars(table_data, "LASSO_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = SHAP_Var, rows = rows_with_vars(table_data, "SHAP_Var", blood_cell_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = SIR1_Var, rows = rows_with_vars(table_data, "SIR1_Var", blood_chem_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = SIR2_Var, rows = rows_with_vars(table_data, "SIR2_Var", blood_chem_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = PC1_Var, rows = rows_with_vars(table_data, "PC1_Var", blood_chem_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = PC2_Var, rows = rows_with_vars(table_data, "PC2_Var", blood_chem_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = LASSO_Var, rows = rows_with_vars(table_data, "LASSO_Var", blood_chem_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightsalmon"),
    locations = cells_body(columns = SHAP_Var, rows = rows_with_vars(table_data, "SHAP_Var", blood_chem_vars))
  ) %>% 
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = SIR1_Var, rows = rows_with_vars(table_data, "SIR1_Var", anthro_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = SIR2_Var, rows = rows_with_vars(table_data, "SIR2_Var", anthro_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = PC1_Var, rows = rows_with_vars(table_data, "PC1_Var", anthro_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = PC2_Var, rows = rows_with_vars(table_data, "PC2_Var", anthro_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = LASSO_Var, rows = rows_with_vars(table_data, "LASSO_Var", anthro_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "khaki"),
    locations = cells_body(columns = SHAP_Var, rows = rows_with_vars(table_data, "SHAP_Var", anthro_vars))
  ) %>% 
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = SIR1_Var, rows = rows_with_vars(table_data, "SIR1_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = SIR2_Var, rows = rows_with_vars(table_data, "SIR2_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = PC1_Var, rows = rows_with_vars(table_data, "PC1_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = PC2_Var, rows = rows_with_vars(table_data, "PC2_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = LASSO_Var, rows = rows_with_vars(table_data, "LASSO_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = SHAP_Var, rows = rows_with_vars(table_data, "SHAP_Var", demo_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = SIR1_Var, rows = rows_with_vars(table_data, "SIR1_Var", survey_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = SIR2_Var, rows = rows_with_vars(table_data, "SIR2_Var", survey_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = PC1_Var, rows = rows_with_vars(table_data, "PC1_Var", survey_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = PC2_Var, rows = rows_with_vars(table_data, "PC2_Var", survey_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = LASSO_Var, rows = rows_with_vars(table_data, "LASSO_Var", survey_vars))
  ) %>%
  tab_style(
    style = cell_fill(color = "violet"),
    locations = cells_body(columns = SHAP_Var, rows = rows_with_vars(table_data, "SHAP_Var", survey_vars))
  ) %>% 
  tab_source_note(
    source_note = md(
      "<span style='background-color:lightblue'>&nbsp;&nbsp;&nbsp;</span> Blood Cell Indices&nbsp;&nbsp;&nbsp;
      <span style='background-color:lightsalmon'>&nbsp;&nbsp;&nbsp;</span> Blood Chemistry/Metabolic&nbsp;&nbsp;&nbsp;
      <span style='background-color:khaki'>&nbsp;&nbsp;&nbsp;</span> Anthropometric&nbsp;&nbsp;&nbsp;
      <span style='background-color:lightgreen'>&nbsp;&nbsp;&nbsp;</span> Demographic/Socioeconomic&nbsp;&nbsp;&nbsp;
      <span style='background-color:violet'>&nbsp;&nbsp;&nbsp;</span> Survey/Medical History"
    )
  )

tab

gtsave(tab, "sum_table.png", expand=10)

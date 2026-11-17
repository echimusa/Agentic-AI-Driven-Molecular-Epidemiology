import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.decomposition import PCA

np.random.seed(99)
n_ind = 1000

# 1. Simulate ancestral split
# 600 individuals are Ancestry A; 400 are Ancestry B
ancestry = np.array(['A'] * 600 + ['B'] * 400)

# Simulate SNP risk allele. In Ancestry B, this SNP is more common.
# Minor allele frequency: Ancestry A = 10%, Ancestry B = 45%
snp = np.zeros(n_ind)
for i in range(n_ind):
    p = 0.10 if ancestry[i] == 'A' else 0.45
    snp[i] = np.random.binomial(n=2, p=p) # Count of risk alleles (0, 1, 2)

# Unmeasured environmental confounder in Ancestry B leads to higher disease rate
base_risk = np.where(ancestry == 'A', -2.5, -0.5)

# The SNP has ZERO causal effect on disease (true beta = 0.0)
log_odds = base_risk + 0.0 * snp 
disease_prob = 1 / (1 + np.exp(-log_odds))
disease = np.random.uniform(0, 1, n_ind) < disease_prob

# Create genotype background Matrix (100 unrelated SNPs to calculate PCA)
genotype_matrix = np.random.binomial(n=2, p=0.2, size=(n_ind, 100))
# Let ancestral differences slide into genotyping background data
for i in range(n_ind):
    if ancestry[i] == 'B':
        genotype_matrix[i, 0:30] = np.random.binomial(n=2, p=0.6, size=30)

df = pd.DataFrame({
    'disease': disease.astype(int),
    'test_snp': snp,
    'ancestry_label': ancestry
})

# 2. Naive Association Test (No Ancestry / PCA Controls)
X_naive = sm.add_constant(df['test_snp'])
model_naive = sm.Logit(df['disease'], X_naive).fit()

# 3. Derive PCA adjustment covariates
pca = PCA(n_components=2)
pcs = pca.fit_transform(genotype_matrix)
df['PC1'] = pcs[:, 0]
df['PC2'] = pcs[:, 1]

# 4. Stratification-Adjusted Model
X_adjusted = sm.add_constant(df[['test_snp', 'PC1', 'PC2']])
model_adjusted = sm.Logit(df['disease'], X_adjusted).fit()

# Display Results (Observe false positive OR in naive model -> resolved in adjusted)
print("\n--- Model Output Summary ---")
print("Naive Model Odds Ratio (SNP):", np.round(np.exp(model_naive.params['test_snp']), 3))
print("Naive Model P-Value:", np.round(model_naive.pvalues['test_snp'], 5))
print("\nAdjusted (PC) Model Odds Ratio (SNP):", np.round(np.exp(model_adjusted.params['test_snp']), 3))
print("Adjusted Model P-Value:", np.round(model_adjusted.pvalues['test_snp'], 5))
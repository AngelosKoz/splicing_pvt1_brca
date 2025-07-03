# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import shapiro
from statsmodels.stats.multitest import multipletests
from tabulate import tabulate

# %%
##--------- Exploration of Data ---------##
## Variance of features (splice sites)
def get_var(directory, data, length=100, primeSS='3', top_var=6, top_ind=0.1, plots=True, gene_naming=''):
    plots_dir = os.path.join(directory, f'SS{primeSS}/{gene_naming}_top_{top_var}/plots_L{length}/fix_variance_plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    variances = data.loc[:, data.columns != 'label'].var().sort_values(ascending=False)
    feature_names = variances.index.tolist()
    total_variance = np.sum(variances)
    var_explained = (variances / total_variance) * 100
    
    top_indices = np.argsort(var_explained)[-top_var:]
    top_features = {feature_names[i]: var_explained[i] for i in top_indices}
    top_features = dict(sorted(top_features.items(), key=lambda item: item[1], reverse=True))
    
    lowest_indices = np.argsort(var_explained)[:10]
    lowest_features = {feature_names[i]: var_explained[i] for i in lowest_indices}
    
    dict_filename = os.path.join(directory, f'SS{primeSS}/{gene_naming}_top_{top_var}/{gene_naming}_Top_{top_var}_Feature_Variances_SE_ss{primeSS}_L{length}.csv')
    dict_filename_low = os.path.join(directory, f'SS{primeSS}/{gene_naming}_top_{top_var}/{gene_naming}_Lowest_10_Feature_Variances_SE_ss{primeSS}_L{length}.csv')

    pd.DataFrame(top_features.items(), columns=['Splice_Site', 'Variance']).to_csv(dict_filename, sep='\t', index=False)
    pd.DataFrame(lowest_features.items(), columns=['Splice_Site', 'Variance']).to_csv(dict_filename_low, sep='\t', index=False)
    
    if plots:
        plot_filename = os.path.join(plots_dir, f'{gene_naming}_EDA_RawVar_SE_ss{primeSS}_L{length}.pdf')
        plt.figure(figsize=(20, 15))
        bars = plt.bar(range(len(variances)), variances, color='skyblue')
        raw_min = np.min(np.sort(variances)[-top_var:])
        plt.axhline(y=raw_min, color='r', linestyle='-', label=f'Min Top Variance: {raw_min:.4f}')
        plt.title('Raw Variance of Splice Sites')
        plt.xlabel('Splice Sites')
        plt.ylabel('Variance')
        plt.xticks(ticks=np.arange(len(variances)), labels=feature_names, rotation=90)
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close()

    return list(top_features.keys()), list(lowest_features.keys())

# %%
# Example usage for pvt1
length = 50
ss = '3'
top_var = 34
main_dir = '/data/chromatin_associated_genes/pvt1/condition_pvt1/'
data_path = os.path.join(main_dir, f'condition_altSS_concat_Control_Tumor_SE_SS{ss}.altSS_L{length}.csv')
data_df = pd.read_csv(data_path, sep='\t', index_col=0)# Adjust as needed
top_features, lowest_features = get_var(main_dir, data_df, length=length, primeSS=ss, top_var=top_var, plots=True)
data_top_SS = data_df[top_features + ['label']]
data_low_SS = data_df[lowest_features + ['label']]
data_label = data_df[['label']].copy()


# %%
### Statistical Test ###
data_top_label = data_top_SS.copy()
data_top_label['label'] = data_label['original_labels']
data_top = data_top_label.drop(columns='label')
data_top_control = data_top_label[data_top_label['label'] == 0].drop(columns='label')
data_top_tumor = data_top_label[data_top_label['label'] == 1].drop(columns='label')

normality_results = {}
stat_test_results = {}

for feature in data_top.columns:
    _, p_shapiro = shapiro(data_top[feature])
    normality_results[feature] = not (p_shapiro < 0.05)
    
    if p_shapiro < 0.05:
        _, p_val = stats.mannwhitneyu(data_top_control[feature], data_top_tumor[feature], alternative='two-sided')
    else:
        _, p_val = stats.ttest_ind(data_top_control[feature], data_top_tumor[feature], equal_var=False)
        
    stat_test_results[feature] = p_val

sorted_results = dict(sorted(stat_test_results.items(), key=lambda item: item[1]))


### Adjusted p-values
# 1) Using FDR
p_values = list(sorted_results.values())
adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
sorted_adjusted_fdr = dict(sorted(zip(sorted_results.keys(), adjusted_p_values), key=lambda item: item[1]))

# 2) Using Bonferroni
adjusted_results = multipletests(p_values, alpha=0.05, method='bonferroni')
pvals_corrected = adjusted_results[1]  # Adjusted p-values
reject_array = adjusted_results[0]  # Boolean array indicating which hypotheses are rejected
sorted_adjusted_bonf = dict(sorted(zip(sorted_results.keys(), pvals_corrected), key=lambda item: item[1]))

print(tabulate(normality_results.items(), headers=['Splice_Site', 'p-value']))
print(tabulate(sorted_results.items(), headers=['Splice_Site', 'p-value']))
print(tabulate(sorted_adjusted_fdr.items(), headers=['Splice_Site', 'p-value (fdr_adj)']))
print(tabulate(sorted_adjusted_bonf.items(), headers=['Splice_Site', 'p-value (bonf_adj)']))

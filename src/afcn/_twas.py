"""Transcritome wide associations

By: 
    *Michael J. Betti
        - Vanderbilt University Medical Center
        - Eric Gamazon Lab

contributors:
    * Robert Vogel
        - Seattle Children's Research Institute;
        - Genomic Data Modeling Lab of Pejman Mohammadi

"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression, LinearRegression
import statsmodels.api as sm
from scipy.stats import norm


def read_pheno(pheno_file):
    pheno = pd.read_csv(pheno_file, header=None, sep="\t")
    if pheno.iloc[0, 0] == "ID":
        pheno.columns = pheno.iloc[0]
        pheno = pheno.iloc[1:]
    return pheno

def reduce_pheno(pheno, pheno_name=None):
    if pheno_name:
        pheno = pheno[['ID', pheno_name]]
    else:
        pheno = pheno[['ID', pheno.columns[-1]]]
    pheno.columns = ["ID", "phenotype"]
    pheno["phenotype"] = pd.to_numeric(pheno["phenotype"])
    return pheno

def read_filter(filter_file, filter_column=3):
    fil = pd.read_csv(filter_file, header=None, sep="\t")
    if fil.iloc[0, 0] == "ID":
        fil.columns = fil.iloc[0]
        fil = fil.iloc[1:]
    fil = fil[['ID', fil.columns[filter_column-1]]]
    fil.columns = ["ID", "fil_val"]
    return fil

def read_predicted(pred_exp_file):
    pred_exp = pd.read_csv(pred_exp_file, sep="\t")
    return pred_exp

# The predicted expression was written in the following format:
# fout.write_line_record(v[fpars.idx("chrom")],
#                                    v[fpars.idx("gene_start")],
#                                    v[fpars.idx("gene_end")],
#                                    gene_id,
#                                    gene_expr)

def merge_and_filter(pheno, pred_exp, fil=None, filter_val=1):
    merged = pd.merge(pheno, pred_exp, on=["ID"], how="inner")
    if fil is not None:
        merged = pd.merge(merged, fil, on=["ID"], how="inner")
        merged = merged[merged["fil_val"] == filter_val]
    return merged

import statsmodels.api as sm

import statsmodels.api as sm
from scipy.stats import norm

def association(merged, genes, test_type="linear"):
    assoc_df = pd.DataFrame(columns=["gene", "beta", "statistic", "p", "se(beta)"])  

    for gene in genes:
        X = sm.add_constant(merged[[gene]])  # Add a constant term for the intercept
        y = merged["phenotype"]

        if test_type == "logistic":
            model = sm.Logit(y, X)
            try:
                result = model.fit_regularized(method='l1')
            except Exception as e:
                print(f"Error fitting gene {gene}: {e}")
                continue
            beta = result.params[gene]
            stat = result.tvalues[gene] if hasattr(result, 'tvalues') else None
            p = result.pvalues[gene] if hasattr(result, 'pvalues') else None
            se = result.bse[gene] if hasattr(result, 'bse') else None
        else:
            model = sm.OLS(y, X)
            result = model.fit()
            beta = result.params[gene]
            stat = result.tvalues[gene] if hasattr(result, 'tvalues') else None
            p = result.pvalues[gene] if hasattr(result, 'pvalues') else None
            se = result.bse[gene] if hasattr(result, 'bse') else None

        assoc_df = pd.concat([assoc_df,
                              pd.DataFrame([{"gene": gene,
                                             "beta": beta,
                                             "statistic": stat,
                                             "p": p,
                                             "se(beta)": se}])],
                             ignore_index=True)

    return assoc_df

def write_association(assoc_df, output_file):
    assoc_df.to_csv(output_file, index=False)

def run(pheno_file, pheno_name,
        filter_file, filter_column, filter_val,
        pred_exp_file,
        test_type,
        missing_phenotype,
        drop_nans,
        file_out):
    # Get Arguments
    
    # Run functions
    pheno = read_pheno(pheno_file)
    pheno = reduce_pheno(pheno, pheno_name)
    
    fil_df = None

    if filter_file is not None:
        fil_df = read_filter(filter_file,
                             filter_column)
    
    pred_exp = read_predicted(pred_exp_file)
    genes = pred_exp.columns[1:]  # First column is ID
    
    merged = merge_and_filter(pheno, pred_exp, fil_df,
                              filter_val)
    
    if drop_nans:
        merged = merged.dropna(subset=['phenotype'])

    if missing_phenotype is not None:
        merged = merged[merged['phenotype'] != missing_phenotype]

    
    if merged.shape[0] == 0:
        raise ValueError("Filtered out all rows of phenotype")
    
    if test_type == 'logistic' and len(merged['phenotype'].unique()) > 2:
        raise ValueError("Logsitic regression requires binary phenotypes")
    
    assoc_df = association(merged, genes, test_type)
    write_association(assoc_df, file_out)
    print("Done. Results saved in", file_out)


def run():
    print(f"Successfully found {__file__}.")
    raise NotImplementedError

import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression, LinearRegression

def read_pheno(pheno_file):
    pheno = pd.read_csv(pheno_file, header=None)
    if pheno.iloc[0, 0] == "ID":
        pheno.columns = pheno.iloc[0]
        pheno = pheno.iloc[1:]
    return pheno

def reduce_pheno(pheno, pheno_name=None):
    if pheno_name:
        pheno = pheno[['ID', pheno_name]]
    else:
        pheno = pheno[['ID', pheno.columns[-1]]]
    pheno.columns = ["id", "phenotype"]
    pheno["phenotype"] = pd.to_numeric(pheno["phenotype"])
    return pheno

def read_filter(filter_file, filter_column=3):
    fil = pd.read_csv(filter_file, header=None)
    if fil.iloc[0, 0] == "ID":
        fil.columns = fil.iloc[0]
        fil = fil.iloc[1:]
    fil = fil[['ID', fil.columns[filter_column-1]]]
    fil.columns = ["id", "fil_val"]
    return fil

def read_predicted(pred_exp_file):
    pred_exp = pd.read_csv(pred_exp_file)
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

def association(merged, genes, test_type="linear"):
    assoc_df = pd.DataFrame(columns=["gene", "beta", "t", "p", "se(beta)"])

    for gene in genes:
        X = merged[[gene]].values.reshape(-1, 1)
        y = merged["phenotype"]

        if test_type == "logistic":
            model = LogisticRegression()
        elif test_type == "linear":
            model = LinearRegression()
        model.fit(X, y)
        if test_type == "logistic":
            beta = model.coef_[0][0]
            t, p, se = None, None, None  # Logistic Regression doesn't have t, p, se(beta)
        else:
            beta = model.coef_[0]
            t = model.coef_[0] / model.coef_std_
            p = model.p_values_
            se = model.coef_std_
        assoc_df = assoc_df.append({"gene": gene, "beta": beta, "t": t, "p": p, "se(beta)": se}, ignore_index=True)

    return assoc_df

def write_association(assoc_df, output_file):
    assoc_df.to_csv(output_file, index=False)

def main():
    # Get Arguments
    args = sys.argv[1:]

    # Set default values for arguments and set to correct data types
    arg_dict = dict(zip(args[::2], args[1::2]))
    arg_dict.setdefault('PHENO_COLUMN', 'None')
    arg_dict.setdefault('PHENO_NAME', 'None')
    arg_dict.setdefault('FILTER_COLUMN', '2')
    arg_dict.setdefault('FILTER_VAL', '1')
    arg_dict.setdefault('TEST_TYPE', 'linear')
    arg_dict.setdefault('MISSING_PHENOTYPE', 'NA')
    
    # Run functions
    pheno = read_pheno(arg_dict['PHENO_FILE'])
    pheno = reduce_pheno(pheno, arg_dict['PHENO_NAME'])
    
    if arg_dict['FILTER_FILE'] == 'None':
        fil_df = None
    else:
        fil_df = read_filter(arg_dict['FILTER_FILE'], int(arg_dict['FILTER_COLUMN']))
    
    pred_exp = read_predicted(arg_dict['PRED_EXP_FILE'])
    genes = pred_exp.columns[1:]  # First column is ID
    
    merged = merge_and_filter(pheno, pred_exp, fil_df, int(arg_dict['FILTER_VAL']))
    
    if arg_dict['MISSING_PHENOTYPE'].upper() == 'NA':
        merged = merged.dropna(subset=['phenotype'])
    else:
        merged = merged[merged['phenotype'] != float(arg_dict['MISSING_PHENOTYPE'])]
    
    if merged.shape[0] == 0:
        print("ERROR: Filtered out all rows of phenotype")
        return
    
    if arg_dict['TEST_TYPE'] == 'logistic' and len(merged['phenotype'].unique()) > 2:
        print("ERROR: For logistic tests, phenotype column can only have 2 values: 0 - unaffected, 1 - affected")
        return
    
    assoc_df = association(merged, genes, arg_dict['TEST_TYPE'])
    write_association(assoc_df, arg_dict['OUT'])
    print("Done. Results saved in", arg_dict['OUT'])

if __name__ == "__main__":
    main()


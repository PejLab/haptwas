from unittest import TestCase, main
import numpy as np
from afcn import model, _twas
from scipy.stats import truncnorm
import pandas as pd

np.random.seed(42)

class TestTWAS(TestCase):

    def setUp(self):
        # Define parameters
        n_genes = 100  # Number of genes
        n_samples = 100  # Number of individuals
        num_regulators_per_gene = 5  # Number of regulators for each gene

        # Generate sample names for individuals
        sample_names = ['Sample' + str(i).zfill(3) for i in range(n_samples)]

        # Generate gene expression data
        gene_ids = ['GENE' + str(i).zfill(3) for i in range(n_genes)]
        gene_expression = np.random.normal(loc=0, scale=1, size=(n_samples, n_genes))

        # Define regulatory relationships between genes
        regulatory_network = {}

        for gene in range(n_genes):
            # Randomly select regulators for the current gene
            regulators = np.random.choice(np.arange(n_genes), size=num_regulators_per_gene, replace=False)
            regulatory_network[gene] = regulators

        # Calculate gene expression levels influenced by regulatory relationships
        regulated_expression = np.zeros_like(gene_expression)

        for gene, regulators in regulatory_network.items():
            # Aggregate expression levels of regulator genes
            regulator_expression = np.sum(gene_expression[:, regulators], axis=1)
            # Calculate regulated expression as a function of regulator expression
            regulated_expression[:, gene] = gene_expression[:, gene] + regulator_expression

        # Normalize regulated expression to have a mean of 0
        regulated_expression -= np.mean(regulated_expression, axis=0)

        # Define logistic function for binary phenotype simulation
        def logistic(x):
            return 1 / (1 + np.exp(-x))

        # Calculate probability of being affected using logistic function
        logistic_alpha = 0
        logistic_beta = 1  # Adjust as needed
        phenotype_prob = logistic(logistic_alpha + logistic_beta * regulated_expression)

        # Reshape phenotype_prob into a 1-dimensional array
        phenotype_prob_1d = phenotype_prob.ravel()

        # Generate binary phenotypes based on the probability for each individual
        binary_phenotypes = np.random.binomial(1, p=phenotype_prob_1d[:n_samples], size=n_samples)

        # Define parameters for simulating continuous phenotype
        mean_continuous_phenotype = 0
        std_continuous_phenotype = 1

        # Generate continuous phenotype using normal distribution
        continuous_phenotypes = np.random.normal(loc=mean_continuous_phenotype, scale=std_continuous_phenotype, size=n_samples)

        # Create DataFrame for gene expression
        self.expression_df = pd.DataFrame(gene_expression, columns=gene_ids, index=sample_names)

        # Reset index and rename the first column to "ID"
        self.expression_df = self.expression_df.reset_index().rename(columns={"index": "ID"})

        # Create DataFrame for binary phenotypes
        self.binary_phenotype_df = pd.DataFrame({'ID': sample_names, 'phenotype': binary_phenotypes})

        # Create DataFrame for continuous phenotypes
        self.continuous_phenotype_df = pd.DataFrame({'ID': sample_names, 'phenotype': continuous_phenotypes})

    def tearDown(self):
        pass

    def test_binary_twas(self):
        #Test that simulated expression and phenotypes have the expected size and format
        self.assertEqual(len(self.expression_df.index), 100)
        assert self.binary_phenotype_df.columns.tolist() == ["ID", "phenotype"]

        merged = _twas.merge_and_filter(self.binary_phenotype_df, self.expression_df, fil=None, filter_val=1)
        assoc = _twas.association(merged, self.expression_df.columns[1:], test_type="logistic")

        #Test that all column names are correct in the association results
        assert assoc.columns.tolist() == ["gene", "beta", "statistic", "p", "se(beta)"]

        #Test that all p-values are between 0 and 1
        assert assoc["p"].between(0, 1).all()
        
    def test_continuous_twas(self):
        #Test that simulated expression and phenotypes have the expected size and format
        self.assertEqual(len(self.expression_df.index), 100)
        assert self.continuous_phenotype_df.columns.tolist() == ["ID", "phenotype"]

        merged = _twas.merge_and_filter(self.continuous_phenotype_df, self.expression_df, fil=None, filter_val=1)
        assoc = _twas.association(merged, self.expression_df.columns[1:], test_type="linear")
       
        #Test that all column names are correct in the association results
        assert assoc.columns.tolist() == ["gene", "beta", "statistic", "p", "se(beta)"]

        #Test that all p-values are between 0 and 1
        assert assoc["p"].between(0, 1).all()

if __name__ == "__main__":
    main()
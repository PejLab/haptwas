import unittest
import numpy as np
from afcn import model


class TestLinearized(unittest.TestCase):
    j_vars = 5
    n_samples = 1000
    q_values = j_vars * n_samples

    def setUp(self):
        if not hasattr(self, "rng"):
            self.rng = np.random.default_rng()

        self.hap_one = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.hap_two = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.effect_sizes = self.rng.normal(0,2,size=self.j_vars)

        self.y = (np.exp(self.hap_one @ self.effect_sizes)
                  + np.exp(self.hap_one @ self.effect_sizes))


    def test_output_types_properties(self):

        lsq = model._linear_expansion_model(self.hap_one,
                                            self.hap_two,
                                            self.y)

        self.assertIsInstance(lsq, dict)
        self.assertEqual(len(lsq), 2)
        for key in ["pars", "rank"]:
            self.assertTrue(key in lsq)


        self.assertEqual(lsq["pars"].size, self.j_vars+1)
        self.assertEqual(lsq["rank"], self.j_vars+1)

    # TODO
    def test_out_parameter_value(self):
        pass


class TestObj(unittest.TestCase):
    j_vars = 5
    n_samples = 1000
    q_values = j_vars * n_samples

    def setUp(self):
        if not hasattr(self, "rng"):
            self.rng = np.random.default_rng()

        self.hap_one = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.hap_two = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.effect_sizes = self.rng.normal(0,2,size=self.j_vars)

        self.y = (np.exp(self.hap_one @ self.effect_sizes)
                  + np.exp(self.hap_one @ self.effect_sizes))

    def test_variables_and_output_type(self):

        reg = "l1"
        regconst = 2

        input_vals = (self.hap_one, self.hap_two, self.y,
                      "l1", 2)
        f = model._obj(*input_vals)
        self.assertTrue(callable(f))

        # test to make sure that the variables associated with
        # the closure are those that I intended.  Moreover,
        # the first elment of the closure is the lambda function
        # as reg and reg_const are not used in _g, they are not
        # stored in the closure.  Consequently, I am only testing
        # the matching of haplotypes and y

        for v, c in zip(input_vals, f.__closure__[1:]):

            self.assertTrue((v == c.cell_contents).all())
        

    def test_l1_penalty(self):
        reg = "l1"
        regconst = 2

        # recall the number of parameters is j_vars +1, where
        # the 1 accounts for the log reference abundance
        pars = np.ones(self.j_vars+1)

        # when all the haplotypes are match the reference
        # allele, then I expect the model prediction to produce
        # 2* exp(log_ref_const)  = 2 * exp2(1) = 4
        y = np.zeros(self.n_samples) + np.log(1 + 4)

        haplotypes = np.zeros(shape=self.hap_one.shape)

        # given these inputs I expect that the only term 
        # of the sum squared residuals that is non-zero
        # is the penalty term.  As each parameter = 1, then the 
        # penalty term under l1 is simply regconst * (j_vars + 1)
        f = model._obj(haplotypes, haplotypes, y, "l1", regconst)

        self.assertEqual(f.__closure__[0].cell_contents.__closure__[0].cell_contents,
                         regconst)

        self.assertEqual(f(pars), regconst*(self.j_vars + 1))
        

        # 2* exp(log_ref_const)  = 2 * exp2(-1) = 1
        pars = -np.ones(self.j_vars+1)
        y = np.zeros(self.n_samples) + np.log(1 + 1)

        f = model._obj(haplotypes, haplotypes, y, "l1", regconst)
        self.assertEqual(f(pars), regconst*(self.j_vars + 1))


    def test_no_penalty(self):
        pass

    def test_l2_penalty(self):
        pass

if __name__ == "__main__":
    unittest.main()

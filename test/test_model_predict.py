"""Test the accuracy and robustness of expression prediction."""

from unittest import TestCase, main
import numpy as np
from afcn import model



class Test_predict(TestCase):

    def test_model_with_one_sample(self):
        # under this genotype the expression should always be
        # e^(3.5) + e^(1.5)
        h1 = np.array([0,1, 1])
        h2 = np.array([1, 0, 1])
        beta = np.array([1, 3, 0.5])
        alpha = 0
        y_true = np.exp2(3.5) + np.exp2(1.5)

        self.assertAlmostEqual(y_true,
                               model._predict(h1, h2, alpha, beta))

        h1 = h1.reshape(1,3)
        self.assertAlmostEqual(y_true,
                               model._predict(h1,
                                              h2,
                                              alpha, beta))

        # haplotype transposing results in a (3,1) dot with a (3,)
        # and should result in value error.  Note that rows should
        # always be number of samples, not SNV alleles

        with self.assertRaises(ValueError):
            self.assertAlmostEqual(y_true,
                               model._predict(h1.transpose(), 
                                              h2, 
                                              alpha, beta))


    def test_with_one_variant(self):
        # under this genotype, the expression should always be
        # an array of length 3 and the corresponding
        # values below
        h1 = np.array([0,1, 1]).reshape(3,1)
        h2 = np.array([1, 0, 1]).reshape(3,1)
        beta = np.array([3])
        alpha = 0

        y_true = np.array([1 + np.exp2(3),
                           np.exp2(3) + 1,
                           2*np.exp2(3)])

        y_predict = model._predict(h1,h2,alpha,beta)

        self.assertTupleEqual(y_true.shape,
                              y_predict.shape)

        for true_val, predict_val in zip(y_true, y_predict):
            self.assertAlmostEqual(true_val,
                                   predict_val)


    def test_array_shape_many_variants(self):
        h1 = np.array([[0,1,1], [1,0,0]])
        h2 = np.array([[1,0,1], [0,1,0]])
        beta = np.array([3, 0.5, 1])
        alpha = 1

        y_predict = model._predict(h1,h2,alpha,beta)

        self.assertTupleEqual((2,), y_predict.shape)


    def test_model_equation(self):

        # y = (e^(0 + 0 * 1) + e^(0 + 0 * 1)) = 2
        self.assertAlmostEqual(2,
                               model._predict(0, 0,
                                              0,
                                              1)
                               )
        
        # y = (e^(1 + 0 * 1) + e^(1 + 0 * 1)) = 2e
        self.assertAlmostEqual(2 * np.exp2(1),
                               model._predict(0, 0,
                                              1,
                                              1)
                               )

        # y = (e^(0 + 1 * 2) + e^(0 + 0 * 2)) = e^2 + 1
        self.assertAlmostEqual(np.exp2(2) + 1,
                               model._predict(1, 0,
                                              0,
                                              2)
                               )

        # y = e^(0 + dot([1,1],[1,2])) + e^(0 + dot([1,0][1,2])
        # = e^3 + e^1
        self.assertAlmostEqual(np.exp2(3) + np.exp2(1),
                               model._predict(np.array([1,1]), 
                                              np.array([1,0]),
                                              0,
                                              np.array([1,2]))
                               )


        # y = e^(0 + dot([1,1],[1,2])) + e^(0 + dot([0,1][1,2])
        # = e^3 + e^2
        self.assertAlmostEqual(np.exp2(3) + np.exp2(2),
                               model._predict(np.array([1,1]), 
                                              np.array([0,1]),
                                              0,
                                              np.array([1,2]))
                               )

        # haplotype consists of 100 variants, all ref allele
        self.assertAlmostEqual(2,
                               model._predict(np.zeros(100), 
                                              np.zeros(100),
                                              0,
                                              np.ones(100))
                               )
        # haplotype consists of 100 variants, on haplotype
        # is only ref, the other only alt

        self.assertAlmostEqual(np.exp2(100) + 1,
                               model._predict(np.zeros(100), 
                                              np.ones(100),
                                              0,
                                              np.ones(100))
                               )
        


class TestPredict(TestCase):
    rng = np.random.default_rng()
    h1 = rng.choice([0,1], size=100, replace=True)
    h2 = rng.choice([0,1], size=100, replace=True)

    alpha = rng.normal()
    beta = rng.normal(size=100)

    def test_biallelic_data_check(self):
        """Test whether biallelic data check is called returned values.

        Unit tests for is_biallelic func more rigorously checks the
        function.  Here, we are just checking that some method is called 
        """
        tmp = self.h1.copy()
        tmp[24] = 4
        with self.assertRaises(ValueError):
            model.predict(self.h1, tmp, self.alpha, self.beta)

        with self.assertRaises(ValueError):
            model.predict(tmp, self.h2, self.alpha, self.beta)


    def test_beta_check(self):
        with self.assertRaises(ValueError):
            model.predict(self.h1, 
                          self.h2, 
                          self.alpha, 
                          np.hstack([self.beta, self.beta]))

        with self.assertRaises(ValueError):
            model.predict(self.h1, 
                          self.h2, 
                          self.alpha, 
                          self.beta[:-2])

        with self.assertRaises(ValueError):
            tmp = self.beta.copy()
            tmp[2] = "the"
            model.predict(self.h1, 
                          self.h2, 
                          self.alpha, 
                          tmp)

        with self.assertRaises(AttributeError):
            model.predict(self.h1, 
                          self.h2, 
                          self.alpha, 
                          self.beta.tolist())

        with self.assertRaises(AttributeError):
            model.predict(self.h1, 
                          self.h2, 
                          self.alpha, 
                          1.3)


    def test_alpha_check(self):
        with self.assertRaises(ValueError):
            model.predict(self.h1,
                          self.h2,
                          np.array([self.alpha]),
                          self.beta)
            
        with self.assertRaises(ValueError):
            model.predict(self.h1,
                          self.h2,
                          "cat",
                          self.beta)

        with self.assertRaises(ValueError):
            model.predict(self.h1,
                          self.h2,
                          [self.alpha],
                          self.beta)


    def test_prediction_vals(self):
        """Test whether the outputs of model._predict are returned."""
        alpha = 1.5
        self.assertAlmostEqual(np.exp2(alpha)*(np.exp2(3) + np.exp2(2)),
                               model.predict(np.array([1,1]), 
                                             np.array([0,1]),
                                             alpha,
                                             np.array([1,2]))
                               )


# TODO
class TestSimulate(TestCase):
    pass


if __name__ == "__main__":
    main()

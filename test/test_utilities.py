from unittest import TestCase, main, mock
import io
import numpy as np
from afcn import utils


class TestIsNumericNparray(TestCase):
    def test_true_numeric(self):
        self.assertTrue(utils.is_numeric_nparray(np.array([1,2])))
        self.assertTrue(utils.is_numeric_nparray(np.array([1,False])))

    def test_False_numeric(self):
        self.assertFalse(utils.is_numeric_nparray(np.array([False])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",2])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",2.1])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",np.pi])))
        self.assertFalse(utils.is_numeric_nparray(np.array([None,2])))
        self.assertFalse(utils.is_numeric_nparray(np.array([None,False])))
        self.assertFalse(utils.is_numeric_nparray(np.array([True,False])))



class TestIsBiallelic(TestCase):
    def test_true_cases(self):
        self.assertTrue(utils.is_biallelic(np.array([1,1,1,1])))
        self.assertTrue(utils.is_biallelic(np.array([1.,1.,1.,1.])))
        self.assertTrue(utils.is_biallelic(np.array([0,0,0,0])))
        self.assertTrue(utils.is_biallelic(np.array([1,0,0,0])))
        self.assertTrue(utils.is_biallelic(np.array([1,0,0,1])))
        self.assertTrue(utils.is_biallelic([1,0,0,1]))
        self.assertTrue(utils.is_biallelic(np.array([1])))
        self.assertTrue(utils.is_biallelic(np.array([0])))
        self.assertTrue(utils.is_biallelic(0))
        self.assertTrue(utils.is_biallelic(1))
        self.assertTrue(utils.is_biallelic(np.nan))
        self.assertTrue(utils.is_biallelic(np.array([0,1,True])))
        self.assertTrue(utils.is_biallelic(np.array([0,1,1, np.nan])))
        self.assertTrue(utils.is_biallelic(np.array([np.nan,1, np.nan])))

    def test_false_cases(self):
        self.assertFalse(utils.is_biallelic([0, 1, 2]))
        self.assertFalse(utils.is_biallelic(np.array([0,1,2])))
        self.assertFalse(utils.is_biallelic(np.array([0,1,1.1])))
        self.assertFalse(utils.is_biallelic(2))
        self.assertFalse(utils.is_biallelic(None))
        self.assertFalse(utils.is_biallelic(True))
        self.assertFalse(utils.is_biallelic(False))
        self.assertFalse(utils.is_biallelic(np.array(["1", "0"])))

    def test_exception(self):
        with self.assertRaises(TypeError):
            utils.is_biallelic(np.array([0,None,2]))


class Test_is_int(TestCase):
    def test_correct(self):
        self.assertTrue(utils.is_int("-2323"))
        self.assertTrue(utils.is_int("3"))
        self.assertTrue(utils.is_int("+363464"))
        self.assertTrue(utils.is_int("3464"))
        self.assertTrue(utils.is_int("0"))

    def test_not_int(self):
        self.assertFalse(utils.is_int("423.32"))
        self.assertFalse(utils.is_int("-423.32"))
        self.assertFalse(utils.is_int("0.0"))
        self.assertFalse(utils.is_int("0."))
        self.assertFalse(utils.is_int("+0."))
        self.assertFalse(utils.is_int("+.5343"))
        self.assertFalse(utils.is_int("-.5343"))
        self.assertFalse(utils.is_int(".5343"))
        self.assertFalse(utils.is_int("a"))
        self.assertFalse(utils.is_int("3a"))
        self.assertFalse(utils.is_int("3.a"))
        self.assertFalse(utils.is_int("3-a"))
        self.assertFalse(utils.is_int("-a3"))
        self.assertFalse(utils.is_int("-\.3"))
        self.assertFalse(utils.is_int("-?3"))

    def test_exceptions(self):
        with self.assertRaises(TypeError):
            utils.is_int(5)

        with self.assertRaises(TypeError):
            utils.is_int(5.3434)

        with self.assertRaises(TypeError):
            utils.is_int(True)

        with self.assertRaises(TypeError):
            utils.is_int(["the", "cat"])

        with self.assertRaises(TypeError):
            utils.is_int(["1", "2"])

        with self.assertRaises(TypeError):
            utils.is_int(("1", "2"))


class Test_is_float(TestCase):
    def test_correct(self):
        self.assertTrue(utils.is_float("-23.23"))
        self.assertTrue(utils.is_float("3."))
        self.assertTrue(utils.is_float("+363464."))
        self.assertTrue(utils.is_float("3.464"))
        self.assertTrue(utils.is_float("0."))

    def test_not_float(self):
        self.assertFalse(utils.is_float(".42332"))
        self.assertFalse(utils.is_float("-42332"))
        self.assertFalse(utils.is_float(".0"))
        self.assertFalse(utils.is_float("0"))
        self.assertFalse(utils.is_float("+0"))
        self.assertFalse(utils.is_float("+5343"))
        self.assertFalse(utils.is_float("-5343"))
        self.assertFalse(utils.is_float(".5343"))
        self.assertFalse(utils.is_float("a"))
        self.assertFalse(utils.is_float("3a"))
        self.assertFalse(utils.is_float("3.a"))
        self.assertFalse(utils.is_float("3-a"))
        self.assertFalse(utils.is_float("-a3"))
        self.assertFalse(utils.is_float("-\.3"))
        self.assertFalse(utils.is_float("-?3"))

    def test_exceptions(self):
        with self.assertRaises(TypeError):
            utils.is_float(5)

        with self.assertRaises(TypeError):
            utils.is_float(5.3434)

        with self.assertRaises(TypeError):
            utils.is_float(True)

        with self.assertRaises(TypeError):
            utils.is_float(["the", "cat"])

        with self.assertRaises(TypeError):
            utils.is_float(["1", "2"])

        with self.assertRaises(TypeError):
            utils.is_float(("1", "2"))



if __name__ == "__main__":
    main()

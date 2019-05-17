from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from numpy.random import permutation
from numpy import array_split, concatenate
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np


class MushroomProblem(object):

    def __init__(self, data_file):
        self.data_frame = pd.read_csv(data_file)
        for k in self.data_frame.columns[1:]:
            self.data_frame[k], _ = pd.factorize(self.data_frame[k])
        categories = sorted(pd.Categorical(self.data_frame['class']).categories)
        self.classes = np.array(categories)
        self.features = self.data_frame.columns[self.data_frame.columns != 'class']
    @staticmethod
    def __factorize(data):
        y, _ = pd.factorize(pd.Categorical(data['class']), sort=True)
        return y


if __name__ == '__main__':

    f = r"C:\Users\Letissier\Desktop\Master\TFM\msd-tfm\Source\agaricus-lepiota.data"

    mp = MushroomProblem(f)
    print(mp.data_frame.to_string())


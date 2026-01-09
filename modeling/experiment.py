from datasets import Data
from models import ml_models
import pandas as pd

ds = Data(filename='example_csv')

test_set = []

count = 0
for x_train, x_test, y_train, y_test in ds.kFolds():
    models = ml_models(x_train, x_test, y_train, y_test)
    rsquared, rf = models.RandomForest()
    if count == 0:
        test_set = x_test.copy()
        test_set.loc[:, 'LogP Prediction'] = rf.round(2)
        test_set.loc[:, 'R^2 Value'] = rsquared.__round__(3)
        count += 1

test_set.to_csv("results.csv", index=True)
import pandas as pd
from sklearn.model_selection import KFold

kf = KFold(n_splits=5, shuffle=True, random_state=42)

class Data:

    def __init__(self, filename: str):
        self.dataset = pd.read_csv(filename + ".csv")
    

    def kFolds(self):
        self.X = self.dataset.drop(columns=['LogP', 'smiles'])
        self.Y = self.dataset['LogP']
        for train_index, test_index in kf.split(self.X):
            x_train, x_test = self.X.iloc[train_index], self.X.iloc[test_index]
            y_train, y_test = self.Y.iloc[train_index], self.Y.iloc[test_index]
            yield x_train, x_test, y_train, y_test



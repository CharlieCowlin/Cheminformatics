from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestRegressor
import xgboost
from xgboost import XGBRegressor

class ml_models:

    def __init__(self, xtrain, xtest, ytrain, ytest):
        self.xtrain = xtrain
        self.xtest = xtest
        self.ytrain = ytrain
        self.ytest = ytest

    def ridge(self):
        scaler = StandardScaler()
        ridge_regression = linear_model.Ridge(alpha=0.1, fit_intercept=True, tol=0.0001, random_state=42)
        x_trained_scaled = scaler.fit_transform(self.xtrain)
        x_test_scaled = scaler.transform(self.xtest)
        ridge_regression.fit(x_trained_scaled, self.ytrain)
        prediction = ridge_regression.predict(x_test_scaled)
        return r2_score(self.ytest, prediction)
    
    def RandomForest(self):
        randomForest = RandomForestRegressor(n_estimators=100, random_state=42, max_depth=5)
        randomForest.fit(self.xtrain, self.ytrain)
        prediction = randomForest.predict(self.xtest)
        return r2_score(self.ytest, prediction), prediction
    
    def XGBoost(self):
        xgboost = XGBRegressor(n_estimators=100, learning_rate=0.1, max_depth=5)
        xgboost.fit(self.xtrain, self.ytrain)
        prediction = xgboost.predict(self.xtest)
        return r2_score(self.ytest, prediction)

    
The modeling folder is used to train and test the datasets, create the machine learning algorithms, and evaluate the models.
datasets.py is used to split the test and train data to 80/20 and to fold/randomize the data 5 times to create more arbitrary data for the models.
models.py is used to craft three different machine learning models including a ridge regressor, random forest, and XGBoost regressor.
experiment.py is the test file to show users how the machine learning models could be implemented with the datasets, while also showing different evaluation functions that can be used to examine the model's accuracy.

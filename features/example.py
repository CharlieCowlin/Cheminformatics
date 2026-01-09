from featurize import create_features

example = create_features('clean_smiles.csv')

example.featurize_with_descriptors()
example.create_csv('example_csv')
The features folder contains files that generate descriptors (in descriptor_sets.py) and Morgan fingerprints (in fingerprints.py) from the molecular smiles.
descriptor_sets.py is used to import a file with clean smiles generated from the standardize.py file and generate different descriptors of those molecules including polarity, molecular weight, and hydrogen bonding.
fingerprints.py is used to generate Morgan fingerprints from a set of clean smiles from the standardize.py file.
featurize.py is used to import a file, generate descriptors, and then exporting the dataset as a csv file that will be used in the future.
example.py is an example of how to run the code with a small dataset (~500 molecules).

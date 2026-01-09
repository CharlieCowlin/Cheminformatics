import pandas as pd
from descriptor_sets import Physiochemical

file = []

class create_features:
    
    
    def __init__(self, csv: str):
        df = pd.read_csv(csv)
        self.file = df.dropna(subset=['clean molecules'])
        self.smiles = self.file['clean molecules']
    
    def featurize_with_descriptors(self):
        for smile in self.smiles:
            root = Physiochemical(smile)

            logp = root.logP()
            weight = root.weight()
            polarity = root.TPSA()
            h_bonding = root.hydrogen_bondingHBD()
            rbonds = root.rotatable_bonds()
            rings = root.num_rings()
            hba = root.hydrogen_bondingHBA()

            file.append({
                "smiles": smile,
                "LogP": round(logp, 2),
                "molecular weight": round(weight, 2),
                "TPSA": round(polarity, 2),
                "HBD": h_bonding,
                "HBA": hba,
                "RotB": rbonds,
                "Ring Count": rings
            })
    
    def create_csv(self, filename: str):
        dataframe = pd.DataFrame(file)
        dataframe.to_csv(filename + ".csv", index=False)

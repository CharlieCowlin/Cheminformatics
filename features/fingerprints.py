from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np

class fingerprint:

    def __init__(self, smile):
        self.mol = Chem.MolFromSmiles(smile)


    def MorganFingerprint(self):
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol=self.mol, radius=3, nBits=2048)
        return fp
    
    def fingerprintBits(self, fp):
        array = np.array(fp)
        shape = np.shape(array)
        return shape



from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf
from rdkit.Chem import rdMolDescriptors

class Physiochemical:

    def __init__(self, smile):
        self.mol = Chem.MolFromSmiles(smile)
        if self.mol is None:
            self.mol = Chem.MolFromSmiles("C")


    def logP(self):
        logp = Crippen._pyMolLogP(self.mol)
        return logp
    
    def weight(self):
        molweight = Descriptors.HeavyAtomMolWt(self.mol)
        return molweight
    
    def polarity(self):
        pol = MolSurf._pyTPSA(self.mol)
        return pol
    
    def TPSA(self):
        top_pol = rdMolDescriptors.CalcTPSA(self.mol)
        return top_pol
    
    def hydrogen_bondingHBD(self):
        hydrogen = rdMolDescriptors.CalcNumHBD(self.mol)
        return hydrogen
    
    def hydrogen_bondingHBA(self):
        hydrogen2 = rdMolDescriptors.CalcNumHBA(self.mol)
        return hydrogen2
    
    def rotatable_bonds(self):
        rbonds = rdMolDescriptors.CalcNumRotatableBonds(self.mol)
        return rbonds
    
    def num_rings(self):
        ring_count = rdMolDescriptors.CalcNumRings(self.mol)
        return ring_count
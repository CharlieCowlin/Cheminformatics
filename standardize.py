import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Load dataset
csv = pd.read_csv('clean.csv') # Ensure correct filename extension

enumerator = rdMolStandardize.TautomerEnumerator()
uncharger = rdMolStandardize.Uncharger()

results = {
    "id": [],
    "original smiles": [],
    "status": [],
    "clean molecules": [],
    "salt removed": [],
    "reason": []
}

# Set to track unique SMILES
seen_smiles = set()

for i in range(len(csv)):
    row = csv.iloc[i]
    smile = str(row['smiles'])
    mol = Chem.MolFromSmiles(smile)

    if mol is None:
        results["id"].append(row['id'])
        results["original smiles"].append(smile)
        results["status"].append("invalid")
        results["clean molecules"].append("NaN")
        results["salt removed"].append("none")
        results["reason"].append("parsing failed")
        continue

    frags = Chem.GetMolFrags(mol, asMols=True)
    organic_frags = [f for f in frags if any(a.GetSymbol() == 'C' for a in f.GetAtoms())]
    inorganic_frags = [f for f in frags if not any(a.GetSymbol() == 'C' for a in f.GetAtoms())]

    if len(organic_frags) == 0:
        results["id"].append(row['id'])
        results["original smiles"].append(smile)
        results["status"].append("inorganic")
        results["clean molecules"].append("NaN")
        results["salt removed"].append(", ".join([Chem.MolToSmiles(x) for x in inorganic_frags]) if inorganic_frags else "none")
        results["reason"].append("no organic fragment")
        continue

    elif len(organic_frags) > 1:
        results["id"].append(row['id'])
        results["original smiles"].append(smile)
        results["status"].append("mixture")
        results["clean molecules"].append("NaN")
        results["salt removed"].append("none")
        results["reason"].append("multiple organic fragments")
        continue

    else:
        parent_mol = organic_frags[0]
        salt_str = ", ".join([Chem.MolToSmiles(x) for x in inorganic_frags]) if inorganic_frags else "none"
        
        try:
            canonical_mol = enumerator.Canonicalize(parent_mol)
            neutral_mol = uncharger.uncharge(canonical_mol)
            final_smile = Chem.MolToSmiles(neutral_mol)
            
            # --- DUPLICATE CHECK ---
            if final_smile in seen_smiles:
                # Option: Skip entirely, or label as duplicate. 
                # Here we label it so you don't lose the ID mapping, but mark clean molecules as NaN
                status = "duplicate"
                reason_str = f"duplicate of previous entry: {final_smile}"
                clean_val = "NaN"
            else:
                seen_smiles.add(final_smile)
                status = "salt_removed" if len(inorganic_frags) > 0 else "ok"
                clean_val = final_smile
                reason_str = "salt removed" if status == "salt_removed" else ("standardized" if final_smile != smile else "no changes")

            results["id"].append(row['id'])
            results["original smiles"].append(smile)
            results["status"].append(status)
            results["clean molecules"].append(clean_val)
            results["salt removed"].append(salt_str)
            results["reason"].append(reason_str)

        except Exception:
            results["id"].append(row['id'])
            results["original smiles"].append(smile)
            results["status"].append("invalid")
            results["clean molecules"].append("NaN")
            results["salt removed"].append("error")
            results["reason"].append("standardization failed")

output_df = pd.DataFrame(results)

# Final cleanup: If you want to strictly remove rows that are duplicates from the CSV:
# output_df = output_df[output_df["status"] != "duplicate"]

output_df.to_csv("clean_smiles.csv", index=False)
print(f"Processing complete. Unique molecules found: {len(seen_smiles)}")

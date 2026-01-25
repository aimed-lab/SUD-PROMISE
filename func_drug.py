"""
Drug Database Utility Functions
Handles searching and extracting drug information from the database
"""

import pandas as pd
import ast
from typing import Optional, List, Tuple


def find_drug_in_database(drug_name: str, drugs_df: pd.DataFrame) -> Tuple[Optional[str], Optional[List[str]]]:
    """
    Extract SMILES and targets for a drug by EXACT name match first, then partial
    
    Args:
        drug_name: Name of the drug to search for
        drugs_df: DataFrame containing drug information with columns:
                 - DrugName: name of the drug
                 - DrugSmile: SMILES structure
                 - DrugTarget: protein targets (as string representation of list)
    
    Returns:
        Tuple of (smiles, targets):
        - smiles: SMILES string representation of molecular structure, or None if not found
        - targets: List of protein target names, or empty list if none found
        
    Examples:
        >>> smiles, targets = find_drug_in_database("Naltrexone", drugs_df)
        >>> print(f"SMILES: {smiles}")
        >>> print(f"Targets: {targets}")
    """
    if drugs_df is None:
        return None, None
    
    # Try EXACT match first (case-insensitive)
    mask = drugs_df['DrugName'].str.lower() == drug_name.lower()
    matches = drugs_df[mask]
    
    # If no exact match, try partial match
    if len(matches) == 0:
        mask = drugs_df['DrugName'].str.lower().str.contains(drug_name.lower(), na=False)
        matches = drugs_df[mask]
    
    if len(matches) == 0:
        return None, None
    
    # Get the first match
    row = matches.iloc[0]
    
    # Extract SMILES - try multiple possible column names
    smiles = None
    smiles_columns = ['DrugSmile', 'SMILES', 'Smiles', 'smiles', 'DrugSMILES', 'Structure']
    for col in smiles_columns:
        if col in row.index and pd.notna(row[col]) and str(row[col]).strip():
            smiles = str(row[col]).strip()
            break
    
    # Extract Targets - try multiple possible column names
    targets = []
    target_columns = ['DrugTarget', 'Targets', 'Target', 'targets', 'ProteinTargets']
    for col in target_columns:
        if col in row.index and pd.notna(row[col]) and str(row[col]).strip():
            raw_targets = str(row[col]).strip()
            
            # Parse if it's a string representation of a list
            if raw_targets.startswith('[') and raw_targets.endswith(']'):
                try:
                    targets = ast.literal_eval(raw_targets)
                except:
                    targets = [raw_targets]
            elif ',' in raw_targets:
                # Handle comma-separated targets
                targets = [t.strip() for t in raw_targets.split(',')]
            else:
                # Single target
                targets = [raw_targets]
            break
    
    return smiles, targets if targets else []


def search_drugs_by_target(target_name: str, drugs_df: pd.DataFrame, 
                           exact_match: bool = False) -> pd.DataFrame:
    """
    Search for drugs that target a specific protein
    
    Args:
        target_name: Name of the protein target to search for
        drugs_df: DataFrame containing drug information
        exact_match: If True, only return exact matches; if False, use partial matching
    
    Returns:
        DataFrame containing all drugs that target the specified protein
    """
    if drugs_df is None:
        return pd.DataFrame()
    
    def has_target(targets_str):
        if pd.isna(targets_str):
            return False
        
        # Parse targets
        targets = []
        if str(targets_str).startswith('['):
            try:
                targets = ast.literal_eval(str(targets_str))
            except:
                targets = [str(targets_str)]
        else:
            targets = [str(targets_str)]
        
        # Check if target matches
        for target in targets:
            if exact_match:
                if target.lower() == target_name.lower():
                    return True
            else:
                if target_name.lower() in target.lower():
                    return True
        return False
    
    mask = drugs_df['DrugTarget'].apply(has_target)
    return drugs_df[mask]


def get_drug_info(drug_name: str, drugs_df: pd.DataFrame) -> dict:
    """
    Get complete information about a drug
    
    Args:
        drug_name: Name of the drug
        drugs_df: DataFrame containing drug information
    
    Returns:
        Dictionary with drug information including:
        - name: Drug name
        - smiles: SMILES structure
        - targets: List of protein targets
        - found: Boolean indicating if drug was found
    """
    smiles, targets = find_drug_in_database(drug_name, drugs_df)
    
    return {
        'name': drug_name,
        'smiles': smiles,
        'targets': targets if targets else [],
        'found': smiles is not None
    }


def validate_smiles(smiles: str) -> bool:
    """
    Validate if a SMILES string is valid (requires RDKit)
    
    Args:
        smiles: SMILES string to validate
    
    Returns:
        True if valid, False otherwise
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


# Example usage and testing
if __name__ == "__main__":
    # Load example data
    try:
        drugs_df = pd.read_csv("data1/drugsInfo.csv")
        print(f" Loaded {len(drugs_df)} drugs from database\n")
        
        # Test 1: Find a specific drug
        print("="*60)
        print("TEST 1: Find Naltrexone")
        print("="*60)
        smiles, targets = find_drug_in_database("Naltrexone", drugs_df)
        if smiles:
            print(f" Found!")
            print(f"   SMILES: {smiles[:60]}...")
            print(f"   Targets: {targets[:3]}{'...' if len(targets) > 3 else ''}")
        else:
            print(f"❌ Not found")
        
        # Test 2: Get complete info
        print("\n" + "="*60)
        print("TEST 2: Get complete drug info")
        print("="*60)
        info = get_drug_info("Buprenorphine", drugs_df)
        print(f"Drug: {info['name']}")
        print(f"Found: {info['found']}")
        print(f"Targets: {len(info['targets'])} targets")
        
        # Test 3: Search by target
        print("\n" + "="*60)
        print("TEST 3: Search drugs targeting opioid receptors")
        print("="*60)
        results = search_drugs_by_target("opioid", drugs_df, exact_match=False)
        print(f"Found {len(results)} drugs")
        if len(results) > 0:
            print("Examples:")
            for i, (idx, row) in enumerate(results.head(3).iterrows(), 1):
                print(f"   {i}. {row['DrugName']}")
        
    except FileNotFoundError:
        print("❌ Error: data1/drugsInfo.csv not found")
        print("   Make sure the database file exists in the correct location")
"""
ML/DL Prediction Utilities
Functions for making predictions with trained ML/DL models

VERSION: 3.5.0 - Single initialization function for easy integration
- Added initialize_ml_system() function
- Consolidates all loading logic
- Returns complete system state
- Uses const_config for default paths
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, List, Optional
import joblib

# Import configuration constants
from const_config import MODEL_DIR, DRUGS_FILE, DISEASES_FILE

def initialize_ml_system(model_dir: str = MODEL_DIR,
                         drugs_file: str = DRUGS_FILE,
                         diseases_file: str = DISEASES_FILE) -> Dict:
    """
    Initialize the complete ML/DL prediction system
    
    This single function replaces all the initialization code in the main script.
    It loads CSV datasets and ML/DL models, handling all errors gracefully.
    
    Args:
        model_dir: Path to directory containing model files (default from const_config.MODEL_DIR)
        drugs_file: Path to drugs CSV file (default from const_config.DRUGS_FILE)
        diseases_file: Path to diseases CSV file (default from const_config.DISEASES_FILE)
        
    Returns:
        Dictionary containing:
        - 'models_available': bool - Whether ML/DL models loaded successfully
        - 'ml_components': dict - Loaded models and preprocessors
        - 'drugs_df': DataFrame or None - Drugs database
        - 'diseases_df': DataFrame or None - Diseases database
    """
    
    result = {
        'models_available': False,
        'ml_components': {},
        'drugs_df': None,
        'diseases_df': None
    }
    
    try:
        # Suppress RDKit warnings
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
        # Suppress TensorFlow warnings
        try:
            import tensorflow as tf
            tf.get_logger().setLevel('ERROR')
        except:
            pass
        
        # Load CSV datasets
        print(" Loading drug and disease datasets...")
        try:
            result['drugs_df'] = pd.read_csv(drugs_file)
            result['diseases_df'] = pd.read_csv(diseases_file)
            print(f" Loaded {len(result['drugs_df'])} drugs and {len(result['diseases_df'])} diseases")
        except Exception as e:
            print(f"  Could not load CSV files: {e}")
            print("   Continuing without drug/disease databases")
        
        # Load ML/DL models
        model_path = Path(model_dir)
        
        if not model_path.exists():
            print(f"  ML/DL models directory not found: {model_path}")
            print("   Using synthetic scores only.")
            return result
        
        print(f"\n Loading ML/DL models from: {model_path}")
        success, ml_components = load_ml_models(model_path)
        
        if success:
            result['models_available'] = True
            result['ml_components'] = ml_components
            print("\n ML/DL Models loaded successfully!")
            print(f"   Available diseases in model: {len(ml_components['le_disease'].classes_)}")
            print(f"   Available targets in model: {len(ml_components['mlb'].classes_)}")
        else:
            print("\n  ML/DL models not found. Using synthetic scores only.")
        
    except Exception as e:
        print(f"\n  Error during ML system initialization: {e}")
        print("   Using synthetic scores only.")
        
        # Try to at least load CSV files if models failed
        if result['drugs_df'] is None:
            try:
                result['drugs_df'] = pd.read_csv(drugs_file)
                result['diseases_df'] = pd.read_csv(diseases_file)
                print(f" Loaded {len(result['drugs_df'])} drugs and {len(result['diseases_df'])} diseases (without ML models)")
            except:
                print(" Could not load CSV files either. Using fully synthetic data.")
    
    return result


def load_ml_models(model_dir: Path) -> Tuple[bool, Dict]:
    """
    Load all ML/DL models and preprocessors with detailed error reporting
    
    Args:
        model_dir: Path to directory containing model files
        
    Returns:
        Tuple of (success: bool, ml_components: dict)
    """
    try:
        from tensorflow import keras
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        ml_components = {}
        
        print(f"\n Loading models from: {model_dir.absolute()}")
        
        # Define all model files
        model_files = {
            'lr_model': ('logistic_regression.pkl', 'joblib'),
            'rf_model': ('random_forest.pkl', 'joblib'),
            'dnn_model': ('mm_dnn_model.keras', 'keras'),
            'mlb': ('target_binarizer.pkl', 'joblib'),
            'le_disease': ('disease_encoder.pkl', 'joblib'),
            'disease_ohe_df': ('disease_ohe_df.pkl', 'pandas'),  #  Special handling
            'scaler': ('scaler.pkl', 'joblib'),
        }
        
        # Load each file with individual error handling
        for component_name, (filename, file_type) in model_files.items():
            file_path = model_dir / filename
            
            if not file_path.exists():
                print(f"     {filename}: File not found, skipping...")
                continue
            
            try:
                if file_type == 'joblib':
                    #  FIXED: Use joblib for sklearn models
                    ml_components[component_name] = joblib.load(file_path)
                    print(f"    {filename}: Loaded successfully")
                    
                elif file_type == 'pandas':
                    #  FIXED: Use pd.read_pickle for pandas DataFrames
                    ml_components[component_name] = pd.read_pickle(file_path)
                    print(f"    {filename}: Loaded successfully")
                    
                elif file_type == 'keras':
                    ml_components[component_name] = keras.models.load_model(str(file_path))
                    print(f"    {filename}: Loaded successfully")
                    
            except Exception as e:
                print(f"    {filename}: Failed to load")
                print(f"      Error: {e}")
                # Continue loading other files
                continue
        
        # Check all required components are loaded
        required = ['lr_model', 'rf_model', 'dnn_model', 'mlb', 'le_disease', 'scaler']
        missing = [k for k in required if k not in ml_components]
        
        if missing:
            print(f"\n  Missing required components: {missing}")
            return False, {}
        else:
            print(f"\n All required components loaded successfully!")
            return True, ml_components
            
    except ImportError as e:
        print(f" Import error: {e}")
        print("   Make sure tensorflow and rdkit are installed:")
        print("   pip install tensorflow rdkit --break-system-packages")
        return False, {}
    except Exception as e:
        print(f" Unexpected error loading ML models: {e}")
        import traceback
        traceback.print_exc()
        return False, {}


def prepare_drug_features(drug_smiles: str, drug_targets: List[str], mlb) -> Optional[np.ndarray]:
    """
    Prepare drug features from SMILES and targets
    
    Args:
        drug_smiles: SMILES string of the drug
        drug_targets: List of protein target names
        mlb: MultiLabelBinarizer for targets
        
    Returns:
        Feature vector or None if error
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # Generate Morgan fingerprint
        mol = Chem.MolFromSmiles(drug_smiles)
        if mol is None:
            return None
            
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        fp_array = np.array(fp)
        
        # Encode targets
        if drug_targets:
            target_encoded = mlb.transform([drug_targets])
        else:
            target_encoded = mlb.transform([[]])
        
        # Concatenate features
        features = np.concatenate([fp_array, target_encoded[0]])
        
        return features
        
    except Exception as e:
        print(f"Error preparing drug features: {e}")
        return None


def predict_with_ml_models(drug_smiles: str, 
                           drug_targets: List[str], 
                           disease_id: str, 
                           ml_components: Dict) -> Tuple[Optional[Dict], str, str]:
    """
    Make predictions using ML/DL models
    
    Args:
        drug_smiles: SMILES string of drug
        drug_targets: List of protein targets
        disease_id: Disease ID (e.g., 'MESH:D009293')
        ml_components: Dictionary containing loaded models and preprocessors
        
    Returns:
        Tuple of (results_dict, message, score_type)
        - results_dict: {'LR': score, 'RF': score, 'DNN': score, 'Ensemble': score} or None
        - message: Status message
        - score_type: 'Real' or 'Synthetic'
    """
    try:
        # Validate inputs
        if not drug_smiles or not disease_id:
            return None, "Missing drug SMILES or disease ID", "Synthetic"
        
        if not ml_components:
            return None, "ML models not available", "Synthetic"
        
        # Convert disease classes to list and check membership
        disease_classes = ml_components['le_disease'].classes_
        if isinstance(disease_classes, pd.Series):
            disease_classes = disease_classes.tolist()
        elif isinstance(disease_classes, np.ndarray):
            disease_classes = disease_classes.tolist()
        else:
            disease_classes = list(disease_classes)
        
        # Check membership in Python list
        if disease_id not in disease_classes:
            return None, f"Disease {disease_id} not in training data", "Synthetic"
        
        # Prepare drug features
        drug_features = prepare_drug_features(drug_smiles, drug_targets, ml_components['mlb'])
        if drug_features is None:
            return None, "Failed to generate drug features", "Synthetic"
        
        # Encode disease (for DNN input)
        disease_encoded = ml_components['le_disease'].transform([disease_id])
        
        #  FIXED: Get disease OHE and drop 'DiseaseID' column (matching working code)
        if 'disease_ohe_df' in ml_components:
            disease_ohe_df = ml_components['disease_ohe_df']
            
            # Use the exact approach from sud_promise_uab_theme.py line 259
            disease_matches = disease_ohe_df[disease_ohe_df['DiseaseID'] == disease_id]
            if len(disease_matches) == 0:
                return None, f"Disease ID '{disease_id}' not found in disease_ohe_df", "Synthetic"
            
            #  CRITICAL: Drop DiseaseID column before getting values
            disease_ohe = disease_matches.drop('DiseaseID', axis=1).values
        else:
            return None, "disease_ohe_df not loaded", "Synthetic"
        
        #  FIXED: Prepare full feature vector for ML models (LR, RF)
        # Drug features (1024 + target_dim) + disease OHE (should give 4615 total)
        X_combined = np.concatenate([
            drug_features.reshape(1, -1),
            disease_ohe  #  This is the full one-hot vector WITHOUT DiseaseID column
        ], axis=1)
        
        # Scale features for LR (RF might not need it, check your training)
        if 'scaler' in ml_components:
            X_scaled = ml_components['scaler'].transform(X_combined)
        else:
            X_scaled = X_combined
        
        # Make predictions with each model
        results = {}
        
        # Logistic Regression (uses scaled features)
        if 'lr_model' in ml_components:
            try:
                lr_pred = ml_components['lr_model'].predict_proba(X_scaled)[0][1]
                results['Logistic Regression'] = float(lr_pred)
            except Exception as e:
                print(f"   LR prediction failed: {e}")
        
        # Random Forest (check if it needs scaled or unscaled)
        # Based on working code, RF uses unscaled ml_input
        if 'rf_model' in ml_components:
            try:
                rf_pred = ml_components['rf_model'].predict_proba(X_combined)[0][1]
                results['Random Forest'] = float(rf_pred)
            except Exception as e:
                print(f"   RF prediction failed: {e}")
        
        #  FIXED: Deep Neural Network uses different input format
        # DNN takes [drug_features, disease_index] as separate inputs
        if 'dnn_model' in ml_components:
            try:
                disease_idx_arr = np.array([disease_encoded[0]], dtype=np.int32)
                dnn_pred = ml_components['dnn_model'].predict(
                    [drug_features.reshape(1, -1), disease_idx_arr], 
                    verbose=0
                )[0][0]
                results['MM-DNN'] = float(dnn_pred)
            except Exception as e:
                print(f"   DNN prediction failed: {e}")
        
        # Calculate ensemble (average of all models)
        if results:
            ensemble_score = np.mean(list(results.values()))
            results['Ensemble'] = float(ensemble_score)
            
            return results, " Prediction successful", "Real"
        else:
            return None, "No models available for prediction", "Synthetic"
            
    except Exception as e:
        import traceback
        error_msg = f"Prediction error: {str(e)}"
        # Uncomment for debugging:
        print(f"\n🔍 DETAILED ERROR:")
        print(traceback.format_exc())
        return None, error_msg, "Synthetic"


def get_ensemble_prediction(results: Dict) -> float:
    """Get ensemble prediction from model results"""
    if not results:
        return 0.5
    
    if 'Ensemble' in results:
        return results['Ensemble']
    
    scores = [v for k, v in results.items() if k != 'Ensemble']
    if scores:
        return np.mean(scores)
    else:
        return 0.5


def interpret_prediction_score(score: float) -> Tuple[str, str, str]:
    """Interpret prediction score"""
    if score >= 0.7:
        return "HIGH", "🟢", "Strong therapeutic potential"
    elif score >= 0.5:
        return "MODERATE", "🟡", "Mixed evidence, further investigation needed"
    else:
        return "LOW", "🔴", "Limited evidence for this indication"


def batch_predict(drug_disease_pairs: List[Tuple[str, List[str], str]], 
                  ml_components: Dict) -> List[Dict]:
    """Make predictions for multiple drug-disease pairs"""
    results = []
    
    for drug_smiles, drug_targets, disease_id in drug_disease_pairs:
        pred_results, message, score_type = predict_with_ml_models(
            drug_smiles, drug_targets, disease_id, ml_components
        )
        
        results.append({
            'drug_smiles': drug_smiles,
            'disease_id': disease_id,
            'predictions': pred_results,
            'message': message,
            'score_type': score_type
        })
    
    return results


if __name__ == "__main__":
    print("="*70)
    print("ML PREDICTION UTILITIES - TESTING")
    print("="*70)
    
    #  NEW: Single function call uses defaults from const_config
    # No need to specify paths - they're already in const_config!
    system = initialize_ml_system()
    
    print("\n" + "="*70)
    print("SYSTEM STATUS")
    print("="*70)
    print(f"Models Available: {system['models_available']}")
    print(f"Drugs Database: {' Loaded' if system['drugs_df'] is not None else ' Not loaded'}")
    print(f"Diseases Database: {' Loaded' if system['diseases_df'] is not None else ' Not loaded'}")
    
    if system['models_available']:
        print("\n ALL TESTS PASSED")
        
        # Test feature dimensions
        print("\n🔍 Testing feature dimensions...")
        test_smiles = "CC(C)NCC(COc1ccccc1)O"  # Example SMILES
        test_targets = []
        test_disease = list(system['ml_components']['le_disease'].classes_)[0]
        
        results, msg, score_type = predict_with_ml_models(
            test_smiles, test_targets, test_disease, system['ml_components']
        )
        
        if results:
            print(f" Prediction successful!")
            print(f"   Score Type: {score_type}")
            print(f"   Results: {results}")
        else:
            print(f"  Prediction failed: {msg}")
    else:
        print("\n  MODEL LOADING FAILED - System will use synthetic scores")
    
    print("="*70)
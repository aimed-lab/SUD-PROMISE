"""
Data Generator Module
Extracted from sud_promise_uab_theme.py
Handles disease ID mapping and synthetic data generation with real ML/DL scores
"""

import random
from datetime import datetime, timedelta
from dataclasses import dataclass, field
from typing import List

# Import dependencies
from const_ui import (
    UAB_GREEN, UAB_DARK_GREEN, UAB_LIGHT_GREEN, UAB_ACCENT_TEAL,
    STAGE_MAPPING
)
from const_data import (
    SUD_CATEGORY_DESCRIPTIONS, DISEASE_SEARCH_CONFIG,
    DRUG_TEMPLATES, PROJECT_TEMPLATES
)
from func_drug import find_drug_in_database

# ========================================
# DATA MODELS
# ========================================

@dataclass
class Project:
    """Research project providing evidence"""
    id: str
    name: str
    project_type: str
    added_date: datetime
    sample_size: int
    impact_score: float
    status: str
    summary: str

@dataclass
class DrugCandidate:
    """Drug being evaluated for repositioning"""
    id: str
    drug_name: str
    current_indication: str
    target_sud_subtype: str
    mechanism: str
    stage: str
    evidence_score: float
    baseline_score: float
    smiles: str
    attached_projects: List[Project]
    last_updated: datetime
    cohort_count: int
    has_market_analysis: bool
    has_validation_plan: bool
    team_members: int
    data_produced: int
    publications: int
    tools_used: int
    data_governance: int
    training_participated: int
    stage_entry_date: datetime = None
    stage_history: List[tuple] = field(default_factory=list)
    score_type: str = "Synthetic"
    model_scores: dict = field(default_factory=dict)
    protein_targets: List[str] = field(default_factory=list)
    disease_id: str = ""

@dataclass
class SUDCategory:
    """SUD disease category"""
    name: str
    color: str
    hex_color: str
    icon: str
    candidate_count: int
    description: str
    disease_id: str = ""

# ========================================
# DISEASE ID MAPPING
# ========================================

def setup_disease_mapping(diseases_df):
    """Setup disease ID mapping from database"""
    DISEASE_ID_MAPPING = {}
    
    if diseases_df is not None:
        print("\n🔍 Searching for best disease name matches in database...")
        
        for category, config in DISEASE_SEARCH_CONFIG.items():
            print(f"\n Searching for: {category}")
            found = False
            
            for search_term in config["search_terms"]:
                # Try exact match first
                exact_mask = diseases_df['DiseaseName'].str.lower() == search_term.lower()
                exact_matches = diseases_df[exact_mask]
                
                if len(exact_matches) > 0:
                    disease_id = exact_matches.iloc[0]['DiseaseID']
                    disease_name = exact_matches.iloc[0]['DiseaseName']
                    DISEASE_ID_MAPPING[category] = disease_id
                    print(f"    EXACT MATCH: '{disease_name}'")
                    print(f"      → {disease_id}")
                    found = True
                    break
                
                # Try partial match if exact fails
                partial_mask = diseases_df['DiseaseName'].str.lower().str.contains(search_term.lower(), na=False)
                partial_matches = diseases_df[partial_mask]
                
                if len(partial_matches) > 0:
                    disease_id = partial_matches.iloc[0]['DiseaseID']
                    disease_name = partial_matches.iloc[0]['DiseaseName']
                    DISEASE_ID_MAPPING[category] = disease_id
                    print(f"    Partial match for '{search_term}': '{disease_name}'")
                    print(f"      → {disease_id}")
                    found = True
                    break
            
            if not found:
                print(f"    No matches found for {category}")
        
        print(f"\n{'='*80}")
        print(f"📋 FINAL DISEASE MAPPING:")
        print(f"{'='*80}")
        for category, disease_id in DISEASE_ID_MAPPING.items():
            disease_row = diseases_df[diseases_df['DiseaseID'] == disease_id]
            if len(disease_row) > 0:
                disease_name = disease_row.iloc[0]['DiseaseName']
                print(f"{category}:")
                print(f"   → {disease_id} ({disease_name})")
        
        print(f"\n Successfully mapped {len(DISEASE_ID_MAPPING)}/6 SUD categories\n")
    
    return DISEASE_ID_MAPPING

# ========================================
# SYNTHETIC DATA GENERATION WITH REAL ML/DL SCORES
# ========================================

def generate_synthetic_data(DISEASE_ID_MAPPING, MODELS_AVAILABLE, ml_components, drugs_df):
    """Generate realistic synthetic SUD repositioning data with REAL ML/DL scores where possible"""
    
    # Import here to avoid circular dependency
    from func_models import predict_with_ml_models
    
    print("\n" + "="*70)
    print("🔬 GENERATING CANDIDATE DATA WITH ML/DL PREDICTIONS")
    print("="*70)
    
    # Use corrected category names
    sud_categories = [
        SUDCategory("Opioid-Related Disorders", "", UAB_GREEN, "", 0, 
                   SUD_CATEGORY_DESCRIPTIONS["Opioid-Related Disorders"],
                   disease_id=DISEASE_ID_MAPPING.get("Opioid-Related Disorders", "")),
        SUDCategory("Alcohol Use Disorder", "", UAB_ACCENT_TEAL, "", 0,
                   SUD_CATEGORY_DESCRIPTIONS["Alcohol Use Disorder"],
                   disease_id=DISEASE_ID_MAPPING.get("Alcohol Use Disorder", "")),
        SUDCategory("Stimulant Use Disorder", "", UAB_LIGHT_GREEN, "", 0,
                   SUD_CATEGORY_DESCRIPTIONS["Stimulant Use Disorder"],
                   disease_id=DISEASE_ID_MAPPING.get("Stimulant Use Disorder", "")),
        SUDCategory("Cannabis Use Disorder", "", UAB_DARK_GREEN, "", 0,
                   SUD_CATEGORY_DESCRIPTIONS["Cannabis Use Disorder"],
                   disease_id=DISEASE_ID_MAPPING.get("Cannabis Use Disorder", "")),
        SUDCategory("Sedative/Hypnotic Disorder", "", "#4A8B7A", "", 0,
                   SUD_CATEGORY_DESCRIPTIONS["Sedative/Hypnotic Disorder"],
                   disease_id=DISEASE_ID_MAPPING.get("Sedative/Hypnotic Disorder", "")),
        SUDCategory("Nicotine Use Disorder", "", "#2C7A64", "", 0,
                   SUD_CATEGORY_DESCRIPTIONS["Nicotine Use Disorder"],
                   disease_id=DISEASE_ID_MAPPING.get("Nicotine Use Disorder", "")),
    ]
    
    candidates = []
    candidate_id = 1
    
    total_real_scores = 0
    total_synthetic_scores = 0
    
    for category in sud_categories:
        if category.name not in DRUG_TEMPLATES:
            continue
        
        category_real = 0
        category_synthetic = 0
        
        print(f"\n🔬 Processing {category.name}...")
        if category.disease_id:
            print(f"   Disease ID: {category.disease_id}")
        else:
            print(f"     No disease ID mapped - using synthetic scores only")
        
        for drug_info in DRUG_TEMPLATES[category.name]:
            drug_name, current_use, mechanism, stage = drug_info
            
            # Try to find drug in database using corrected function
            drug_smiles, drug_targets = find_drug_in_database(drug_name, drugs_df)
            
            if not drug_smiles:
                # Generate random SMILES if not found
                drug_smiles = f"C{'C' * random.randint(5, 15)}N"
                drug_targets = []
                print(f"     {drug_name}: Not in database, using synthetic SMILES")
            else:
                print(f"    {drug_name}: Found in database")
            
            # Generate stage history
            stage_num = int(stage[1])
            stage_history = []
            
            days_in_pipeline = random.randint(360, 900)
            start_date_candidate = datetime.now() - timedelta(days=days_in_pipeline)
            
            current_date = start_date_candidate
            for s in range(stage_num + 1):
                stage_name = f"S{s}"
                days_in_stage = random.randint(60, 120)
                stage_entry = current_date
                stage_history.append((stage_name, stage_entry))
                current_date = current_date + timedelta(days=days_in_stage)
            
            stage_entry_date = stage_history[-1][1] if stage_history else datetime.now()
            
            score_min, score_max = STAGE_MAPPING[stage]
            baseline = random.uniform(0.40, 0.55)
            
            # Try to get REAL ML/DL scores
            score_type = "Synthetic"
            model_scores = {}
            evidence_score = baseline
            baseline_score = baseline  # Will be updated for Real predictions
            
            if category.disease_id and drug_smiles and MODELS_AVAILABLE:
                print(f"      Predicting {drug_name} → {category.disease_id}...", end=" ")
                ml_results, message, ml_score_type = predict_with_ml_models(
                    drug_smiles, drug_targets, category.disease_id, ml_components
                )
                
                if ml_results is not None and ml_score_type == "Real":
                    score_type = "Real"
                    model_scores = ml_results
                    
                    # FIXED: Use ensemble as BASELINE
                    baseline_score = ml_results.get('Ensemble', baseline)
                    
                    # Calculate target evidence score
                    num_projects = random.randint(3, 5)
                    total_impact_needed = random.uniform(0.15, 0.35)  # Total impact from projects
                    evidence_score = baseline_score + total_impact_needed
                    evidence_score = max(0.20, min(0.95, evidence_score))
                    
                    # Recalculate actual impact needed after clamping
                    total_impact_needed = evidence_score - baseline_score
                    
                    print(f" Real ML/DL: {baseline_score:.3f} → {evidence_score:.3f}")
                    category_real += 1
                    total_real_scores += 1
                else:
                    # Generate synthetic score with projects
                    baseline_score = baseline
                    num_projects = random.randint(3, 5)
                    total_impact_needed = random.uniform(0.20, 0.40)
                    evidence_score = baseline_score + total_impact_needed
                    evidence_score = max(0.20, min(0.95, evidence_score))
                    total_impact_needed = evidence_score - baseline_score
                    
                    print(f"  Synthetic score: {evidence_score:.3f} ({message})")
                    category_synthetic += 1
                    total_synthetic_scores += 1
            else:
                # Generate synthetic score
                baseline_score = baseline
                num_projects = random.randint(3, 5)
                total_impact_needed = random.uniform(0.20, 0.40)
                evidence_score = baseline_score + total_impact_needed
                evidence_score = max(0.20, min(0.95, evidence_score))
                total_impact_needed = evidence_score - baseline_score
                
                print(f"      {drug_name}: Synthetic score: {evidence_score:.3f}")
                category_synthetic += 1
                total_synthetic_scores += 1
            
            # NOW Generate projects with impacts that sum EXACTLY to total_impact_needed
            projects = []
            
            # Distribute the total impact across projects
            remaining_impact = total_impact_needed
            for i in range(num_projects):
                template = random.choice(PROJECT_TEMPLATES)
                days_ago = random.randint(30 + i*50, days_in_pipeline - (num_projects - i)*30)
                
                if i == num_projects - 1:
                    # Last project gets exactly the remaining impact
                    impact = remaining_impact
                else:
                    # Random portion of remaining, but ensure we don't exhaust it
                    max_this_impact = remaining_impact * 0.6  # Use at most 60% of remaining
                    impact = random.uniform(-0.05, max_this_impact)
                    remaining_impact -= impact
                
                project = Project(
                    id=f"proj_{candidate_id}_{i}",
                    name=template["name"],
                    project_type=template["type"],
                    added_date=datetime.now() - timedelta(days=days_ago),
                    sample_size=random.randint(*template["size_range"]),
                    impact_score=impact,  # Use calculated impact
                    status="Completed" if random.random() > 0.2 else "Active",
                    summary=template["summary"]
                )
                projects.append(project)
            
            # Sort by date
            projects.sort(key=lambda p: p.added_date)
            
            stage_num = int(stage[1])
            
            candidate = DrugCandidate(
                id=f"cand_{candidate_id}",
                drug_name=drug_name,
                current_indication=current_use,
                target_sud_subtype=category.name,
                mechanism=mechanism,
                stage=stage,
                evidence_score=evidence_score,
                baseline_score=baseline_score,  # Now uses ensemble for Real scores
                smiles=drug_smiles,
                attached_projects=projects,
                last_updated=datetime.now() - timedelta(days=random.randint(1, 30)),
                cohort_count=random.randint(1 + stage_num, 4 + stage_num * 2),
                has_market_analysis=random.random() > 0.3,
                has_validation_plan=random.random() > 0.4,
                team_members=random.randint(2 + stage_num, 5 + stage_num * 2),
                data_produced=random.randint(1 + stage_num, 3 + stage_num * 2),
                publications=random.randint(0 + stage_num // 2, 2 + stage_num),
                tools_used=random.randint(1 + stage_num // 2, 3 + stage_num),
                data_governance=random.randint(1, 2 + stage_num // 2),
                training_participated=random.randint(stage_num, 2 + stage_num),
                stage_entry_date=stage_entry_date,
                stage_history=stage_history,
                score_type=score_type,
                model_scores=model_scores,
                protein_targets=drug_targets,
                disease_id=category.disease_id
            )
            candidates.append(candidate)
            candidate_id += 1
        
        # Update category count
        category.candidate_count = category_real + category_synthetic
        print(f"   Summary: {category_real} real, {category_synthetic} synthetic")
    
    print(f"\n{'='*70}")
    print(f" FINAL STATISTICS")
    print(f"{'='*70}")
    print(f"Total candidates: {len(candidates)}")
    print(f"Real ML/DL scores: {total_real_scores}")
    print(f"Synthetic scores: {total_synthetic_scores}")
    print(f"Real score percentage: {total_real_scores / len(candidates) * 100:.1f}%")
    print(f"{'='*70}\n")
    
    return sud_categories, candidates
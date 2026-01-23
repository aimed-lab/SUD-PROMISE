"""
🧬 SUD-PROMISE Drug Repositioning Assessment Platform 
UAB Color Theme: Forest Green/Blue (#1E6B52) & White
Font: Times New Roman
Optimized for HuggingFace Spaces - WITH WORKING BUTTONS
"""

import gradio as gr
import pandas as pd
import plotly.graph_objects as go
from datetime import datetime, timedelta
import random
from dataclasses import dataclass
from typing import List, Optional
import json

# ========================================
# UAB COLOR PALETTE (Blue/Green & White)
# ========================================
UAB_GREEN = "#1E6B52"  # Primary UAB Green
UAB_DARK_GREEN = "#16533E"  # Darker shade
UAB_LIGHT_GREEN = "#5A9B7F"  # Lighter shade
UAB_ACCENT_TEAL = "#008C95"  # Teal accent
UAB_PALE_GREEN = "#E8F4F0"  # Very light green for backgrounds
UAB_WHITE = "#FFFFFF"  # White

# ========================================
# SYNTHETIC DATA MODELS
# ========================================

@dataclass
class Project:
    """Research project providing evidence"""
    id: str
    name: str
    project_type: str  # Clinical Trial, Meta-Analysis, RWE, Biomarker
    added_date: datetime
    sample_size: int
    impact_score: float  # How much it changes prediction
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
    status: str  # Discovery, Preclinical, Phase I, II, III
    evidence_score: float  # 0-1
    baseline_score: float  # Score without projects
    smiles: str
    attached_projects: List[Project]
    last_updated: datetime
    cohort_count: int
    has_market_analysis: bool
    has_validation_plan: bool

@dataclass
class SUDCategory:
    """SUD disease category"""
    name: str
    color: str
    hex_color: str
    icon: str
    candidate_count: int
    description: str

# ========================================
# SYNTHETIC DATA GENERATION
# ========================================

def generate_synthetic_data():
    """Generate realistic synthetic SUD repositioning data"""
    
    # SUD Categories - UAB Green/Blue theme only
    sud_categories = [
        SUDCategory("Opioid Use Disorder", "", UAB_GREEN, "", 12, 
                   "Addiction to opioids including prescription painkillers, heroin, and fentanyl"),
        SUDCategory("Alcohol Use Disorder", "", UAB_ACCENT_TEAL, "", 8,
                   "Problematic pattern of alcohol use leading to clinically significant impairment"),
        SUDCategory("Stimulant Use Disorder", "", UAB_LIGHT_GREEN, "", 6,
                   "Addiction to cocaine, methamphetamine, or prescription stimulants"),
        SUDCategory("Cannabis Use Disorder", "", UAB_DARK_GREEN, "", 5,
                   "Problematic cannabis use with withdrawal and tolerance symptoms"),
        SUDCategory("Sedative/Hypnotic Disorder", "", "#4A8B7A", "", 3,
                   "Dependence on benzodiazepines or other sedative medications"),
        SUDCategory("Nicotine Use Disorder", "", "#2C7A64", "", 4,
                   "Tobacco/nicotine dependence and addiction"),
    ]
    
    # Sample projects
    project_templates = [
        {
            "name": "NIDA Clinical Trial NCT02892123",
            "type": "Clinical Trial",
            "size_range": (150, 500),
            "impact_range": (0.10, 0.20),
            "summary": "Randomized controlled trial showing significant reduction in relapse rates"
        },
        {
            "name": "Meta-Analysis Lee et al. 2024",
            "type": "Meta-Analysis",
            "size_range": (800, 3000),
            "impact_range": (0.08, 0.15),
            "summary": "Pooled analysis of multiple RCTs demonstrating efficacy across populations"
        },
        {
            "name": "RWE Claims Database Analysis",
            "type": "Real-World Evidence",
            "size_range": (2000, 10000),
            "impact_range": (0.05, 0.12),
            "summary": "Insurance claims data showing real-world treatment adherence and outcomes"
        },
        {
            "name": "Biomarker Validation Study",
            "type": "Biomarker Study",
            "size_range": (50, 200),
            "impact_range": (0.06, 0.10),
            "summary": "Identification of predictive biomarkers for treatment response"
        },
        {
            "name": "Phase II Safety Trial",
            "type": "Clinical Trial",
            "size_range": (100, 300),
            "impact_range": (0.07, 0.13),
            "summary": "Safety and tolerability data in target SUD population"
        }
    ]
    
    # Sample drugs for each category
    drug_templates = {
        "Opioid Use Disorder": [
            ("Naltrexone", "Alcohol/opioid dependence", "μ-opioid receptor antagonist", "Phase III"),
            ("Buprenorphine", "Opioid dependence", "Partial μ-opioid agonist", "Phase III"),
            ("Lofexidine", "Hypertension", "α2-adrenergic agonist", "Phase II"),
            ("Gabapentin", "Epilepsy, neuropathic pain", "GABA analog", "Phase II"),
            ("Topiramate", "Epilepsy, migraine", "GABA modulator", "Phase I"),
            ("Ondansetron", "Nausea/vomiting", "5-HT3 antagonist", "Preclinical"),
        ],
        "Alcohol Use Disorder": [
            ("Acamprosate", "Alcohol dependence", "NMDA antagonist", "Phase III"),
            ("Naltrexone", "Opioid dependence", "μ-opioid antagonist", "Phase III"),
            ("Baclofen", "Muscle spasticity", "GABA-B agonist", "Phase II"),
            ("Topiramate", "Epilepsy", "GABA modulator", "Phase II"),
            ("Varenicline", "Smoking cessation", "Nicotinic receptor agonist", "Phase I"),
        ],
        "Stimulant Use Disorder": [
            ("Modafinil", "Narcolepsy", "Dopamine reuptake inhibitor", "Phase II"),
            ("Bupropion", "Depression", "NDRI", "Phase II"),
            ("Topiramate", "Epilepsy", "GABA modulator", "Phase I"),
            ("N-acetylcysteine", "Acetaminophen overdose", "Antioxidant", "Phase I"),
        ],
        "Cannabis Use Disorder": [
            ("N-acetylcysteine", "Acetaminophen overdose", "Antioxidant", "Phase II"),
            ("Gabapentin", "Epilepsy", "GABA analog", "Phase I"),
            ("Zolpidem", "Insomnia", "GABA-A modulator", "Preclinical"),
        ],
        "Sedative/Hypnotic Disorder": [
            ("Flumazenil", "Benzodiazepine overdose", "GABA-A antagonist", "Phase I"),
            ("Gabapentin", "Epilepsy", "GABA analog", "Phase I"),
        ],
        "Nicotine Use Disorder": [
            ("Varenicline", "Smoking cessation", "Nicotinic receptor agonist", "Phase III"),
            ("Cytisine", "Smoking cessation (EU)", "Nicotinic receptor agonist", "Phase II"),
            ("Bupropion", "Depression", "NDRI", "Phase III"),
        ],
    }
    
    # Generate candidates
    candidates = []
    candidate_id = 1
    
    for category in sud_categories:
        if category.name not in drug_templates:
            continue
            
        for drug_info in drug_templates[category.name]:
            drug_name, current_use, mechanism, status = drug_info
            
            # Generate baseline score
            baseline = random.uniform(0.45, 0.65)
            
            # Generate 1-4 projects
            num_projects = random.randint(1, 4)
            projects = []
            
            for i in range(num_projects):
                template = random.choice(project_templates)
                days_ago = random.randint(30, 365)
                
                project = Project(
                    id=f"proj_{candidate_id}_{i}",
                    name=template["name"],
                    project_type=template["type"],
                    added_date=datetime.now() - timedelta(days=days_ago),
                    sample_size=random.randint(*template["size_range"]),
                    impact_score=random.uniform(*template["impact_range"]),
                    status="Completed" if random.random() > 0.2 else "Active",
                    summary=template["summary"]
                )
                projects.append(project)
            
            # Calculate evidence score (baseline + project impacts)
            total_impact = sum(p.impact_score for p in projects)
            evidence_score = min(0.95, baseline + total_impact)
            
            # Generate SMILES (simplified)
            smiles = f"C{'C' * random.randint(5, 15)}{'N' if random.random() > 0.5 else 'O'}"
            
            candidate = DrugCandidate(
                id=f"cand_{candidate_id}",
                drug_name=drug_name,
                current_indication=current_use,
                target_sud_subtype=category.name,
                mechanism=mechanism,
                status=status,
                evidence_score=evidence_score,
                baseline_score=baseline,
                smiles=smiles,
                attached_projects=sorted(projects, key=lambda p: p.added_date),
                last_updated=datetime.now() - timedelta(days=random.randint(1, 30)),
                cohort_count=random.randint(0, 3),
                has_market_analysis=random.random() > 0.3,
                has_validation_plan=random.random() > 0.4
            )
            candidates.append(candidate)
            candidate_id += 1
    
    return sud_categories, candidates

# Generate data
SUD_CATEGORIES, CANDIDATES = generate_synthetic_data()

# ========================================
# HELPER FUNCTIONS
# ========================================

def get_status_badge(status):
    """Return colored badge HTML for status - UAB blue/green theme"""
    colors = {
        "Discovery": "#CBD5E0",
        "Preclinical": UAB_GREEN,
        "Phase I": UAB_GREEN,
        "Phase II": UAB_GREEN,
        "Phase III": UAB_DARK_GREEN,
        "FDA Review": UAB_DARK_GREEN
    }
    color = colors.get(status, "#CBD5E0")
    return f'<span style="background: {color}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white; font-family: \'Times New Roman\', Times, serif;">{status}</span>'

def get_evidence_stars(score):
    """Convert score to star rating"""
    stars = int(score * 5)
    return "⭐" * stars

def format_date_ago(date):
    """Format datetime as 'X days ago'"""
    delta = datetime.now() - date
    if delta.days == 0:
        return "Today"
    elif delta.days == 1:
        return "Yesterday"
    elif delta.days < 7:
        return f"{delta.days} days ago"
    elif delta.days < 30:
        return f"{delta.days // 7} weeks ago"
    else:
        return f"{delta.days // 30} months ago"

def get_impact_badge(impact):
    """Return badge for project impact - UAB theme"""
    if impact > 0.10:
        return f'<span style="background: {UAB_DARK_GREEN}; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">High +{impact:.2f}</span>'
    elif impact > 0.05:
        return f'<span style="background: {UAB_GREEN}; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Moderate +{impact:.2f}</span>'
    else:
        return f'<span style="background: #A0AEC0; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Low +{impact:.2f}</span>'

def render_candidate_card(candidate, rank, border_color=None):
    """Render a single candidate card HTML - compact design with evidence score inline"""
    stars = get_evidence_stars(candidate.evidence_score)
    status_badge = get_status_badge(candidate.status)
    
    # Always use UAB green for consistency
    border = UAB_GREEN
    
    # Market analysis as number (1 or 0)
    market_analysis_num = 1 if candidate.has_market_analysis else 0
    
    html = f"""
    <div style="background: white; border: 2px solid {border}; padding: 15px 20px; margin: 10px 20px; border-radius: 12px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); font-family: 'Times New Roman', Times, serif;">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px;">
            <div style="display: flex; align-items: center; gap: 12px;">
                <span style="background: #EDF2F7; padding: 4px 10px; border-radius: 20px; font-size: 14px; color: {UAB_DARK_GREEN}; font-weight: bold; font-family: 'Times New Roman', Times, serif;">{rank}</span>
                <span style="font-size: 24px; font-weight: bold; color: {UAB_DARK_GREEN}; font-family: 'Times New Roman', Times, serif;">{candidate.drug_name}</span>
                <span style="font-size: 16px; color: {UAB_GREEN}; font-weight: bold; font-family: 'Times New Roman', Times, serif;">{stars} {candidate.evidence_score:.2f}</span>
            </div>
            {status_badge}
        </div>
        
        <p style="margin: 8px 0; color: #4A5568; font-family: 'Times New Roman', Times, serif;"><b>Current Use:</b> {candidate.current_indication}</p>
        <p style="margin: 8px 0; color: #4A5568; font-family: 'Times New Roman', Times, serif;"><b>Mechanism:</b> {candidate.mechanism}</p>
        
        <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin: 12px 0;">
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px; font-family: 'Times New Roman', Times, serif;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{len(candidate.attached_projects)}</div>
                <div style="font-size: 12px; color: white;">Projects</div>
            </div>
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px; font-family: 'Times New Roman', Times, serif;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{candidate.cohort_count}</div>
                <div style="font-size: 12px; color: white;">Cohorts</div>
            </div>
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px; font-family: 'Times New Roman', Times, serif;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{market_analysis_num}</div>
                <div style="font-size: 12px; color: white;">Market Analysis</div>
            </div>
        </div>
        
        <p style="margin: 8px 0 0 0; font-size: 13px; color: #A0AEC0; font-family: 'Times New Roman', Times, serif;">Last updated: {format_date_ago(candidate.last_updated)}</p>
    </div>
    """
    return html

# ========================================
# VIEW RENDERERS
# ========================================

def render_landing_dashboard():
    """Step 1: Landing dashboard showing SUD landscape - UAB theme"""
    
    total_candidates = len(CANDIDATES)
    active_projects = sum(len(c.attached_projects) for c in CANDIDATES)
    total_cohorts = sum(c.cohort_count for c in CANDIDATES)
    
    # Top candidates
    top_candidates = sorted(CANDIDATES, key=lambda c: c.evidence_score, reverse=True)[:5]
    
    html = f"""
    <div style="padding: 20px; font-family: 'Times New Roman', Times, serif;">
        <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 30px 0;">
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 36px; font-family: 'Times New Roman', Times, serif; color: white;">{total_candidates}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Drug Candidates</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 36px; font-family: 'Times New Roman', Times, serif; color: white;">{active_projects}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Evidence Projects</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 36px; font-family: 'Times New Roman', Times, serif; color: white;">{total_cohorts}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Patient Cohorts</p>
            </div>
        </div>
    </div>
    """
    
    # Category buttons
    category_choices = [f"{cat.name} ({cat.candidate_count} candidates)" for cat in SUD_CATEGORIES]
    
    return html, category_choices, top_candidates

def render_category_view(category_selection, sort_by="Evidence Score"):
    """Step 2: Show candidates for selected SUD category - UAB theme"""
    
    if not category_selection:
        return "<p>Please select a category</p>", [], []
    
    # Parse category name
    category_name = category_selection.split(" (")[0].strip()
    
    # Get category
    category = next((c for c in SUD_CATEGORIES if c.name == category_name), None)
    if not category:
        return "<p>Category not found</p>", [], []
    
    # Filter candidates
    filtered = [c for c in CANDIDATES if c.target_sud_subtype == category_name]
    
    # Sort - map display names to sort keys
    if sort_by == "Evidence Score" or sort_by == "evidence_score":
        filtered.sort(key=lambda c: c.evidence_score, reverse=True)
    elif sort_by == "Recent" or sort_by == "recent":
        filtered.sort(key=lambda c: c.last_updated, reverse=True)
    elif sort_by == "Name" or sort_by == "name":
        filtered.sort(key=lambda c: c.drug_name)
    
    html = f"""
    <div style="padding: 20px; font-family: 'Times New Roman', Times, serif;">
        <h1 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">{category.name}</h1>
        <p style="color: #718096; margin-bottom: 20px; font-family: 'Times New Roman', Times, serif;">{category.description}</p>
        <p style="background: #F7FAFC; padding: 15px; border-radius: 10px; border-left: 4px solid {category.hex_color}; font-family: 'Times New Roman', Times, serif;">
            <b style="color: {UAB_DARK_GREEN};">{len(filtered)} drug candidates</b> being evaluated for repositioning
        </p>
        
        <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; 
                    margin: 20px 0; text-align: center; color: white; box-shadow: 0 4px 6px rgba(30,107,82,0.3); font-family: 'Times New Roman', Times, serif;">
            <h3 style="margin: 0 0 10px 0; font-family: 'Times New Roman', Times, serif; color: white;">View Detailed Dashboard</h3>
            <p style="margin: 0; font-size: 14px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Click the "📊 View Details" button below any candidate card to view complete analysis with evidence timeline</p>
        </div>
    </div>
    """
    
    # Return candidate choices for dropdown
    candidate_choices = [f"{c.drug_name} (Score: {c.evidence_score:.2f})" for c in filtered]
    
    return html, candidate_choices, filtered

def render_candidate_dashboard(candidate_selection, category_selection):
    """Step 4-5: Detailed dashboard for selected candidate - UAB theme"""
    
    if not candidate_selection:
        return "<p>Please select a candidate</p>", None, ""
    
    # Parse candidate name
    candidate_name = candidate_selection.split(" (Score:")[0].strip()
    
    # Get candidate
    candidate = next((c for c in CANDIDATES if c.drug_name == candidate_name), None)
    if not candidate:
        return "<p>Candidate not found</p>", None, ""
    
    stars = get_evidence_stars(candidate.evidence_score)
    status_badge = get_status_badge(candidate.status)
    
    # Create timeline figure
    timeline_fig = create_evidence_timeline(candidate)
    
    # PART 1: Everything BEFORE the timeline plot
    html_before_plot = f"""
    <div style="padding: 20px; font-family: 'Times New Roman', Times, serif;">
        <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 30px; border-radius: 15px; color: white; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
            <h1 style="margin: 0; font-family: 'Times New Roman', Times, serif; color: white;">{candidate.drug_name}</h1>
            <p style="margin: 10px 0 0 0; font-size: 18px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Repositioning Candidate for {candidate.target_sud_subtype}</p>
            <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.8; font-family: 'Times New Roman', Times, serif; color: white;">University of Alabama at Birmingham - SUD-PROMISE</p>
        </div>
        
        <div style="display: flex; gap: 10px; margin-bottom: 30px;">
            <div>{status_badge}</div>
            <div style="background: {UAB_GREEN}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white; font-family: 'Times New Roman', Times, serif;">
                {stars} {candidate.evidence_score:.2f}
            </div>
        </div>
        
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Drug Information</h2>
        <div style="background: #F0F9F6; padding: 20px; border-radius: 10px; margin-bottom: 30px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <p style="margin: 10px 0;"><b>Current Indication:</b> {candidate.current_indication}</p>
            <p style="margin: 10px 0;"><b>Target SUD:</b> {candidate.target_sud_subtype}</p>
            <p style="margin: 10px 0;"><b>Mechanism of Action:</b> {candidate.mechanism}</p>
            <p style="margin: 10px 0;"><b>SMILES:</b> <code>{candidate.smiles}</code></p>
        </div>
        
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Evidence Evolution</h2>
        <div style="background: {UAB_PALE_GREEN}; padding: 20px; border-radius: 10px; margin-bottom: 20px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <p><b>Baseline Score:</b> {candidate.baseline_score:.2f} (no evidence)</p>
            <p><b>Current Score:</b> <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.2f}</span></p>
            <p style="color: {UAB_GREEN}; font-weight: bold;">Total Improvement: +{candidate.evidence_score - candidate.baseline_score:.2f} ({((candidate.evidence_score - candidate.baseline_score) / candidate.baseline_score * 100):.1f}%)</p>
        </div>
    </div>
    """
    
    # PART 2: Everything AFTER the timeline plot
    html_after_plot = f"""
    <div style="padding: 0 20px 20px 20px; font-family: 'Times New Roman', Times, serif;">
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Attached Evidence Projects ({len(candidate.attached_projects)})</h2>
        <p style="color: #718096; margin-bottom: 15px; font-family: 'Times New Roman', Times, serif;">Projects contributing to prediction confidence</p>
    """
    
    # Render projects
    for project in candidate.attached_projects:
        impact_badge = get_impact_badge(project.impact_score)
        
        html_after_plot += f"""
        <div style="background: white; border-left: 5px solid {UAB_GREEN}; padding: 20px; margin: 15px 0; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); font-family: 'Times New Roman', Times, serif;">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <h3 style="margin: 0; color: {UAB_DARK_GREEN}; font-family: 'Times New Roman', Times, serif;">{project.name}</h3>
                {impact_badge}
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px; margin: 15px 0;">
                <div>
                    <p style="margin: 5px 0; color: #4A5568;"><b>Type:</b> {project.project_type}</p>
                    <p style="margin: 5px 0; color: #4A5568;"><b>Sample Size:</b> {project.sample_size:,} patients</p>
                </div>
                <div>
                    <p style="margin: 5px 0; color: #4A5568;"><b>Status:</b> {project.status}</p>
                    <p style="margin: 5px 0; color: #4A5568;"><b>Added:</b> {format_date_ago(project.added_date)}</p>
                </div>
            </div>
            
            <div style="background: #F0F9F6; padding: 15px; border-radius: 8px; margin-top: 10px;">
                <p style="margin: 0; font-style: italic; color: #4A5568;">"{project.summary}"</p>
            </div>
        </div>
        """
    
    html_after_plot += f"""
        <div style="margin-top: 30px; padding: 20px; background: #F0F9F6; border-radius: 10px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <h3 style="color: {UAB_DARK_GREEN}; font-family: 'Times New Roman', Times, serif;">Summary Statistics</h3>
            <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px;">
                <div style="text-align: center; font-family: 'Times New Roman', Times, serif;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{len(candidate.attached_projects)}</div>
                    <div style="color: #4A5568;">Evidence Projects</div>
                </div>
                <div style="text-align: center; font-family: 'Times New Roman', Times, serif;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.cohort_count}</div>
                    <div style="color: #4A5568;">Patient Cohorts</div>
                </div>
                <div style="text-align: center; font-family: 'Times New Roman', Times, serif;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{"✅" if candidate.has_validation_plan else "—"}</div>
                    <div style="color: #4A5568;">Validation Plan</div>
                </div>
                <div style="text-align: center; font-family: 'Times New Roman', Times, serif;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{"✅" if candidate.has_market_analysis else "—"}</div>
                    <div style="color: #4A5568;">Market Analysis</div>
                </div>
            </div>
        </div>
    </div>
    """
    
    return html_before_plot, timeline_fig, html_after_plot

def create_evidence_timeline(candidate):
    """Create plotly timeline showing evidence evolution - UAB colors with Times New Roman"""
    
    # Build timeline data
    timeline_data = []
    
    # Baseline
    timeline_data.append({
        'date': candidate.attached_projects[0].added_date - timedelta(days=30) if candidate.attached_projects else datetime.now() - timedelta(days=365),
        'score': candidate.baseline_score,
        'label': 'Baseline (No Evidence)',
        'color': '#A0AEC0'
    })
    
    # Add each project cumulatively
    cumulative_score = candidate.baseline_score
    for project in candidate.attached_projects:
        cumulative_score += project.impact_score
        cumulative_score = min(0.95, cumulative_score)  # Cap at 0.95
        
        timeline_data.append({
            'date': project.added_date,
            'score': cumulative_score,
            'label': f'+ {project.name}',
            'color': UAB_GREEN,
            'impact': f'+{project.impact_score:.2f}'
        })
    
    # Create figure
    fig = go.Figure()
    
    # Add line - UAB green
    fig.add_trace(go.Scatter(
        x=[d['date'] for d in timeline_data],
        y=[d['score'] for d in timeline_data],
        mode='lines+markers',
        line=dict(color=UAB_GREEN, width=3),
        marker=dict(
            size=12,
            color=[d['color'] for d in timeline_data],
            line=dict(width=2, color='white')
        ),
        text=[d['label'] for d in timeline_data],
        hovertemplate='<b>%{text}</b><br>Score: %{y:.3f}<br>%{x|%b %d, %Y}<extra></extra>'
    ))
    
    # Add annotations for significant changes - UAB green
    for i in range(1, len(timeline_data)):
        if 'impact' in timeline_data[i]:
            fig.add_annotation(
                x=timeline_data[i]['date'],
                y=timeline_data[i]['score'],
                text=timeline_data[i]['impact'],
                showarrow=True,
                arrowhead=2,
                arrowcolor=UAB_GREEN,
                bgcolor=UAB_GREEN,
                bordercolor=UAB_GREEN,
                font=dict(color='white', size=10, family="Times New Roman, serif"),
                yshift=10
            )
    
    fig.update_layout(
        title="Evidence Score Evolution Over Time",
        xaxis_title="Date",
        yaxis_title="Prediction Score",
        hovermode='closest',
        height=400,
        plot_bgcolor='#F0F9F6',
        paper_bgcolor='white',
        font=dict(family="Times New Roman, serif", color=UAB_DARK_GREEN),
        yaxis=dict(range=[0, 1])
    )
    
    return fig

# ========================================
# GRADIO INTERFACE
# ========================================

def create_interface():
    """Create the Gradio interface with UAB theme and Times New Roman font"""
    
    # Create custom UAB theme
    uab_theme = gr.themes.Soft(
        primary_hue="emerald",
        secondary_hue="slate",
    ).set(
        button_primary_background_fill=UAB_GREEN,
        button_primary_background_fill_hover=UAB_DARK_GREEN,
        button_primary_text_color="white",
        button_secondary_border_color=UAB_GREEN,
        button_secondary_text_color=UAB_GREEN,
        body_text_color="#2D3748",
        input_border_color="#CBD5E0",
        input_border_color_focus=UAB_GREEN,
    )
    
    with gr.Blocks(title="SUD-PROMISE | UAB", theme=uab_theme) as demo:
        
        # Add global CSS for Times New Roman font
        gr.HTML("""
        <style>
            * {
                font-family: 'Times New Roman', Times, serif !important;
            }
            
            /* Font family for buttons */
            .gr-button {
                font-family: 'Times New Roman', Times, serif !important;
            }
            
            /* Selected dropdown items */
            .choices__item--selectable.is-highlighted {
                background-color: #1E6B52 !important;
            }
        </style>
        """)
        
        # State management
        current_view = gr.State("dashboard")
        selected_category = gr.State(None)
        selected_candidate = gr.State(None)
        
        with gr.Column():
            # Header - UAB colors with icon
            gr.Markdown(f"""
            <div style="text-align: center; padding: 20px; background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); border-radius: 15px; margin-bottom: 20px; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <div style="font-size: 50px; margin-bottom: 10px;">🧬</div>
                <h1 style="color: white; margin: 0; font-family: 'Times New Roman', Times, serif;">SUD-PROMISE Drug Repositioning Assessment Platform</h1>
                <p style="color: white; opacity: 0.9; margin: 10px 0 0 0; font-family: 'Times New Roman', Times, serif;">University of Alabama at Birmingham - Evidence-Based Drug Discovery for Substance Use Disorders</p>
            </div>
            """)
            
            # Navigation breadcrumb
            breadcrumb = gr.Markdown("")
            
            # Main content area
            with gr.Column() as main_content:
                
                # Dashboard view
                with gr.Column(visible=True) as dashboard_view:
                    dashboard_html = gr.HTML()

                    gr.Markdown("---")
                    gr.Markdown("### Top Candidates by Evidence Score")
                    
                    # Store top candidate buttons
                    top_candidate_components = []
                    for i in range(5):
                        with gr.Column():
                            top_html = gr.HTML()
                            top_btn = gr.Button(f"📊 View Details", variant="primary", size="sm")
                            top_candidate_components.append((top_html, top_btn))
                    
                    gr.Markdown("---")
                    gr.Markdown("### Select a SUD Category to Explore")
                    category_selector = gr.Dropdown(
                        choices=[],
                        label="",
                        show_label=False,  # Completely hide label container
                        interactive=True
                    )
                    view_category_btn = gr.Button("View Candidates", variant="primary", size="lg")                
                # Category view
                with gr.Column(visible=False) as category_view:
                    category_html = gr.HTML()
                    
                    gr.Markdown("**Sort by**")
                    with gr.Row():
                        sort_dropdown = gr.Dropdown(
                            choices=["Evidence Score", "Recent", "Name"],
                            value="Evidence Score",
                            label="",
                            show_label=False,  # Completely hide label container
                            interactive=True
                        )
                    
                    # Hidden candidate selector (kept for compatibility but not displayed)
                    candidate_selector = gr.Dropdown(
                        choices=[],
                        label="",
                        show_label=False,  # Completely hide label container
                        interactive=True,
                        visible=False
                    )
                    
                    # Store category candidate buttons - will be created dynamically
                    category_candidate_components = []
                    for i in range(12):  # Max candidates per category
                        with gr.Column(visible=False) as cat_card_col:
                            cat_html = gr.HTML()
                            cat_btn = gr.Button(f"📊 View Details", variant="primary", size="sm")
                            category_candidate_components.append((cat_card_col, cat_html, cat_btn))
                    
                    with gr.Row():
                        back_to_dashboard_btn = gr.Button("← Back to Dashboard", variant="secondary")
                
                # Candidate detail view
                with gr.Column(visible=False) as candidate_view:
                    candidate_html_before = gr.HTML()
                    
                    # Timeline plot positioned between evidence evolution and projects
                    gr.Markdown("### Evidence Score Evolution Timeline")
                    timeline_plot = gr.Plot()
                    
                    gr.Markdown("---")
                    
                    candidate_html_after = gr.HTML()
                    
                    back_to_category_btn = gr.Button("← Back to Category", variant="secondary", size="lg")
        
        # ========================================
        # EVENT HANDLERS
        # ========================================
        
        def show_dashboard():
            """Show landing dashboard with top 5 candidates"""
            dash_html, categories, top_5 = render_landing_dashboard()
            
            # Render top 5 candidate cards
            top_outputs = [dash_html, gr.update(choices=categories, value=None)]
            
            for i, (html_component, btn_component) in enumerate(top_candidate_components):
                if i < len(top_5):
                    candidate = top_5[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    top_outputs.append(card_html)
                else:
                    top_outputs.append("")
            
            # Hide ALL category cards (col, html, btn)
            for col, html, btn in category_candidate_components:
                top_outputs.extend([
                    gr.update(visible=False),  # column
                    "",  # html
                    gr.update(visible=False)  # button
                ])
            
            top_outputs.extend([
                gr.update(visible=True),  # dashboard_view
                gr.update(visible=False),  # category_view
                gr.update(visible=False),  # candidate_view
                "dashboard",  # current_view
                "📊 Dashboard"  # breadcrumb
            ])
            
            return top_outputs
        
        def show_category(category_selection):
            """Show category view with candidate cards"""
            if not category_selection:
                return [gr.update()] * (8 + len(category_candidate_components) * 3)
            
            cat_html, candidates, filtered = render_category_view(category_selection)
            category_name = category_selection.split(" (")[0]
            
            outputs = [
                cat_html,  # category_html
                gr.update(choices=candidates, value=None),  # candidate_selector
            ]
            
            # Render candidate cards - explicitly control all three components (col, html, btn)
            for i, (col, html, btn) in enumerate(category_candidate_components):
                if i < len(filtered):
                    candidate = filtered[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    outputs.extend([
                        gr.update(visible=True),  # column
                        card_html,  # html
                        gr.update(visible=True)  # button
                    ])
                else:
                    # Explicitly hide ALL three components for unused slots
                    outputs.extend([
                        gr.update(visible=False),  # column
                        "",  # html
                        gr.update(visible=False)  # button
                    ])
            
            outputs.extend([
                gr.update(visible=False),  # dashboard_view
                gr.update(visible=True),  # category_view
                gr.update(visible=False),  # candidate_view
                category_selection,  # selected_category
                "category",  # current_view
                f"📊 Dashboard > {category_name}",  # breadcrumb
            ])
            
            return outputs
        
        def show_candidate(candidate_selection, category_selection):
            """Show candidate detail view"""
            if not candidate_selection:
                return [gr.update()] * 8
            
            html_before, timeline_fig, html_after = render_candidate_dashboard(candidate_selection, category_selection)
            candidate_name = candidate_selection.split(" (Score:")[0]
            category_name = category_selection.split(" (")[0] if category_selection else "Category"
            
            return (
                f"📊 Dashboard > {category_name} > {candidate_name}",  # breadcrumb
                html_before,  # candidate_html_before
                timeline_fig,  # timeline_plot
                html_after,   # candidate_html_after
                gr.update(visible=False),  # dashboard_view
                gr.update(visible=False),  # category_view
                gr.update(visible=True),  # candidate_view
                category_selection,  # selected_category - IMPORTANT!
            )
        
        def update_category_sort(category_selection, sort_by):
            """Update category view when sorting changes"""
            if not category_selection:
                return [gr.update()] * (2 + len(category_candidate_components) * 3)
                
            cat_html, candidates, filtered = render_category_view(category_selection, sort_by)
            
            outputs = [cat_html, gr.update(choices=candidates)]
            
            # Re-render candidate cards with explicit button visibility
            for i, (col, html, btn) in enumerate(category_candidate_components):
                if i < len(filtered):
                    candidate = filtered[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    outputs.extend([
                        gr.update(visible=True),  # column
                        card_html,  # html
                        gr.update(visible=True)  # button
                    ])
                else:
                    outputs.extend([
                        gr.update(visible=False),  # column
                        "",  # html
                        gr.update(visible=False)  # button
                    ])
            
            return outputs
        
        def make_view_candidate_handler(candidate_index, is_top=False):
            """Factory function to create handlers for View Details buttons"""
            def handler():
                if is_top:
                    # Top 5 candidates
                    top_5 = sorted(CANDIDATES, key=lambda c: c.evidence_score, reverse=True)[:5]
                    if candidate_index >= len(top_5):
                        return [gr.update()] * 8
                    candidate = top_5[candidate_index]
                else:
                    # Category candidates - need current category
                    return [gr.update()] * 8  # Will be overridden by actual handler
                
                category_str = f"{candidate.target_sud_subtype} ({sum(1 for c in CANDIDATES if c.target_sud_subtype == candidate.target_sud_subtype)} candidates)"
                candidate_str = f"{candidate.drug_name} (Score: {candidate.evidence_score:.2f})"
                
                html_before, timeline_fig, html_after = render_candidate_dashboard(candidate_str, category_str)
                candidate_name = candidate.drug_name
                category_name = candidate.target_sud_subtype
                
                return (
                    f"📊 Dashboard > {category_name} > {candidate_name}",  # breadcrumb
                    html_before,
                    timeline_fig,
                    html_after,
                    gr.update(visible=False),  # dashboard_view
                    gr.update(visible=False),  # category_view
                    gr.update(visible=True),  # candidate_view
                    category_str,  # selected_category - IMPORTANT!
                )
            return handler
        
        def make_category_candidate_handler(candidate_index):
            """Create handler for category candidate buttons"""
            def handler(category_selection):
                if not category_selection:
                    return [gr.update()] * 8
                
                category_name = category_selection.split(" (")[0].strip()
                filtered = [c for c in CANDIDATES if c.target_sud_subtype == category_name]
                filtered.sort(key=lambda c: c.evidence_score, reverse=True)
                
                if candidate_index >= len(filtered):
                    return [gr.update()] * 8
                
                candidate = filtered[candidate_index]
                candidate_str = f"{candidate.drug_name} (Score: {candidate.evidence_score:.2f})"
                
                html_before, timeline_fig, html_after = render_candidate_dashboard(candidate_str, category_selection)
                
                return (
                    f"📊 Dashboard > {category_name} > {candidate.drug_name}",
                    html_before,
                    timeline_fig,
                    html_after,
                    gr.update(visible=False),
                    gr.update(visible=False),
                    gr.update(visible=True),
                    category_selection,  # selected_category - IMPORTANT!
                )
            return handler
        
        # Connect events
        demo.load(
            fn=show_dashboard,
            outputs=[dashboard_html, category_selector] + 
                    [comp for pair in top_candidate_components for comp in [pair[0]]] +
                    [comp for triple in category_candidate_components for comp in [triple[0], triple[1], triple[2]]] +
                    [dashboard_view, category_view, candidate_view, current_view, breadcrumb]
        )
        
        view_category_btn.click(
            fn=show_category,
            inputs=[category_selector],
            outputs=[category_html, candidate_selector] +
                    [comp for triple in category_candidate_components for comp in [triple[0], triple[1], triple[2]]] +
                    [dashboard_view, category_view, candidate_view, selected_category, current_view, breadcrumb]
        )
        
        sort_dropdown.change(
            fn=update_category_sort,
            inputs=[selected_category, sort_dropdown],
            outputs=[category_html, candidate_selector] +
                    [comp for triple in category_candidate_components for comp in [triple[0], triple[1], triple[2]]]
        )
        
        # Wire top candidate buttons
        for i, (html, btn) in enumerate(top_candidate_components):
            btn.click(
                fn=make_view_candidate_handler(i, is_top=True),
                outputs=[breadcrumb, candidate_html_before, timeline_plot, candidate_html_after,
                        dashboard_view, category_view, candidate_view, selected_category]
            )
        
        # Wire category candidate buttons
        for i, (col, html, btn) in enumerate(category_candidate_components):
            btn.click(
                fn=make_category_candidate_handler(i),
                inputs=[selected_category],
                outputs=[breadcrumb, candidate_html_before, timeline_plot, candidate_html_after,
                        dashboard_view, category_view, candidate_view, selected_category]
            )
        
        # Dropdown auto-open
        candidate_selector.change(
            fn=show_candidate,
            inputs=[candidate_selector, selected_category],
            outputs=[breadcrumb, candidate_html_before, timeline_plot, candidate_html_after,
                    dashboard_view, category_view, candidate_view, selected_category]
        )
        
        back_to_dashboard_btn.click(
            fn=show_dashboard,
            outputs=[dashboard_html, category_selector] +
                    [comp for pair in top_candidate_components for comp in [pair[0]]] +
                    [comp for triple in category_candidate_components for comp in [triple[0], triple[1], triple[2]]] +
                    [dashboard_view, category_view, candidate_view, current_view, breadcrumb]
        )
        
        back_to_category_btn.click(
            fn=show_category,
            inputs=[selected_category],
            outputs=[category_html, candidate_selector] +
                    [comp for triple in category_candidate_components for comp in [triple[0], triple[1], triple[2]]] +
                    [dashboard_view, category_view, candidate_view, selected_category, current_view, breadcrumb]
        )
    
    return demo

# ========================================
# MAIN
# ========================================

if __name__ == "__main__":
    demo = create_interface()
    demo.launch()
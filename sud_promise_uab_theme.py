"""
🧬 SUD-PROMISE Drug Repositioning Assessment Platform 
UAB Color Theme: Forest Green (#1E6B52) & Gold (#FFB81C)
Optimized for HuggingFace Spaces
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
# UAB COLOR PALETTE
# ========================================
UAB_GREEN = "#1E6B52"
UAB_DARK_GREEN = "#16533E"
UAB_GOLD = "#FFB81C"
UAB_DARK_GOLD = "#E69D00"
UAB_LIGHT_GREEN = "#5A9B7F"
UAB_ACCENT_TEAL = "#008C95"

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
    
    # SUD Categories - UAB Green/Gold theme
    sud_categories = [
        SUDCategory("Opioid Use Disorder", "🟩", UAB_GREEN, "💊", 12, 
                   "Addiction to opioids including prescription painkillers, heroin, and fentanyl"),
        SUDCategory("Alcohol Use Disorder", "🟨", UAB_GOLD, "🍺", 8,
                   "Problematic pattern of alcohol use leading to clinically significant impairment"),
        SUDCategory("Stimulant Use Disorder", "🟢", UAB_LIGHT_GREEN, "⚡", 6,
                   "Addiction to cocaine, methamphetamine, or prescription stimulants"),
        SUDCategory("Cannabis Use Disorder", "🟩", UAB_ACCENT_TEAL, "🌿", 5,
                   "Problematic cannabis use with withdrawal and tolerance symptoms"),
        SUDCategory("Sedative/Hypnotic Disorder", "🟫", "#8B6914", "😴", 3,
                   "Dependence on benzodiazepines or other sedative medications"),
        SUDCategory("Nicotine Use Disorder", "🟧", UAB_DARK_GOLD, "🚬", 4,
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
    """Return colored badge HTML for status - UAB theme"""
    colors = {
        "Discovery": "#CBD5E0",
        "Preclinical": UAB_ACCENT_TEAL,
        "Phase I": UAB_GOLD,
        "Phase II": UAB_DARK_GOLD,
        "Phase III": UAB_GREEN,
        "FDA Review": UAB_DARK_GREEN
    }
    color = colors.get(status, "#CBD5E0")
    return f'<span style="background: {color}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white;">{status}</span>'

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
        return f'<span style="background: {UAB_GREEN}; padding: 3px 10px; border-radius: 10px; color: white;">🔥 High +{impact:.2f}</span>'
    elif impact > 0.05:
        return f'<span style="background: {UAB_GOLD}; padding: 3px 10px; border-radius: 10px; color: white;">⚡ Moderate +{impact:.2f}</span>'
    else:
        return f'<span style="background: #A0AEC0; padding: 3px 10px; border-radius: 10px; color: white;">✨ Low +{impact:.2f}</span>'

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
    <div style="padding: 20px;">
        <h1 style="color: {UAB_GREEN};">🧬 SUD-PROMISE Drug Repositioning Assessment Platform</h1>
        <p style="font-size: 16px; color: {UAB_DARK_GREEN};">University of Alabama at Birmingham - Evidence-Based Drug Discovery for Substance Use Disorders</p>
        
        <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 30px 0;">
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 36px;">{total_candidates}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9;">Drug Candidates</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GOLD} 0%, {UAB_DARK_GOLD} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(255,184,28,0.3);">
                <h2 style="margin: 0; font-size: 36px;">{active_projects}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9;">Evidence Projects</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_GOLD} 100%); padding: 25px; border-radius: 15px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 36px;">{total_cohorts}</h2>
                <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.9;">Patient Cohorts</p>
            </div>
        </div>
        
        <h2 style="color: {UAB_GREEN};">🎯 SUD Categories</h2>
        <p style="color: #718096; margin-bottom: 20px;">Select a category to view repositioning candidates</p>
    </div>
    """
    
    # Category buttons
    category_choices = [f"{cat.color} {cat.name} ({cat.candidate_count} candidates)" for cat in SUD_CATEGORIES]
    
    # Top candidates preview
    top_html = '<div style="padding: 0 20px;"><h2 style="color: ' + UAB_GREEN + ';">⭐ Top Candidates by Evidence Score</h2>'
    for i, candidate in enumerate(top_candidates, 1):
        stars = get_evidence_stars(candidate.evidence_score)
        status_badge = get_status_badge(candidate.status)
        top_html += f"""
        <div style="background: white; border: 2px solid #E2E8F0; padding: 15px; margin: 10px 0; border-radius: 10px; border-left: 5px solid {UAB_GREEN};">
            <div style="display: flex; justify-content: space-between; align-items: center;">
                <div>
                    <h3 style="margin: 0; color: {UAB_DARK_GREEN};">#{i}. {candidate.drug_name}</h3>
                    <p style="margin: 5px 0; color: #718096; font-size: 14px;">{candidate.target_sud_subtype}</p>
                </div>
                {status_badge}
            </div>
            <p style="margin: 10px 0; font-size: 16px;"><b>Evidence:</b> {stars} <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.2f}</span></p>
        </div>
        """
    top_html += '</div>'
    
    return html, top_html, category_choices

def render_category_view(category_selection, sort_by="evidence_score"):
    """Step 2: Show candidates for selected SUD category - UAB theme"""
    
    if not category_selection:
        return "<p>Please select a category</p>"
    
    # Parse category name
    category_name = category_selection.split(" (")[0].replace("🟩 ", "").replace("🟨 ", "").replace("🟢 ", "").replace("🟫 ", "").replace("🟧 ", "").strip()
    
    # Get category
    category = next((c for c in SUD_CATEGORIES if c.name == category_name), None)
    if not category:
        return "<p>Category not found</p>"
    
    # Filter candidates
    filtered = [c for c in CANDIDATES if c.target_sud_subtype == category_name]
    
    # Sort
    if sort_by == "evidence_score":
        filtered.sort(key=lambda c: c.evidence_score, reverse=True)
    elif sort_by == "recent":
        filtered.sort(key=lambda c: c.last_updated, reverse=True)
    elif sort_by == "name":
        filtered.sort(key=lambda c: c.drug_name)
    
    html = f"""
    <div style="padding: 20px;">
        <h1 style="color: {UAB_GREEN};">{category.color} {category.name}</h1>
        <p style="color: #718096; margin-bottom: 20px;">{category.description}</p>
        <p style="background: #F7FAFC; padding: 15px; border-radius: 10px; border-left: 4px solid {category.hex_color};">
            <b style="color: {UAB_DARK_GREEN};">{len(filtered)} drug candidates</b> being evaluated for repositioning
        </p>
        
        <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; 
                    margin: 20px 0; text-align: center; color: white; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
            <h3 style="margin: 0 0 10px 0;">📊 View Detailed Dashboard</h3>
            <p style="margin: 0; font-size: 14px; opacity: 0.9;">Select a candidate from the dropdown below to view complete analysis with evidence timeline</p>
        </div>
    </div>
    """
    
    # Render candidate cards
    for i, candidate in enumerate(filtered, 1):
        stars = get_evidence_stars(candidate.evidence_score)
        status_badge = get_status_badge(candidate.status)
        
        # Rank badge - UAB colors
        if i == 1:
            rank_badge = "🥇 #1"
            border_color = UAB_GOLD
        elif i == 2:
            rank_badge = "🥈 #2"
            border_color = UAB_LIGHT_GREEN
        elif i == 3:
            rank_badge = "🥉 #3"
            border_color = UAB_ACCENT_TEAL
        else:
            rank_badge = f"#{i}"
            border_color = "#E2E8F0"
        
        html += f"""
        <div style="background: white; border: 2px solid {border_color}; padding: 20px; margin: 15px 20px; border-radius: 12px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <div>
                    <span style="background: #EDF2F7; padding: 5px 10px; border-radius: 20px; font-size: 14px; margin-right: 10px; color: {UAB_DARK_GREEN}; font-weight: bold;">{rank_badge}</span>
                    <span style="font-size: 24px; font-weight: bold; color: {UAB_DARK_GREEN};">💊 {candidate.drug_name}</span>
                </div>
                {status_badge}
            </div>
            
            <p style="margin: 10px 0; color: #4A5568;"><b>Current Use:</b> {candidate.current_indication}</p>
            <p style="margin: 10px 0; color: #4A5568;"><b>Mechanism:</b> {candidate.mechanism}</p>
            
            <div style="margin: 15px 0; padding: 15px; background: #F0F9F6; border-radius: 8px; border: 1px solid {UAB_LIGHT_GREEN};">
                <p style="margin: 5px 0; font-size: 18px;"><b>Evidence Score:</b> {stars} <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.2f}</span></p>
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; margin: 15px 0;">
                <div style="text-align: center; padding: 10px; background: #F0F9F6; border-radius: 8px; border: 1px solid {UAB_LIGHT_GREEN};">
                    <div style="font-size: 24px;">📁</div>
                    <div style="font-weight: bold; color: {UAB_DARK_GREEN};">{len(candidate.attached_projects)}</div>
                    <div style="font-size: 12px; color: #4A5568;">Projects</div>
                </div>
                <div style="text-align: center; padding: 10px; background: #FFFBF0; border-radius: 8px; border: 1px solid {UAB_GOLD};">
                    <div style="font-size: 24px;">👥</div>
                    <div style="font-weight: bold; color: {UAB_DARK_GOLD};">{candidate.cohort_count}</div>
                    <div style="font-size: 12px; color: #4A5568;">Cohorts</div>
                </div>
                <div style="text-align: center; padding: 10px; background: #F0F9F6; border-radius: 8px; border: 1px solid {UAB_LIGHT_GREEN};">
                    <div style="font-size: 24px;">💰</div>
                    <div style="font-weight: bold; color: {UAB_GREEN};">{"✅" if candidate.has_market_analysis else "❌"}</div>
                    <div style="font-size: 12px; color: #4A5568;">Market</div>
                </div>
            </div>
            
            <p style="margin: 10px 0; font-size: 13px; color: #A0AEC0;">Last updated: {format_date_ago(candidate.last_updated)}</p>
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
    <div style="padding: 20px;">
        <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 30px; border-radius: 15px; color: white; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
            <h1 style="margin: 0;">💊 {candidate.drug_name}</h1>
            <p style="margin: 10px 0 0 0; font-size: 18px; opacity: 0.9;">Repositioning Candidate for {candidate.target_sud_subtype}</p>
            <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.8;">University of Alabama at Birmingham - SUD-PROMISE</p>
        </div>
        
        <div style="display: flex; gap: 10px; margin-bottom: 30px;">
            <div>{status_badge}</div>
            <div style="background: {UAB_GREEN}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white;">
                {stars} {candidate.evidence_score:.2f}
            </div>
        </div>
        
        <h2 style="color: {UAB_GREEN};">📋 Drug Information</h2>
        <div style="background: #F0F9F6; padding: 20px; border-radius: 10px; margin-bottom: 30px; border-left: 4px solid {UAB_GREEN};">
            <p style="margin: 10px 0;"><b>Current Indication:</b> {candidate.current_indication}</p>
            <p style="margin: 10px 0;"><b>Target SUD:</b> {candidate.target_sud_subtype}</p>
            <p style="margin: 10px 0;"><b>Mechanism of Action:</b> {candidate.mechanism}</p>
            <p style="margin: 10px 0;"><b>SMILES:</b> <code>{candidate.smiles}</code></p>
        </div>
        
        <h2 style="color: {UAB_GREEN};">📊 Evidence Evolution</h2>
        <div style="background: #FFFBF0; padding: 20px; border-radius: 10px; margin-bottom: 20px; border-left: 4px solid {UAB_GOLD};">
            <p><b>Baseline Score:</b> {candidate.baseline_score:.2f} (no evidence)</p>
            <p><b>Current Score:</b> <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.2f}</span></p>
            <p style="color: {UAB_GREEN}; font-weight: bold;">Total Improvement: +{candidate.evidence_score - candidate.baseline_score:.2f} ({((candidate.evidence_score - candidate.baseline_score) / candidate.baseline_score * 100):.1f}%)</p>
        </div>
    </div>
    """
    
    # PART 2: Everything AFTER the timeline plot
    html_after_plot = f"""
    <div style="padding: 0 20px 20px 20px;">
        <h2 style="color: {UAB_GREEN};">📁 Attached Evidence Projects ({len(candidate.attached_projects)})</h2>
        <p style="color: #718096; margin-bottom: 15px;">Projects contributing to prediction confidence</p>
    """
    
    # Render projects
    for project in candidate.attached_projects:
        impact_badge = get_impact_badge(project.impact_score)
        
        html_after_plot += f"""
        <div style="background: white; border-left: 5px solid {UAB_GREEN}; padding: 20px; margin: 15px 0; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <h3 style="margin: 0; color: {UAB_DARK_GREEN};">📁 {project.name}</h3>
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
        <div style="margin-top: 30px; padding: 20px; background: #F0F9F6; border-radius: 10px; border-left: 4px solid {UAB_GREEN};">
            <h3 style="color: {UAB_DARK_GREEN};">💡 Quick Stats</h3>
            <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px;">
    """
    
    html_after_plot += f"""
                <div style="text-align: center;">
                    <div style="font-size: 32px;">📁</div>
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{len(candidate.attached_projects)}</div>
                    <div style="color: #4A5568;">Projects</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-size: 32px;">👥</div>
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GOLD};">{candidate.cohort_count}</div>
                    <div style="color: #4A5568;">Cohorts</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-size: 32px;">🧪</div>
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{"✅" if candidate.has_validation_plan else "❌"}</div>
                    <div style="color: #4A5568;">Validation Plan</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-size: 32px;">💰</div>
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GOLD};">{"✅" if candidate.has_market_analysis else "❌"}</div>
                    <div style="color: #4A5568;">Market Analysis</div>
                </div>
    """
    
    html_after_plot += """
            </div>
        </div>
    </div>
    """
    
    return html_before_plot, timeline_fig, html_after_plot

def create_evidence_timeline(candidate):
    """Create plotly timeline showing evidence evolution - UAB colors"""
    
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
    
    # Add annotations for significant changes - UAB gold
    for i in range(1, len(timeline_data)):
        if 'impact' in timeline_data[i]:
            fig.add_annotation(
                x=timeline_data[i]['date'],
                y=timeline_data[i]['score'],
                text=timeline_data[i]['impact'],
                showarrow=True,
                arrowhead=2,
                arrowcolor=UAB_GOLD,
                bgcolor=UAB_GOLD,
                bordercolor=UAB_GOLD,
                font=dict(color='white', size=10),
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
        font=dict(family="Arial, sans-serif", color=UAB_DARK_GREEN),
        yaxis=dict(range=[0, 1])
    )
    
    return fig

# ========================================
# GRADIO INTERFACE
# ========================================

def create_interface():
    """Create the Gradio interface with UAB theme"""
    
    with gr.Blocks(title="SUD-PROMISE | UAB", theme=gr.themes.Soft()) as demo:
        
        # State management
        current_view = gr.State("dashboard")
        selected_category = gr.State(None)
        selected_candidate = gr.State(None)
        
        with gr.Column():
            # Header - UAB colors
            gr.Markdown(f"""
            <div style="text-align: center; padding: 20px; background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); border-radius: 15px; margin-bottom: 20px; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <h1 style="color: white; margin: 0;">🧬 SUD-PROMISE Drug Repositioning Assessment Platform</h1>
                <p style="color: white; opacity: 0.9; margin: 10px 0 0 0;">University of Alabama at Birmingham - Evidence-Based Drug Discovery for Substance Use Disorders</p>
            </div>
            """)
            
            # Navigation breadcrumb
            breadcrumb = gr.Markdown("")
            
            # Main content area
            with gr.Column() as main_content:
                
                # Dashboard view
                with gr.Column(visible=True) as dashboard_view:
                    dashboard_html = gr.HTML()
                    top_candidates_html = gr.HTML()
                    
                    gr.Markdown("### Select a SUD Category to Explore")
                    category_selector = gr.Dropdown(
                        choices=[],
                        label="SUD Categories",
                        interactive=True
                    )
                    view_category_btn = gr.Button("🔍 View Candidates", variant="primary", size="lg")
                
                # Category view
                with gr.Column(visible=False) as category_view:
                    category_html = gr.HTML()
                    
                    with gr.Row():
                        sort_dropdown = gr.Dropdown(
                            choices=["evidence_score", "recent", "name"],
                            value="evidence_score",
                            label="Sort by",
                            interactive=True
                        )
                    
                    gr.Markdown("### 🔍 Click a candidate below to view detailed dashboard")
                    candidate_selector = gr.Dropdown(
                        choices=[],
                        label="Select Candidate",
                        interactive=True,
                        info="Auto-opens detailed dashboard when selected"
                    )
                    
                    with gr.Row():
                        back_to_dashboard_btn = gr.Button("← Back to Dashboard", variant="secondary")
                        view_candidate_btn = gr.Button("📊 View Selected Dashboard", variant="primary", size="lg", visible=False)
                
                # Candidate detail view
                with gr.Column(visible=False) as candidate_view:
                    candidate_html_before = gr.HTML()
                    
                    # Timeline plot positioned between evidence evolution and projects
                    gr.Markdown("### 📈 Evidence Score Evolution Timeline")
                    timeline_plot = gr.Plot()
                    
                    gr.Markdown("---")
                    
                    candidate_html_after = gr.HTML()
                    
                    back_to_category_btn = gr.Button("← Back to Category", variant="secondary", size="lg")
        
        # ========================================
        # EVENT HANDLERS
        # ========================================
        
        def show_dashboard():
            """Show landing dashboard"""
            dash_html, top_html, categories = render_landing_dashboard()
            return (
                "🏠 Dashboard",  # breadcrumb
                dash_html,  # dashboard_html
                top_html,  # top_candidates_html
                gr.update(choices=categories, value=None),  # category_selector
                gr.update(visible=True),  # dashboard_view
                gr.update(visible=False),  # category_view
                gr.update(visible=False),  # candidate_view
                "dashboard"  # current_view
            )
        
        def show_category(category_selection):
            """Show category view"""
            if not category_selection:
                return [gr.update()] * 8
            
            cat_html, candidates, filtered_candidates = render_category_view(category_selection)
            category_name = category_selection.split(" (")[0]
            
            return (
                f"🏠 Dashboard > 🎯 {category_name}",  # breadcrumb
                cat_html,  # category_html
                gr.update(choices=candidates, value=None),  # candidate_selector
                gr.update(visible=False),  # dashboard_view
                gr.update(visible=True),  # category_view
                gr.update(visible=False),  # candidate_view
                category_selection,  # selected_category
                "category"  # current_view
            )
        
        def show_candidate(candidate_selection, category_selection):
            """Show candidate detail view"""
            if not candidate_selection:
                return [gr.update()] * 7
            
            html_before, timeline_fig, html_after = render_candidate_dashboard(candidate_selection, category_selection)
            candidate_name = candidate_selection.split(" (Score:")[0]
            category_name = category_selection.split(" (")[0] if category_selection else "Category"
            
            return (
                f"🏠 Dashboard > 🎯 {category_name} > 💊 {candidate_name}",  # breadcrumb
                html_before,  # candidate_html_before
                timeline_fig,  # timeline_plot
                html_after,   # candidate_html_after
                gr.update(visible=False),  # dashboard_view
                gr.update(visible=False),  # category_view
                gr.update(visible=True),  # candidate_view
            )
        
        def update_category_sort(category_selection, sort_by):
            """Update category view when sorting changes"""
            cat_html, candidates, filtered = render_category_view(category_selection, sort_by)
            return cat_html, gr.update(choices=candidates)
        
        # Connect events
        demo.load(
            fn=show_dashboard,
            outputs=[breadcrumb, dashboard_html, top_candidates_html, category_selector,
                    dashboard_view, category_view, candidate_view, current_view]
        )
        
        view_category_btn.click(
            fn=show_category,
            inputs=[category_selector],
            outputs=[breadcrumb, category_html, candidate_selector,
                    dashboard_view, category_view, candidate_view, selected_category, current_view]
        )
        
        sort_dropdown.change(
            fn=update_category_sort,
            inputs=[selected_category, sort_dropdown],
            outputs=[category_html, candidate_selector]
        )
        
        # Auto-open dashboard when candidate selected from dropdown
        candidate_selector.change(
            fn=show_candidate,
            inputs=[candidate_selector, selected_category],
            outputs=[breadcrumb, candidate_html_before, timeline_plot, candidate_html_after,
                    dashboard_view, category_view, candidate_view]
        )
        
        view_candidate_btn.click(
            fn=show_candidate,
            inputs=[candidate_selector, selected_category],
            outputs=[breadcrumb, candidate_html_before, timeline_plot, candidate_html_after,
                    dashboard_view, category_view, candidate_view]
        )
        
        back_to_dashboard_btn.click(
            fn=show_dashboard,
            outputs=[breadcrumb, dashboard_html, top_candidates_html, category_selector,
                    dashboard_view, category_view, candidate_view, current_view]
        )
        
        back_to_category_btn.click(
            fn=show_category,
            inputs=[selected_category],
            outputs=[breadcrumb, category_html, candidate_selector,
                    dashboard_view, category_view, candidate_view, selected_category, current_view]
        )
    
    return demo

# ========================================
# MAIN
# ========================================

if __name__ == "__main__":
    demo = create_interface()
    demo.launch()

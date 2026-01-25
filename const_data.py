"""
Data Constants - SUD Categories, Drug/Project Templates, Synthetic Data Settings
All data generation and template related constants
"""

# ========================================
# SUD CATEGORIES (Used in generate_synthetic_data)
# ========================================

SUD_CATEGORY_DESCRIPTIONS = {
    "Opioid-Related Disorders": "Addiction to opioids including prescription painkillers, heroin, and fentanyl",
    "Alcohol Use Disorder": "Problematic pattern of alcohol use leading to clinically significant impairment",
    "Stimulant Use Disorder": "Addiction to cocaine, methamphetamine, or prescription stimulants",
    "Cannabis Use Disorder": "Problematic cannabis use with withdrawal and tolerance symptoms",
    "Sedative/Hypnotic Disorder": "Dependence on benzodiazepines or other sedative medications",
    "Nicotine Use Disorder": "Tobacco/nicotine dependence and addiction",
}

SUD_CATEGORY_COLORS = {
    "Opioid-Related Disorders": "#1E6B52",
    "Alcohol Use Disorder": "#008C95",
    "Stimulant Use Disorder": "#5A9B7F",
    "Cannabis Use Disorder": "#16533E",
    "Sedative/Hypnotic Disorder": "#4A8B7A",
    "Nicotine Use Disorder": "#2C7A64",
}

DISEASE_SEARCH_CONFIG = {
    "Opioid-Related Disorders": {
        "search_terms": ["Opioid-Related Disorders"],
    },
    "Alcohol Use Disorder": {
        "search_terms": ["Alcoholism"],
    },
    "Stimulant Use Disorder": {
        "search_terms": ["Stimulant Use Disorder"],
    },
    "Cannabis Use Disorder": {
        "search_terms": ["Cannabis Use Disorder"],
    },
    "Nicotine Use Disorder": {
        "search_terms": ["Tobacco Use Disorder"],
    },
    "Sedative/Hypnotic Disorder": {
        "search_terms": ["Sedative/Hypnotic Disorder"],
    },
}

# ========================================
# DRUG TEMPLATES (Used in generate_synthetic_data)
# ========================================

DRUG_TEMPLATES = {
    "Opioid-Related Disorders": [
        ("Naltrexone", "Alcohol/opioid dependence", "μ-opioid receptor antagonist", "S5"),
        ("Buprenorphine", "Opioid dependence", "Partial μ-opioid agonist", "S6"),
        ("Lofexidine", "Hypertension", "α2-adrenergic agonist", "S4"),
        ("Gabapentin", "Epilepsy, neuropathic pain", "GABA analog", "S4"),
        ("Topiramate", "Epilepsy, migraine", "GABA modulator", "S3"),
        ("Ondansetron", "Nausea/vomiting", "5-HT3 antagonist", "S2"),
        ("Clonidine", "Hypertension", "α2-adrenergic agonist", "S3"),
        ("Memantine", "Alzheimer's disease", "NMDA antagonist", "S2"),
        ("Prazosin", "Hypertension", "α1-adrenergic antagonist", "S1"),
        ("Mirtazapine", "Depression", "α2-adrenergic antagonist", "S2"),
    ],
    "Alcohol Use Disorder": [
        ("Acamprosate", "Alcohol dependence", "NMDA antagonist", "S5"),
        ("Naltrexone", "Opioid dependence", "μ-opioid antagonist", "S5"),
        ("Baclofen", "Muscle spasticity", "GABA-B agonist", "S4"),
        ("Topiramate", "Epilepsy", "GABA modulator", "S4"),
        ("Varenicline", "Smoking cessation", "Nicotinic receptor agonist", "S3"),
        ("Ondansetron", "Nausea/vomiting", "5-HT3 antagonist", "S3"),
        ("Gabapentin", "Epilepsy", "GABA analog", "S2"),
        ("Zonisamide", "Epilepsy", "Carbonic anhydrase inhibitor", "S1"),
    ],
    "Stimulant Use Disorder": [
        ("Modafinil", "Narcolepsy", "Dopamine reuptake inhibitor", "S4"),
        ("Bupropion", "Depression", "NDRI", "S4"),
        ("Topiramate", "Epilepsy", "GABA modulator", "S3"),
        ("N-acetylcysteine", "Acetaminophen overdose", "Antioxidant", "S3"),
        ("Mirtazapine", "Depression", "α2-antagonist", "S2"),
        ("Atomoxetine", "ADHD", "NRI", "S1"),
    ],
    "Cannabis Use Disorder": [
        ("N-acetylcysteine", "Acetaminophen overdose", "Antioxidant", "S3"),
        ("Gabapentin", "Epilepsy", "GABA analog", "S3"),
        ("Zolpidem", "Insomnia", "GABA-A modulator", "S2"),
        ("Dronabinol", "Nausea/vomiting", "CB1 agonist", "S2"),
        ("Buspirone", "Anxiety", "5-HT1A agonist", "S1"),
    ],
    "Sedative/Hypnotic Disorder": [
        ("Flumazenil", "Benzodiazepine overdose", "GABA-A antagonist", "S3"),
        ("Gabapentin", "Epilepsy", "GABA analog", "S3"),
        ("Pregabalin", "Neuropathic pain", "GABA analog", "S2"),
    ],
    "Nicotine Use Disorder": [
        ("Varenicline", "Smoking cessation", "Nicotinic receptor agonist", "S6"),
        ("Cytisine", "Smoking cessation (EU)", "Nicotinic receptor agonist", "S5"),
        ("Bupropion", "Depression", "NDRI", "S5"),
        ("Naltrexone", "Opioid dependence", "μ-opioid antagonist", "S4"),
    ],
}

# ========================================
# PROJECT TEMPLATES (Used in generate_synthetic_data)
# ========================================

PROJECT_TEMPLATES = [
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
    },
    {
        "name": "Safety Concern Report FDA-2024",
        "type": "Safety Analysis",
        "size_range": (80, 250),
        "impact_range": (-0.15, -0.08),
        "summary": "Adverse event reports indicating potential safety concerns in target population"
    },
    {
        "name": "Negative RCT Results Study",
        "type": "Clinical Trial",
        "size_range": (120, 400),
        "impact_range": (-0.12, -0.06),
        "summary": "Randomized trial showing no significant efficacy compared to placebo"
    },
]

# ========================================
# METRICS & ANALYTICS (Used in analytics functions)
# ========================================

AVERAGE_TIME_TO_IND = 12.7  # months
CURRENT_IND_READY_COUNT = 2

TIME_BINS_LABELS = ['6-9', '9-12', '12-15', '15-18', '18-21', '21-24', '24+']
TIME_BINS_COUNTS = [10, 16, 25, 31, 39, 42, 49]

QUARTERLY_DELIVERIES = [
    {'quarter': 'Sep 2025', 'year': 2025, 'count': 0, 'cumulative': 0},
    {'quarter': 'Jan 2026', 'year': 2026, 'count': 0.5, 'cumulative': 0.5},
    {'quarter': 'Apr 2026', 'year': 2026, 'count': 0.3, 'cumulative': 0.8},
    {'quarter': 'Jun 2026', 'year': 2026, 'count': 1.2, 'cumulative': 2.0},
    {'quarter': 'Sep 2026', 'year': 2026, 'count': 0, 'cumulative': 2.0},
]
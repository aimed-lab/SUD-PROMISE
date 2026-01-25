"""
UI Constants - Colors, Styling, Display Settings
All visual/interface related constants
"""

# ========================================
# COLOR PALETTE (Used throughout pipeline)
# ========================================

UAB_GREEN = "#1E6B52"
UAB_DARK_GREEN = "#16533E"
UAB_LIGHT_GREEN = "#5A9B7F"
UAB_ACCENT_TEAL = "#008C95"
UAB_PALE_GREEN = "#E8F4F0"
UAB_WHITE = "#FFFFFF"

# ========================================
# STAGE CONFIGURATION (Used in generate_synthetic_data, badges, etc.)
# ========================================

STAGE_MAPPING = {
    "S0": (0.35, 0.50),
    "S1": (0.45, 0.58),
    "S2": (0.55, 0.68),
    "S3": (0.65, 0.78),
    "S4": (0.73, 0.85),
    "S5": (0.80, 0.92),
    "S6": (0.88, 0.95),
}

STAGE_NAMES = {
    "S0": "Ideation",
    "S1": "Intake",
    "S2": "Feasibility",
    "S3": "In Silico",
    "S4": "Wet Lab",
    "S5": "IND-Ready",
    "S6": "Handoff"
}

STAGE_COLORS = {
    "S0": "#94A3B8",
    "S1": "#A8CABD",
    "S2": UAB_LIGHT_GREEN,
    "S3": UAB_GREEN,
    "S4": UAB_GREEN,
    "S5": UAB_DARK_GREEN,
    "S6": UAB_DARK_GREEN,
}

# ========================================
# SCORE THRESHOLDS (Used in charts and interpretation)
# ========================================

SCORE_THRESHOLD_HIGH = 0.7
SCORE_THRESHOLD_MODERATE = 0.5

# ========================================
# DISPLAY SETTINGS (Used in UI rendering)
# ========================================

MAX_CATEGORY_CANDIDATES_DISPLAY = 12
MAX_TOP_CANDIDATES_DISPLAY = 5
SMILES_DISPLAY_LENGTH = 60
MAX_TARGET_DISPLAY = 5

# Chart heights
CHART_HEIGHT_SMALL = 190
CHART_HEIGHT_MEDIUM = 300
CHART_HEIGHT_LARGE = 400

# ========================================
# HELPER FUNCTIONS
# ========================================

def get_status_badge(stage):
    """Generate HTML badge for stage status"""
    colors = STAGE_COLORS
    names = STAGE_NAMES
    
    color = colors.get(stage, "#CBD5E0")
    name = names.get(stage, stage)
    
    return f'<span style="background: {color}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white; font-family: \'Times New Roman\', Times, serif;">{stage} - {name}</span>'


def get_score_type_badge(score_type):
    """Get badge for score type (Real ML/DL or Synthetic/AI)"""
    if score_type == "Real":
        return f'<span style="background: {UAB_DARK_GREEN}; padding: 4px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif; font-size: 11px;">🤖 Real ML/DL</span>'
    else:
        return f'<span style="background: #FFA500; padding: 4px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif; font-size: 11px;">🔮 Synthetic/AI</span>'


def get_evidence_stars(score):
    """Convert score to star rating"""
    stars = int(score * 5)
    return "⭐" * stars


def get_impact_badge(impact):
    """Generate HTML badge for project impact"""
    if impact >= 0.10:
        return f'<span style="background: {UAB_DARK_GREEN}; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">High +{impact:.2f}</span>'
    elif impact >= 0.05:
        return f'<span style="background: {UAB_GREEN}; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Moderate +{impact:.2f}</span>'
    elif impact > 0:
        return f'<span style="background: #A0AEC0; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Low +{impact:.2f}</span>'
    elif impact <= -0.10:
        return f'<span style="background: #DC2626; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">High Risk {impact:.2f}</span>'
    elif impact <= -0.05:
        return f'<span style="background: #EF4444; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Moderate Risk {impact:.2f}</span>'
    else:
        return f'<span style="background: #F87171; padding: 3px 10px; border-radius: 10px; color: white; font-family: \'Times New Roman\', Times, serif;">Low Risk {impact:.2f}</span>'
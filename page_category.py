"""
Category Page Module
Extracted from sud_promise_uab_theme.py
"""

# Import UI constants
from const_ui import (
    UAB_GREEN, UAB_DARK_GREEN,
    MAX_TOP_CANDIDATES_DISPLAY
)

from const_data import SUD_CATEGORY_DESCRIPTIONS


def render_landing_dashboard():
    from sud_promise_uab_theme import CANDIDATES, SUD_CATEGORIES
    
    top_candidates = sorted(CANDIDATES, key=lambda c: c.evidence_score, reverse=True)[:MAX_TOP_CANDIDATES_DISPLAY]
    
    html = ""
    
    category_choices = [f"{cat.name} ({cat.candidate_count} candidates)" for cat in SUD_CATEGORIES]
    
    return html, category_choices, top_candidates

def render_category_view(category_selection, sort_by="Evidence Score"):
    from sud_promise_uab_theme import CANDIDATES, SUD_CATEGORIES
    
    if not category_selection:
        return "<p>Please select a category</p>", [], []
    
    category_name = category_selection.split(" (")[0].strip()
    category = next((c for c in SUD_CATEGORIES if c.name == category_name), None)
    if not category:
        return "<p>Category not found</p>", [], []
    
    filtered = [c for c in CANDIDATES if c.target_sud_subtype == category_name]
    
    if sort_by == "Evidence Score" or sort_by == "evidence_score":
        filtered.sort(key=lambda c: c.evidence_score, reverse=True)
    elif sort_by == "Recent" or sort_by == "recent":
        filtered.sort(key=lambda c: c.last_updated, reverse=True)
    elif sort_by == "Name" or sort_by == "name":
        filtered.sort(key=lambda c: c.drug_name)
    elif sort_by == "Stage":
        filtered.sort(key=lambda c: c.stage, reverse=True)
    
    
    candidate_choices = [f"{c.drug_name} (Score: {c.evidence_score:.2f})" for c in filtered]

    html = ""
    
    return html, candidate_choices, filtered
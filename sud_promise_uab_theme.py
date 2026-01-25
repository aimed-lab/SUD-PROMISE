"""
🧬 SUD-PROMISE Drug Repositioning Assessment Platform 
UAB Color Theme: Forest Green/Blue (#1E6B52) & White
Font: Times New Roman
Stage System: S0 (Ideation) → S1 (Intake) → S2 (Feasibility) → S3 (In Silico) → S4 (Wet Lab) → S5 (IND-Ready) → S6 (Handoff)

VERSION: 3.4.0 - Clickable Card Overlay
- NEW: Entire card is clickable with invisible overlay button
- CLEAN: No visible "View Details" button - card itself is the button
- FIXED: Timeline shows ensemble baseline + all project impacts
- FIXED: Real ML/DL scores start from ensemble, then show evidence progression

UPDATED THEME:
- Changed Gradio default primary button color ("View Candidates") to UAB GREEN
"""

import gradio as gr
import pandas as pd
import plotly.graph_objects as go
from datetime import datetime, timedelta
import random
import json
import numpy as np
from pathlib import Path
import warnings

# Import custom utilities
from func_drug import find_drug_in_database
from func_models import initialize_ml_system, predict_with_ml_models, load_ml_models
from module_dashboard import render_executive_dashboard
from page_candidate import render_candidate_dashboard
from page_category import render_landing_dashboard, render_category_view

# Import UI constants and helpers
from const_ui import (
    # Colors
    UAB_GREEN, UAB_DARK_GREEN, UAB_LIGHT_GREEN, UAB_ACCENT_TEAL, UAB_PALE_GREEN, UAB_WHITE,
    # Stage configuration
    STAGE_MAPPING, STAGE_NAMES, STAGE_COLORS,
    # Thresholds
    SCORE_THRESHOLD_HIGH, SCORE_THRESHOLD_MODERATE,
    # Display limits
    MAX_CATEGORY_CANDIDATES_DISPLAY, MAX_TOP_CANDIDATES_DISPLAY,
    SMILES_DISPLAY_LENGTH, MAX_TARGET_DISPLAY,
    CHART_HEIGHT_SMALL, CHART_HEIGHT_MEDIUM, CHART_HEIGHT_LARGE,
    # Helper functions
    get_status_badge, get_score_type_badge, get_evidence_stars, get_impact_badge
)

# Import data constants
from const_data import (
    # SUD Categories
    SUD_CATEGORY_DESCRIPTIONS, SUD_CATEGORY_COLORS, DISEASE_SEARCH_CONFIG,
    # Templates
    DRUG_TEMPLATES, PROJECT_TEMPLATES,
    # Metrics
    AVERAGE_TIME_TO_IND, CURRENT_IND_READY_COUNT,
    TIME_BINS_LABELS, TIME_BINS_COUNTS, QUARTERLY_DELIVERIES,
)

# Import config constants
from const_config import (
    APP_NAME, APP_VERSION, INSTITUTION,
    DATA_DIR, MODEL_DIR, DRUGS_FILE, DISEASES_FILE,
)

# Import data generation module
from data_generator import (
    Project, DrugCandidate, SUDCategory,
    setup_disease_mapping,
    generate_synthetic_data
)

warnings.filterwarnings('ignore')

uab_blue_theme = gr.themes.Soft(
    primary_hue="blue",
    secondary_hue="slate",
).set(
    # Primary buttons
    button_primary_background_fill=UAB_GREEN,
    button_primary_background_fill_hover=UAB_DARK_GREEN,
    button_primary_text_color="white",

    # Secondary buttons
    button_secondary_border_color=UAB_GREEN,
    button_secondary_text_color=UAB_GREEN,

    # Inputs focus ring / highlights
    input_border_color_focus=UAB_GREEN,
    block_label_background_fill=UAB_GREEN,
    block_label_text_color="white",
    link_text_color=UAB_GREEN,
)

# ========================================
# ML/DL MODEL LOADING
# ========================================

system = initialize_ml_system()

# Extract for backwards compatibility
MODELS_AVAILABLE = system['models_available']
ml_components = system['ml_components']
drugs_df = system['drugs_df']
diseases_df = system['diseases_df']

# ========================================
# DISEASE ID MAPPING
# ========================================

DISEASE_ID_MAPPING = setup_disease_mapping(diseases_df)

# ========================================
# GENERATE DATA
# ========================================

SUD_CATEGORIES, CANDIDATES = generate_synthetic_data(
    DISEASE_ID_MAPPING, 
    MODELS_AVAILABLE, 
    ml_components, 
    drugs_df
)

# ========================================
# HELPER FUNCTIONS
# ========================================

def format_date_ago(date):
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

def render_candidate_card(candidate, rank, border_color=None):
    stars = get_evidence_stars(candidate.evidence_score)
    status_badge = get_status_badge(candidate.stage)
    score_type_badge = get_score_type_badge(candidate.score_type)
    border = UAB_GREEN
    market_analysis_num = 1 if hasattr(candidate, 'has_market_analysis') and candidate.has_market_analysis else 0

    html = f"""
    <div style="background: white; border: 2px solid {border}; padding: 15px 20px; margin: 10px 20px;
                border-radius: 12px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                font-family: 'Times New Roman', Times, serif; cursor: pointer; transition: all 0.3s ease;">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px;">
            <div style="display: flex; align-items: center; gap: 12px;">
                <span style="background: #EDF2F7; padding: 4px 10px; border-radius: 20px;
                             font-size: 14px; color: {UAB_DARK_GREEN}; font-weight: bold;">
                    {rank}
                </span>
                <span style="font-size: 24px; font-weight: bold; color: {UAB_DARK_GREEN};">
                    {candidate.drug_name}
                </span>
                <span style="font-size: 16px; color: {UAB_GREEN}; font-weight: bold;">
                    {stars} {candidate.evidence_score:.2f}
                </span>
                {score_type_badge}
            </div>
            {status_badge}
        </div>

        <p style="margin: 8px 0; color: #4A5568;"><b>Current Use:</b> {candidate.current_indication}</p>
        <p style="margin: 8px 0; color: #4A5568;"><b>Mechanism:</b> {candidate.mechanism}</p>

        <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin: 12px 0;">
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{len(candidate.attached_projects)}</div>
                <div style="font-size: 12px; color: white;">Projects</div>
            </div>
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{candidate.cohort_count}</div>
                <div style="font-size: 12px; color: white;">Cohorts</div>
            </div>
            <div style="text-align: center; padding: 10px; background: {UAB_GREEN}; border-radius: 8px;">
                <div style="font-weight: bold; color: white; font-size: 20px;">{market_analysis_num}</div>
                <div style="font-size: 12px; color: white;">Market Analysis</div>
            </div>
        </div>

        <div style="background: {UAB_PALE_GREEN}; padding: 12px; border-radius: 8px; margin: 12px 0;
                    border-left: 3px solid {UAB_GREEN};">
            <p style="margin: 5px 0; color: #4A5568; font-size: 13px;">
                <b>Team Members:</b> {candidate.team_members} |
                <b>Data Produced:</b> {candidate.data_produced} |
                <b>Publications:</b> {candidate.publications}
            </p>
            <p style="margin: 5px 0; color: #4A5568; font-size: 13px;">
                <b>Tools Used:</b> {candidate.tools_used} |
                <b>Data Governance:</b> {candidate.data_governance} |
                <b>Training:</b> {candidate.training_participated}
            </p>
        </div>

        <p style="margin: 8px 0 0 0; font-size: 13px; color: #A0AEC0;">
            Last updated: {format_date_ago(candidate.last_updated)}
        </p>
    </div>
    """
    return html

# ========================================
# GRADIO INTERFACE
# ========================================

def create_interface():
    CSS = f"""
    * {{
        font-family: 'Times New Roman', Times, serif !important;
    }}

    .gr-button {{
        font-family: 'Times New Roman', Times, serif !important;
    }}

    /* Dropdown highlight: UAB GREEN */
    .choices__item--selectable.is-highlighted {{
        background-color: rgba(30, 58, 138, 0.15) !important;
        color: {UAB_DARK_GREEN} !important;
    }}

    /* Force Gradio primary button color to UAB GREEN (HF Spaces safe) */
    .gr-button-primary,
    .gr-button-primary button {{
        background: {UAB_GREEN} !important;
        border-color: {UAB_GREEN} !important;
        color: white !important;
        font-weight: bold !important;
    }}

    .gr-button-primary:hover,
    .gr-button-primary button:hover {{
        background: {UAB_DARK_GREEN} !important;
        border-color: {UAB_DARK_GREEN} !important;
    }}

    /* Card container with relative positioning */
    .card-wrap {{
        position: relative;
        width: 100%;
        margin: 10px 0;
    }}

    /* Invisible overlay button covering entire card */
    .overlay-btn {{
        position: absolute !important;
        inset: 0 !important;
        z-index: 10 !important;
        opacity: 0 !important;
        cursor: pointer !important;
    }}

    .overlay-btn button {{
        width: 100% !important;
        height: 100% !important;
        border-radius: 12px !important;
    }}

    /* Hover effect on card when hovering over invisible button */
    .card-wrap:hover > div:first-child > div {{
        box-shadow: 0 4px 12px rgba(30, 58, 138, 0.22) !important;
        transform: translateY(-2px);
    }}
    """

    with gr.Blocks(
        title="SUD-PROMISE | UAB",
        css=CSS,
        theme=uab_blue_theme,  #  APPLY THEME HERE
    ) as demo:

        current_view = gr.State("dashboard")
        selected_category = gr.State(None)
        selected_candidate = gr.State(None)

        with gr.Column():
            gr.HTML(f"""
            <div style="text-align: center; padding: 40px 20px;
                        background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%);
                        border-radius: 15px; margin-bottom: 20px;
                        box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
                <div style="font-size: 60px; margin-bottom: 20px;">🧬</div>
                <h1 style="color: white; margin: 0; font-size: 42px; font-weight: bold;">
                    {APP_NAME}
                </h1>
                <p style="color: white; opacity: 0.95; margin: 15px 0 0 0; font-size: 18px;">
                    {INSTITUTION} - Evidence-Based Drug Discovery for Substance Use Disorders
                </p>
                <p style="color: white; opacity: 0.85; margin: 10px 0 0 0; font-size: 14px;">
                    🤖 v8: Clickable Card Overlay + UAB Green Primary Buttons
                </p>
            </div>
            """)

            breadcrumb = gr.Markdown("")

            with gr.Column() as main_content:

                # DASHBOARD VIEW
                with gr.Column(visible=True) as dashboard_view:
                    exec_html = gr.HTML()

                    gr.Markdown("---")

                    with gr.Row():
                        with gr.Column():
                            gr.Markdown("### Average Time to IND-Ready")
                            with gr.Row():
                                with gr.Column(scale=1):
                                    avg_time_display = gr.Markdown("")
                                with gr.Column(scale=2):
                                    time_dist_plot = gr.Plot()

                        with gr.Column():
                            gr.Markdown("### Quarterly Deliveries")
                            with gr.Row():
                                with gr.Column(scale=1):
                                    quarterly_display = gr.Markdown("")
                                with gr.Column(scale=2):
                                    quarterly_plot = gr.Plot()               

                    with gr.Row():
                        with gr.Column():
                            gr.Markdown("### Pipeline Progression - Stage Distribution Over Time")
                            pipeline_plot = gr.Plot()

                        with gr.Column():
                            gr.Markdown("### Active Portfolio Breakdown")
                            portfolio_plot = gr.Plot()

                    gr.Markdown("---")

                    dashboard_html = gr.HTML()
                    gr.Markdown("### Top Candidates by Evidence Score")
                    gr.Markdown("*Click on any card to view detailed information*")

                    top_candidate_components = []
                    for i in range(MAX_TOP_CANDIDATES_DISPLAY):
                        with gr.Column(elem_classes=["card-wrap"]):
                            top_html = gr.HTML()
                            top_btn = gr.Button("View Details", elem_classes=["overlay-btn"])
                            top_candidate_components.append((top_html, top_btn))

                    gr.Markdown("---")
                    gr.Markdown("### Select a SUD Category to Explore")

                    category_selector = gr.Dropdown(
                        choices=[],
                        label="",
                        show_label=False,
                        interactive=True
                    )

                    #  This is the button in your screenshot: now it becomes UAB GREEN
                    view_category_btn = gr.Button("View Candidates", variant="primary", size="lg")

                # CATEGORY VIEW
                with gr.Column(visible=False) as category_view:
                    category_html = gr.HTML()

                    gr.Markdown("**Sort by**")
                    with gr.Row():
                        sort_dropdown = gr.Dropdown(
                            choices=["Evidence Score", "Recent", "Name", "Stage"],
                            value="Evidence Score",
                            label="",
                            show_label=False,
                            interactive=True
                        )

                    candidate_selector = gr.Dropdown(
                        choices=[],
                        label="",
                        show_label=False,
                        interactive=True,
                        visible=False
                    )

                    gr.Markdown("*Click on any card to view detailed information*")

                    category_candidate_components = []
                    category_candidate_buttons = []
                    for i in range(MAX_CATEGORY_CANDIDATES_DISPLAY):
                        with gr.Column(visible=False, elem_classes=["card-wrap"]) as cat_card_col:
                            cat_html = gr.HTML()
                            cat_btn = gr.Button("View Details", elem_classes=["overlay-btn"])
                            category_candidate_components.append((cat_card_col, cat_html))
                            category_candidate_buttons.append(cat_btn)

                    with gr.Row():
                        back_to_dashboard_btn = gr.Button("← Back to Dashboard", variant="secondary")

                # CANDIDATE DETAIL VIEW
                with gr.Column(visible=False) as candidate_view:
                    candidate_html_before = gr.HTML()

                    with gr.Row():
                        with gr.Column():
                            gr.Markdown("### Evidence Score Gauge")
                            gauge_plot = gr.Plot()
                        with gr.Column():
                            gr.Markdown("### Model Score Comparison")
                            model_comparison_plot = gr.Plot()

                    gr.Markdown("### Evidence Score Evolution Timeline")
                    timeline_plot = gr.Plot()
                    gr.Markdown("---")
                    candidate_html_after = gr.HTML()
                    back_to_category_btn = gr.Button("← Back to Category", variant="secondary", size="lg")

        # ==========================================================
        # EVENT HANDLERS
        # ==========================================================

        def show_dashboard():
            exec_html_content, time_fig, quart_fig, pipe_fig, port_fig, avg_time, ind_ready = render_executive_dashboard(
                candidates=CANDIDATES,
                categories=SUD_CATEGORIES,
                models_available=MODELS_AVAILABLE
            )

            avg_time_md = f"""
            <div style="text-align: center; padding: 15px; background: white;
                        border: 2px solid {UAB_GREEN}; border-radius: 10px;">
                <div style="font-size: 36px; font-weight: bold; color: {UAB_GREEN};">{avg_time:.1f}</div>
                <div style="font-size: 14px; color: {UAB_DARK_GREEN};">months</div>
            </div>
            """

            quarterly_md = f"""
            <div style="text-align: center; padding: 15px; background: white;
                        border: 2px solid {UAB_GREEN}; border-radius: 10px;">
                <div style="font-size: 36px; font-weight: bold; color: {UAB_GREEN};">{ind_ready}</div>
                <div style="font-size: 14px; color: {UAB_DARK_GREEN};">IND-Ready</div>
            </div>
            """

            dash_html, categories, top_5 = render_landing_dashboard()

            top_outputs = [
                exec_html_content,
                avg_time_md,
                time_fig,
                quarterly_md,
                quart_fig,
                pipe_fig,
                port_fig,
                dash_html,
                gr.update(choices=categories, value=None)
            ]

            for i, (html_component, btn_component) in enumerate(top_candidate_components):
                if i < len(top_5):
                    candidate = top_5[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    top_outputs.append(card_html)
                else:
                    top_outputs.append("")

            for col, html in category_candidate_components:
                top_outputs.extend([gr.update(visible=False), ""])

            top_outputs.extend([
                gr.update(visible=True),
                gr.update(visible=False),
                gr.update(visible=False),
                "dashboard",
                "📊 Dashboard"
            ])

            return top_outputs

        def show_category(category_selection):
            if not category_selection:
                return [gr.update()] * (2 + len(category_candidate_components) * 2 + 6)

            cat_html, candidates, filtered = render_category_view(category_selection)
            category_name = category_selection.split(" (")[0]

            outputs = [
                cat_html,
                gr.update(choices=candidates, value=None),
            ]

            for i, (col, html) in enumerate(category_candidate_components):
                if i < len(filtered):
                    candidate = filtered[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    outputs.extend([gr.update(visible=True), card_html])
                else:
                    outputs.extend([gr.update(visible=False), ""])

            outputs.extend([
                gr.update(visible=False),
                gr.update(visible=True),
                gr.update(visible=False),
                category_selection,
                "category",
                f"📊 Dashboard > {category_name}",
            ])
            return outputs

        def show_candidate(candidate_selection, category_selection):
            if not candidate_selection:
                return [gr.update()] * 10

            html_before, gauge_fig, model_fig, timeline_fig, html_after = render_candidate_dashboard(
                candidate_selection, category_selection
            )
            candidate_name = candidate_selection.split(" (Score:")[0]
            category_name = category_selection.split(" (")[0] if category_selection else "Category"

            return (
                f"📊 Dashboard > {category_name} > {candidate_name}",
                html_before,
                gauge_fig,
                model_fig,
                timeline_fig,
                html_after,
                gr.update(visible=False),
                gr.update(visible=False),
                gr.update(visible=True),
                category_selection,
            )

        def update_category_sort(category_selection, sort_by):
            if not category_selection:
                return [gr.update()] * (2 + len(category_candidate_components) * 2)

            cat_html, candidates, filtered = render_category_view(category_selection, sort_by)
            outputs = [cat_html, gr.update(choices=candidates)]

            for i, (col, html) in enumerate(category_candidate_components):
                if i < len(filtered):
                    candidate = filtered[i]
                    card_html = render_candidate_card(candidate, f"#{i+1}")
                    outputs.extend([gr.update(visible=True), card_html])
                else:
                    outputs.extend([gr.update(visible=False), ""])

            return outputs

        def make_view_candidate_handler(candidate_index, is_top=False):
            def handler():
                if is_top:
                    top_5 = sorted(CANDIDATES, key=lambda c: c.evidence_score, reverse=True)[:MAX_TOP_CANDIDATES_DISPLAY]
                    if candidate_index >= len(top_5):
                        return [gr.update()] * 10
                    candidate = top_5[candidate_index]
                else:
                    return [gr.update()] * 10

                category_str = f"{candidate.target_sud_subtype} ({sum(1 for c in CANDIDATES if c.target_sud_subtype == candidate.target_sud_subtype)} candidates)"
                candidate_str = f"{candidate.drug_name} (Score: {candidate.evidence_score:.2f})"

                html_before, gauge_fig, model_fig, timeline_fig, html_after = render_candidate_dashboard(candidate_str, category_str)
                return (
                    f"📊 Dashboard > {candidate.target_sud_subtype} > {candidate.drug_name}",
                    html_before,
                    gauge_fig,
                    model_fig,
                    timeline_fig,
                    html_after,
                    gr.update(visible=False),
                    gr.update(visible=False),
                    gr.update(visible=True),
                    category_str,
                )
            return handler

        def make_category_candidate_handler(candidate_index):
            def handler(category_selection):
                if not category_selection:
                    return [gr.update()] * 10

                category_name = category_selection.split(" (")[0].strip()
                filtered = [c for c in CANDIDATES if c.target_sud_subtype == category_name]
                filtered.sort(key=lambda c: c.evidence_score, reverse=True)

                if candidate_index >= len(filtered):
                    return [gr.update()] * 10

                candidate = filtered[candidate_index]
                candidate_str = f"{candidate.drug_name} (Score: {candidate.evidence_score:.2f})"

                html_before, gauge_fig, model_fig, timeline_fig, html_after = render_candidate_dashboard(candidate_str, category_selection)
                return (
                    f"📊 Dashboard > {category_name} > {candidate.drug_name}",
                    html_before,
                    gauge_fig,
                    model_fig,
                    timeline_fig,
                    html_after,
                    gr.update(visible=False),
                    gr.update(visible=False),
                    gr.update(visible=True),
                    category_selection,
                )
            return handler

        # ==========================================================
        # CONNECT EVENTS
        # ==========================================================

        demo.load(
            fn=show_dashboard,
            outputs=[
                exec_html, avg_time_display, time_dist_plot, quarterly_display, quarterly_plot, pipeline_plot, portfolio_plot,
                dashboard_html, category_selector
            ] +
            [comp for pair in top_candidate_components for comp in [pair[0]]] +
            [comp for pair in category_candidate_components for comp in [pair[0], pair[1]]] +
            [dashboard_view, category_view, candidate_view, current_view, breadcrumb]
        )

        view_category_btn.click(
            fn=show_category,
            inputs=[category_selector],
            outputs=[category_html, candidate_selector] +
                    [comp for pair in category_candidate_components for comp in [pair[0], pair[1]]] +
                    [dashboard_view, category_view, candidate_view, selected_category, current_view, breadcrumb]
        )

        sort_dropdown.change(
            fn=update_category_sort,
            inputs=[selected_category, sort_dropdown],
            outputs=[category_html, candidate_selector] +
                    [comp for pair in category_candidate_components for comp in [pair[0], pair[1]]]
        )

        for i, (html, btn) in enumerate(top_candidate_components):
            btn.click(
                fn=make_view_candidate_handler(i, is_top=True),
                outputs=[
                    breadcrumb, candidate_html_before, gauge_plot, model_comparison_plot, timeline_plot, candidate_html_after,
                    dashboard_view, category_view, candidate_view, selected_category
                ]
            )

        for i, btn in enumerate(category_candidate_buttons):
            btn.click(
                fn=make_category_candidate_handler(i),
                inputs=[selected_category],
                outputs=[
                    breadcrumb, candidate_html_before, gauge_plot, model_comparison_plot, timeline_plot, candidate_html_after,
                    dashboard_view, category_view, candidate_view, selected_category
                ]
            )

        candidate_selector.change(
            fn=show_candidate,
            inputs=[candidate_selector, selected_category],
            outputs=[
                breadcrumb, candidate_html_before, gauge_plot, model_comparison_plot, timeline_plot, candidate_html_after,
                dashboard_view, category_view, candidate_view, selected_category
            ]
        )

        back_to_dashboard_btn.click(
            fn=show_dashboard,
            outputs=[
                exec_html, avg_time_display, time_dist_plot, quarterly_display, quarterly_plot, pipeline_plot, portfolio_plot,
                dashboard_html, category_selector
            ] +
            [comp for pair in top_candidate_components for comp in [pair[0]]] +
            [comp for pair in category_candidate_components for comp in [pair[0], pair[1]]] +
            [dashboard_view, category_view, candidate_view, current_view, breadcrumb]
        )

        back_to_category_btn.click(
            fn=show_category,
            inputs=[selected_category],
            outputs=[category_html, candidate_selector] +
                    [comp for pair in category_candidate_components for comp in [pair[0], pair[1]]] +
                    [dashboard_view, category_view, candidate_view, selected_category, current_view, breadcrumb]
        )

    return demo

# ========================================
# MAIN
# ========================================

if __name__ == "__main__":
    print("\n" + "="*60)
    print(f"🧬 {APP_NAME}")
    print(f"Version: 3.4.0 - Clickable Card Overlay + UAB Green Primary Buttons")
    print("="*60)
    print(f"ML/DL Models Status: {' LOADED' if MODELS_AVAILABLE else '  NOT AVAILABLE'}")
    if MODELS_AVAILABLE:
        print(f"Available Diseases: {len(ml_components['le_disease'].classes_)}")
        print(f"Available Targets: {len(ml_components['mlb'].classes_)}")
    if drugs_df is not None:
        print(f"Drugs Database:  {len(drugs_df)} drugs loaded")
    if diseases_df is not None:
        print(f"Diseases Database:  {len(diseases_df)} diseases loaded")
    print(f"Disease ID Mappings: {len(DISEASE_ID_MAPPING)} SUD categories mapped")
    print("="*60 + "\n")

    demo = create_interface()
    demo.launch()
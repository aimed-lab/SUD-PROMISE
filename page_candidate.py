"""
Candidate Detail Page Module
Extracted from sud_promise_uab_theme.py
"""

import plotly.graph_objects as go
from datetime import datetime, timedelta

# Import UI constants
from const_ui import (
    UAB_GREEN, UAB_DARK_GREEN, UAB_LIGHT_GREEN, UAB_ACCENT_TEAL, UAB_PALE_GREEN,
    SCORE_THRESHOLD_HIGH, SCORE_THRESHOLD_MODERATE,
    SMILES_DISPLAY_LENGTH, MAX_TARGET_DISPLAY,
    get_status_badge, get_score_type_badge, get_evidence_stars, get_impact_badge
)

from const_config import INSTITUTION


def create_score_gauge(score: float, score_type: str, model_scores: dict = None) -> go.Figure:
    """Create a gauge chart for the evidence score"""
    
    # Determine color based on score
    if score >= SCORE_THRESHOLD_HIGH:
        color = UAB_DARK_GREEN
    elif score >= SCORE_THRESHOLD_MODERATE:
        color = UAB_GREEN
    else:
        color = UAB_LIGHT_GREEN
    
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=score,
        domain={'x': [0, 1], 'y': [0, 1]},
        title={'text': f"Evidence Score ({score_type})", 
               'font': {'size': 16, 'family': "Times New Roman, serif", 'color': UAB_DARK_GREEN}},
        number={'font': {'size': 40, 'family': "Times New Roman, serif"}},
        gauge={
            'axis': {'range': [0, 1], 'tickwidth': 1, 'tickcolor': UAB_DARK_GREEN},
            'bar': {'color': color},
            'bgcolor': "white",
            'borderwidth': 2,
            'bordercolor': UAB_DARK_GREEN,
            'steps': [
                {'range': [0, 0.5], 'color': '#FFE5E5'},
                {'range': [0.5, 0.7], 'color': '#FFF4E5'},
                {'range': [0.7, 1], 'color': UAB_PALE_GREEN}
            ],
            'threshold': {
                'line': {'color': "red", 'width': 4},
                'thickness': 0.75,
                'value': 0.7
            }
        }
    ))
    
    fig.update_layout(
        height=300,
        margin=dict(l=20, r=20, t=60, b=20),
        paper_bgcolor='white',
        font={'family': 'Times New Roman, serif', 'color': UAB_DARK_GREEN}
    )
    
    return fig

def create_model_comparison_chart(model_scores: dict, score_type: str) -> go.Figure:
    """Create a bar chart comparing different model scores"""
    
    if not model_scores or score_type != "Real":
        # Return empty figure for synthetic scores
        fig = go.Figure()
        fig.update_layout(
            title="Model Scores (Not Available for Synthetic Predictions)",
            height=300,
            paper_bgcolor='white',
            font={'family': 'Times New Roman, serif', 'color': UAB_DARK_GREEN}
        )
        return fig
    
    models = list(model_scores.keys())
    scores = list(model_scores.values())
    
    colors = [UAB_GREEN if s >= SCORE_THRESHOLD_HIGH else UAB_LIGHT_GREEN if s >= SCORE_THRESHOLD_MODERATE else '#FFA500' for s in scores]
    
    fig = go.Figure(data=[
        go.Bar(
            x=models,
            y=scores,
            marker_color=colors,
            text=[f'{s:.3f}' for s in scores],
            textposition='outside',
            textfont=dict(size=12, family="Times New Roman, serif", color=UAB_DARK_GREEN)
        )
    ])
    
    fig.update_layout(
        title="Individual Model Predictions",
        xaxis_title="",
        yaxis_title="Prediction Score",
        height=300,
        margin=dict(l=40, r=20, t=50, b=80),
        plot_bgcolor=UAB_PALE_GREEN,
        paper_bgcolor='white',
        font={'family': 'Times New Roman, serif', 'color': UAB_DARK_GREEN},
        yaxis=dict(range=[0, 1]),
        showlegend=False
    )
    
    fig.add_hline(y=0.5, line_dash="dash", line_color="orange", annotation_text="Threshold 0.5")
    fig.add_hline(y=SCORE_THRESHOLD_HIGH, line_dash="dash", line_color=UAB_DARK_GREEN, annotation_text="Threshold 0.7")
    
    return fig

def create_evidence_timeline(candidate):
    """Create timeline showing baseline + all project impacts"""
    timeline_data = []
    
    # Start with baseline (either ML/DL ensemble or random baseline)
    start_date = candidate.attached_projects[0].added_date - timedelta(days=30) if candidate.attached_projects else datetime.now() - timedelta(days=365)
    
    if candidate.score_type == "Real":
        timeline_data.append({
            'date': start_date,
            'score': candidate.baseline_score,
            'label': f'ML/DL Ensemble Baseline: {candidate.baseline_score:.3f}',
            'color': UAB_DARK_GREEN,
            'impact': 0
        })
    else:
        timeline_data.append({
            'date': start_date,
            'score': candidate.baseline_score,
            'label': 'Baseline (No Evidence)',
            'color': '#A0AEC0',
            'impact': 0
        })
    
    # Add each project's impact
    cumulative_score = candidate.baseline_score
    for project in candidate.attached_projects:
        cumulative_score += project.impact_score
        cumulative_score = max(0.20, min(0.95, cumulative_score))
        
        marker_color = UAB_GREEN if project.impact_score > 0 else '#DC2626'
        
        timeline_data.append({
            'date': project.added_date,
            'score': cumulative_score,
            'label': f'+ {project.name}' if project.impact_score > 0 else f'- {project.name}',
            'color': marker_color,
            'impact': project.impact_score,
            'impact_text': f'+{project.impact_score:.2f}' if project.impact_score > 0 else f'{project.impact_score:.2f}'
        })
    
    fig = go.Figure()
    
    # Draw line segments
    for i in range(len(timeline_data) - 1):
        segment_color = timeline_data[i+1]['color']
        
        fig.add_trace(go.Scatter(
            x=[timeline_data[i]['date'], timeline_data[i+1]['date']],
            y=[timeline_data[i]['score'], timeline_data[i+1]['score']],
            mode='lines',
            line=dict(color=segment_color, width=3),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Draw markers
    fig.add_trace(go.Scatter(
        x=[d['date'] for d in timeline_data],
        y=[d['score'] for d in timeline_data],
        mode='markers',
        marker=dict(
            size=12,
            color=[d['color'] for d in timeline_data],
            line=dict(width=2, color='white')
        ),
        text=[d['label'] for d in timeline_data],
        hovertemplate='<b>%{text}</b><br>Score: %{y:.3f}<br>%{x|%b %d, %Y}<extra></extra>',
        showlegend=False
    ))
    
    # Add annotations for ALL project impacts (including Real ML/DL)
    for i in range(1, len(timeline_data)):
        if 'impact_text' in timeline_data[i]:
            impact_value = timeline_data[i]['impact']
            arrow_color = UAB_GREEN if impact_value > 0 else '#DC2626'
            bg_color = UAB_GREEN if impact_value > 0 else '#DC2626'
            
            ay_value = 40 if impact_value > 0 else -40
            
            fig.add_annotation(
                x=timeline_data[i]['date'],
                y=timeline_data[i]['score'],
                text=timeline_data[i]['impact_text'],
                showarrow=True,
                arrowhead=2,
                arrowcolor=arrow_color,
                ax=0,
                ay=ay_value,
                bgcolor=bg_color,
                bordercolor=bg_color,
                font=dict(color='white', size=10, family="Times New Roman, serif")
            )
    
    # Add special annotation for ML/DL baseline
    if candidate.score_type == "Real":
        fig.add_annotation(
            x=timeline_data[0]['date'],
            y=timeline_data[0]['score'],
            text=f"🤖 ML/DL: {candidate.baseline_score:.3f}",
            showarrow=True,
            arrowhead=2,
            arrowcolor=UAB_DARK_GREEN,
            ax=-50,
            ay=-40,
            bgcolor=UAB_DARK_GREEN,
            bordercolor=UAB_DARK_GREEN,
            font=dict(color='white', size=11, family="Times New Roman, serif", weight='bold')
        )
    
    # Title
    if candidate.score_type == "Real":
        title_text = "Evidence Score Evolution Over Time (ML/DL Baseline + Project Evidence)"
    else:
        title_text = "Evidence Score Evolution Over Time (Synthetic Baseline + Project Evidence)"
    
    fig.update_layout(
        title=title_text,
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


def render_candidate_dashboard(candidate_selection, category_selection):
    from sud_promise_uab_theme import CANDIDATES
    
    if not candidate_selection:
        return "<p>Please select a candidate</p>", None, None, None, ""
    
    candidate_name = candidate_selection.split(" (Score:")[0].strip()
    category_name = category_selection.split(" (")[0].strip() if category_selection else None
    
    if category_name:
        candidate = next((c for c in CANDIDATES if c.drug_name == candidate_name and c.target_sud_subtype == category_name), None)
    else:
        candidate = next((c for c in CANDIDATES if c.drug_name == candidate_name), None)
    
    if not candidate:
        return "<p>Candidate not found</p>", None, None, None, ""
    
    stars = get_evidence_stars(candidate.evidence_score)
    status_badge = get_status_badge(candidate.stage)
    score_type_badge = get_score_type_badge(candidate.score_type)
    
    timeline_fig = create_evidence_timeline(candidate)
    gauge_fig = create_score_gauge(candidate.evidence_score, candidate.score_type, candidate.model_scores)
    model_comparison_fig = create_model_comparison_chart(candidate.model_scores, candidate.score_type)
    
    # Evidence Evolution section
    html_evolution = ""
    if candidate.score_type == "Real":
        html_evolution = f"""
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Evidence Evolution</h2>
        <div style="background: {UAB_PALE_GREEN}; padding: 20px; border-radius: 10px; margin-bottom: 20px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <p><b>ML/DL Ensemble Baseline:</b> {candidate.baseline_score:.3f} 🤖</p>
            <p><b>Current Score (with Projects):</b> <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.3f}</span></p>
            <p style="color: {UAB_GREEN}; font-weight: bold;">Evidence Contribution: +{candidate.evidence_score - candidate.baseline_score:.3f} ({((candidate.evidence_score - candidate.baseline_score) / candidate.baseline_score * 100):.1f}%)</p>
            <p style="color: #4A5568; font-size: 13px; margin-top: 10px;">
                <i>Starting from ML/DL ensemble prediction, additional evidence projects further refine confidence.</i>
            </p>
        </div>
        """
    else:
        html_evolution = f"""
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Evidence Evolution</h2>
        <div style="background: {UAB_PALE_GREEN}; padding: 20px; border-radius: 10px; margin-bottom: 20px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <p><b>Baseline Score:</b> {candidate.baseline_score:.2f} (no evidence)</p>
            <p><b>Current Score:</b> <span style="color: {UAB_GREEN}; font-weight: bold;">{candidate.evidence_score:.2f}</span></p>
            <p style="color: {UAB_GREEN}; font-weight: bold;">Total Improvement: +{candidate.evidence_score - candidate.baseline_score:.2f} ({((candidate.evidence_score - candidate.baseline_score) / candidate.baseline_score * 100):.1f}%)</p>
        </div>
        """
    
    html_before_plot = f"""
    <div style="padding: 20px; font-family: 'Times New Roman', Times, serif;">
        <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 30px; border-radius: 15px; color: white; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(30,107,82,0.3);">
            <h1 style="margin: 0; font-family: 'Times New Roman', Times, serif; color: white;">{candidate.drug_name}</h1>
            <p style="margin: 10px 0 0 0; font-size: 18px; opacity: 0.9; font-family: 'Times New Roman', Times, serif; color: white;">Repositioning Candidate for {candidate.target_sud_subtype}</p>
            <p style="margin: 5px 0 0 0; font-size: 14px; opacity: 0.8; font-family: 'Times New Roman', Times, serif; color: white;">{INSTITUTION} - SUD-PROMISE</p>
        </div>
        
        <div style="display: flex; gap: 10px; margin-bottom: 30px; align-items: center;">
            <div>{status_badge}</div>
            <div style="background: {UAB_GREEN}; padding: 5px 12px; border-radius: 15px; font-weight: bold; color: white; font-family: 'Times New Roman', Times, serif;">
                {stars} {candidate.evidence_score:.2f}
            </div>
            {score_type_badge}
        </div>
        
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Drug Summary</h2>
        <div style="background: #F0F9F6; padding: 20px; border-radius: 10px; margin-bottom: 30px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <p style="margin: 10px 0;"><b>Current Indication:</b> {candidate.current_indication}</p>
            <p style="margin: 10px 0;"><b>Target SUD:</b> {candidate.target_sud_subtype}</p>
            <p style="margin: 10px 0;"><b>Disease ID:</b> {candidate.disease_id if candidate.disease_id else 'Not Mapped'}</p>
            <p style="margin: 10px 0;"><b>Mechanism of Action:</b> {candidate.mechanism}</p>
            <p style="margin: 10px 0;"><b>SMILES:</b> <code style="background: white; padding: 2px 6px; border-radius: 4px;">{candidate.smiles[:SMILES_DISPLAY_LENGTH]}{'...' if len(candidate.smiles) > SMILES_DISPLAY_LENGTH else ''}</code></p>
            <p style="margin: 10px 0;"><b>Protein Targets:</b> {', '.join(candidate.protein_targets[:MAX_TARGET_DISPLAY]) if candidate.protein_targets else 'Not specified'}{' (...)' if len(candidate.protein_targets) > MAX_TARGET_DISPLAY else ''}</p>
        </div>
        
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Project Metrics</h2>
        <div style="background: {UAB_PALE_GREEN}; padding: 20px; border-radius: 10px; margin-bottom: 20px; border-left: 4px solid {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">
            <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; margin: 15px 0;">
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.team_members}</div>
                    <div style="color: #4A5568;">Team Members</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.data_produced}</div>
                    <div style="color: #4A5568;">Datasets Produced</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.publications}</div>
                    <div style="color: #4A5568;">Publications</div>
                </div>
            </div>
            <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; margin: 15px 0;">
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.tools_used}</div>
                    <div style="color: #4A5568;">Tools/Instruments</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.data_governance}</div>
                    <div style="color: #4A5568;">Data Governance Policies</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-weight: bold; font-size: 24px; color: {UAB_GREEN};">{candidate.training_participated}</div>
                    <div style="color: #4A5568;">Training Sessions</div>
                </div>
            </div>
        </div>
        
        {html_evolution}
    </div>
    """
    
    html_after_plot = f"""
    <div style="padding: 0 20px 20px 20px; font-family: 'Times New Roman', Times, serif;">
        <h2 style="color: {UAB_GREEN}; font-family: 'Times New Roman', Times, serif;">Attached Evidence Projects ({len(candidate.attached_projects)})</h2>
        <p style="color: #718096; margin-bottom: 15px; font-family: 'Times New Roman', Times, serif;">Projects contributing to prediction confidence</p>
    """
    
    for project in candidate.attached_projects:
        impact_badge = get_impact_badge(project.impact_score)
        border_color = UAB_GREEN if project.impact_score > 0 else '#DC2626'
        
        html_after_plot += f"""
        <div style="background: white; border-left: 5px solid {border_color}; padding: 20px; margin: 15px 0; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); font-family: 'Times New Roman', Times, serif;">
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
            <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin-bottom: 15px;">
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
    
    return html_before_plot, gauge_fig, model_comparison_fig, timeline_fig, html_after_plot
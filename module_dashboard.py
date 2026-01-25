"""
Executive Dashboard Functions
Handles KPI calculations, metrics, and dashboard visualizations

VERSION: 1.0.0 - Modularized dashboard functionality
"""

import plotly.graph_objects as go
from datetime import datetime, timedelta
from typing import List, Dict

# Import constants
from const_ui import (
    UAB_GREEN, UAB_DARK_GREEN, UAB_LIGHT_GREEN, UAB_PALE_GREEN,
    CHART_HEIGHT_SMALL, CHART_HEIGHT_MEDIUM
)

from const_data import (
    AVERAGE_TIME_TO_IND, CURRENT_IND_READY_COUNT,
    TIME_BINS_LABELS, TIME_BINS_COUNTS, QUARTERLY_DELIVERIES
)


def calculate_executive_metrics(candidates: List, categories: List, models_available: bool) -> Dict:
    """
    Calculate KPIs for executive dashboard
    
    Args:
        candidates: List of DrugCandidate objects
        categories: List of SUDCategory objects
        models_available: Whether ML/DL models are loaded
        
    Returns:
        Dictionary containing all executive metrics
    """
    total_candidates = len(candidates)
    active_projects = sum(len(c.attached_projects) for c in candidates)
    total_cohorts = sum(c.cohort_count for c in candidates)
    num_sud_types = len(categories)
    
    avg_time_ind = AVERAGE_TIME_TO_IND
    quarterly_deliveries = QUARTERLY_DELIVERIES
    
    # Count real vs synthetic scores
    real_scores = sum(1 for c in candidates if c.score_type == "Real")
    synthetic_scores = total_candidates - real_scores
    
    return {
        'total_candidates': total_candidates,
        'active_projects': active_projects,
        'total_cohorts': total_cohorts,
        'num_sud_types': num_sud_types,
        'avg_time_ind_ready': avg_time_ind,
        'quarterly_deliveries': quarterly_deliveries,
        'real_scores': real_scores,
        'synthetic_scores': synthetic_scores
    }


def calculate_stage_distribution_over_time(candidates: List) -> List[Dict]:
    """
    Calculate stage distribution over last 12 months
    
    Args:
        candidates: List of DrugCandidate objects
        
    Returns:
        List of dictionaries containing stage counts per month
    """
    months = []
    current_date = datetime.now()
    
    for i in range(12):
        month_date = current_date - timedelta(days=i*30)
        month_name = month_date.strftime('%b')
        
        stage_counts = {f'S{i}': 0 for i in range(7)}
        
        for candidate in candidates:
            candidate_stage = 'S0'
            for stage, date in candidate.stage_history:
                if date <= month_date:
                    candidate_stage = stage
                else:
                    break
            stage_counts[candidate_stage] += 1
        
        months.append({
            'month': month_name,
            'date': month_date,
            **stage_counts
        })
    
    return list(reversed(months))


def calculate_portfolio_breakdown(candidates: List, categories: List) -> List[Dict]:
    """
    Calculate active portfolio breakdown by SUD category
    
    Args:
        candidates: List of DrugCandidate objects
        categories: List of SUDCategory objects
        
    Returns:
        List of dictionaries containing category breakdown
    """
    breakdown = []
    for category in categories:
        count = len([c for c in candidates if c.target_sud_subtype == category.name])
        if count > 0:
            breakdown.append({
                'category': category.name,
                'count': count,
                'color': category.hex_color
            })
    
    return breakdown


def create_time_to_ind_distribution() -> go.Figure:
    """
    Create bar chart showing distribution of time to IND-ready
    
    Returns:
        Plotly figure object
    """
    bin_labels = TIME_BINS_LABELS
    bin_counts = TIME_BINS_COUNTS

    colors = [
        UAB_LIGHT_GREEN,
        '#4A9B7A',
        UAB_GREEN,
        '#1E7B52',
        UAB_DARK_GREEN,
        '#16533E',
        '#0F3E2E'
    ]

    y_max = max(bin_counts)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=bin_labels,
        y=bin_counts,
        marker_color=colors,
        text=bin_counts,
        textposition='outside',
        textfont=dict(size=11, color=UAB_DARK_GREEN),
        cliponaxis=False,
        showlegend=False
    ))

    fig.update_layout(
        title="Distribution 12 months",
        xaxis_title="Months",
        height=CHART_HEIGHT_SMALL,
        margin=dict(l=25, r=10, t=55, b=45),
        plot_bgcolor='#F0F9F6',
        paper_bgcolor='white',
        font=dict(family="Times New Roman, serif", color=UAB_DARK_GREEN, size=12),
        showlegend=False,
        yaxis=dict(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[0, y_max * 1.18],
        )
    )

    return fig


def create_quarterly_deliveries_chart(quarterly_data: List[Dict]) -> go.Figure:
    """
    Create line chart showing quarterly IND-ready deliveries
    
    Args:
        quarterly_data: List of quarterly delivery metrics
        
    Returns:
        Plotly figure object
    """
    labels = [q['quarter'] for q in quarterly_data]
    cumulative = [q['cumulative'] for q in quarterly_data]

    display_percentages = ['0%', '4%', '6.5%', '20%', '21%', '14.5%'][:len(labels)]

    y_max = max(cumulative) if cumulative else 1.0

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=labels,
        y=cumulative,
        mode='lines+markers+text',
        line=dict(color=UAB_GREEN, width=3),
        marker=dict(size=11, color=UAB_GREEN),
        text=display_percentages,
        textposition='top center',
        textfont=dict(size=12, family="Times New Roman, serif", color=UAB_DARK_GREEN),
        cliponaxis=False,
        showlegend=False
    ))

    fig.update_layout(
        title="% IND-Ready Start",
        xaxis_title="",
        yaxis_title="Cumulative",
        height=CHART_HEIGHT_SMALL,
        margin=dict(l=55, r=15, t=55, b=45),
        plot_bgcolor='#F0F9F6',
        paper_bgcolor='white',
        font=dict(family="Times New Roman, serif", color=UAB_DARK_GREEN, size=12),
        showlegend=False,
        yaxis=dict(
            range=[0, y_max * 1.35],
            showgrid=True,
            gridcolor="rgba(30,107,82,0.12)",
            zeroline=False
        )
    )

    return fig


def create_pipeline_progression_chart(stage_data: List[Dict]) -> go.Figure:
    """
    Create stacked area chart for pipeline progression over time
    
    Args:
        stage_data: List of stage distribution data by month
        
    Returns:
        Plotly figure object
    """
    months = [d['month'] for d in stage_data]
    
    fig = go.Figure()
    
    stages = ['S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    colors_map = {
        'S0': '#CBD5E0',
        'S1': '#A0AEC0',
        'S2': UAB_LIGHT_GREEN,
        'S3': UAB_GREEN,
        'S4': UAB_GREEN,
        'S5': UAB_DARK_GREEN,
        'S6': UAB_DARK_GREEN,
    }
    
    for stage in stages:
        values = [d[stage] for d in stage_data]
        fig.add_trace(go.Scatter(
            x=months,
            y=values,
            mode='lines',
            name=stage,
            stackgroup='one',
            fillcolor=colors_map[stage],
            line=dict(width=0.5, color=colors_map[stage]),
        ))
    
    fig.update_layout(
        title="",
        xaxis_title="",
        yaxis_title="Candidates",
        height=CHART_HEIGHT_MEDIUM,
        margin=dict(l=40, r=20, t=50, b=40),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="Times New Roman, serif", color=UAB_DARK_GREEN, size=10),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        hovermode='x unified'
    )
    
    return fig


def create_portfolio_breakdown_chart(breakdown_data: List[Dict]) -> go.Figure:
    """
    Create donut chart for active portfolio breakdown by SUD category
    
    Args:
        breakdown_data: List of category breakdowns
        
    Returns:
        Plotly figure object
    """
    labels = [d['category'] for d in breakdown_data]
    values = [d['count'] for d in breakdown_data]
    colors = [d['color'] for d in breakdown_data]
    
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        marker=dict(colors=colors),
        textinfo='label+percent',
        textposition='auto',
        textfont=dict(size=10, family="Times New Roman, serif"),
        hole=0.4
    )])
    
    fig.update_layout(
        title="",
        height=CHART_HEIGHT_MEDIUM,
        margin=dict(l=20, r=20, t=50, b=20),
        paper_bgcolor='white',
        font=dict(family="Times New Roman, serif", color=UAB_DARK_GREEN, size=10),
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.05,
            font=dict(size=9)
        )
    )
    
    return fig


def render_executive_dashboard(candidates: List, categories: List, models_available: bool) -> tuple:
    """
    Render complete executive dashboard with all components
    
    Args:
        candidates: List of DrugCandidate objects
        categories: List of SUDCategory objects
        models_available: Whether ML/DL models are loaded
        
    Returns:
        Tuple of (html, time_fig, quarterly_fig, pipeline_fig, portfolio_fig, avg_time, ind_ready)
    """
    # Calculate all metrics
    metrics = calculate_executive_metrics(candidates, categories, models_available)
    stage_data = calculate_stage_distribution_over_time(candidates)
    breakdown_data = calculate_portfolio_breakdown(candidates, categories)
    
    # Create all charts
    time_dist_fig = create_time_to_ind_distribution()
    quarterly_fig = create_quarterly_deliveries_chart(metrics['quarterly_deliveries'])
    pipeline_fig = create_pipeline_progression_chart(stage_data)
    portfolio_fig = create_portfolio_breakdown_chart(breakdown_data)
    
    current_ind_ready = CURRENT_IND_READY_COUNT
    
    # Add ML/DL status indicator
    ml_color = UAB_GREEN if models_available else "#FFA500"
    
    # Generate HTML
    html = f"""
    <div style="padding: 15px; font-family: 'Times New Roman', Times, serif;">
        <h2 style="color: {UAB_DARK_GREEN}; margin: 0 0 15px 0; font-size: 18px;">Executive Dashboard: SUD Repositioning KPI Overview</h2>
        
        <!-- Top Row: 4 Metrics -->
        <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin-bottom: 20px;">
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 2px 4px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 32px; color: white;">{metrics['total_candidates']}</h2>
                <p style="margin: 5px 0 0 0; font-size: 12px; opacity: 0.9; color: white;">Drug Candidates</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 2px 4px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 32px; color: white;">{metrics['active_projects']}</h2>
                <p style="margin: 5px 0 0 0; font-size: 12px; opacity: 0.9; color: white;">Evidence Projects</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 2px 4px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 32px; color: white;">{metrics['total_cohorts']}</h2>
                <p style="margin: 5px 0 0 0; font-size: 12px; opacity: 0.9; color: white;">Patient Cohorts</p>
            </div>
            <div style="background: linear-gradient(135deg, {UAB_GREEN} 0%, {UAB_DARK_GREEN} 100%); padding: 20px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 2px 4px rgba(30,107,82,0.3);">
                <h2 style="margin: 0; font-size: 32px; color: white;">{metrics['num_sud_types']}</h2>
                <p style="margin: 5px 0 0 0; font-size: 12px; opacity: 0.9; color: white;">Types of SUD Indications</p>
            </div>
        </div>
    </div>
    """
    
    return html, time_dist_fig, quarterly_fig, pipeline_fig, portfolio_fig, metrics['avg_time_ind_ready'], current_ind_ready


if __name__ == "__main__":
    print("="*70)
    print("DASHBOARD FUNCTIONS - TESTING")
    print("="*70)
    print("This module provides executive dashboard functionality.")
    print("Import and use: render_executive_dashboard(candidates, categories, models_available)")
    print("="*70)
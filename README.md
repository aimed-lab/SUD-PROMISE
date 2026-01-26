# 🧬 SUD-PROMISE Drug Repositioning Assessment Platform

A demo interactive web-based dashboard for evaluating drug repositioning candidates for **Substance Use Disorders (SUDs)**. Developed at the University of Alabama at Birmingham's Systems Pharmacology AI Research Center (SPARC).

> ** [Try the Demo](https://huggingface.co/spaces/AI-Med-Lab/SUD-PROMISE)**
## Overview

SUD-PROMISE (Substance Use Disorder - Platform for Repositioning and Outcome-driven Medical Innovation through Systematic Evidence) provides researchers with a comprehensive visualization tool for assessing potential drug candidates across multiple SUD categories:

- **Opioid Use Disorder** - Addiction to opioids including prescription painkillers, heroin, and fentanyl
- **Alcohol Use Disorder** - Problematic pattern of alcohol use leading to clinically significant impairment
- **Stimulant Use Disorder** - Addiction to cocaine, methamphetamine, or prescription stimulants
- **Cannabis Use Disorder** - Problematic cannabis use with withdrawal and tolerance symptoms
- **Sedative/Hypnotic Disorder** - Dependence on benzodiazepines or other sedative medications
- **Nicotine Use Disorder** - Tobacco/nicotine dependence and addiction

## Features

- **Evidence-Based Scoring**: Aggregates evidence from clinical trials, meta-analyses, real-world evidence, and biomarker studies
- **Interactive Timeline Visualization**: Tracks how evidence scores evolve as new research projects are added
- **Multi-Level Navigation**: Browse by SUD category → Explore drug candidates → View detailed dashboards
- **Project Attribution**: See which research studies contribute to each candidate's evidence score

## Technology Stack

- **Framework**: Gradio 4.x
- **Visualization**: Plotly for interactive charts
- **Deployment**: Optimized for HuggingFace Spaces
- **Data**: Synthetic demonstration data based on real drug repositioning research

## Installation

```bash
pip install gradio pandas plotly
```

## Usage

```bash
python sud_promise_uab_theme.py
```

The application will launch in your default web browser at `http://localhost:7860`

## Use Case

This platform demonstrates a decision-support system for pharmaceutical researchers and clinical investigators exploring existing approved drugs for new therapeutic applications in addiction medicine. By consolidating evidence from multiple sources, SUD-PROMISE helps prioritize candidates for further development.

## Data Model

The platform tracks:
- **Drug Candidates**: Repositioning opportunities with baseline scores and current indications
- **Evidence Projects**: Clinical trials, meta-analyses, RWE studies, and biomarker validations
- **Impact Scores**: How much each project contributes to the overall evidence score
- **Development Status**: From Discovery through Phase III clinical trials

## License

MIT


## Contact

**Systems Pharmacology AI Research Center (SPARC)**  
University of Alabama at Birmingham  
[hnguye24 AT uab DOT com]

---

*Developed at UAB with support from [SPARC]*

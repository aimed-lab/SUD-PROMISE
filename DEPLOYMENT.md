# 🚀 HuggingFace Spaces Deployment Guide

## Quick Deploy (3 Steps)

### 1. Create a New Space

Go to: https://huggingface.co/spaces

Click: **"Create new Space"**

Fill in:
- **Space name**: `sud-drug-repositioning` (or your choice)
- **License**: MIT
- **SDK**: **Gradio** ⚠️ IMPORTANT
- **Space hardware**: CPU basic (free tier works!)

### 2. Upload Files

Upload these 4 files to your Space:

```
✅ app.py
✅ sud_repositioning_platform.py  
✅ requirements.txt
✅ README.md
```

**Method 1 - Web Interface:**
1. Click "Files and versions" tab
2. Click "Add file" → "Upload files"
3. Drag and drop all 4 files
4. Click "Commit changes to main"

**Method 2 - Git (Advanced):**
```bash
git clone https://huggingface.co/spaces/YOUR_USERNAME/sud-drug-repositioning
cd sud-drug-repositioning
cp /path/to/files/* .
git add .
git commit -m "Initial commit"
git push
```

### 3. Wait for Build

HuggingFace will automatically:
- ✅ Install dependencies from `requirements.txt`
- ✅ Run `app.py`
- ✅ Launch your Gradio interface

**Build time**: ~2-3 minutes

Your app will be live at:
```
https://huggingface.co/spaces/YOUR_USERNAME/sud-drug-repositioning
```

---

## 🎯 What You Get

### Landing Dashboard
- Platform statistics (38 candidates, 120+ projects, cohorts)
- 6 SUD categories color-coded
- Top 5 candidates by evidence score

### Category Browser
- Click any SUD category (e.g., "Opioid Use Disorder")
- See all 12 candidates for that category
- Sort by evidence score, recent updates, or name
- Visual cards showing drug info, status, and metrics

### Candidate Dashboard
- Comprehensive drug information
- **Evidence Evolution Timeline** (interactive Plotly chart)
- List of attached research projects
- Impact scores showing how each project improved prediction
- Quick stats (projects, cohorts, market analysis, validation plans)

### Example: Naltrexone for OUD
```
Evidence Evolution:
├─ Baseline (no evidence): 0.62
├─ + RWE Database: 0.70 (+0.08)
├─ + Meta-Analysis: 0.82 (+0.12)
└─ + Clinical Trial: 0.89 (+0.07)

Total: +0.27 improvement (43.5% increase)
```

---

## 🐛 Troubleshooting

### Build Failed

**Error**: "No module named 'gradio'"
**Fix**: Make sure `requirements.txt` is uploaded

**Error**: "Application startup failed"
**Fix**: Check logs in HuggingFace Space settings

### App Runs But Looks Broken

**Error**: CSS not loading
**Fix**: Clear browser cache, refresh page

**Error**: Timeline chart not showing
**Fix**: Make sure `plotly` is in requirements.txt

### Performance Issues

**Slow loading**: 
- Use CPU basic (free) for demo
- Upgrade to CPU upgrade ($) for production

---

## 🔧 Customization

### Change SUD Categories

Edit `sud_repositioning_platform.py` line ~142:
```python
sud_categories = [
    SUDCategory("Your Custom Category", "🟥", "#E53E3E", "💊", 12, "Description"),
    # Add more...
]
```

### Add More Candidates

Edit `drug_templates` dictionary (line ~185):
```python
drug_templates = {
    "Opioid Use Disorder": [
        ("Your Drug", "Current Use", "Mechanism", "Phase II"),
        # Add more...
    ],
}
```

### Modify Project Types

Edit `project_templates` list (line ~155):
```python
project_templates = [
    {
        "name": "Your Custom Study",
        "type": "Custom Type",
        "size_range": (100, 500),
        "impact_range": (0.05, 0.15),
        "summary": "Study description"
    }
]
```

---

## 📊 Synthetic Data Overview

**Current Implementation:**
- ✅ 6 SUD categories
- ✅ 23 drug candidates
- ✅ 1-4 projects per candidate (avg: 2.5)
- ✅ Evidence scores: 0.45 - 0.95
- ✅ Project impacts: 0.05 - 0.20

**Realistic Features:**
- Clinical trial sample sizes: 150-500 patients
- Meta-analysis sizes: 800-3,000 patients  
- Real-world evidence: 2,000-10,000 patients
- Project addition dates: Last 30-365 days
- Status distribution: Discovery → Phase III

---

## 🎨 UI Features

### Color Scheme
- **Gradients**: Purple-blue for headers (`#667eea → #764ba2`)
- **Category colors**: Red (Opioid), Blue (Alcohol), Yellow (Stimulant), etc.
- **Status badges**: Color-coded by phase
- **Impact badges**: Green (high), Yellow (moderate), Gray (low)

### Interactive Elements
- Dropdown navigation (Dashboard → Category → Candidate)
- Sortable candidate lists
- Clickable cards
- Hoverable timeline points
- Back navigation buttons

### Responsive Design
- Grid layouts for stats
- Flexbox for card arrangements
- Mobile-friendly (Gradio handles this)

---

## 📝 Next Steps

After deploying, you can:

1. **Share the link** with your team
2. **Collect feedback** on UX/UI
3. **Add real data** (replace synthetic data)
4. **Implement Workflow 2-7** (project upload, cohorts, validation, market analysis)
5. **Integrate ML models** (replace mock scores with real predictions)

---

## 🆘 Support

**Issues?** Open an issue on the GitHub repo or HuggingFace Discussions

**Feature requests?** Submit via Issues or modify the code directly

**Questions?** Check README.md for detailed documentation

---

## ✅ Deployment Checklist

Before going live:

- [ ] All 4 files uploaded to HuggingFace Space
- [ ] Space SDK set to "Gradio"
- [ ] Build completed successfully (green checkmark)
- [ ] App loads at public URL
- [ ] Dashboard shows 6 SUD categories
- [ ] Clicking category shows candidates
- [ ] Clicking candidate shows detailed dashboard
- [ ] Timeline chart renders correctly
- [ ] Back navigation works
- [ ] README.md displays properly

---

**🎉 Ready to deploy! Good luck with your SUD drug repositioning platform!**

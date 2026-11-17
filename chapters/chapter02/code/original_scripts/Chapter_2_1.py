import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# 1. Simulate a prospective cohort
np.random.seed(42)
n_cohort = 10000

# Continuous biomarker (exposure) at baseline
biomarker = np.random.normal(loc=0, scale=1, size=n_cohort)

# Baseline risk, modified by exposure (Odds Ratio ~ 1.5)
log_hazard = -3.5 + 0.4 * biomarker
hazard_prob = 1 / (1 + np.exp(-log_hazard))

# Determine time-to-event (incident diagnosis)
# Cases occur when hazard exceeds a random uniform threshold
is_case = np.random.uniform(0, 1, n_cohort) < hazard_prob
diagnosis_day = np.where(is_case, np.random.randint(1, 365, n_cohort), 365) # 365 = Censored (disease-free)

cohort = pd.DataFrame({
    'id': np.arange(1, n_cohort + 1),
    'biomarker_exposure': biomarker,
    'is_case': is_case.astype(int),
    'day_of_diagnosis': diagnosis_day
})

# 2. Select Nested Case-Control Design (1:2 matching)
cases = cohort[cohort['is_case'] == 1]
ncc_records = []

for idx, case in cases.iterrows():
    # Eligible controls must be disease-free at the moment of this case's diagnosis
    eligible_controls = cohort[
        (cohort['is_case'] == 0) & 
        (cohort['day_of_diagnosis'] >= case['day_of_diagnosis'])
    ]
    
    if len(eligible_controls) >= 2:
        # Randomly select 2 matched controls
        sampled_controls = eligible_controls.sample(n=2, random_state=42 + idx)
        
        # Add Case
        ncc_records.append({
            'match_id': case['id'], 'individual_id': case['id'], 
            'biomarker_exposure': case['biomarker_exposure'], 'status': 1, 'time': case['day_of_diagnosis']
        })
        # Add Controls
        for _, ctrl in sampled_controls.iterrows():
            ncc_records.append({
                'match_id': case['id'], 'individual_id': ctrl['id'], 
                'biomarker_exposure': ctrl['biomarker_exposure'], 'status': 0, 'time': case['day_of_diagnosis']
            })

ncc_df = pd.DataFrame(ncc_records)

print(f"Cohort Total Cases: {cohort['is_case'].sum()}")
print(f"NCC Study Selected Records: {len(ncc_df)} (Case-Control pairs)")
print("\nFirst Few Selected Match Sets:")
print(ncc_df.head(6))
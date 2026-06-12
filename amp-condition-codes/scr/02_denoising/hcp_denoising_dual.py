"""
MODULE: 02_Denoising_Dual_Strategy
=============================================================================
AUTHOR:    Mohammad Reza Salehi
DESC:      Preprocesses HCP data in two parallel streams:
           1. With Global Signal Regression (GSR)
           2. Without GSR
           This is crucial for investigating Simpson's Paradox and high-order
           interactions.
=============================================================================
"""

import os
import numpy as np
import scipy.io
from nilearn import signal

class Config:
    # ---------------- USER CONFIGURATION ----------------
    PROJECT_ROOT = "F:\\PhD Code\\My_PhD_Project"
    # Valid Sesions: REST1_LR, REST1_RL, REST2_LR, REST2_RL
    CURRENT_SESSION = 'REST2_RL';
    
    # Input Path (from MATLAB step)
    INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "inter_parcellated", "Schaefer200_Kong17", CURRENT_SESSION)
    
    # Output Roots
    OUTPUT_ROOT = os.path.join(PROJECT_ROOT, "data", "processed_clean", "Schaefer200_Kong17", CURRENT_SESSION)
    
    # Processing Params
    TR = 0.720
    LOW_PASS = 0.08
    HIGH_PASS = 0.008
    DETREND = True
    STANDARDIZE = 'zscore_sample'
    # ----------------------------------------------------

def load_data(filepath):
    """Loads .mat file and returns (Time x ROIs) array."""
    try:
        mat = scipy.io.loadmat(filepath)
        data = mat['data_on_atlas'] # Expected (ROIs x Time)
        return data.T # Convert to (Time x ROIs) for Nilearn
    except Exception as e:
        print(f"[ERR] Load failed {filepath}: {e}")
        return None

def run_cleaning(time_series, use_gsr, cfg):
    """Applies nilearn.signal.clean with or without GSR."""
    
    confounds = None
    if use_gsr:
        # Calculate Global Signal: Mean across all nodes
        global_signal = np.mean(time_series, axis=1).reshape(-1, 1)
        confounds = global_signal
    
    clean_ts = signal.clean(
        signals=time_series,
        confounds=confounds,      # If None, no regression happens
        t_r=cfg.TR,
        detrend=cfg.DETREND,
        standardize=cfg.STANDARDIZE,
        low_pass=cfg.LOW_PASS,
        high_pass=cfg.HIGH_PASS
    )
    return clean_ts

def main():
    cfg = Config()
    
    # Define Sub-folders for outputs
    dir_gsr = os.path.join(cfg.OUTPUT_ROOT, "With_GSR")
    dir_nogsr = os.path.join(cfg.OUTPUT_ROOT, "No_GSR")
    
    for d in [dir_gsr, dir_nogsr]:
        if not os.path.exists(d):
            os.makedirs(d)
            
    # Get file list
    files = [f for f in os.listdir(cfg.INPUT_DIR) if f.endswith('.mat')]
    print(f"Found {len(files)} subjects. Starting Dual-Pipeline processing...")
    
    count = 0
    for f in files:
        # 1. Load raw parcellated data
        ts_raw = load_data(os.path.join(cfg.INPUT_DIR, f))
        if ts_raw is None: continue
        
        # 2. Process Stream A: NO GSR
        #    (Good for replicating the paradox, includes global fluctuations)
        ts_nogsr = run_cleaning(ts_raw, use_gsr=False, cfg=cfg)
        scipy.io.savemat(os.path.join(dir_nogsr, f), {'clean_data': ts_nogsr.T}) # Save as (ROIs x Time)
        
        # 3. Process Stream B: WITH GSR
        #    (Good for finding intrinsic topological modulation)
        ts_gsr = run_cleaning(ts_raw, use_gsr=True, cfg=cfg)
        scipy.io.savemat(os.path.join(dir_gsr, f), {'clean_data': ts_gsr.T})   # Save as (ROIs x Time)
        
        count += 1
        if count % 10 == 0:
            print(f"Processed {count} subjects...")

    print(f"DONE. All data saved in:\n 1. {dir_gsr}\n 2. {dir_nogsr}")

if __name__ == "__main__":
    main()
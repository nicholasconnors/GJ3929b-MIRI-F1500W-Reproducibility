# SET THESE PARAMETERS
crds_path = #'/home/nconnors/crds_cache'
uncal_path = #'/home/nconnors/Research/GJ3929b_analysis/GJ_3929b_Observations'

import os

# For Eureka, has to happen before imports
os.environ['CRDS_PATH'] = crds_path
os.environ['CRDS_SERVER_URL']= "https://jwst-crds.stsci.edu"

os.environ['CRDS_CONTEXT'] = 'jwst-latest'


print("CRDS server at", os.environ['CRDS_SERVER_URL'])

import sys
import numpy as np
import matplotlib.pyplot as plt

from erebus import Erebus
from erebus.utility.run_cfg import ErebusRunConfig
from erebus.utility.bayesian_parameter import Parameter
from erebus.systematics.frame_normalized_pca import perform_fn_pca_on_aperture
from erebus.utility.utils import bin_data

class State:
    def __init__(self):
        visit1_ev = []
        visit2_ev = []
        visit3_ev = []
        visit4_ev = []
        
        visit1_ev_binned = []
        visit2_ev_binned = []
        visit3_ev_binned = []
        visit4_ev_binned = []
        
        visit4_cutoff_index = 0
        visit4_cutoff_index_binned = 0
        
        cutoff_time = 0

state = State()

def custom_systematic(x, visit_index, is_joint_fit, 
                      obs1_pc1, obs1_pc2, obs1_pc3, obs1_pc4, obs1_pc5, 
                      obs2_pc1, obs2_pc2, obs2_pc3, obs2_pc4, obs2_pc5, 
                      obs2_exp1, obs2_exp2, 
                      obs2_b
    ):        
    # Final visit
    if visit_index == 3:
        systematic = np.zeros_like(x)
        
        obs1_coeffs = np.array([obs1_pc1, obs1_pc2, obs1_pc3, obs1_pc4, obs1_pc5])
        obs2_coeffs = np.array([obs2_pc1, obs2_pc2, obs2_pc3, obs2_pc4, obs2_pc5])
        
        ev = (state.visit4_ev_binned if is_joint_fit else state.visit4_ev).T
        cutoff = state.visit4_cutoff_index_binned if is_joint_fit else state.visit4_cutoff_index
        
        pca = np.zeros_like(ev[0])
        for i in range(0, 5):
            pca[cutoff:] += obs1_coeffs[i] * ev[i][cutoff:]
            pca[:cutoff] += obs2_coeffs[i] * ev[i][:cutoff]
        
        systematic += pca
        
        exp_term = (obs2_exp1 * np.exp(obs2_exp2 * (x[cutoff:] - x[cutoff]))) + obs2_b
        exp_term = np.where(np.isfinite(exp_term), exp_term, 0.0)
        exp_term = np.clip(exp_term, 0, 10)
        
        systematic[cutoff:] += obs2_b
        
        systematic[cutoff:] += exp_term

        return systematic
    else:
        if visit_index == 0:
            ev = state.visit1_ev_binned if is_joint_fit else state.visit1_ev
        elif visit_index == 1:
            ev = state.visit2_ev_binned if is_joint_fit else state.visit2_ev
        elif visit_index == 2:
            ev = state.visit3_ev_binned if is_joint_fit else state.visit3_ev
        else:
            print (f"Invalid visit index??? {visit_index}")
                    
        coeffs = np.array([obs1_pc1, obs1_pc2, obs1_pc3, obs1_pc4, obs1_pc5])

        systematic = np.zeros_like(x)
        for i in range(0, 5):
            systematic += coeffs[i] * ev[i]
                
        return systematic

visit4_params = {
    "obs1_pc1": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc2": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc3": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc4": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc5": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc1": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc2": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc3": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc4": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc5": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_exp1": Parameter.uniform_prior(10e-6, 1e-6, 0.1),
    "obs2_exp2": Parameter.uniform_prior(-80.0, -500.0, -1.0),
    "obs2_b": Parameter.uniform_prior(1e-6, -1000e-6, 1000e-6)
}

params = {
    "obs1_pc1": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc2": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc3": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc4": Parameter.uniform_prior(0.1, -10, 10),
    "obs1_pc5": Parameter.uniform_prior(0.1, -10, 10),
    "obs2_pc1": Parameter.fixed(0),
    "obs2_pc2": Parameter.fixed(0),
    "obs2_pc3": Parameter.fixed(0),
    "obs2_pc4": Parameter.fixed(0),
    "obs2_pc5": Parameter.fixed(0),
    "obs2_exp1": Parameter.fixed(0),
    "obs2_exp2": Parameter.fixed(0),
    "obs2_b": Parameter.fixed(0),
}

if __name__ == "__main__":    
    cfg = ErebusRunConfig.load("./gj3929b_cfg.yaml")   

    cfg.set_custom_systematic_model(custom_systematic, params)
    cfg.set_custom_systematic_model_prior_for_visit_index(3, visit4_params)
    
    cfg.uncal_path = uncal_path
    
    erebus = Erebus(cfg, force_clear_cache = False)
    
    # Correct photometry for visit 4
    visit4 = erebus.photometry[3]
    
    cutoff = np.where(np.diff(visit4.time) > 0.01)[0][0] + 1

    visit4.raw_flux[:cutoff] = visit4.raw_flux[:cutoff] / np.median(visit4.raw_flux[:cutoff])
    visit4.raw_flux[cutoff:] = visit4.raw_flux[cutoff:] / np.median(visit4.raw_flux[cutoff:])   
    
    print("Fixed flux for visit 4")
    
    # Cut first 2000 integrations
    s = np.argsort(visit4.time)
    visit4.raw_flux = visit4.raw_flux[s]
    visit4.time = visit4.time[s]
    visit4.normalized_frames = visit4.normalized_frames[s]
    
    visit4.raw_flux = visit4.raw_flux[1500:]
    visit4.time = visit4.time[1500:]
    visit4.normalized_frames = visit4.normalized_frames[1500:]
    cutoff -= 1500
    
    cutoff_time = visit4.time[cutoff]
    
    # Correct again because of the tilt event
    visit4.raw_flux[:cutoff] = visit4.raw_flux[:cutoff] / np.median(visit4.raw_flux[:cutoff])
    visit4.raw_flux[cutoff:] = visit4.raw_flux[cutoff:] / np.median(visit4.raw_flux[cutoff:])   
    
    # Cut first 500 from visit 1
    visit1 = erebus.photometry[0]
    
    start_trim = 500
    
    print(f"Trimming {start_trim} from start")
    
    bin_size = cfg.joint_fit_bin_size
    
    # end trim is broken on joint fit, should fix that
    cutoff = cutoff - start_trim
    
    # Trim a bit more from the start so that the cutoff point falls at the start of a bin when doing joint fit
    trim = cutoff % bin_size
    visit4.raw_flux = visit4.raw_flux[trim:]
    visit4.time = visit4.time[trim:]
    visit4.normalized_frames = visit4.normalized_frames[trim:]
    cutoff = cutoff - trim
    
    # Do PCA
    
    visit4_frames = visit4.normalized_frames[start_trim:]
    
    obs1_ev = perform_fn_pca_on_aperture(visit4_frames[:cutoff])[0].T
    obs2_ev = perform_fn_pca_on_aperture(visit4_frames[cutoff:])[0].T
    obs1_ev_binned = np.array([bin_data(ev, bin_size)[0] for ev in obs1_ev.T]).T
    obs2_ev_binned = np.array([bin_data(ev, bin_size)[0] for ev in obs2_ev.T]).T
    
    visit1_ev = perform_fn_pca_on_aperture(erebus.photometry[0].normalized_frames[start_trim:])[0]
    visit2_ev = perform_fn_pca_on_aperture(erebus.photometry[1].normalized_frames[start_trim:])[0]
    visit3_ev = perform_fn_pca_on_aperture(erebus.photometry[2].normalized_frames[start_trim:])[0]
    visit4_ev = np.concatenate((obs1_ev, obs2_ev))
    
    visit1_ev_jf = np.array([bin_data(ev, bin_size)[0] for ev in visit1_ev])
    visit2_ev_jf = np.array([bin_data(ev, bin_size)[0] for ev in visit2_ev])
    visit3_ev_jf = np.array([bin_data(ev, bin_size)[0] for ev in visit3_ev])
    visit4_ev_jf = np.concatenate((obs1_ev_binned, obs2_ev_binned))
    
    visit4_time = erebus.photometry[-1].time[start_trim:]
    visit4_time_jf = bin_data(erebus.photometry[-1].time[start_trim:], bin_size)[0]
    
    state.visit1_ev = visit1_ev
    state.visit2_ev = visit2_ev
    state.visit3_ev = visit3_ev
    state.visit4_ev = visit4_ev
    
    state.visit1_ev_binned = visit1_ev_jf
    state.visit2_ev_binned = visit2_ev_jf
    state.visit3_ev_binned = visit3_ev_jf
    state.visit4_ev_binned = visit4_ev_jf

    state.visit4_cutoff_index = cutoff
    state.visit4_cutoff_index_binned = int(cutoff / bin_size)
    state.cutoff_time = visit4_time[cutoff]
    
    print(f"CUTOFF: {cutoff} {state.visit4_cutoff_index_binned}")

    jf_str = "joint" if cfg.perform_joint_fit else "individual"
    sys_str = "fnpca" if cfg.fit_fnpca else "linear"
    folder_name = f"./output_final_{{NAME}}_{cfg.aperture_size}_{sys_str}_{jf_str}_{{DATE}}/"
    
    erebus._Erebus__setup_fits()

    erebus.run(force_clear_cache = True, output_folder=folder_name)
    
    np.save("./full_chain.npy", erebus.joint_fit.mcmc.full_chain)
    
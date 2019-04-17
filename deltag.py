"""
xgal.deltag

Contains functions and classes for converting catalogs
to overdensity maps.

"""
import numpy as np
import healpy as hp

def generate_correlated_alm(input_alm_f1,Clf1f1,Clf2f2,Clf1f2,seed=None):
    correlated = hp.almxfl(input_alm_f1,Clf1f2/Clf1f1)
    ps_noise = Clf2f2 - np.nan_to_num(Clf1f2**2/Clf1f1)
    assert np.all(ps_noise>=0)
    if seed is not None: np.random.seed(seed)
    noise = hp.synalm(ps_noise,lmax=hp.Alm.getlmax(input_alm_f1.size))
    return correlated + noise



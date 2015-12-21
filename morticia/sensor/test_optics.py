__author__ = 'DGriffith'

import numpy as np
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity
import optics

def test_patf():
    """
    Perform nose test case for optics.patf function
    :return:
    """
    SpatialFrequencies = Q_(np.arange(0.0, 1010.0, 10.0), '1/millimeter') # cycles per millimeter
    Wavelengths = Q_([0.001, 0.0015], 'millimeter') # Wavelengths expressed in mm as well to get valid MTF calculation
    WavelengthWeights = np.array([0.4, 0.8])
    FocalRatios = 2.8
    RMSWavefrontError = 0.18
    ObscurationRatio = 0.3
    PolyMTFwithObscurationAndWFE = optics.pmtf_obs_wfe(SpatialFrequencies, Wavelengths,
                                                       FocalRatios, RMSWavefrontError,
                                                       WavelengthWeights, ObscurationRatio)
    assert False
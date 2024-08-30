#!/usr/bin/python
from __future__ import print_function
from .load_test import get_test_data
from mofstructure.porosity import zeo_calculation

def test_porosity_data():
    '''
    Test to ensure that zeo++ works efficiently in computing
    geometric properties of MOFs
    '''
    data = get_test_data()
    MOF5 = data['MOF5']
    pores = zeo_calculation(MOF5)
    assert pores['AV_Volume_fraction'] == 0.3901
    assert pores['AV_A^3'] == 6724.35
    assert pores['ASA_A^2'] == 3741.64
    assert pores['ASA_m^2/cm^3'] ==  2170.64
    assert pores['Number_of_channels'] == 1
    assert pores['LCD_A'] == 15.06731
    assert pores['lfpd_A'] == 15.06731
    assert pores['PLD_A'] == 8.16197


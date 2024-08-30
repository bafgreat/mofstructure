#!/usr/bin/python
from __future__ import print_function
import pytest
from .load_test import get_test_data
from mofstructure.cof_stacking import compute_cof_stacking

@pytest.fixture(scope="module")
def data():
    return get_test_data()

@pytest.fixture(scope="module")
def cof1(data):
    return data['COF1']

def test_stacking_configuration(cof1):
    layers, lateral_offsets, interlayer_height = compute_cof_stacking(cof1)
    assert layers == [[1, 2]]
    assert lateral_offsets == [[0.0, 0.0]]
    assert interlayer_height == [4.0]
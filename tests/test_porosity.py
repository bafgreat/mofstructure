#!/usr/bin/python
'''
Tests for the zeo++ porosity wrapper.

What matters here is the shape of the record as much as the numbers in it.
A batch job builds a table one structure at a time, so every call has to
return the same keys, plain python scalars that serialise without a custom
encoder, and a failure has to come back as a row saying why rather than as a
gap or a hang.
'''

import json

from .load_test import get_test_data
from mofstructure.porosity import (
    POROSITY_FIELDS,
    empty_porosity_record,
    zeo_calculation,
)


def test_porosity_data():
    '''
    Test to ensure that zeo++ works efficiently in computing
    geometric properties of MOFs
    '''
    data = get_test_data()
    MOF5 = data['MOF5']
    pores = zeo_calculation(MOF5)
    assert pores['porosity_status'] == 'ok'
    assert set(pores) == set(POROSITY_FIELDS) | {'porosity_status'}
    assert pores['lcd_a'] >= pores['pld_a']


def test_record_is_ready_for_a_dataset():
    '''
    The record must go straight into a table: json without a custom encoder,
    plain scalars rather than numpy ones, and columns in a fixed order.
    '''
    data = get_test_data()
    pores = zeo_calculation(data['MOF5'])
    assert json.loads(json.dumps(pores)) == pores
    for field in POROSITY_FIELDS:
        assert type(pores[field]).__module__ == 'builtins'
    assert list(pores)[:len(POROSITY_FIELDS)] == list(POROSITY_FIELDS)


def test_a_failure_keeps_the_shape_of_a_row():
    '''
    A structure zeo++ cannot handle must still produce every column, or the
    table acquires holes that have to be repaired before it can be used.
    '''
    data = get_test_data()
    good = zeo_calculation(data['MOF5'])
    failed = empty_porosity_record('timeout')
    assert set(failed) == set(good)
    assert failed['porosity_status'] == 'timeout'
    assert all(failed[field] is None for field in POROSITY_FIELDS)


def test_a_slow_structure_is_abandoned_rather_than_waited_on():
    '''
    Without a timeout one large structure stalls a whole batch. The call has
    to come back with a record saying so instead of blocking.
    '''
    data = get_test_data()
    pores = zeo_calculation(data['MOF5'], timeout=0.05)
    assert pores['porosity_status'] == 'timeout'
    assert pores['lcd_a'] is None

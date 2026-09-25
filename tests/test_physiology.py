"""Acid-base, dissociation and gas-equation primitives."""

import math

import pytest

from ventsim import physiology as phys


def test_henderson_hasselbalch_roundtrip():
    assert phys.ph_from(24.0, 40.0) == pytest.approx(7.40, abs=0.005)
    assert phys.hco3_from(7.40, 40.0) == pytest.approx(24.0, abs=0.05)
    # a fixed bicarbonate with a rising CO2 must acidise, and not linearly
    low = phys.ph_from(24.0, 40.0)
    high = phys.ph_from(24.0, 80.0)
    assert high < 7.25 and high > low - 0.35


def test_base_excess_anchors():
    # standard base excess at 15 g/dL is ~0 when everything is normal
    assert phys.base_excess(24.0, 7.40, 15.0) == pytest.approx(0.0, abs=0.2)
    # HCO3 14 with the pH that goes with it: a metabolic acidosis of about -13
    assert phys.base_excess(14.0, phys.ph_from(14.0, 40.0), 15.0) == pytest.approx(-13.3, abs=0.6)
    # base excess is haemoglobin corrected: buffering comes from blood, not plasma
    be_low_hb = phys.base_excess(14.0, phys.ph_from(14.0, 40.0), 7.0)
    be_hi_hb = phys.base_excess(14.0, phys.ph_from(14.0, 40.0), 20.0)
    assert be_low_hb > be_hi_hb


def test_p50_shifts_the_way_it_should():
    normal = phys.p50()
    assert normal == pytest.approx(26.6, abs=0.3)
    assert phys.p50(ph=7.2) > normal                 # acidosis: right shift
    assert phys.p50(ph=7.6) < normal                 # alkalosis: left shift
    assert phys.p50(temp_c=40.0) > normal            # fever: right shift
    assert phys.p50(pco2=70.0) > normal              # hypercapnia (Bohr)


def test_o2_saturation_curve():
    assert phys.so2(26.6) == pytest.approx(0.50, abs=0.01)
    assert phys.so2(60.0) == pytest.approx(0.90, abs=0.01)
    assert phys.so2(100.0) == pytest.approx(0.97, abs=0.01)
    assert phys.so2(0.0) == pytest.approx(0.0, abs=1e-6)
    # acidotic, hypercapnic blood holds on to oxygen: a right shift costs
    # saturation at any given PO2
    assert phys.so2(40.0, ph=7.2, pco2=60.0) < phys.so2(40.0)
    assert phys.so2(40.0, p50_mm=phys.p50(7.2, 60.0)) == pytest.approx(
        phys.so2(40.0, ph=7.2, pco2=60.0), abs=1e-12)


def test_contents_and_the_curve_inverse():
    content = phys.o2_content(100.0, 15.0)
    assert content == pytest.approx(0.199, abs=0.003)      # mL O2 per mL blood
    back, sat = phys.content_to_o2(content, 15.0)
    assert back == pytest.approx(100.0, abs=0.5)
    assert sat == pytest.approx(phys.so2(100.0), abs=1e-3)


def test_co2_content_is_linear_in_pco2():
    a = phys.co2_content(40.0)
    b = phys.co2_content(80.0)
    assert b == pytest.approx(2 * a, rel=1e-9)
    assert a > 0.0


def test_alveolar_gas_equation():
    pio2 = phys.pio2(0.21)
    assert pio2 == pytest.approx(149.7, abs=0.5)
    assert phys.pio2(1.0) == pytest.approx(713.0, abs=0.5)
    pao2 = phys.alveolar_gas(pio2, 40.0, 0.8, 0.21)
    assert pao2 == pytest.approx(101.8, abs=1.0)
    # on pure oxygen the equation collapses towards PiO2 - PaCO2
    assert phys.alveolar_gas(phys.pio2(1.0), 40.0, 0.8, 1.0) == pytest.approx(673.0, abs=1.0)


def test_alveolar_ventilation_pco2():
    # PaCO2 = 0.863 * V'CO2 / V'A  for the classic 200 mL/min at 4 L/min
    assert phys.minute_to_alveolar_pco2(4.0, 200.0) == pytest.approx(43.2, abs=0.5)


def test_anion_gap():
    assert phys.anion_gap(140.0, 104.0, 24.0) == pytest.approx(12.0, abs=0.1)
    # albumin correction: a low albumin hides an anion gap
    assert phys.anion_gap(140.0, 104.0, 24.0, albumin=2.0) > phys.anion_gap(140.0, 104.0, 24.0)

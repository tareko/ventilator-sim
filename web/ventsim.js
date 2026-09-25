/* ventsim.js - JavaScript port of the Python ventsim package.
 *
 * Faithful port of ventsim/physiology.py, ventsim/lung.py, ventsim/model.py and
 * ventsim/dynamics.py.  Same constants, same iteration counts, same damping,
 * same order of arithmetic, same clamps.  Contents are always mL of gas per
 * mL of whole blood (0.20 == 20 mL/dL).
 *
 * Usable as a browser <script src="ventsim.js"> (defines window.Ventsim) and
 * under Node (module.exports).  No dependencies, no bundler.
 */
(function (globalRoot) {
  'use strict';

  var M = {};

  /* ------------------------------------------------------------------ utils */
  function clamp(x, lo, hi) { return x < lo ? lo : (x > hi ? hi : x); }
  /* Python int(round(x)): banker's rounding on exact halves, truncation toward
     zero for int().  Step counts must match Python bit for bit. */
  function pyRound(x) {
    var f = Math.floor(x), d = x - f, r;
    if (d > 0.5) r = f + 1;
    else if (d < 0.5) r = f;
    else r = (f % 2 === 0) ? f : f + 1;
    return r;
  }
  function pyInt(x) { return x < 0 ? Math.ceil(x) : Math.floor(x); }
  function sum(arr) { var s = 0, i; for (i = 0; i < arr.length; i++) s += arr[i]; return s; }

  /* =========================================================== physiology */
  var PHYS = {
    PB: 760.0,
    PH2O: 47.0,
    DRY: 760.0 - 47.0,                    /* 713 mmHg */
    ALPHA_O2: 0.00003,
    HUNTER: 1.34,
    BETA_CO2: 0.0121,
    PK_CO2: 6.101,
    ALPHA_CO2: 0.0301,
    HILL_N: 2.7,
    P50_NORMAL: 26.6
  };
  /* BTPS -> STPD, computed exactly as the Python module-level constant does */
  PHYS.BTPS_TO_STPD = (PHYS.DRY / PHYS.PB) * (273.15 / 310.35);

  function phFrom(hco3, pco2) {
    hco3 = Math.max(hco3, 1e-4);
    pco2 = Math.max(pco2, 1e-4);
    return PHYS.PK_CO2 + Math.log10(hco3 / (PHYS.ALPHA_CO2 * pco2));
  }

  function hco3From(ph, pco2) {
    return Math.pow(10.0, ph - PHYS.PK_CO2) * PHYS.ALPHA_CO2 * pco2;
  }

  function baseExcess(hco3, ph, hgb) {
    if (hgb === undefined || hgb === null) hgb = 15.0;
    return (1.0 - 0.014 * hgb) * (hco3 - 24.0 + (1.43 * hgb + 7.7) * (ph - 7.4));
  }

  function anionGap(na, cl, hco3, albumin) {
    if (albumin === undefined || albumin === null) albumin = 4.0;
    return (na - cl - hco3) + 2.5 * (4.0 - albumin);
  }

  function p50(ph, pco2, tempC) {
    if (ph === undefined || ph === null) ph = 7.4;
    if (pco2 === undefined || pco2 === null) pco2 = 40.0;
    if (tempC === undefined || tempC === null) tempC = 37.0;
    return PHYS.P50_NORMAL * Math.pow(
      10.0,
      -0.44 * (ph - 7.4) + 0.0010 * (pco2 - 40.0) + 0.024 * (tempC - 37.0)
    );
  }

  function so2(po2, p50mm, ph, pco2, tempC) {
    /* No explicit P50: shift the curve for pH/CO2/temperature, as Python does. */
    po2 = Math.max(po2, 0.0);
    var p50_ = (p50mm === undefined || p50mm === null) ? p50(ph, pco2, tempC) : p50mm;
    if (po2 <= 0.0) return 0.0;
    var x = Math.pow(po2, PHYS.HILL_N);
    return x / (Math.pow(p50_, PHYS.HILL_N) + x);
  }

  function o2Content(po2, hgb, ph, pco2, tempC) {
    if (hgb === undefined || hgb === null) hgb = 15.0;
    var s = so2(po2, null, ph, pco2, tempC);
    return PHYS.HUNTER * (hgb / 100.0) * s + PHYS.ALPHA_O2 * po2;
  }

  /* Invert the O2-Hb curve: returns [po2, so2] for a content value. */
  function contentToO2(content, hgb, ph, pco2, tempC, po2Hi) {
    if (hgb === undefined || hgb === null) hgb = 15.0;
    if (po2Hi === undefined || po2Hi === null) po2Hi = 3000.0;
    if (content <= 0.0) return [0.0, 0.0];
    var c = function (po2) { return o2Content(po2, hgb, ph, pco2, tempC); };
    var lo = 0.0, hi = po2Hi;
    if (c(hi) < content) {                       /* supersaturated */
      var boundCap = PHYS.HUNTER * (hgb / 100.0);
      return [Math.max((content - boundCap) / PHYS.ALPHA_O2, hi), 1.0];
    }
    for (var _i = 0; _i < 60; _i++) {
      var mid = 0.5 * (lo + hi);
      if (c(mid) < content) lo = mid; else hi = mid;
    }
    var po2 = 0.5 * (lo + hi);
    return [po2, so2(po2, null, ph, pco2, tempC)];
  }

  function oxygenDelivery(coLMin, cao2) { return coLMin * 1000.0 * cao2; }

  function co2Content(pco2) { return PHYS.BETA_CO2 * Math.max(pco2, 0.0); }
  function pco2FromContent(content) { return Math.max(content, 0.0) / PHYS.BETA_CO2; }

  function pio2(fio2, pb) {
    if (pb === undefined || pb === null) pb = PHYS.PB;
    return fio2 * (pb - PHYS.PH2O);
  }

  function alveolarGas(pio2Val, pco2, rq, fio2) {
    if (rq === undefined || rq === null) rq = 0.8;
    if (fio2 === undefined || fio2 === null) fio2 = 0.21;
    var coef = (1.0 - fio2 * (1.0 - rq)) / rq;
    return Math.max(pio2Val - pco2 * coef, 0.0);
  }

  function minuteToAlveolarPco2(vaLMin, vco2MlMin) {
    if (vaLMin <= 1e-6) return 200.0;
    return 0.863 * vco2MlMin / vaLMin;
  }

  /* ================================================================ lung */
  var MAX_PACO2 = 160.0;

  function makeLung(opts) {
    var o = opts || {};
    var g = function (name, dflt) {
      return (o[name] === undefined) ? dflt : o[name];
    };
    var L = {
      n_units: pyInt(Math.max(4, pyInt(g('n_units', 24)))),
      severity: g('severity', 0.0),
      affected_fraction: (o.affected_fraction === undefined) ? null : o.affected_fraction,
      heterogeneity: g('heterogeneity', 0.30),
      perfusion_gradient: g('perfusion_gradient', 1.0),
      vascular_loss: g('vascular_loss', 0.0),
      healthy_closing: g('healthy_closing', 3.0),
      affected_closing_lo: g('affected_closing_lo', 7.0),
      affected_closing_hi: g('affected_closing_hi', 30.0),
      od_lo_healthy: g('od_lo_healthy', 17.0),
      od_span_healthy: g('od_span_healthy', 14.0),
      od_lo_affected: g('od_lo_affected', 30.0),
      od_span_affected: g('od_span_affected', 25.0),
      stiff_affected: g('stiff_affected', 0.7),
      recruit_width: g('recruit_width', 5.0),
      base_shunt: g('base_shunt', 0.02),
      max_closing: (o.max_closing === undefined) ? null : o.max_closing
    };
    var sev = clamp(L.severity, 0.0, 1.0);
    L.affected_frac = (L.affected_fraction === null)
      ? (0.55 * sev)
      : clamp(L.affected_fraction, 0.0, 1.0);
    L.build = function () { lungBuild(L); };
    L.build();
    return L;
  }

  function lungBuild(L) {
    var n = L.n_units, af = L.affected_frac;
    var edge = Math.max(1.5 / n, 0.02);          /* soft healthy/diseased edge */
    var start = 1.0 - af;
    L.z = new Array(n); L.perf = new Array(n); L.p_open = new Array(n);
    L.od_limit = new Array(n); L.stiff = new Array(n);
    L.vent_het = new Array(n); L.aff = new Array(n);
    var raw = new Array(n);
    L.perf_scale = new Array(n);
    var tot = 0.0;
    for (var i = 0; i < n; i++) {
      var z = (i + 0.5) / n;
      L.z[i] = z;
      var mask = 0.5 + 0.5 * Math.sin((i + 1) * 1.7138750);
      var scale = clamp(1.0 - 2.0 * L.vascular_loss * mask, 0.02, 1.0);
      L.perf_scale[i] = scale;
      var p = (0.55 + L.perfusion_gradient * z) * scale;
      raw[i] = p;
      tot += p;
      var aff = clamp((z - start) / edge, 0.0, 1.0);
      L.aff[i] = aff;

      var healthyClose = L.healthy_closing * (0.25 + 0.75 * z);
      var span = Math.max(af, 1e-6);
      var disClose = (L.affected_closing_lo
        + (L.affected_closing_hi - L.affected_closing_lo)
          * clamp((z - start) / span, 0.0, 1.0));
      if (L.max_closing !== null) disClose = L.max_closing * z;
      L.p_open[i] = healthyClose + aff * (disClose - healthyClose);

      var odH = L.od_lo_healthy + L.od_span_healthy * z;
      var odA = L.od_lo_affected + L.od_span_affected * z;
      L.od_limit[i] = odH + aff * (odA - odH);
      L.stiff[i] = 1.0 + aff * (L.stiff_affected - 1.0);
      /* golden-angle scatter: deterministic, identical formula in Python */
      L.vent_het[i] = 1.0 + L.heterogeneity * Math.sin((i + 1) * 2.3999632);
    }
    for (var j = 0; j < n; j++) L.perf[j] = raw[j] / tot;
  }

  /* PO2 of one unit from its O2 mass balance, by bisection. */
  function solveUnitPo2(a, cvo2, fio2, bCap, p50mm, dry, alpha) {
    var k = alpha + a / dry;
    var c0 = cvo2 + a * fio2;
    var hi = fio2 * dry;
    if (hi <= 0.0) return 0.0;
    var g = function (p) {
      var sRhs = (c0 - k * p) / bCap;
      if (p <= 0.0) return -sRhs;
      var lnX = PHYS.HILL_N * Math.log(p / p50mm);
      var sH;
      if (lnX > 40.0) sH = 1.0;
      else if (lnX < -40.0) sH = Math.exp(lnX);
      else { var e = Math.exp(lnX); sH = e / (1.0 + e); }
      return sH - sRhs;
    };
    var loP = 0.0, hiP = hi;
    if (g(hiP) < 0.0) return hiP;              /* pathological parameters */
    for (var _i = 0; _i < 26; _i++) {
      var mid = 0.5 * (loP + hiP);
      if (g(mid) < 0.0) loP = mid; else hiP = mid;
    }
    return 0.5 * (loP + hiP);
  }

  /* vaMl: alveolar ventilation per unit, mL/min BTPS. perf: any scale. */
  function exchange(lung, vaMl, perf, opts) {
    opts = opts || {};
    var vco2MlMin = opts.vco2_ml_min;
    var coMlMin = opts.co_ml_min;
    var fio2 = opts.fio2;
    var hgb = (opts.hgb === undefined) ? 15.0 : opts.hgb;
    var ph = (opts.ph === undefined) ? 7.40 : opts.ph;
    var rq = (opts.rq === undefined) ? 0.8 : opts.rq;
    var tempC = (opts.temp_c === undefined) ? 37.0 : opts.temp_c;
    var recruited = (opts.recruited === undefined) ? null : opts.recruited;
    var overdist = (opts.overdist === undefined) ? null : opts.overdist;
    var extraShunt = (opts.extra_shunt === undefined) ? 0.0 : opts.extra_shunt;

    var n = lung.n_units;
    if (vaMl.length !== n || perf.length !== n) {
      throw new Error('per-unit vectors must match lung.n_units');
    }
    var totPerf = sum(perf);
    if (totPerf <= 0) throw new Error('perfusion distribution is empty');
    var w = perf.map(function (p) { return p / totPerf; });

    var fs = clamp(lung.base_shunt + extraShunt, 0.0, 0.5);
    var wEff = w.map(function (p) { return p * (1.0 - fs); });

    var co = Math.max(coMlMin, 500.0);
    var vo2MlMin = vco2MlMin / rq;
    var kO2 = vo2MlMin / co;
    var beta = PHYS.BETA_CO2;
    var dry = PHYS.DRY;
    var denom = beta * dry;

    var vst = vaMl.map(function (v) { return Math.max(v, 0.0) * PHYS.BTPS_TO_STPD; });

    /* ---- CO2, closed form ---- */
    var aList = new Array(n), qList = new Array(n), shuntI = new Array(n);
    var sSum = fs, shunt = fs;
    for (var i = 0; i < n; i++) {
      var q = wEff[i] * co;
      qList[i] = q;
      var a = (q > 1e-9) ? (vst[i] / (q * denom)) : Infinity;
      aList[i] = a;
      sSum += wEff[i] / (1.0 + a);
      if (q > 1e-9 && vaMl[i] / q < 0.05) { shuntI[i] = true; shunt += wEff[i]; }
      else shuntI[i] = false;
    }
    var kCo2 = vco2MlMin / (co * beta);
    var paco2;
    if (sSum >= 0.999) paco2 = MAX_PACO2;
    else paco2 = Math.min(sSum * kCo2 / (1.0 - sSum), MAX_PACO2);
    var pvco2 = paco2 + kCo2;

    /* ---- O2: mass balance per unit, venous content iterated ---- */
    var bCap = PHYS.HUNTER * (hgb / 100.0);
    var alpha = PHYS.ALPHA_O2;
    var p50mm = p50(ph, pvco2, tempC);
    var vented = new Array(n);
    for (var v0 = 0; v0 < n; v0++) {
      vented[v0] = qList[v0] > 1e-9 && vaMl[v0] > 1e-9 && !shuntI[v0];
    }
    var cvo2 = Math.max(Math.min(0.13, bCap * 0.65), 0.02);
    var pList = new Array(n);
    for (var v1 = 0; v1 < n; v1++) pList[v1] = 0.0;
    var cao2 = 0.15;
    for (var it = 0; it < 12; it++) {
      cao2 = shunt * cvo2;
      for (var i2 = 0; i2 < n; i2++) {
        if (!vented[i2]) continue;
        pList[i2] = solveUnitPo2(vst[i2] / qList[i2], cvo2, fio2,
          bCap, p50mm, dry, alpha);
        cao2 += wEff[i2] * o2Content(pList[i2], hgb, ph, pvco2, tempC);
      }
      var cvo2Next = Math.max(cao2 - kO2, 0.01);
      if (Math.abs(cvo2Next - cvo2) < 1e-7) { cvo2 = cvo2Next; break; }
      cvo2 = cvo2 + 0.55 * (cvo2Next - cvo2);
    }

    var pvo2 = contentToO2(cvo2, hgb, ph, pvco2, tempC)[0];
    for (var i3 = 0; i3 < n; i3++) {
      if (!vented[i3]) {
        pList[i3] = qList[i3] > 1e-9 ? pvo2 : fio2 * dry;
      }
    }

    var art = contentToO2(cao2, hgb, ph, paco2, tempC);
    var pao2 = art[0], sao2 = art[1];

    /* ---- per-unit report rows ---- */
    var ventTotal = sum(vaMl);
    var vqHigh = 0.0, vqLow = 0.0;
    var rows = [];
    for (var i4 = 0; i4 < n; i4++) {
      var q4 = qList[i4];
      var vq = q4 > 1e-9 ? vaMl[i4] / q4 : Infinity;
      if (vaMl[i4] > 1e-9 && ventTotal > 0) {
        var frac = vaMl[i4] / ventTotal;
        if (vq > 2.0) vqHigh += frac;
        else if (vq < 0.6) vqLow += frac;
      }
      rows.push({
        z: lung.z[i4], perf: wEff[i4], va: vaMl[i4], vq: vq,
        pc_co2: (aList[i4] === Infinity) ? 0.0 : pvco2 / (1.0 + aList[i4]),
        po2: pList[i4],
        cc_o2: o2Content(pList[i4], hgb, ph, pvco2, tempC),
        shunt: shuntI[i4],
        diseased: lung.aff[i4],
        recruited: (recruited === null) ? 1.0 : recruited[i4],
        overdist: (overdist === null) ? 0.0 : overdist[i4]
      });
    }

    return {
      pao2: pao2, paco2: paco2, sao2: sao2, cao2: cao2, cvo2: cvo2,
      pvco2: pvco2, pvo2: pvo2, vo2_ml_min: vo2MlMin, shunt_frac: shunt,
      vq_high: vqHigh, vq_low: vqLow, units: rows
    };
  }

  /* =============================================================== model */
  var MODES = ['vc', 'pc', 'psv'];

  var ANATOMIC_VD_PER_KG = 2.2;
  var FRC_PER_KG = 0.035;
  var OD_WIDTH = 12.0;
  var OD_STIFFEN = 0.5;
  var OD_CAPILLARY = 0.75;
  var NONDEP_TRANSMISSION = 0.5;
  var OD_CO = 0.30;
  var OD_MEAN_BASE = 14.0;
  var OD_MEAN_CO = 0.25;
  var RECRUIT_TIDAL = 0.35;
  var DRIVE_SLOPE = 1.05;
  var HYPOXIC_GAIN = 0.30;
  var PMUS_MAX = 12.0;
  var PMUS_CAP = 18.0;
  var ACUTE_HCO3_PER_MMHG = 0.10;
  var CHRONIC_HCO3_PER_MMHG = 0.35;
  var CHRONIC_HYPO_PER_MMHG = 0.10;
  var TAU_ACUTE_MIN = 0.6;
  var TAU_RENAL_MIN = 720.0;
  var BUFFER_SPACE_FRACTION = 0.5;

  function ibwKg(heightCm, sex) {
    if (sex === undefined || sex === null) sex = 'm';
    if (String(sex).toLowerCase().indexOf('f') === 0) {
      return 45.5 + 0.91 * (heightCm - 152.4);
    }
    return 50.0 + 0.91 * (heightCm - 152.4);
  }

  /* ---------------------------------------------------- vent / patient */
  var VENT_FIELDS = ['mode', 'fio2', 'peep', 'vt_ml_kg', 'rr', 'ie', 'dp', 'ti_s'];

  function Vent(opts) {
    var o = opts || {};
    this.mode = (o.mode === undefined) ? 'vc' : o.mode;
    this.fio2 = (o.fio2 === undefined) ? 0.5 : o.fio2;
    this.peep = (o.peep === undefined) ? 5.0 : o.peep;
    this.vt_ml_kg = (o.vt_ml_kg === undefined) ? 8.0 : o.vt_ml_kg;
    this.rr = (o.rr === undefined) ? 12.0 : o.rr;
    this.ie = (o.ie === undefined) ? 2.0 : o.ie;
    this.dp = (o.dp === undefined) ? 12.0 : o.dp;
    this.ti_s = (o.ti_s === undefined) ? null : o.ti_s;
  }
  Vent.prototype.clamped = function () {
    var v = new Vent(this);
    v.mode = String(v.mode).toLowerCase().trim();
    if (MODES.indexOf(v.mode) < 0) {
      throw new Error('mode must be one of ' + MODES + ', got ' + JSON.stringify(this.mode));
    }
    v.fio2 = clamp(v.fio2, 0.21, 1.0);
    v.peep = clamp(v.peep, 0.0, 45.0);
    v.vt_ml_kg = clamp(v.vt_ml_kg, 2.0, 15.0);
    v.rr = clamp(v.rr, 1.0, 60.0);
    v.ie = clamp(v.ie, 0.35, 10.0);
    v.dp = clamp(v.dp, 0.0, 60.0);
    if (v.ti_s !== null && v.ti_s !== undefined) v.ti_s = clamp(v.ti_s, 0.3, 4.0);
    return v;
  };
  Vent.prototype.with = function (kw) {
    var v = new Vent(this);
    if (kw) for (var k in kw) if (Object.prototype.hasOwnProperty.call(kw, k)) v[k] = kw[k];
    return v;
  };
  Object.defineProperty(Vent.prototype, 'label', {
    get: function () {
      return { vc: 'VC-AC', pc: 'PC-AC', psv: 'PSV' }[this.mode];
    }
  });
  Vent.prototype.summary = function (ibw) {
    var mid;
    if (this.mode === 'vc') {
      mid = 'Vt ' + this.vt_ml_kg.toFixed(1) + ' ml/kg = ' + Math.round(this.vt_ml_kg * ibw) +
        ' ml   RR ' + Math.round(this.rr) + '   I:E 1:' + this.ie;
    } else {
      var tag = this.mode === 'psv' ? 'PS' : 'dP';
      mid = tag + ' ' + Math.round(this.dp) + ' cmH2O   RR ' + Math.round(this.rr) +
        '   I:E 1:' + this.ie;
    }
    return this.label + '   FiO2 ' + this.fio2.toFixed(2) + '   PEEP ' + Math.round(this.peep) +
      '   ' + mid;
  };

  var PATIENT_FIELDS = ['name', 'sex', 'height_cm', 'weight_kg', 'compliance', 'resistance',
    'exp_limit', 'severity', 'affected_fraction', 'heterogeneity', 'perfusion_gradient',
    'vascular_loss', 'base_shunt', 'pfo', 'lung_overrides', 'hpv', 'n_units', 'hgb',
    'co_l_min', 'rv_sensitivity', 'vco2', 'rq', 'temp_c', 'hco3', 'renal', 'acid_mmol_h',
    'bicarb_mmol_h', 'sedation', 'set_point', 'muscle', 've_rest', 'rr_rest', 'frc_l',
    'notes'];

  function Patient(opts) {
    var o = opts || {};
    var g = function (name, dflt) {
      return (o[name] === undefined) ? dflt : o[name];
    };
    this.name = g('name', 'adult');
    this.sex = g('sex', 'm');
    this.height_cm = g('height_cm', 172.0);
    this.weight_kg = g('weight_kg', 75.0);
    this.compliance = g('compliance', 60.0);
    this.resistance = g('resistance', 5.0);
    this.exp_limit = g('exp_limit', 1.0);
    this.severity = g('severity', 0.0);
    this.affected_fraction = (o.affected_fraction === undefined) ? null : o.affected_fraction;
    this.heterogeneity = g('heterogeneity', 0.30);
    this.perfusion_gradient = g('perfusion_gradient', 1.0);
    this.vascular_loss = g('vascular_loss', 0.0);
    this.base_shunt = g('base_shunt', 0.02);
    this.pfo = g('pfo', 0.0);
    this.lung_overrides = g('lung_overrides', {});
    this.hpv = g('hpv', 0.35);
    this.n_units = g('n_units', 24);
    this.hgb = g('hgb', 15.0);
    this.co_l_min = g('co_l_min', 5.0);
    this.rv_sensitivity = g('rv_sensitivity', 1.0);
    this.vco2 = (o.vco2 === undefined) ? null : o.vco2;
    this.rq = g('rq', 0.8);
    this.temp_c = g('temp_c', 37.0);
    this.hco3 = g('hco3', 24.0);
    this.renal = g('renal', 1.0);
    this.acid_mmol_h = g('acid_mmol_h', 0.0);
    this.bicarb_mmol_h = g('bicarb_mmol_h', 0.0);
    this.sedation = g('sedation', 1.0);
    this.set_point = g('set_point', 34.0);
    this.muscle = g('muscle', 1.0);
    this.ve_rest = (o.ve_rest === undefined) ? null : o.ve_rest;
    this.rr_rest = g('rr_rest', 15.0);
    this.frc_l = (o.frc_l === undefined) ? null : o.frc_l;
    this.notes = g('notes', '');
    this._init = {};
    for (var i = 0; i < PATIENT_FIELDS.length; i++) {
      var f = PATIENT_FIELDS[i];
      this._init[f] = this[f];
    }
    this.postInit();
  }
  Patient.prototype.postInit = function () {
    this.ibw = ibwKg(this.height_cm, this.sex);
    if (this.vco2 === null) this.vco2 = 3.3 * this.ibw;
    if (this.frc_l === null) this.frc_l = FRC_PER_KG * this.ibw;
    if (this.ve_rest === null) this.ve_rest = 0.105 * this.ibw;
    var lo = this.lung_overrides || {};
    var lungOpts = {
      n_units: this.n_units, severity: this.severity,
      affected_fraction: this.affected_fraction,
      heterogeneity: this.heterogeneity,
      perfusion_gradient: this.perfusion_gradient,
      vascular_loss: this.vascular_loss,
      base_shunt: this.base_shunt
    };
    for (var k in lo) {
      if (Object.prototype.hasOwnProperty.call(lo, k)) lungOpts[k] = lo[k];
    }
    this.lung = makeLung(lungOpts);
  };
  Object.defineProperty(Patient.prototype, 'weight', {
    get: function () { return this.weight_kg > 0 ? this.weight_kg : this.ibw; }
  });
  Object.defineProperty(Patient.prototype, 'buffer_space_l', {
    get: function () { return BUFFER_SPACE_FRACTION * Math.max(this.weight_kg, 1.0); }
  });
  /* equivalent of dataclasses.replace(patient, **kw) */
  Patient.prototype.with = function (kw) {
    var copy = {};
    for (var k in this._init) {
      if (Object.prototype.hasOwnProperty.call(this._init, k)) {
        var v = this._init[k];
        copy[k] = (v && typeof v === 'object') ? JSON.parse(JSON.stringify(v)) : v;
      }
    }
    if (kw) for (var f in kw) {
      if (Object.prototype.hasOwnProperty.call(kw, f)) copy[f] = kw[f];
    }
    return new Patient(copy);
  };

  function makePatient(obj) { return new Patient(obj || {}); }
  function makeVent(obj) { return new Vent(obj || {}); }

  /* ------------------------------------------------------ respiratory drive */
  function minuteTarget(paco2, sao2v, patient) {
    if (patient.sedation >= 1.0) return 0.0;
    var ve = DRIVE_SLOPE * Math.max(0.0, paco2 - patient.set_point);
    var hyp = HYPOXIC_GAIN * Math.max(0.0, 0.90 - sao2v) / 0.08;
    return Math.max(0.0, (1.0 - patient.sedation) * ve * (1.0 + hyp));
  }

  function driveIndex(paco2, sao2v, patient) {
    return minuteTarget(paco2, sao2v, patient) / Math.max(patient.ve_rest, 0.5);
  }

  /* --------------------------------------------------- the operating point */
  function meanAirwayPressure(peepT, elastic, ti, te, tt, vtL, tauExp, resistance) {
    if (tt <= 1e-9) return peepT;
    var flow = ti > 1e-9 ? vtL / ti : 0.0;
    var resistive = resistance * flow;
    var meanInsp = peepT + resistive + 0.5 * elastic;
    var frac;
    if (te > 1e-9 && tauExp > 1e-9) {
      frac = (tauExp / te) * (1.0 - Math.exp(-Math.min(te / tauExp, 40.0)));
    } else {
      frac = 1.0;
    }
    var meanExp = peepT + elastic * frac;
    return (ti * meanInsp + te * meanExp) / tt;
  }

  function evaluate(patient, vent, opts) {
    opts = opts || {};
    var hco3 = opts.hco3;
    var phForO2 = opts.phForO2;
    var drive = opts.drive;
    var tMin = (opts.t_min === undefined) ? 0.0 : opts.t_min;
    var keepUnits = !!(opts.keep_units || opts.keepUnits);   /* both spellings */
    var state = (opts.state === undefined) ? null : opts.state;

    vent = vent.clamped();
    var lung = patient.lung;
    var n = lung.n_units;
    var ibw = patient.ibw;
    var vco2 = patient.vco2;

    var stiffSum = sum(lung.stiff);
    if (!stiffSum) stiffSum = 1.0;
    var c0 = new Array(n);
    for (var i0 = 0; i0 < n; i0++) {
      c0[i0] = patient.compliance / 1000.0 * lung.stiff[i0] / stiffSum;
    }

    /* ---- breath timing ---- */
    var veTarget = drive * patient.ve_rest;
    var rrEff, vtCap, dpEff, pmus;
    if (vent.mode === 'psv' && veTarget > 1e-6) {
      var rrPat = clamp(patient.rr_rest * Math.pow(veTarget / patient.ve_rest, 0.6), 4.0, 45.0);
      rrEff = Math.max(rrPat, vent.rr);
      vtCap = rrPat > 1e-6 ? veTarget / rrPat : null;
      pmus = clamp(PMUS_MAX * (drive - 1.0), 0.0, PMUS_CAP) * patient.muscle;
      dpEff = vent.dp + pmus;
    } else {
      rrEff = vent.rr;
      vtCap = null;
      dpEff = vent.dp;
    }
    var tt = 60.0 / rrEff;

    var opened = new Array(n), od = new Array(n), stress = new Array(n);
    for (var i1 = 0; i1 < n; i1++) { opened[i1] = 1.0; od[i1] = 0.0; stress[i1] = 0.0; }
    var shares = c0.slice();
    var sumShares = sum(c0) || 1e-12;
    var vtL = vent.vt_ml_kg * ibw / 1000.0;
    var ti = Math.min(Math.max(vent.ti_s || tt / (1.0 + vent.ie), 0.25),
      Math.max(tt - 0.15, 0.25));
    var te = tt - ti;
    var cTot = patient.compliance / 1000.0;
    var auto = 0.0;
    var peepT = vent.peep;
    var plat = 0.0, elastic = 0.0, dpRecruit = 0.0;
    var tau = 0.0, tauExp = 0.0;

    for (var loop = 0; loop < 24; loop++) {
      var cs = 0.0;
      for (var a0 = 0; a0 < n; a0++) cs += c0[a0] * opened[a0] * (1.0 - OD_STIFFEN * od[a0]);
      cTot = Math.max(cs, 1e-5);
      tau = patient.resistance * cTot;
      tauExp = tau * Math.max(patient.exp_limit, 1e-3);
      if (vent.mode === 'psv') {
        ti = clamp(0.45 + 1.4 * vtL, 0.45, Math.max(0.85 * tt, 0.3));
      } else {
        ti = clamp(vent.ti_s || tt / (1.0 + vent.ie), 0.25, Math.max(tt - 0.15, 0.25));
      }
      te = Math.max(tt - ti, 1e-4);
      if (vtL > 0) {
        auto = (vtL / cTot) / Math.max(Math.exp(Math.min(te / tauExp, 40.0)) - 1.0, 1e-9);
      } else {
        auto = 0.0;
      }
      auto = Math.min(auto, 40.0);
      peepT = vent.peep + auto;

      if (vent.mode === 'vc') {
        vtL = vent.vt_ml_kg * ibw / 1000.0;
      } else {
        vtL = dpEff * cTot * (tau > 1e-9 ? (1.0 - Math.exp(-ti / tau)) : 1.0);
        if (vtCap !== null) vtL = Math.min(vtL, vtCap);
      }

      elastic = vtL / cTot;
      plat = peepT + elastic;
      dpRecruit = peepT + RECRUIT_TIDAL * elastic;

      for (var i2 = 0; i2 < n; i2++) {
        var tgtOpen = clamp((dpRecruit - lung.p_open[i2]) / lung.recruit_width, 0.0, 1.0);
        opened[i2] += 0.4 * (tgtOpen - opened[i2]);
      }
      shares = new Array(n);
      var sSum2 = 0.0;
      for (var i3 = 0; i3 < n; i3++) {
        shares[i3] = c0[i3] * opened[i3] * (1.0 - OD_STIFFEN * od[i3]) * lung.vent_het[i3];
        sSum2 += shares[i3];
      }
      sumShares = Math.max(sSum2, 1e-12);
      for (var i4 = 0; i4 < n; i4++) {
        stress[i4] = vtL * shares[i4] / sumShares / Math.max(c0[i4], 1e-12);
        var dist = peepT * (1.0 - NONDEP_TRANSMISSION * lung.z[i4]);
        var target = clamp((stress[i4] + dist - lung.od_limit[i4]) / OD_WIDTH, 0.0, 1.0);
        if (opened[i4] <= 1e-9) target = 0.0;
        od[i4] += 0.5 * (target - od[i4]);
      }
    }

    /* ---- perfusion ---- */
    var meanOd = 0.0;
    for (var i5 = 0; i5 < n; i5++) meanOd += lung.perf[i5] * od[i5];
    var perf = new Array(n), perfSum = 0.0;
    for (var i6 = 0; i6 < n; i6++) {
      perf[i6] = lung.perf[i6] * (1.0 - OD_CAPILLARY * od[i6])
        * (1.0 - patient.hpv * (1.0 - opened[i6]));
      perfSum += perf[i6];
    }
    if (perfSum <= 1e-9) perf = lung.perf.slice();

    /* ---- ventilation distribution ---- */
    var vtDl = ANATOMIC_VD_PER_KG * ibw / 1000.0;
    var veLMin = rrEff * vtL;
    var vaLMin = Math.max(rrEff * Math.max(vtL - vtDl, 0.0), 1e-4);
    var vaMl = new Array(n);
    for (var i7 = 0; i7 < n; i7++) vaMl[i7] = vaLMin * 1000.0 * shares[i7] / sumShares;

    var pmean = meanAirwayPressure(peepT, elastic, ti, te, tt, vtL, tauExp, patient.resistance);
    var peakFlow = ti > 1e-9 ? (vtL / ti * 60.0) : 0.0;
    var resistive = patient.resistance * peakFlow / 60.0;
    var pip = plat + resistive;
    var coLMin = patient.co_l_min * (
      1.0 - patient.rv_sensitivity * (OD_CO * meanOd
        + OD_MEAN_CO * Math.max(pmean - OD_MEAN_BASE, 0.0) / 10.0));
    coLMin = Math.max(coLMin, 0.40 * patient.co_l_min);

    var extraShunt = 0.0;
    if (patient.pfo > 0.0) {
      extraShunt = patient.pfo * clamp(1.0 + (pmean - 12.0) / 24.0, 0.0, 2.5);
    }

    var ex = exchange(lung, vaMl, perf, {
      vco2_ml_min: vco2, co_ml_min: coLMin * 1000.0, fio2: vent.fio2,
      hgb: patient.hgb, ph: phForO2, rq: patient.rq, temp_c: patient.temp_c,
      recruited: opened, overdist: od, extra_shunt: extraShunt
    });

    var paco2 = ex.paco2, pao2 = ex.pao2, sao2v = ex.sao2;
    if (state !== null) { paco2 = state[0]; pao2 = state[1]; sao2v = state[2]; }
    var ph = phFrom(hco3, paco2);
    var be = baseExcess(hco3, ph, patient.hgb);

    /* ---- indices ---- */
    var cDyn = elastic > 1e-9 ? (vtL / elastic * 1000.0) : 0.0;
    var dpStat = elastic;
    var peco2;
    if (veLMin > 1e-6) {
      var faco2 = vco2 / (veLMin * 1000.0 * PHYS.BTPS_TO_STPD);
      peco2 = faco2 * PHYS.PB;
    } else peco2 = 0.0;
    var vdVt = clamp(paco2 > 1e-6 ? (1.0 - peco2 / paco2) : 1.0, 0.0, 0.99);
    var pf = vent.fio2 > 0.05 ? (pao2 / vent.fio2) : 0.0;
    var oi = pao2 > 1.0 ? (vent.fio2 * pmean * 100.0 / pao2) : 9999.0;
    var recruitPct = 0.0;
    for (var i8 = 0; i8 < n; i8++) recruitPct += lung.perf[i8] * opened[i8];
    recruitPct = 100.0 * recruitPct;
    var maxStress = null;
    for (var i9 = 0; i9 < n; i9++) {
      var sv = stress[i9] * opened[i9] + peepT * (1.0 - NONDEP_TRANSMISSION * lung.z[i9]);
      if (maxStress === null || sv > maxStress) maxStress = sv;
    }
    if (maxStress === null) maxStress = 0.0;
    var mpElas = rrEff * 0.098 * 0.5 * elastic * vtL;
    var mpRes = rrEff * 0.098 * patient.resistance * vtL * vtL / Math.max(ti, 1e-6);
    var mpJMin = mpElas + mpRes;
    var cao2v = o2Content(pao2, patient.hgb, ph, paco2, patient.temp_c);
    var cvo2v = Math.max(cao2v - (vco2 / patient.rq) / Math.max(coLMin * 1000.0, 1.0), 0.02);
    var do2 = coLMin * 1000.0 * cao2v;
    var pao2Ideal = alveolarGas(pio2(vent.fio2), paco2, patient.rq, vent.fio2);
    var aa = Math.max(pao2Ideal - pao2, 0.0);

    return {
      t_min: tMin,
      ph: ph, paco2: paco2, pao2: pao2, sao2: sao2v, hco3: hco3, be: be,
      pvco2: ex.pvco2, pvo2: ex.pvo2,
      mode: vent.mode, vt_ml: vtL * 1000.0, vt_ml_kg: vtL * 1000.0 / ibw, rr: rrEff,
      ti_s: ti, ve_l_min: veLMin, va_l_min: vaLMin, vd_vt: vdVt, peco2: peco2,
      peep: vent.peep, auto_peep: auto, peep_total: peepT, pplat: plat, pip: pip,
      dp_stat: dpStat, pmean: pmean, peak_flow: peakFlow, c_dyn: cDyn,
      local_stress: maxStress, mp_j_min: mpJMin, mp_elas: mpElas, mp_res: mpRes,
      fio2: vent.fio2, pf: pf, oi: oi, shunt_pct: 100.0 * ex.shunt_frac,
      vq_low_pct: 100.0 * ex.vq_low, vq_high_pct: 100.0 * ex.vq_high,
      recruit_pct: recruitPct, aa_gradient: aa, cao2: cao2v,
      cvao2_diff: cao2v - cvo2v, do2: do2, co_l_min: coLMin, drive: drive,
      effort_ratio: (vtCap && vtCap > 1e-6) ? (vtL / vtCap) : 1.0,
      flags: [], units: keepUnits ? ex.units : []
    };
  }

  function solveSteady(patient, vent, opts) {
    opts = opts || {};
    vent = vent.clamped();
    var hco3 = (opts.hco3 === undefined || opts.hco3 === null) ? patient.hco3 : opts.hco3;
    var paco2State = (opts.paco2_state === undefined || opts.paco2_state === null)
      ? 40.0 : opts.paco2_state;
    var sao2State = (opts.sao2_state === undefined || opts.sao2_state === null)
      ? 0.97 : opts.sao2_state;
    var drive = (opts.drive === undefined) ? null : opts.drive;
    var tMin = (opts.t_min === undefined) ? 0.0 : opts.t_min;
    var keepUnits = !!(opts.keep_units || opts.keepUnits);   /* both spellings */
    var state = (opts.state === undefined) ? null : opts.state;
    var phForO2 = phFrom(hco3, Math.max(paco2State, 5.0));

    var evOpts = {
      hco3: hco3, phForO2: phForO2, t_min: tMin, keep_units: keepUnits, state: state
    };

    if (vent.mode !== 'psv' || drive !== null) {
      var dFixed = drive !== null ? drive : driveIndex(paco2State, sao2State, patient);
      evOpts.drive = dFixed;
      var resFixed = evaluate(patient, vent, evOpts);
      resFixed.flags = resultFlags(resFixed, patient);
      return resFixed;
    }

    /* spontaneous mode: solve the drive <-> PaCO2 equilibrium exactly */
    var res = driveEquilibrium(patient, vent, hco3, phForO2,
      driveIndex(paco2State, sao2State, patient),
      { t_min: tMin, keep_units: keepUnits, state: state });
    res.flags = resultFlags(res, patient);
    return res;
  }

  /* Solve drive <-> PaCO2 self-consistently by bracketing + bisection on the
     monotone gap(d) = d - driveIndex(PaCO2(d), SaO2(d)).  Port of
     model._drive_equilibrium, iteration counts included. */
  function driveEquilibrium(patient, vent, hco3, phForO2, seed, opts) {
    var tMin = (opts.t_min === undefined) ? 0.0 : opts.t_min;
    var keepUnits = !!(opts.keep_units || opts.keepUnits);   /* both spellings */
    var state = (opts.state === undefined) ? null : opts.state;
    var gap = function (d) {
      var o = {
        hco3: hco3, phForO2: phForO2, drive: d, t_min: tMin,
        keep_units: keepUnits, state: state
      };
      var r = evaluate(patient, vent, o);
      return [d - driveIndex(r.paco2, r.sao2, patient), r];
    };

    var first = gap(0.0), g0 = first[0], res0 = first[1];
    if (g0 >= 0.0) { res0.drive = 0.0; return res0; }

    var lo = 0.0, resLo = res0, hi = null, resHi = null;
    var d = Math.max(seed, 0.5), rd = null;
    for (var k = 0; k < 8; k++) {
      var gdR = gap(d), gd = gdR[0]; rd = gdR[1];
      if (gd >= 0.0) { hi = d; resHi = rd; break; }
      lo = d; resLo = rd;
      d *= 2.0;
    }
    if (hi === null) {           /* never crosses: patient's own limit */
      resHi = rd;
      resHi.drive = d;
      return resHi;
    }

    for (var m = 0; m < 20; m++) {
      var mid = 0.5 * (lo + hi);
      var gmR = gap(mid), gm = gmR[0], rm = gmR[1];
      if (gm >= 0.0) { hi = mid; resHi = rm; } else { lo = mid; resLo = rm; }
    }
    resHi.drive = hi;
    return resHi;
  }

  /* ------------------------------------------------------------- flags */
  function fmt(x, digits) {
    /* round-half-even formatting, close enough to Python's f-string */
    var m = Math.pow(10, digits || 0), v = x * m;
    var fl = Math.floor(v), diff = v - fl, r;
    if (diff > 0.5) r = fl + 1;
    else if (diff < 0.5) r = fl;
    else r = (fl % 2 === 0) ? fl : fl + 1;
    return (r / m).toFixed(digits || 0);
  }

  function resultFlags(res, patient) {
    var f = [];
    if (res.pplat > 30) f.push(['alert', 'Pplat ' + fmt(res.pplat, 0) + ' > 30 cmH2O']);
    if (res.dp_stat > 17) f.push(['alert', 'driving pressure ' + fmt(res.dp_stat, 0) + ' > 17 cmH2O']);
    else if (res.dp_stat > 15) f.push(['warn', 'driving pressure ' + fmt(res.dp_stat, 0) + ' > 15 cmH2O']);
    if (res.vt_ml_kg > 9.5) f.push(['warn', 'Vt ' + fmt(res.vt_ml_kg, 1) + ' ml/kg IBW is large']);
    if (res.local_stress > 28) {
      f.push(['warn', 'peak local stress ' + fmt(res.local_stress, 0) +
        ' cmH2O (regional overstretch)']);
    }
    if (res.auto_peep > 8) {
      f.push(['alert', 'auto-PEEP ' + fmt(res.auto_peep, 0) + ' cmH2O: breath stacking']);
    } else if (res.auto_peep > 4) {
      f.push(['warn', 'auto-PEEP ' + fmt(res.auto_peep, 0) + ' cmH2O']);
    }
    if (res.mp_j_min > 10) f.push(['warn', 'tidal mechanical power ' + fmt(res.mp_j_min, 1) + ' J/min']);
    if (res.sao2 < 0.90) {
      f.push(['alert', 'SaO2 ' + fmt(100 * res.sao2, 0) + '% - hypoxaemia']);
    }
    if (res.pao2 > 130 && res.fio2 > 0.35) {
      f.push(['warn', 'PaO2 ' + fmt(res.pao2, 0) + ' on FiO2 ' + res.fio2.toFixed(2) +
        ': hyperoxia, wean FiO2']);
    }
    if (res.ph < 7.30) f.push(['alert', 'pH ' + res.ph.toFixed(2) + ' acidemia']);
    else if (res.ph > 7.50) f.push(['alert', 'pH ' + res.ph.toFixed(2) + ' alkalemia']);
    if (res.paco2 > 55 && res.hco3 < 30) {
      f.push(['warn', 'CO2 ' + fmt(res.paco2, 0) + ' with little bicarbonate buffer']);
    }
    if (res.vd_vt > 0.65) f.push(['warn', 'Vd/Vt ' + res.vd_vt.toFixed(2) + ': high dead space']);
    if (res.co_l_min < 0.75 * patient.co_l_min) {
      f.push(['warn', 'cardiac output ' + fmt(res.co_l_min, 1) + ' L/min, down ' +
        fmt(100 * (1 - res.co_l_min / patient.co_l_min), 0) + '% from pressure effects']);
    }
    if (res.drive > 1.8) {
      f.push(['warn', 'respiratory drive ' + fmt(res.drive, 1) +
        'x resting: high work of breathing']);
    }
    if (res.mode === 'psv' && res.effort_ratio < 0.97) {
      f.push(['warn', 'support-limited: only ' + fmt(100 * res.effort_ratio, 0) +
        '% of the targeted tidal volume is delivered']);
    }
    if (patient.severity >= 0.25 && res.vt_ml_kg > 8.0) {
      f.push(['warn', 'stiff lung at ' + fmt(res.vt_ml_kg, 1) +
        ' ml/kg - protective ceiling is 8']);
    }
    if (res.sao2 >= 0.90 && res.pao2 >= 62 && 7.35 <= res.ph && res.ph <= 7.48
      && res.dp_stat <= 15 && res.auto_peep <= 4) {
      f.push(['ok', 'oxygenation, ventilation and stress targets all met']);
    }
    return f;
  }

  /* =========================================================== dynamics */
  var MAX_STEPS = 400000;
  var TAU_DRIVE_MIN = 0.4;
  var MAX_DRIVE_STEP_MIN = 0.1;

  function relax(x, target, tau, dt) {
    if (tau <= 1e-9) return target;
    return target + (x - target) * Math.exp(-dt / tau);
  }

  function Simulation(patient, vent, opts) {
    opts = opts || {};
    this.patient = patient;
    this.vent = (vent || new Vent()).clamped();
    this._hco3Base = parseFloat((opts.hco3 === undefined || opts.hco3 === null)
      ? patient.hco3 : opts.hco3);
    this._hco3Buf = 0.0;
    this.t_min = (opts.t_min === undefined) ? 0.0 : opts.t_min;
    var base = solveSteady(patient, this.vent, { hco3: this.hco3 });
    this.paco2 = (opts.paco2 === undefined || opts.paco2 === null) ? base.paco2 : opts.paco2;
    this.pao2 = (opts.pao2 === undefined || opts.pao2 === null) ? base.pao2 : opts.pao2;
    this._fastCo2 = this.paco2;
    this._slowCo2 = this.paco2;
    this._fastO2 = this.pao2;
    this._slowO2 = this.pao2;
    this.drive = driveIndex(this.paco2, this.sao2Current(), patient);
  }
  Object.defineProperty(Simulation.prototype, 'hco3', {
    get: function () { return this._hco3Base + this._hco3Buf; },
    set: function (value) { this._hco3Base = parseFloat(value); this._hco3Buf = 0.0; }
  });
  Simulation.prototype.sao2Current = function () {
    var ph = phFrom(this.hco3, Math.max(this.paco2, 4.0));
    return so2(this.pao2, null, ph, Math.max(this.paco2, 4.0), this.patient.temp_c);
  };
  Simulation.prototype.target = function (keepUnits) {
    return solveSteady(this.patient, this.vent, {
      hco3: this.hco3, paco2_state: this.paco2, sao2_state: this.sao2Current(),
      drive: this.drive, keep_units: !!keepUnits
    });
  };
  Simulation.prototype.current = function (keepUnits) {
    return solveSteady(this.patient, this.vent, {
      hco3: this.hco3,
      paco2_state: this.paco2, sao2_state: this.sao2Current(),
      drive: this.drive, t_min: this.t_min, keep_units: !!keepUnits,
      state: [this.paco2, this.pao2, this.sao2Current()]
    });
  };
  Simulation.prototype.step = function (dtMin) {
    var dtTotal = Math.max(Math.min(dtMin, 30.0), 0.0);
    if (dtTotal <= 0) return this.current();
    /* Pressure support closes a feedback loop - drive responds to the blood
       gases, the gases to the ventilation that drive produces.  Steps longer
       than a tenth of the drive time constant make that loop overshoot, so
       sub-step whenever the loop is closed, exactly as Python does. */
    var n = 1;
    if (this.vent.mode === 'psv' && this.patient.sedation < 1.0) {
      n = Math.ceil(dtTotal / MAX_DRIVE_STEP_MIN);
    }
    var h = dtTotal / n;
    for (var i = 0; i < n; i++) this._stepOnce(h);
    return this.current();
  };
  Simulation.prototype._stepOnce = function (dt) {
    var p = this.patient;
    var tgt = this.target();
    var va = Math.max(tgt.va_l_min, 0.15);
    var co = Math.max(tgt.co_l_min, 0.4);

    var tauFastCo2 = Math.max(0.05, 0.6 * p.frc_l / va);
    var tauSlowCo2 = 30.0 / co;
    var tauFastO2 = Math.max(0.03, 0.5 * p.frc_l / va);
    var tauSlowO2 = 10.0 / co;

    this._fastCo2 = relax(this._fastCo2, tgt.paco2, tauFastCo2, dt);
    this._slowCo2 = relax(this._slowCo2, this._fastCo2, tauSlowCo2, dt);
    this.paco2 = Math.max(0.5 * (this._fastCo2 + this._slowCo2), 3.0);

    this._fastO2 = relax(this._fastO2, tgt.pao2, tauFastO2, dt);
    this._slowO2 = relax(this._slowO2, this._fastO2, tauSlowO2, dt);
    this.pao2 = Math.max(0.5 * (this._fastO2 + this._slowO2), 4.0);

    this._hco3Buf = relax(this._hco3Buf,
      ACUTE_HCO3_PER_MMHG * (this.paco2 - 40.0), TAU_ACUTE_MIN, dt);
    if (p.renal > 1e-6) {
      var gain = this.paco2 >= 40.0 ? CHRONIC_HCO3_PER_MMHG : CHRONIC_HYPO_PER_MMHG;
      var baseTarget = p.hco3 + gain * (this.paco2 - 40.0);
      this._hco3Base = relax(this._hco3Base, baseTarget,
        TAU_RENAL_MIN / Math.max(p.renal, 1e-6), dt);
    }
    var load = (p.bicarb_mmol_h - p.acid_mmol_h) / 60.0 / p.buffer_space_l;
    this._hco3Base = clamp(this._hco3Base + load * dt, 5.0, 60.0);

    var driveNow = driveIndex(this.paco2, this.sao2Current(), p);
    this.drive = relax(this.drive, driveNow, TAU_DRIVE_MIN, dt);
    this.t_min += dt;
  };
  Simulation.prototype.run = function (minutes, sampleEvery, dt) {
    if (sampleEvery === undefined) sampleEvery = 1.0;
    if (minutes <= 0) return [];
    if (dt === undefined || dt === null) dt = Math.min(sampleEvery / 4.0, 0.25);
    dt = Math.max(dt, 1e-4);
    var nSteps = Math.max(pyRound(minutes / dt), 1);
    if (nSteps > MAX_STEPS) {
      throw new Error('run too long: ' + nSteps + ' steps (max ' + MAX_STEPS +
        '); pass a coarser dt');
    }
    var nSamp = Math.max(pyRound(sampleEvery / dt), 1);
    var out = [];
    for (var k = 1; k <= nSteps; k++) {
      this.step(dt);
      if (k % nSamp === 0) out.push(this.current());
    }
    if (out.length === 0 || nSteps % nSamp) out.push(this.current());
    return out;
  };
  Simulation.prototype.adapt = function (hours) {
    var minutes = hours * 60.0;
    return this.run(minutes, Math.max(minutes / 4.0, 5.0), 2.0);
  };
  Simulation.prototype.setVent = function (kw) {
    this.vent = this.vent.with(kw).clamped();
    return this.vent;
  };
  Simulation.prototype.setPatient = function (kw) {
    this.patient = this.patient.with(kw);
    return this.patient;
  };
  Simulation.prototype.equilibrate = function () {
    var tgt = this.target();
    this.paco2 = tgt.paco2; this.pao2 = tgt.pao2;
    this._fastCo2 = this._slowCo2 = tgt.paco2;
    this._fastO2 = this._slowO2 = tgt.pao2;
    this.drive = driveIndex(this.paco2, this.sao2Current(), this.patient);
    return this.current();
  };
  Simulation.prototype.flagsNow = function () {
    return resultFlags(this.current(), this.patient);
  };

  /* ------------------------------------------------------------ exports */
  M.PHYS = PHYS;
  M.MAX_PACO2 = MAX_PACO2;
  M.MODES = MODES;
  M.makeLung = makeLung;
  M.lungBuild = lungBuild;
  M.exchange = exchange;
  M.solveUnitPo2 = solveUnitPo2;
  M.Patient = Patient;
  M.Vent = Vent;
  M.makePatient = makePatient;
  M.makeVent = makeVent;
  M.PATIENT_FIELDS = PATIENT_FIELDS;
  M.VENT_FIELDS = VENT_FIELDS;
  M.ibwKg = ibwKg;
  M.phFrom = phFrom;
  M.hco3From = hco3From;
  M.baseExcess = baseExcess;
  M.anionGap = anionGap;
  M.p50 = p50;
  M.so2 = so2;
  M.o2Content = o2Content;
  M.contentToO2 = contentToO2;
  M.oxygenDelivery = oxygenDelivery;
  M.co2Content = co2Content;
  M.pco2FromContent = pco2FromContent;
  M.pio2 = pio2;
  M.alveolarGas = alveolarGas;
  M.minuteToAlveolarPco2 = minuteToAlveolarPco2;
  M.minuteTarget = minuteTarget;
  M.driveIndex = driveIndex;
  M.solveSteady = solveSteady;
  M.driveEquilibrium = driveEquilibrium;
  M.evaluate = evaluate;
  M.resultFlags = resultFlags;
  M.meanAirwayPressure = meanAirwayPressure;
  M.Simulation = Simulation;
  M.clamp = clamp;
  M.CONST = {
    ANATOMIC_VD_PER_KG: ANATOMIC_VD_PER_KG, FRC_PER_KG: FRC_PER_KG,
    OD_WIDTH: OD_WIDTH, OD_STIFFEN: OD_STIFFEN, OD_CAPILLARY: OD_CAPILLARY,
    NONDEP_TRANSMISSION: NONDEP_TRANSMISSION, OD_CO: OD_CO,
    OD_MEAN_BASE: OD_MEAN_BASE, OD_MEAN_CO: OD_MEAN_CO,
    RECRUIT_TIDAL: RECRUIT_TIDAL, DRIVE_SLOPE: DRIVE_SLOPE,
    HYPOXIC_GAIN: HYPOXIC_GAIN, PMUS_MAX: PMUS_MAX, PMUS_CAP: PMUS_CAP,
    ACUTE_HCO3_PER_MMHG: ACUTE_HCO3_PER_MMHG,
    CHRONIC_HCO3_PER_MMHG: CHRONIC_HCO3_PER_MMHG,
    CHRONIC_HYPO_PER_MMHG: CHRONIC_HYPO_PER_MMHG,
    TAU_ACUTE_MIN: TAU_ACUTE_MIN, TAU_RENAL_MIN: TAU_RENAL_MIN,
    BUFFER_SPACE_FRACTION: BUFFER_SPACE_FRACTION,
    MAX_STEPS: MAX_STEPS, TAU_DRIVE_MIN: TAU_DRIVE_MIN,
    MAX_DRIVE_STEP_MIN: MAX_DRIVE_STEP_MIN
  };

  globalRoot.Ventsim = M;
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = M;
  }
})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));

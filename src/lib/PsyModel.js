/**
 * @fileoverview Psychoacoustic model implementation for LAME MP3 encoder.
 * Ported from psymodel.c. Computes masking thresholds based on FFT analysis.
 * Uses ES Module syntax.
 *
 * Original C Source Header:
 *      psymodel.c
 *
 *      Copyright (c) 1999-2000 Mark Taylor
 *      Copyright (c) 2001-2002 Naoki Shibata
 *      Copyright (c) 2000-2003 Takehiro Tominaga
 *      Copyright (c) 2000-2008 Robert Hegemann
 *      Copyright (c) 2000-2005 Gabriel Bouvigne
 *      Copyright (c) 2000-2005 Alexander Leidinger
 *      ... (License details omitted for brevity) ...
 *
 * $Id: PsyModel.java,v 1.27 2011/05/24 20:48:06 kenchis Exp $
 *
 * @module PsyModel
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import { FFT } from "./FFT.js";
import { Encoder } from "./Encoder.js";
import { MPEGMode } from './MPEGMode.js';

// Type Imports for JSDoc
/** @typedef {import('./LameGlobalFlags.js').default} LameGlobalFlags */
/** @typedef {import('./LameInternalFlags.js').default} LameInternalFlags */
/** @typedef {import('./III_psy_ratio.js').III_psy_ratio} III_psy_ratio */

// Destructure common utilities for easier access
const {
    // System, // Not used directly here
    VbrMode,
    Float,
    ShortBlock,
    Util,
    Arrays,
    // new_array_n, // Used indirectly via new_float_n etc.
    // new_byte, // Not used directly
    // new_double, // Not used directly
    new_float,
    new_float_n,
    new_int,
    new_int_n,
    assert
} = common;

// --- Psychoacoustic Constants ---

/** Natural logarithm of 10. @const */
const LOG10 = 2.30258509299404568402;
/** Log base 10 of e. @const */
const LN_TO_LOG10 = 0.2302585093;

/** Pre-echo control RPELEV value for long blocks. @const */
const rpelev = 2;
/** Pre-echo control RPELEV2 value for long blocks. @const */
const rpelev2 = 16;
/** Pre-echo control RPELEV value for short blocks. @const */
const rpelev_s = 2;
/** Pre-echo control RPELEV2 value for short blocks. @const */
const rpelev2_s = 16;

/** Size of each partition band, in barks. @const */
const DELBARK = .34;

/** Scale factor for loudness approximation, tuned for output level. @const */
const VO_SCALE = (1. / (14752 * 14752) / (Encoder.BLKSIZE / 2));

/** Temporal masking sustain time in seconds. @const */
const temporalmask_sustain_sec = 0.01;

/** Pre-echo attenuation factor 0 (NSPSY). @const */
const NS_PREECHO_ATT0 = 0.8;
/** Pre-echo attenuation factor 1 (NSPSY). @const */
const NS_PREECHO_ATT1 = 0.6;
/** Pre-echo attenuation factor 2 (NSPSY). @const */
const NS_PREECHO_ATT2 = 0.3;

/** Mid/Side fixing factor (NSPSY). @const */
const NS_MSFIX = 3.5;

/** Attack detection threshold for long blocks (NSPSY). @const */
const NSATTACKTHRE = 4.4;
/** Attack detection threshold for short blocks (NSPSY). @const */
const NSATTACKTHRE_S = 25;

/** FIR filter length for attack detection (NSPSY). @const */
const NSFIRLEN = 21;

/**
 * Placeholder for non-linear energy scaling. Currently identity.
 * @param {number} x - Input energy value.
 * @returns {number} Scaled energy value.
 */
function NON_LINEAR_SCALE_ENERGY(x) {
    return x;
}

// --- Masking Addition Optimization Constants ---
/** Optimization limit for mask_add (i > 8). @const */
const I1LIMIT = 8;
/** Optimization limit for mask_add (i > 24, originally 23). @const */
const I2LIMIT = 23;
/** Optimization limit for mask_add (m < 15). @const */
const MLIMIT = 15;

/** Precalculated limit for mask_add optimization (10^((I1LIMIT+1)/16)). */
let ma_max_i1;
/** Precalculated limit for mask_add optimization (10^((I2LIMIT+1)/16)). */
let ma_max_i2;
/** Precalculated limit for mask_add optimization (10^(MLIMIT/10)). */
let ma_max_m;

/**
 * Initializes precalculated maximum values for mask_add optimization.
 */
function init_mask_add_max_values() {
    ma_max_i1 = Math.pow(10, (I1LIMIT + 1) / 16.0);
    ma_max_i2 = Math.pow(10, (I2LIMIT + 1) / 16.0);
    ma_max_m = Math.pow(10, (MLIMIT) / 10.0);
}

// --- Masking Tables ---

/**
 * Masking table based on tonality. Represents pow(10, -0.0..-0.6).
 * @const {number[]}
 */
const tab = [1.0, 0.79433, 0.63096, 0.63096,
    0.63096, 0.63096, 0.63096, 0.25119, 0.11749];

/** Precomputed table 1 for mask_add function (squared values). @const {number[]} */
const table1 = [3.3246 * 3.3246, 3.23837 * 3.23837, 3.15437 * 3.15437, 3.00412 * 3.00412, 2.86103 * 2.86103, 2.65407 * 2.65407, 2.46209 * 2.46209, 2.284 * 2.284, 2.11879 * 2.11879, 1.96552 * 1.96552, 1.82335 * 1.82335, 1.69146 * 1.69146, 1.56911 * 1.56911, 1.46658 * 1.46658, 1.37074 * 1.37074, 1.31036 * 1.31036, 1.25264 * 1.25264, 1.20648 * 1.20648, 1.16203 * 1.16203, 1.12765 * 1.12765, 1.09428 * 1.09428, 1.0659 * 1.0659, 1.03826 * 1.03826, 1.01895 * 1.01895, 1];
/** Precomputed table 2 for mask_add function (squared values). @const {number[]} */
const table2 = [1.33352 * 1.33352, 1.35879 * 1.35879, 1.38454 * 1.38454, 1.39497 * 1.39497, 1.40548 * 1.40548, 1.3537 * 1.3537, 1.30382 * 1.30382, 1.22321 * 1.22321, 1.14758 * 1.14758, 1];
/** Precomputed table 3 for mask_add function (squared values). @const {number[]} */
const table3 = [2.35364 * 2.35364, 2.29259 * 2.29259, 2.23313 * 2.23313, 2.12675 * 2.12675, 2.02545 * 2.02545, 1.87894 * 1.87894, 1.74303 * 1.74303, 1.61695 * 1.61695, 1.49999 * 1.49999, 1.39148 * 1.39148, 1.29083 * 1.29083, 1.19746 * 1.19746, 1.11084 * 1.11084, 1.03826 * 1.03826];
/** Precomputed table 2 (copy) for vbrpsy_mask_add function (squared values). @const {number[]} */
const table2_ = [1.33352 * 1.33352, 1.35879 * 1.35879, 1.38454 * 1.38454, 1.39497 * 1.39497, 1.40548 * 1.40548, 1.3537 * 1.3537, 1.30382 * 1.30382, 1.22321 * 1.22321, 1.14758 * 1.14758, 1];

// --- FIR Filter Coefficients ---
/** Coefficients for the FIR high-pass filter used in attack detection (NSPSY). Stored as `coeff * 2`. @const {number[]} */
const fircoef = [-8.65163e-18 * 2, -0.00851586 * 2, -6.74764e-18 * 2, 0.0209036 * 2, -3.36639e-17 * 2, -0.0438162 * 2, -1.54175e-17 * 2, 0.0931738 * 2, -5.52212e-17 * 2, -0.313819 * 2];
const fircoef_ = fircoef; // Alias used in vbrpsy_attack_detection

// --- Perceptual Entropy Coefficients ---
/** Coefficients for perceptual entropy calculation (short blocks). Tuned for 44.1kHz. @const {number[]} */
const regcoef_s = [11.8, 13.6, 17.2, 32, 46.5, 51.3, 57.5, 67.1, 71.5, 84.6, 97.6, 130, /* 255.8 */];
/** Coefficients for perceptual entropy calculation (long blocks). Tuned for 44.1kHz. @const {number[]} */
const regcoef_l = [6.8, 5.8, 5.8, 6.4, 6.5, 9.9, 12.1, 14.4, 15, 18.9, 21.6, 26.9, 34.2, 40.2, 46.8, 56.5, 60.7, 73.9, 85.7, 93.4, 126.1, /* 241.3 */];


/**
 * @classdesc Implements the psychoacoustic model for LAME MP3 encoding.
 *
 * This class encapsulates the algorithms for computing masking thresholds based on
 * the spectral analysis of the input audio signal. It handles both the standard
 * psychoacoustic model and the NSPSYTUNE variant, as well as VBR-specific logic.
 *
 * The model performs the following main steps:
 * 1. FFT analysis (long and short blocks).
 * 2. Calculation of energy and tonality per partition band.
 * 3. Calculation of masker strength.
 * 4. Convolution with spreading functions to determine masking thresholds.
 * 5. Inter-channel and Mid/Side masking adjustments.
 * 6. Block type switching logic based on attack detection.
 * 7. Perceptual Entropy (PE) calculation.
 *
 * Note: Data returned (maskings, PE, block type) is typically delayed by one granule
 * compared to the input buffer being processed.
 *
 * @constructs PsyModel
 */
class PsyModel {
    /**
     * FFT computation object instance.
     * @type {FFT}
     * @private
     */
    fft;

    constructor() {
        /**
         * @private
         */
        this.fft = new FFT();

        // Initialize s3ind as a 2D array
        this.s3ind = new_int_n([Encoder.CBANDS, 2]);
        
        // Initialize s3 as a 2D array
        this.s3 = new_float_n([Encoder.CBANDS, Encoder.CBANDS]);
    }

    // --- Private Helper Methods ---
    // (JSDoc omitted for brevity as requested, focusing on public methods)

    /** @private */
    _psycho_loudness_approx(energy, gfc) {
        let loudness_power = 0.0;
        for (let i = 0; i < Encoder.BLKSIZE / 2; ++i) {
            loudness_power += energy[i] * gfc.ATH.eql_w[i];
        }
        loudness_power *= VO_SCALE;
        return loudness_power;
    }

    /** @private */
    _compute_ffts(gfp, fftenergy, fftenergy_s, wsamp_l, wsamp_lPos, wsamp_s, wsamp_sPos, gr_out, chn, buffer, bufPos) {
        const gfc = gfp.internal_flags;
        const current_chn_fft_l = wsamp_l[wsamp_lPos];
        const current_chn_fft_s = wsamp_s[wsamp_sPos];

        if (chn < 2) {
            this.fft.fft_long(gfc, current_chn_fft_l, chn, buffer, bufPos);
            this.fft.fft_short(gfc, current_chn_fft_s, chn, buffer, bufPos);
        }
        else if (chn === 2) {
            const left_fft_l = wsamp_l[0];
            const right_fft_l = wsamp_l[1];
            for (let j = Encoder.BLKSIZE - 1; j >= 0; --j) {
                const l = left_fft_l[j];
                const r = right_fft_l[j];
                wsamp_l[0][j] = (l + r) * Util.SQRT2 * 0.5; // Mid
                wsamp_l[1][j] = (l - r) * Util.SQRT2 * 0.5; // Side
            }
            const left_fft_s = wsamp_s[0];
            const right_fft_s = wsamp_s[1];
            for (let b = 2; b >= 0; --b) {
                for (let j = Encoder.BLKSIZE_s - 1; j >= 0; --j) {
                    const l = left_fft_s[b][j];
                    const r = right_fft_s[b][j];
                    wsamp_s[0][b][j] = (l + r) * Util.SQRT2 * 0.5; // Mid
                    wsamp_s[1][b][j] = (l - r) * Util.SQRT2 * 0.5; // Side
                }
            }
        }

        const target_fft_l = wsamp_l[wsamp_lPos];
        fftenergy[0] = NON_LINEAR_SCALE_ENERGY(target_fft_l[0]);
        fftenergy[0] *= fftenergy[0];
        for (let j = Encoder.BLKSIZE / 2 - 1; j >= 0; --j) {
            const re = target_fft_l[Encoder.BLKSIZE / 2 - j];
            const im = target_fft_l[Encoder.BLKSIZE / 2 + j];
            fftenergy[Encoder.BLKSIZE / 2 - j] = NON_LINEAR_SCALE_ENERGY((re * re + im * im) * 0.5);
        }

        const target_fft_s = wsamp_s[wsamp_sPos];
        for (let b = 2; b >= 0; --b) {
            fftenergy_s[b][0] = target_fft_s[b][0];
            fftenergy_s[b][0] *= fftenergy_s[b][0];
            for (let j = Encoder.BLKSIZE_s / 2 - 1; j >= 0; --j) {
                const re = target_fft_s[b][Encoder.BLKSIZE_s / 2 - j];
                const im = target_fft_s[b][Encoder.BLKSIZE_s / 2 + j];
                fftenergy_s[b][Encoder.BLKSIZE_s / 2 - j] = NON_LINEAR_SCALE_ENERGY((re * re + im * im) * 0.5);
            }
        }

        {
            let totalenergy = 0.0;
            for (let j = 11; j < Encoder.HBLKSIZE; j++)
                totalenergy += fftenergy[j];
            gfc.tot_ener[chn] = totalenergy;
        }

        if (gfp.analysis) {
            for (let j = 0; j < Encoder.HBLKSIZE; j++) {
                gfc.pinfo.energy[gr_out][chn][j] = gfc.pinfo.energy_save[chn][j];
                gfc.pinfo.energy_save[chn][j] = fftenergy[j];
            }
             gfc.pinfo.pe[gr_out][chn] = gfc.pe[chn];
        }

        if (gfp.athaa_loudapprox === 2 && chn < 2) {
            gfc.loudness_sq[gr_out][chn] = gfc.loudness_sq_save[chn];
            gfc.loudness_sq_save[chn] = this._psycho_loudness_approx(fftenergy, gfc);
        }
    }

    /** @private */
    _mask_add(m1, m2, kk, b, gfc, shortblock) {
        let ratio;
        if (m2 > m1) {
            if (m2 < (m1 * ma_max_i2)) ratio = m2 / m1;
            else return (m1 + m2);
        } else {
            if (m1 >= (m2 * ma_max_i2)) return (m1 + m2);
            ratio = m1 / m2;
        }
        assert(m1 >= 0); assert(m2 >= 0);
        let combined_masking = m1 + m2;
        if ((b + 3) <= 6) {
            if (ratio >= ma_max_i1) {
                return combined_masking;
            }
            const i = 0 | (Util.FAST_LOG10_X(ratio, 16.0));
            return combined_masking * table2[i];
        }
        const i = 0 | Util.FAST_LOG10_X(ratio, 16.0);
        let ath_cb;
        if (shortblock !== 0) {
            ath_cb = gfc.ATH.cb_s[kk] * gfc.ATH.adjust;
        } else {
            ath_cb = gfc.ATH.cb_l[kk] * gfc.ATH.adjust;
        }
        assert(ath_cb >= 0);
        if (combined_masking < ma_max_m * ath_cb) {
            if (combined_masking > ath_cb) {
                let f = 1.0;
                if (i <= 13) f = table3[i];
                const r_m_ath = Util.FAST_LOG10_X(combined_masking / ath_cb, 10.0 / 15.0);
                return combined_masking * ((table1[i] - f) * r_m_ath + f);
            }
            if (i > 13) return combined_masking;
            return combined_masking * table3[i];
        }
        return combined_masking * table1[i];
    }

    /** @private */
    _vbrpsy_mask_add(m1, m2, b) {
        let ratio;
        if (m1 < 0) m1 = 0;
        if (m2 < 0) m2 = 0;
        if (m1 <= 0) return m2;
        if (m2 <= 0) return m1;
        if (m2 > m1) {
            ratio = m2 / m1;
        } else {
            ratio = m1 / m2;
        }
        if (-2 <= b && b <= 2) {
            if (ratio >= ma_max_i1) {
                return m1 + m2;
            } else {
                const i = 0 | (Util.FAST_LOG10_X(ratio, 16.0));
                return (m1 + m2) * table2_[i];
            }
        }
        if (ratio < ma_max_i2) {
            return m1 + m2;
        }
        if (m1 < m2) return m2;
        return m1;
    }

    /** @private */
    _calc_interchannel_masking(gfp, ratio) {
        const gfc = gfp.internal_flags;
        if (gfc.channels_out > 1) {
            for (let sb = 0; sb < Encoder.SBMAX_l; sb++) {
                const l = gfc.thm[0].l[sb];
                const r = gfc.thm[1].l[sb];
                gfc.thm[0].l[sb] += r * ratio;
                gfc.thm[1].l[sb] += l * ratio;
            }
            for (let sb = 0; sb < Encoder.SBMAX_s; sb++) {
                for (let sblock = 0; sblock < 3; sblock++) {
                    const l = gfc.thm[0].s[sb][sblock];
                    const r = gfc.thm[1].s[sb][sblock];
                    gfc.thm[0].s[sb][sblock] += r * ratio;
                    gfc.thm[1].s[sb][sblock] += l * ratio;
                }
            }
        }
    }

    /** @private */
    _msfix1(gfc) {
        for (let sb = 0; sb < Encoder.SBMAX_l; sb++) {
            const thmL = gfc.thm[0].l[sb]; const thmR = gfc.thm[1].l[sb];
            if (thmL > 1.58 * thmR || thmR > 1.58 * thmL) continue;
            const thmM = gfc.thm[2].l[sb]; const thmS = gfc.thm[3].l[sb];
            const enM = gfc.en[2].l[sb]; const enS = gfc.en[3].l[sb];
            const mld_sb = gfc.mld_l[sb];
            const mld_mid = mld_sb * enS; const rmid = Math.max(thmM, Math.min(thmS, mld_mid));
            const mld_side = mld_sb * enM; const rside = Math.max(thmS, Math.min(thmM, mld_side));
            gfc.thm[2].l[sb] = rmid; gfc.thm[3].l[sb] = rside;
        }
        for (let sb = 0; sb < Encoder.SBMAX_s; sb++) {
            for (let sblock = 0; sblock < 3; sblock++) {
                const thmL = gfc.thm[0].s[sb][sblock]; const thmR = gfc.thm[1].s[sb][sblock];
                if (thmL > 1.58 * thmR || thmR > 1.58 * thmL) continue;
                const thmM = gfc.thm[2].s[sb][sblock]; const thmS = gfc.thm[3].s[sb][sblock];
                const enM = gfc.en[2].s[sb][sblock]; const enS = gfc.en[3].s[sb][sblock];
                const mld_sb = gfc.mld_s[sb];
                const mld_mid = mld_sb * enS; const rmid = Math.max(thmM, Math.min(thmS, mld_mid));
                const mld_side = mld_sb * enM; const rside = Math.max(thmS, Math.min(thmM, mld_side));
                gfc.thm[2].s[sb][sblock] = rmid; gfc.thm[3].s[sb][sblock] = rside;
            }
        }
    }

    /** @private */
    _ns_msfix(gfc, msfix, athadjust) {
        const msfix2 = msfix * 2.0;
        let athlower = Math.pow(10, athadjust); // Expects linear adjustment factor
        for (let sb = 0; sb < Encoder.SBMAX_l; sb++) {
            const ath = gfc.ATH.cb_l[gfc.bm_l[sb]] * athlower;
            const thmL = Math.max(gfc.thm[0].l[sb], ath);
            const thmR = Math.max(gfc.thm[1].l[sb], ath);
            let thmM = Math.max(gfc.thm[2].l[sb], ath);
            let thmS = Math.max(gfc.thm[3].l[sb], ath);
            const thmLR_min = Math.min(thmL, thmR);
            const thmMS_sum = thmM + thmS;
            if (thmMS_sum > 0 && thmLR_min * msfix2 < thmMS_sum) {
                const f = thmLR_min * msfix2 / thmMS_sum;
                thmM *= f; thmS *= f;
                assert(thmM + thmS >= 0); // Allow zero
            }
            gfc.thm[2].l[sb] = Math.min(thmM, gfc.thm[2].l[sb]);
            gfc.thm[3].l[sb] = Math.min(thmS, gfc.thm[3].l[sb]);
        }
        athlower *= (Encoder.BLKSIZE_s / Encoder.BLKSIZE);
        for (let sb = 0; sb < Encoder.SBMAX_s; sb++) {
            for (let sblock = 0; sblock < 3; sblock++) {
                const ath = gfc.ATH.cb_s[gfc.bm_s[sb]] * athlower;
                const thmL = Math.max(gfc.thm[0].s[sb][sblock], ath);
                const thmR = Math.max(gfc.thm[1].s[sb][sblock], ath);
                let thmM = Math.max(gfc.thm[2].s[sb][sblock], ath);
                let thmS = Math.max(gfc.thm[3].s[sb][sblock], ath);
                const thmLR_min = Math.min(thmL, thmR);
                const thmMS_sum = thmM + thmS;
                // Check C code again: NS uses msfix, VBR uses msfix2 here... Let's follow C for NSPSY
                 if (thmMS_sum > 0 && thmLR_min * msfix < thmMS_sum) { // Using msfix here
                    const f = thmLR_min * msfix / thmMS_sum;
                    thmM *= f; thmS *= f;
                    assert(thmM + thmS >= 0); // Allow zero
                }
                gfc.thm[2].s[sb][sblock] = Math.min(gfc.thm[2].s[sb][sblock], thmM);
                gfc.thm[3].s[sb][sblock] = Math.min(gfc.thm[3].s[sb][sblock], thmS);
            }
        }
    }

    /** @private */
    _convert_partition2scalefac_s(gfc, eb, thr, chn, sblock) {
        let sb, b; let enn = 0.0; let thmm = 0.0; const npart_s = gfc.npart_s;
        for (sb = b = 0; sb < Encoder.SBMAX_s; ++sb) {
            const bo_s_sb = gfc.bo_s[sb]; const b_lim = bo_s_sb < npart_s ? bo_s_sb : npart_s;
            while (b < b_lim) { assert(eb[b] >= 0); assert(thr[b] >= 0); enn += eb[b]; thmm += thr[b]; b++; }
            gfc.en[chn].s[sb][sblock] = enn; gfc.thm[chn].s[sb][sblock] = thmm;
            if (b >= npart_s) { ++sb; break; }
            assert(eb[b] >= 0); assert(thr[b] >= 0);
            {
                const w_curr = gfc.PSY.bo_s_weight[sb]; const w_next = 1.0 - w_curr;
                const enn_curr = w_curr * eb[b]; const thmm_curr = w_curr * thr[b];
                gfc.en[chn].s[sb][sblock] += enn_curr; gfc.thm[chn].s[sb][sblock] += thmm_curr;
                enn = w_next * eb[b]; thmm = w_next * thr[b];
            }
            b++;
        }
        for (; sb < Encoder.SBMAX_s; ++sb) { gfc.en[chn].s[sb][sblock] = 0; gfc.thm[chn].s[sb][sblock] = 0; }
    }

    /** @private */
    _convert_partition2scalefac_l(gfc, eb, thr, chn) {
        let sb, b; let enn = 0.0; let thmm = 0.0; const npart_l = gfc.npart_l;
        for (sb = b = 0; sb < Encoder.SBMAX_l; ++sb) {
            const bo_l_sb = gfc.bo_l[sb]; const b_lim = bo_l_sb < npart_l ? bo_l_sb : npart_l;
            while (b < b_lim) { assert(eb[b] >= 0); assert(thr[b] >= 0); enn += eb[b]; thmm += thr[b]; b++; }
            gfc.en[chn].l[sb] = enn; gfc.thm[chn].l[sb] = thmm;
            if (b >= npart_l) { ++sb; break; }
            assert(eb[b] >= 0); assert(thr[b] >= 0);
            {
                const w_curr = gfc.PSY.bo_l_weight[sb]; const w_next = 1.0 - w_curr;
                const enn_curr = w_curr * eb[b]; const thmm_curr = w_curr * thr[b];
                gfc.en[chn].l[sb] += enn_curr; gfc.thm[chn].l[sb] += thmm_curr;
                enn = w_next * eb[b]; thmm = w_next * thr[b];
            }
             b++;
        }
        for (; sb < Encoder.SBMAX_l; ++sb) { gfc.en[chn].l[sb] = 0; gfc.thm[chn].l[sb] = 0; }
    }

    /** @private */
    _compute_masking_s(gfp, fftenergy_s, eb, thr, chn, sblock) {
        const gfc = gfp.internal_flags; let j, b;
        for (b = j = 0; b < gfc.npart_s; ++b) {
            let ebb = 0; const n = gfc.numlines_s[b];
            for (let i = 0; i < n; ++i, ++j) ebb += fftenergy_s[sblock][j];
            eb[b] = ebb;
        }
        assert(b === gfc.npart_s); assert(j === Encoder.HBLKSIZE_s);
        let s3_idx = 0;
        for (b = 0; b < gfc.npart_s; b++) {
            const first_masker_idx = gfc.s3ind_s[b][0]; const last_masker_idx = gfc.s3ind_s[b][1];
            let kk = first_masker_idx; let ecb = 0;
            while (kk <= last_masker_idx) { ecb += gfc.s3_ss[s3_idx] * eb[kk]; s3_idx++; kk++; }
            {
                const limit1 = rpelev_s * gfc.nb_s1[chn][b];
                thr[b] = Math.min(ecb, limit1);
                if (gfc.blocktype_old[chn & 1] === Encoder.SHORT_TYPE) {
                    const limit2 = rpelev2_s * gfc.nb_s2[chn][b];
                    thr[b] = Math.min(thr[b], limit2);
                }
            }
            gfc.nb_s2[chn][b] = gfc.nb_s1[chn][b]; gfc.nb_s1[chn][b] = ecb;
            assert(thr[b] >= 0);
        }
        assert(b === gfc.npart_s);
        for (; b <= Encoder.CBANDS; ++b) { eb[b] = 0; thr[b] = 0; }
        if (b <= Encoder.CBANDS + 1) thr[b] = 0;
    }

    /** @private */
    _block_type_set(gfp, uselongblock, blocktype_d, blocktype) {
        const gfc = gfp.internal_flags;
        if (gfp.short_blocks === ShortBlock.short_block_coupled && !(uselongblock[0] !== 0 && uselongblock[1] !== 0)) {
            uselongblock[0] = uselongblock[1] = 0;
        }
        for (let chn = 0; chn < gfc.channels_out; chn++) {
            if (gfp.short_blocks === ShortBlock.short_block_dispensed) uselongblock[chn] = 1;
            if (gfp.short_blocks === ShortBlock.short_block_forced) uselongblock[chn] = 0;
            blocktype[chn] = Encoder.NORM_TYPE;
            let final_prev_type = gfc.blocktype_old[chn]; // Start assuming previous state persists

            if (uselongblock[chn] !== 0) { // Current is LONG
                 assert(gfc.blocktype_old[chn] !== Encoder.START_TYPE);
                 if (gfc.blocktype_old[chn] === Encoder.SHORT_TYPE) {
                      blocktype[chn] = Encoder.STOP_TYPE; // Current logic determines STOP
                      final_prev_type = Encoder.STOP_TYPE; // Previous becomes STOP
                 }
                  // else: Previous was NORM or STOP, current is NORM. final_prev_type remains NORM or STOP. blocktype[chn] is NORM.
            } else { // Current is SHORT
                 blocktype[chn] = Encoder.SHORT_TYPE;
                 if (gfc.blocktype_old[chn] === Encoder.NORM_TYPE) {
                     // Previous was NORM, becomes START
                     final_prev_type = Encoder.START_TYPE;
                 } else if (gfc.blocktype_old[chn] === Encoder.STOP_TYPE) {
                     // Previous was STOP, becomes SHORT
                     final_prev_type = Encoder.SHORT_TYPE;
                 }
                 // else: Previous was START or SHORT, remains START or SHORT. final_prev_type = gfc.blocktype_old[chn]. blocktype[chn] is SHORT.
            }
            blocktype_d[chn] = final_prev_type;
            gfc.blocktype_old[chn] = blocktype[chn]; // Save current decision as history
        }
    }

    /** @private */
    _NS_INTERP(x, y, r) {
        if (r >= 1.0) return x;
        if (r <= 0.0) return y;
        if (y > 0.0) return (Math.pow(x / y, r) * y);
        return 0.0;
    }

    /** @private */
    _pecalc_s(mr, masking_lower) {
        let pe_s = 309.07;
        for (let sb = 0; sb < Encoder.SBMAX_s - 1; sb++) {
            assert(sb < regcoef_s.length);
            for (let sblock = 0; sblock < 3; sblock++) {
                const thm = mr.thm.s[sb][sblock];
                if (thm > 0.0) {
                    const x = thm * masking_lower; const en = mr.en.s[sb][sblock];
                    if (en > x) {
                        if (en > x * 1e10) { pe_s += regcoef_s[sb] * (10.0 * LOG10); }
                        else { assert(x > 0); pe_s += regcoef_s[sb] * Util.FAST_LOG10(en / x); }
                    }
                }
            }
        }
        return pe_s;
    }

    /** @private */
    _pecalc_l(mr, masking_lower) {
        let pe_l = 281.0575;
        for (let sb = 0; sb < Encoder.SBMAX_l - 1; sb++) {
            assert(sb < regcoef_l.length);
            const thm = mr.thm.l[sb];
            if (thm > 0.0) {
                const x = thm * masking_lower; const en = mr.en.l[sb];
                if (en > x) {
                    if (en > x * 1e10) { pe_l += regcoef_l[sb] * (10.0 * LOG10); }
                    else { assert(x > 0); pe_l += regcoef_l[sb] * Util.FAST_LOG10(en / x); }
                }
            }
        }
        return pe_l;
    }

    /** @private */
    _calc_energy(gfc, fftenergy, eb, max, avg) {
        let b, j;
        for (b = j = 0; b < gfc.npart_l; ++b) {
            let ebb = 0, m = 0; const numlines = gfc.numlines_l[b];
            for (let i = 0; i < numlines; ++i, ++j) { const el = fftenergy[j]; assert(el >= 0); ebb += el; if (m < el) m = el; }
            eb[b] = ebb; max[b] = m; avg[b] = ebb * gfc.rnumlines_l[b];
            assert(gfc.rnumlines_l[b] >= 0); assert(ebb >= 0); assert(eb[b] >= 0); assert(max[b] >= 0); assert(avg[b] >= 0);
        }
         assert(b === gfc.npart_l); assert(j === Encoder.HBLKSIZE);
    }

    /** @private */
    _calc_mask_index_l(gfc, max, avg, mask_idx) {
        const last_tab_entry = tab.length - 1; let b = 0;
        let a = avg[b] + avg[b + 1]; assert(a >= 0);
        if (a > 0.0) {
            let m = max[b]; if (m < max[b + 1]) m = max[b + 1];
            const nl = gfc.numlines_l[b] + gfc.numlines_l[b + 1]; assert(nl - 1 > 0);
            a = 20.0 * (m * 2.0 - a) / (a * (nl - 1));
            let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
        } else mask_idx[b] = 0;
        for (b = 1; b < gfc.npart_l - 1; b++) {
            a = avg[b - 1] + avg[b] + avg[b + 1]; assert(a >= 0);
            if (a > 0.0) {
                let m = max[b - 1]; if (m < max[b]) m = max[b]; if (m < max[b + 1]) m = max[b + 1];
                const nl = gfc.numlines_l[b - 1] + gfc.numlines_l[b] + gfc.numlines_l[b + 1]; assert(nl - 1 > 0);
                a = 20.0 * (m * 3.0 - a) / (a * (nl - 1));
                let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
            } else mask_idx[b] = 0;
        }
        assert(b > 0); assert(b === gfc.npart_l - 1);
        a = avg[b - 1] + avg[b]; assert(a >= 0);
        if (a > 0.0) {
            let m = max[b - 1]; if (m < max[b]) m = max[b];
            const nl = gfc.numlines_l[b - 1] + gfc.numlines_l[b]; assert(nl - 1 > 0);
            a = 20.0 * (m * 2.0 - a) / (a * (nl - 1));
            let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
        } else mask_idx[b] = 0;
        assert(b === (gfc.npart_l - 1));
    }

    /** @private */
    _psyvbr_calc_mask_index_s(gfc, max, avg, mask_idx) {
        const last_tab_entry = tab.length - 1; let b = 0;
        let a = avg[b] + avg[b + 1]; assert(a >= 0);
        if (a > 0.0) {
            let m = max[b]; if (m < max[b + 1]) m = max[b + 1];
            const nl = gfc.numlines_s[b] + gfc.numlines_s[b + 1]; assert(nl - 1 > 0);
            a = 20.0 * (m * 2.0 - a) / (a * (nl - 1));
            let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
        } else mask_idx[b] = 0;
        for (b = 1; b < gfc.npart_s - 1; b++) {
            a = avg[b - 1] + avg[b] + avg[b + 1]; assert(b + 1 < gfc.npart_s); assert(a >= 0);
            if (a > 0.0) {
                let m = max[b - 1]; if (m < max[b]) m = max[b]; if (m < max[b + 1]) m = max[b + 1];
                const nl = gfc.numlines_s[b - 1] + gfc.numlines_s[b] + gfc.numlines_s[b + 1]; assert(nl - 1 > 0);
                a = 20.0 * (m * 3.0 - a) / (a * (nl - 1));
                let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
            } else mask_idx[b] = 0;
        }
        assert(b > 0); assert(b === gfc.npart_s - 1);
        a = avg[b - 1] + avg[b]; assert(a >= 0);
        if (a > 0.0) {
            let m = max[b - 1]; if (m < max[b]) m = max[b];
            const nl = gfc.numlines_s[b - 1] + gfc.numlines_s[b]; assert(nl - 1 > 0);
            a = 20.0 * (m * 2.0 - a) / (a * (nl - 1));
            let k = 0 | a; if(k<0)k=0; if (k > last_tab_entry) k = last_tab_entry; mask_idx[b] = k;
        } else mask_idx[b] = 0;
        assert(b === (gfc.npart_s - 1));
    }

    /** @private */
    _vbrpsy_compute_masking_s(gfp, fftenergy_s, eb, thr, chn, sblock) {
        const gfc = gfp.internal_flags;
        const max = new_float(Encoder.CBANDS); const avg = new_float(Encoder.CBANDS);
        const mask_idx_s = new_int(Encoder.CBANDS);
        let i, j, b;
        for (b = j = 0; b < gfc.npart_s; ++b) {
            let ebb = 0, m = 0; const n = gfc.numlines_s[b];
            for (i = 0; i < n; ++i, ++j) { const el = fftenergy_s[sblock][j]; ebb += el; if (m < el) m = el; }
            eb[b] = ebb; assert(ebb >= 0); max[b] = m; assert(n > 0); avg[b] = ebb / n; assert(avg[b] >= 0);
        }
        assert(b === gfc.npart_s); assert(j === Encoder.HBLKSIZE_s);
        for (; b < Encoder.CBANDS; ++b) { max[b] = 0; avg[b] = 0; }
        this._psyvbr_calc_mask_index_s(gfc, max, avg, mask_idx_s);
        let s3_idx = 0;
        for (b = 0; b < gfc.npart_s; b++) {
            const first_masker_idx = gfc.s3ind_s[b][0]; const last_masker_idx = gfc.s3ind_s[b][1];
            let kk = first_masker_idx; let ecb = 0; let dd = 0; let dd_n = 0;
            if (kk <= last_masker_idx) {
                 dd = mask_idx_s[kk]; dd_n = 1;
                 ecb = gfc.s3_ss[s3_idx] * eb[kk] * tab[mask_idx_s[kk]]; s3_idx++; kk++;
                 while (kk <= last_masker_idx) {
                     dd += mask_idx_s[kk]; dd_n += 1;
                     const masker_contribution = gfc.s3_ss[s3_idx] * eb[kk] * tab[mask_idx_s[kk]];
                     ecb = this._vbrpsy_mask_add(ecb, masker_contribution, kk - b); s3_idx++; kk++;
                 }
                 dd = (1 + 2 * dd) / (2 * dd_n); const avg_mask = tab[dd] * 0.5; ecb *= avg_mask;
            } else { ecb = 0; }
            thr[b] = ecb;
            gfc.nb_s2[chn][b] = gfc.nb_s1[chn][b]; gfc.nb_s1[chn][b] = ecb;
            {
                 let max_energy_limit = max[b] * gfc.minval_s[b] * tab[(1 + 2 * dd) / (2 * dd_n)] * 0.5;
                 if (thr[b] > max_energy_limit) thr[b] = max_energy_limit;
            }
            if (gfc.masking_lower > 1) thr[b] *= gfc.masking_lower;
            if (thr[b] > eb[b]) thr[b] = eb[b];
            if (gfc.masking_lower < 1) thr[b] *= gfc.masking_lower;
            assert(thr[b] >= 0);
        }
         assert(b === gfc.npart_s);
         for (; b < Encoder.CBANDS; ++b) { eb[b] = 0; thr[b] = 0; }
         for (; b <= Encoder.CBANDS + 1; ++b) { thr[b] = 0; }
    }

    /** @private */
    _vbrpsy_compute_masking_l(gfc, fftenergy, eb_l, thr, chn) {
        const max = new_float(Encoder.CBANDS); const avg = new_float(Encoder.CBANDS);
        const mask_idx_l = new_int(Encoder.CBANDS + 2); let b;
        this._calc_energy(gfc, fftenergy, eb_l, max, avg);
        this._calc_mask_index_l(gfc, max, avg, mask_idx_l);
        let s3_idx = 0;
        for (b = 0; b < gfc.npart_l; b++) {
            let ecb = 0; let dd = 0; let dd_n = 0;
            const first_masker_idx = gfc.s3ind[b][0]; const last_masker_idx = gfc.s3ind[b][1];
            let kk = first_masker_idx;
            if (kk <= last_masker_idx) {
                 dd = mask_idx_l[kk]; dd_n = 1;
                 ecb = gfc.s3_ll[s3_idx] * eb_l[kk] * tab[mask_idx_l[kk]]; s3_idx++; kk++;
                 while (kk <= last_masker_idx) {
                     dd += mask_idx_l[kk]; dd_n += 1;
                     const masker_contribution = gfc.s3_ll[s3_idx] * eb_l[kk] * tab[mask_idx_l[kk]];
                     ecb = this._vbrpsy_mask_add(ecb, masker_contribution, kk - b); s3_idx++; kk++;
                 }
                 dd = (1 + 2 * dd) / (2 * dd_n); const avg_mask = tab[dd] * 0.5; ecb *= avg_mask;
             } else { ecb = 0; }

            const prev_blocktype = gfc.blocktype_old[chn & 0x01];
            if (prev_blocktype === Encoder.SHORT_TYPE || prev_blocktype === Encoder.START_TYPE) {
                const ecb_limit_1 = rpelev * gfc.nb_1[chn][b];
                if (ecb_limit_1 > 0) thr[b] = Math.min(ecb, ecb_limit_1);
                else thr[b] = Math.min(ecb, eb_l[b] * NS_PREECHO_ATT2);
            } else {
                let ecb_limit_2 = rpelev2 * gfc.nb_2[chn][b]; let ecb_limit_1 = rpelev * gfc.nb_1[chn][b];
                if (ecb_limit_2 <= 0) ecb_limit_2 = ecb; if (ecb_limit_1 <= 0) ecb_limit_1 = ecb;
                let ecb_limit = (prev_blocktype === Encoder.NORM_TYPE) ? Math.min(ecb_limit_1, ecb_limit_2) : ecb_limit_1;
                thr[b] = Math.min(ecb, ecb_limit);
            }
            gfc.nb_2[chn][b] = gfc.nb_1[chn][b]; gfc.nb_1[chn][b] = ecb;
            {
                 let max_energy_limit = max[b] * gfc.minval_l[b] * tab[(1 + 2 * dd) / (2 * dd_n)] * 0.5;
                 if (thr[b] > max_energy_limit) thr[b] = max_energy_limit;
            }
            if (gfc.masking_lower > 1) thr[b] *= gfc.masking_lower;
            if (thr[b] > eb_l[b]) thr[b] = eb_l[b];
            if (gfc.masking_lower < 1) thr[b] *= gfc.masking_lower;
            assert(thr[b] >= 0);
        }
        assert(b === gfc.npart_l);
        for (; b < Encoder.CBANDS; ++b) { eb_l[b] = 0; thr[b] = 0; }
        for (; b <= Encoder.CBANDS + 1; ++b) { thr[b] = 0; }
    }

    /** @private */
    _vbrpsy_compute_block_type(gfp, uselongblock) {
        const gfc = gfp.internal_flags;
        if (gfp.short_blocks === ShortBlock.short_block_coupled && !(uselongblock[0] !== 0 && uselongblock[1] !== 0)) {
            uselongblock[0] = uselongblock[1] = 0;
        }
        for (let chn = 0; chn < gfc.channels_out; chn++) {
            if (gfp.short_blocks === ShortBlock.short_block_dispensed) uselongblock[chn] = 1;
            if (gfp.short_blocks === ShortBlock.short_block_forced) uselongblock[chn] = 0;
        }
    }

    /** @private */
    _vbrpsy_apply_block_type(gfp, uselongblock, blocktype_d) {
        const gfc = gfp.internal_flags;
        for (let chn = 0; chn < gfc.channels_out; chn++) {
            let current_blocktype = Encoder.NORM_TYPE;
            let final_prev_type = gfc.blocktype_old[chn]; // Start assuming previous state persists
            if (uselongblock[chn] !== 0) { // Current wants LONG
                 assert(gfc.blocktype_old[chn] !== Encoder.START_TYPE);
                 if (gfc.blocktype_old[chn] === Encoder.SHORT_TYPE) {
                      current_blocktype = Encoder.STOP_TYPE;
                      final_prev_type = Encoder.STOP_TYPE;
                 }
            } else { // Current wants SHORT
                 current_blocktype = Encoder.SHORT_TYPE;
                 if (gfc.blocktype_old[chn] === Encoder.NORM_TYPE) final_prev_type = Encoder.START_TYPE;
                 else if (gfc.blocktype_old[chn] === Encoder.STOP_TYPE) final_prev_type = Encoder.SHORT_TYPE;
            }
            blocktype_d[chn] = final_prev_type;
            gfc.blocktype_old[chn] = current_blocktype;
        }
    }

    /** @private */
    _vbrpsy_compute_MS_thresholds(eb, thr, cb_mld, ath_cb, athadjust, msfix, n) {
        const msfix2 = msfix * 2.0;
        const athlower = msfix > 0 ? athadjust : 1.0; // athadjust is already linear
        let rside, rmid;
        for (let b = 0; b < n; ++b) {
            const ebM = eb[2][b]; const ebS = eb[3][b];
            const thmL = thr[0][b]; const thmR = thr[1][b];
            const thmM_orig = thr[2][b]; const thmS_orig = thr[3][b];
            if (thmL <= 1.58 * thmR && thmR <= 1.58 * thmL) {
                const mld_factor = cb_mld[b];
                const mld_mid = mld_factor * ebS; const mld_side = mld_factor * ebM;
                rmid = Math.max(thmM_orig, Math.min(thmS_orig, mld_mid));
                rside = Math.max(thmS_orig, Math.min(thmM_orig, mld_side));
            } else { rmid = thmM_orig; rside = thmS_orig; }
            if (msfix > 0) {
                const ath = ath_cb[b] * athlower;
                const thmL_ath = Math.max(thmL, ath); const thmR_ath = Math.max(thmR, ath);
                let thmM_ath = Math.max(rmid, ath); let thmS_ath = Math.max(rside, ath);
                const thmLR_min_ath = Math.min(thmL_ath, thmR_ath);
                const thmMS_sum_ath = thmM_ath + thmS_ath;
                if (thmMS_sum_ath > 0 && (thmLR_min_ath * msfix2) < thmMS_sum_ath) { // VBR uses msfix2 here
                    const f = thmLR_min_ath * msfix2 / thmMS_sum_ath;
                    thmM_ath *= f; thmS_ath *= f;
                     assert(thmM_ath + thmS_ath >= 0);
                }
                rmid = Math.min(thmM_ath, rmid); rside = Math.min(thmS_ath, rside);
            }
            if (rmid > ebM) rmid = ebM; if (rside > ebS) rside = ebS;
            thr[2][b] = rmid; thr[3][b] = rside;
        }
    }

    /** @private */
    _vbrpsy_compute_fft_l(gfp, buffer, bufPos, chn, gr_out, fftenergy, wsamp_l, wsamp_lPos) {
        const gfc = gfp.internal_flags;
        const target_fft_l = wsamp_l[wsamp_lPos];
        if (chn < 2) {
            this.fft.fft_long(gfc, target_fft_l, chn, buffer, bufPos);
        } else if (chn === 2) {
            const left_fft_l = wsamp_l[0]; const right_fft_l = wsamp_l[1];
            for (let j = Encoder.BLKSIZE - 1; j >= 0; --j) {
                const l = left_fft_l[j]; const r = right_fft_l[j];
                wsamp_l[0][j] = (l + r) * Util.SQRT2 * 0.5; wsamp_l[1][j] = (l - r) * Util.SQRT2 * 0.5;
            }
        }
        fftenergy[0] = NON_LINEAR_SCALE_ENERGY(target_fft_l[0]); fftenergy[0] *= fftenergy[0];
        for (let j = Encoder.BLKSIZE / 2 - 1; j >= 0; --j) {
            const re = target_fft_l[Encoder.BLKSIZE / 2 - j]; const im = target_fft_l[Encoder.BLKSIZE / 2 + j];
            fftenergy[Encoder.BLKSIZE / 2 - j] = NON_LINEAR_SCALE_ENERGY((re * re + im * im) * 0.5);
        }
        { let totalenergy = 0.0; for (let j = 11; j < Encoder.HBLKSIZE; j++) totalenergy += fftenergy[j]; gfc.tot_ener[chn] = totalenergy; }
        if (gfp.analysis) {
            for (let j = 0; j < Encoder.HBLKSIZE; j++) { gfc.pinfo.energy[gr_out][chn][j] = gfc.pinfo.energy_save[chn][j]; gfc.pinfo.energy_save[chn][j] = fftenergy[j]; }
             gfc.pinfo.pe[gr_out][chn] = gfc.pe[chn];
        }
    }

    /** @private */
    _vbrpsy_compute_fft_s(gfp, buffer, bufPos, chn, sblock, fftenergy_s, wsamp_s, wsamp_sPos) {
        const gfc = gfp.internal_flags;
        const target_fft_s = wsamp_s[wsamp_sPos];
        if (sblock === 0 && chn < 2) {
            this.fft.fft_short(gfc, target_fft_s, chn, buffer, bufPos);
        }
        if (chn === 2) {
            const left_fft_s = wsamp_s[0]; const right_fft_s = wsamp_s[1];
            for (let j = Encoder.BLKSIZE_s - 1; j >= 0; --j) {
                const l = left_fft_s[sblock][j]; const r = right_fft_s[sblock][j];
                wsamp_s[0][sblock][j] = (l + r) * Util.SQRT2 * 0.5; wsamp_s[1][sblock][j] = (l - r) * Util.SQRT2 * 0.5;
            }
        }
        const current_sblock_fft = target_fft_s[sblock];
        fftenergy_s[sblock][0] = current_sblock_fft[0]; fftenergy_s[sblock][0] *= fftenergy_s[sblock][0];
        for (let j = Encoder.BLKSIZE_s / 2 - 1; j >= 0; --j) {
            const re = current_sblock_fft[Encoder.BLKSIZE_s / 2 - j]; const im = current_sblock_fft[Encoder.BLKSIZE_s / 2 + j];
            fftenergy_s[sblock][Encoder.BLKSIZE_s / 2 - j] = NON_LINEAR_SCALE_ENERGY((re * re + im * im) * 0.5);
        }
    }

    /** @private */
    _vbrpsy_compute_loudness_approximation_l(gfp, gr_out, chn, fftenergy) {
        const gfc = gfp.internal_flags;
        if (gfp.athaa_loudapprox === 2 && chn < 2) {
            gfc.loudness_sq[gr_out][chn] = gfc.loudness_sq_save[chn];
            gfc.loudness_sq_save[chn] = this._psycho_loudness_approx(fftenergy, gfc);
        }
    }

    /** @private */
    _vbrpsy_attack_detection(gfp, buffer, bufPos, gr_out, masking_ratio, masking_MS_ratio, energy, sub_short_factor, ns_attacks, uselongblock) {
        const ns_hpfsmpl = new_float_n([2, 576]);
        const gfc = gfp.internal_flags;
        const n_chn_out = gfc.channels_out;
        const n_chn_psy = (gfp.mode === MPEGMode.JOINT_STEREO) ? 4 : n_chn_out;
        for (let chn = 0; chn < n_chn_out; chn++) {
            const firbuf = buffer[chn]; const input_start_index = bufPos + 288 - (NSFIRLEN - 1) / 2;
            for (let i = 0; i < 576; i++) {
                let c_sum1 = 0.0, c_sum2 = 0.0; const base_idx = input_start_index + i;
                c_sum1 = firbuf[base_idx + 10];
                for (let j = 0; j < (NSFIRLEN - 1) / 2; j += 2) {
                    c_sum1 += fircoef_[j] * (firbuf[base_idx + j] + firbuf[base_idx + NSFIRLEN - 1 - j]);
                    c_sum2 += fircoef_[j + 1] * (firbuf[base_idx + j + 1] + firbuf[base_idx + NSFIRLEN - 1 - (j + 1)]);
                } ns_hpfsmpl[chn][i] = c_sum1 + c_sum2;
            }
            masking_ratio[gr_out][chn].en.assign(gfc.en[chn]); masking_ratio[gr_out][chn].thm.assign(gfc.thm[chn]);
            if (n_chn_psy > 2) { masking_MS_ratio[gr_out][chn].en.assign(gfc.en[chn + 2]); masking_MS_ratio[gr_out][chn].thm.assign(gfc.thm[chn + 2]); }
        }
        for (let chn = 0; chn < n_chn_psy; chn++) {
            const attack_intensity = new_float(12); const en_subshort = new_float(12); const en_short = [0.0, 0.0, 0.0, 0.0];
            const attackThreshold = (chn === 3) ? NSATTACKTHRE_S : NSATTACKTHRE; let ns_uselongblock = 1;
            if (chn === 2) { for (let i = 0; i < 576; ++i) { const l = ns_hpfsmpl[0][i]; const r = ns_hpfsmpl[1][i]; ns_hpfsmpl[0][i] = l + r; ns_hpfsmpl[1][i] = l - r; } }
             const pf = ns_hpfsmpl[chn & 1];
            let pfPos = 0;
            for (let i = 0; i < 3; i++) { en_subshort[i] = gfc.nsPsy.last_en_subshort[chn][i + 6]; assert(gfc.nsPsy.last_en_subshort[chn][i + 4] > 0); attack_intensity[i] = en_subshort[i] / gfc.nsPsy.last_en_subshort[chn][i + 4]; en_short[0] += en_subshort[i]; }
            for (let i = 0; i < 9; i++) {
                const pfe = pfPos + 64; let p = 1.0; for (; pfPos < pfe; pfPos++) { const abs_sample = Math.abs(pf[pfPos]); if (p < abs_sample) p = abs_sample; }
                gfc.nsPsy.last_en_subshort[chn][i] = en_subshort[i + 3] = p; en_short[1 + Math.floor(i / 3)] += p;
                const prev_en = en_subshort[i + 3 - 2];
                if (p > prev_en) { assert(prev_en > 0); attack_intensity[i + 3] = p / prev_en; }
                else if (prev_en > p * 10.0) { assert(p > 0); attack_intensity[i + 3] = prev_en / (p * 10.0); }
                else attack_intensity[i + 3] = 0.0;
            }
            for (let i = 0; i < 3; ++i) {
                const sub_idx = i * 3 + 3; const enn = en_subshort[sub_idx] + en_subshort[sub_idx + 1] + en_subshort[sub_idx + 2]; let factor = 1.0;
                if(enn > 0) { if (en_subshort[sub_idx + 2] * 6 < enn) { factor *= 0.5; if (en_subshort[sub_idx + 1] * 6 < enn) factor *= 0.5; } }
                sub_short_factor[chn][i] = factor;
            }
            if (gfp.analysis) { let max_attack = attack_intensity[0]; for (let i = 1; i < 12; i++) if (max_attack < attack_intensity[i]) max_attack = attack_intensity[i]; gfc.pinfo.ers[gr_out][chn] = gfc.pinfo.ers_save[chn]; gfc.pinfo.ers_save[chn] = max_attack; }
            Arrays.fill(ns_attacks[chn], 0);
            for (let i = 0; i < 12; i++) { const block_idx = Math.floor(i / 3); if (ns_attacks[chn][block_idx] === 0 && attack_intensity[i] > attackThreshold) ns_attacks[chn][block_idx] = (i % 3) + 1; }
            for (let i = 1; i < 4; i++) {
                const u = en_short[i - 1]; const v = en_short[i]; const m = Math.max(u, v);
                if (m < 40000) { if (u < 1.7 * v && v < 1.7 * u) { if (i === 1 && ns_attacks[chn][0] <= ns_attacks[chn][i]) ns_attacks[chn][0] = 0; ns_attacks[chn][i] = 0; } }
            }
            if (ns_attacks[chn][0] <= gfc.nsPsy.lastAttacks[chn]) ns_attacks[chn][0] = 0;
            if (gfc.nsPsy.lastAttacks[chn] === 3 || (ns_attacks[chn][0] + ns_attacks[chn][1] + ns_attacks[chn][2] + ns_attacks[chn][3]) !== 0) {
                ns_uselongblock = 0;
                if (ns_attacks[chn][1] !== 0 && ns_attacks[chn][0] !== 0) ns_attacks[chn][1] = 0;
                if (ns_attacks[chn][2] !== 0 && ns_attacks[chn][1] !== 0) ns_attacks[chn][2] = 0;
                if (ns_attacks[chn][3] !== 0 && ns_attacks[chn][2] !== 0) ns_attacks[chn][3] = 0;
            }
            if (chn < 2) uselongblock[chn] = ns_uselongblock;
            else if (ns_uselongblock === 0) uselongblock[0] = uselongblock[1] = 0;
            energy[chn] = gfc.tot_ener[chn];
        }
    }

    /** @private */
    _vbrpsy_skip_masking_s(gfc, chn, sblock) {
        if (sblock === 0) { for (let b = 0; b < gfc.npart_s; b++) { gfc.nb_s2[chn][b] = gfc.nb_s1[chn][b]; gfc.nb_s1[chn][b] = 0; } }
    }

    /** @private */
    _vbrpsy_skip_masking_l(gfc, chn) {
        for (let b = 0; b < gfc.npart_l; b++) { gfc.nb_2[chn][b] = gfc.nb_1[chn][b]; gfc.nb_1[chn][b] = 0; }
    }

    /** @private */
    _s3_func_x(bark, hf_slope) {
        let tempx = bark; let tempy;
        if (tempx >= 0) tempy = -tempx * 27.0;
        else tempy = tempx * hf_slope;
        if (tempy <= -72.0) return 0.0;
        return Math.exp(tempy * LN_TO_LOG10);
    }

    /** @private */
    _norm_s3_func_x(hf_slope) {
        let lim_a = 0, lim_b = 0;
        { let x = 0, l, h; for (x = 0; this._s3_func_x(x, hf_slope) > 1e-20; x -= 1) ; l = x; h = 0; while (Math.abs(h - l) > 1e-12) { x = (h + l) / 2.0; if (this._s3_func_x(x, hf_slope) > 0) h = x; else l = x; } lim_a = l; }
        { let x = 0, l, h; for (x = 0; this._s3_func_x(x, hf_slope) > 1e-20; x += 1) ; l = 0; h = x; while (Math.abs(h - l) > 1e-12) { x = (h + l) / 2.0; if (this._s3_func_x(x, hf_slope) > 0) l = x; else h = x; } lim_b = h; }
        { let sum = 0; const m = 1000; const step = (lim_b - lim_a) / m; for (let i = 0; i <= m; ++i) { const x = lim_a + i * step; sum += this._s3_func_x(x, hf_slope); } if (sum === 0 || lim_b === lim_a) return 1.0; const norm = (m + 1) / (sum * (lim_b - lim_a)); return norm; }
    }

    /** @private */
    _s3_func(bark) {
        let tempx, x, tempy, temp; tempx = bark;
        if (tempx >= 0) tempx *= 3.0; else tempx *= 1.5;
        if (tempx >= 0.5 && tempx <= 2.5) { temp = tempx - 0.5; x = 8.0 * (temp * temp - 2.0 * temp); } else x = 0.0;
        tempx += 0.474; tempy = 15.811389 + 7.5 * tempx - 17.5 * Math.sqrt(1.0 + tempx * tempx);
        if (tempy <= -60.0) return 0.0;
        tempx = Math.exp((x + tempy) * LN_TO_LOG10); tempx /= 0.6609193; return tempx;
    }

    /** @private */
    _freq2bark(freq) {
        if (freq < 0) freq = 0; freq = freq * 0.001;
        return 13.0 * Math.atan(0.76 * freq) + 3.5 * Math.atan(freq * freq / (7.5 * 7.5));
    }

    /** @private */
    _init_numline(numlines, bo, bm, bval, bval_width, mld, bo_w, sfreq, blksize, scalepos, deltafreq, sbmax) {
        const b_frq = new_float(Encoder.CBANDS + 1);
        const granule_samples = (sbmax > 15) ? 576 : 192;
        const sample_freq_frac = sfreq / (2.0 * granule_samples);
        const partition = new_int(Encoder.HBLKSIZE);
        let ni = 0;
        const freq_per_line = sfreq / blksize;
        let fft_line_idx = 0;
        for (let part_idx = 0; part_idx < Encoder.CBANDS; part_idx++) {
            const bark_start = this._freq2bark(freq_per_line * fft_line_idx);
            b_frq[part_idx] = freq_per_line * fft_line_idx;
            let fft_line_end = fft_line_idx;
            while (this._freq2bark(freq_per_line * fft_line_end) - bark_start < DELBARK && fft_line_end <= blksize / 2) { fft_line_end++; }
            numlines[part_idx] = fft_line_end - fft_line_idx;
            ni = part_idx + 1;
            while (fft_line_idx < fft_line_end) { assert(fft_line_idx < Encoder.HBLKSIZE); partition[fft_line_idx++] = part_idx; }
            if (fft_line_idx > blksize / 2) { fft_line_idx = blksize / 2; break; }
        }
        assert(ni < Encoder.CBANDS);
        b_frq[ni] = freq_per_line * fft_line_idx;
        for (let sfb = 0; sfb < sbmax; sfb++) {
            const sfb_start_line = scalepos[sfb]; const sfb_end_line = scalepos[sfb + 1];
            let i1 = 0 | Math.floor(0.5 + deltafreq * (sfb_start_line - 0.5)); if (i1 < 0) i1 = 0;
            let i2 = 0 | Math.floor(0.5 + deltafreq * (sfb_end_line - 0.5)); if (i2 >= blksize / 2) i2 = blksize / 2 -1; if(i2 < i1) i2 = i1;
            const part1 = (i1 < Encoder.HBLKSIZE) ? partition[i1] : ni - 1; const part2 = (i2 < Encoder.HBLKSIZE) ? partition[i2] : ni - 1;
            bm[sfb] = Math.floor((part1 + part2) / 2); bo[sfb] = part2;
            const f_sfb_end = sample_freq_frac * sfb_end_line;
            const f_part_start = b_frq[bo[sfb]]; const f_part_end = b_frq[bo[sfb] + 1];
             if (f_part_end > f_part_start) bo_w[sfb] = (f_sfb_end - f_part_start) / (f_part_end - f_part_start);
             else bo_w[sfb] = (f_sfb_end >= f_part_start) ? 1.0 : 0.0;
             if (bo_w[sfb] < 0) bo_w[sfb] = 0.0; else if (bo_w[sfb] > 1) bo_w[sfb] = 1.0;
             const mld_freq = freq_per_line * i1; mld[sfb] = this._stereo_demask(mld_freq);
        }
        let current_fft_line = 0;
        for (let k = 0; k < ni; k++) {
            const w = numlines[k];
            if (w > 0) {
                 const bark1 = this._freq2bark(freq_per_line * current_fft_line); const bark2 = this._freq2bark(freq_per_line * (current_fft_line + w - 1));
                 bval[k] = 0.5 * (bark1 + bark2);
                 const bark_edge1 = this._freq2bark(freq_per_line * (current_fft_line - 0.5)); const bark_edge2 = this._freq2bark(freq_per_line * (current_fft_line + w - 0.5));
                 bval_width[k] = bark_edge2 - bark_edge1;
            } else { bval[k] = (k > 0) ? bval[k-1] : 0; bval_width[k] = 0; }
            current_fft_line += w;
        }
        return ni;
    }

    /** @private */
    _init_s3_values(npart, bval, bval_width, norm, use_old_s3) {
        const s3 = new_float_n([Encoder.CBANDS, Encoder.CBANDS]); let numberOfNoneZero = 0;
        if (use_old_s3) {
             for (let i = 0; i < npart; i++) for (let j = 0; j < npart; j++) s3[i][j] = this._s3_func(bval[i] - bval[j]) * bval_width[j] * norm[i];
        } else {
            for (let j = 0; j < npart; j++) {
                const hf_slope = 15.0 + Math.min(21.0 / (bval[j]>0?bval[j]:1e-6), 12.0); // Avoid div by zero
                const s3_x_norm = this._norm_s3_func_x(hf_slope);
                for (let i = 0; i < npart; i++) s3[i][j] = s3_x_norm * this._s3_func_x(bval[i] - bval[j], hf_slope) * bval_width[j] * norm[i];
            }
        }

        // Initialize s3ind with default values
        for (let i = 0; i < npart; i++) {
            this.s3ind[i][0] = 0;
            this.s3ind[i][1] = npart - 1;
        }

        // Then update with actual masking ranges
        for (let i = 0; i < npart; i++) {
             let start_masker, end_masker;
             for (start_masker = 0; start_masker < npart; start_masker++) if (s3[i][start_masker] > 1e-20) break; // Use threshold
             for (end_masker = npart - 1; end_masker >= start_masker; end_masker--) if (s3[i][end_masker] > 1e-20) break;
             if (start_masker <= end_masker) { this.s3ind[i][0] = start_masker; this.s3ind[i][1] = end_masker; numberOfNoneZero += (end_masker - start_masker + 1); }
             else { this.s3ind[i][0] = 0; this.s3ind[i][1] = -1; }
        }
        const p = new_float(numberOfNoneZero); let k = 0;
        for (let i = 0; i < npart; i++) { const start = this.s3ind[i][0]; const end = this.s3ind[i][1]; if (start <= end) for (let j = start; j <= end; j++) p[k++] = s3[i][j]; }
        assert(k === numberOfNoneZero); return p;
    }

    /** @private */
    _stereo_demask(f) {
        let arg = this._freq2bark(f); arg = (Math.min(arg, 15.5) / 15.5);
        return Math.pow(10.0, 1.25 * (1.0 - Math.cos(Math.PI * arg)) - 2.5);
    }

    /** @private */
    _ATHformula_GB(f, value) {
        if (f < -0.3) f = 3410.0;
        f /= 1000.0; f = Math.max(0.1, f);
        const ath = 3.640 * Math.pow(f, -0.8) - 6.800 * Math.exp(-0.6 * Math.pow(f - 3.4, 2.0)) + 6.000 * Math.exp(-0.15 * Math.pow(f - 8.7, 2.0)) + (0.6 + 0.04 * value) * 0.001 * Math.pow(f, 4.0);
        return ath;
    }


    // --- Public Methods ---

    /**
     * Performs psychoacoustic analysis using the NSPSYTUNE model modifications.
     * Computes masking thresholds, perceptual entropy, and determines block types
     * based on attack detection using a high-pass filtered signal.
     * Returns results delayed by one granule.
     *
     * @public
     * @param {object} gfp - LAME global flags and settings.
     * @param {Array<Float32Array>} buffer - Input PCM buffer [channels][samples]. Contains 1152 samples per channel.
     * @param {number} bufPos - Starting index within the buffer (typically 0).
     * @param {number} gr_out - Granule index (0 or 1) for storing output results.
     * @param {Array<Array<object>>} masking_ratio - Output: Array [2][2] storing masking results (en/thm) for L/R channels of the *previous* granule. Each element has `.en` and `.thm` properties, which are themselves objects with `.l` (Float32Array[SBMAX_l]) and `.s` (Array<Float32Array>[SBMAX_s][3]) properties.
     * @param {Array<Array<object>>} masking_MS_ratio - Output: Array [2][2] storing masking results (en/thm) for M/S channels of the *previous* granule. Structure same as `masking_ratio`.
     * @param {Float32Array} percep_entropy - Output: Array [2] for perceptual entropy (PE) for L/R channels of the *previous* granule.
     * @param {Float32Array} percep_MS_entropy - Output: Array [2] for perceptual entropy (PE) for M/S channels of the *previous* granule.
     * @param {Float32Array} energy - Output: Array [4] for total energy (L, R, M, S) of the *previous* granule.
     * @param {Int32Array} blocktype_d - Output: Array [2] indicating the determined block type (see `Encoder` constants like `NORM_TYPE`, `SHORT_TYPE` etc.) for L/R channels of the *previous* granule.
     * @returns {number} Status code (0 for success).
     */
    L3psycho_anal_ns(gfp, buffer, bufPos, gr_out, masking_ratio, masking_MS_ratio, percep_entropy, percep_MS_entropy, energy, blocktype_d) {
        const gfc = gfp.internal_flags;
        const wsamp_L = new_float_n([2, Encoder.BLKSIZE]);
        const wsamp_S = new_float_n([2, 3, Encoder.BLKSIZE_s]);
        const eb_l = new_float(Encoder.CBANDS + 1);
        const eb_s = new_float(Encoder.CBANDS + 1);
        const thr = new_float(Encoder.CBANDS + 2);
        const blocktype = new_int(2);
        const uselongblock = new_int(2);
        const ns_hpfsmpl = new_float_n([2, 576]);
        let pcfact;
        const mask_idx_l = new_int(Encoder.CBANDS + 2);

        let numchn = gfc.channels_out;
        if (gfp.mode === MPEGMode.JOINT_STEREO) numchn = 4;

        if (gfp.VBR === VbrMode.vbr_off) pcfact = gfc.ResvMax === 0 ? 0 : (gfc.ResvSize / gfc.ResvMax) * 0.5;
        else if (gfp.VBR === VbrMode.vbr_rh || gfp.VBR === VbrMode.vbr_mtrh || gfp.VBR === VbrMode.vbr_mt) pcfact = 0.6;
        else pcfact = 1.0;

        // Apply HPF and copy previous results
        for (let chn = 0; chn < gfc.channels_out; chn++) {
            const firbuf = buffer[chn];
            const input_start_index = bufPos + 288 - (NSFIRLEN - 1) / 2; // Center 576 samples start at 288 in 1152 buffer
            for (let i = 0; i < 576; i++) {
                let c_sum1 = 0.0, c_sum2 = 0.0;
                const base_idx = input_start_index + i;
                c_sum1 = firbuf[base_idx + 10];
                for (let j = 0; j < (NSFIRLEN - 1) / 2; j += 2) {
                    c_sum1 += fircoef[j] * (firbuf[base_idx + j] + firbuf[base_idx + NSFIRLEN - 1 - j]);
                    c_sum2 += fircoef[j + 1] * (firbuf[base_idx + j + 1] + firbuf[base_idx + NSFIRLEN - 1 - (j + 1)]);
                }
                ns_hpfsmpl[chn][i] = c_sum1 + c_sum2;
            }
            masking_ratio[gr_out][chn].en.assign(gfc.en[chn]);
            masking_ratio[gr_out][chn].thm.assign(gfc.thm[chn]);
            if (numchn > 2) {
                masking_MS_ratio[gr_out][chn].en.assign(gfc.en[chn + 2]);
                masking_MS_ratio[gr_out][chn].thm.assign(gfc.thm[chn + 2]);
            }
        }

        // Main processing loop per channel
        for (let chn = 0; chn < numchn; chn++) {
            const current_wsamp_l_idx = chn & 1;
            const current_wsamp_s_idx = chn & 1;
            const en_subshort = new_float(12);
            const en_short = [0.0, 0.0, 0.0, 0.0];
            const attack_intensity = new_float(12);
            let ns_uselongblock = 1;
            const max = new_float(Encoder.CBANDS);
            const avg = new_float(Encoder.CBANDS);
            const ns_attacks = [0, 0, 0, 0]; // Renamed from C's static/global-like array
            const fftenergy = new_float(Encoder.HBLKSIZE);
            const fftenergy_s = new_float_n([3, Encoder.HBLKSIZE_s]);

            // Block Type Determination (Attack Detection)
             for (let i = 0; i < 3; i++) { en_subshort[i] = gfc.nsPsy.last_en_subshort[chn][i + 6]; assert(gfc.nsPsy.last_en_subshort[chn][i + 4] > 0); attack_intensity[i] = en_subshort[i] / gfc.nsPsy.last_en_subshort[chn][i + 4]; en_short[0] += en_subshort[i]; }
             if (chn === 2) { for (let i = 0; i < 576; i++) { const l = ns_hpfsmpl[0][i]; const r = ns_hpfsmpl[1][i]; ns_hpfsmpl[0][i] = l + r; ns_hpfsmpl[1][i] = l - r; } }
             { const pf = ns_hpfsmpl[chn & 1]; let pfPos = 0; for (let i = 0; i < 9; i++) { const pfe = pfPos + 64; let p = 1.0; for (; pfPos < pfe; pfPos++) { const abs_sample = Math.abs(pf[pfPos]); if (p < abs_sample) p = abs_sample; } gfc.nsPsy.last_en_subshort[chn][i] = en_subshort[i + 3] = p; en_short[1 + Math.floor(i / 3)] += p; const prev_en = en_subshort[i + 3 - 2]; if (p > prev_en) { assert(prev_en > 0); attack_intensity[i + 3] = p / prev_en; } else if (prev_en > p * 10.0) { assert(p > 0); attack_intensity[i + 3] = prev_en / (p * 10.0); } else attack_intensity[i + 3] = 0.0; } assert(pfPos === 576);}
             if (gfp.analysis) { let max_attack = attack_intensity[0]; for (let i = 1; i < 12; i++) if (max_attack < attack_intensity[i]) max_attack = attack_intensity[i]; gfc.pinfo.ers[gr_out][chn] = gfc.pinfo.ers_save[chn]; gfc.pinfo.ers_save[chn] = max_attack; }
             const attackThreshold = (chn === 3) ? NSATTACKTHRE_S : NSATTACKTHRE; for (let i = 0; i < 12; i++) { if (ns_attacks[Math.floor(i / 3)] === 0 && attack_intensity[i] > attackThreshold) ns_attacks[Math.floor(i / 3)] = (i % 3) + 1; }
             for (let i = 1; i < 4; i++) { const en_curr = en_short[i]; const en_prev = en_short[i - 1]; let ratio = 1.0; if (en_prev > 0 && en_curr > 0) ratio = (en_prev > en_curr) ? (en_prev / en_curr) : (en_curr / en_prev); else if (en_prev > 0 || en_curr > 0) ratio = 100; if (ratio < 1.7) { ns_attacks[i] = 0; if (i === 1) ns_attacks[0] = 0; } }
             if (ns_attacks[0] !== 0 && gfc.nsPsy.lastAttacks[chn] !== 0) ns_attacks[0] = 0;
             if (gfc.nsPsy.lastAttacks[chn] === 3 || (ns_attacks[0] + ns_attacks[1] + ns_attacks[2] + ns_attacks[3]) !== 0) { ns_uselongblock = 0; if (ns_attacks[1] !== 0 && ns_attacks[0] !== 0) ns_attacks[1] = 0; if (ns_attacks[2] !== 0 && ns_attacks[1] !== 0) ns_attacks[2] = 0; if (ns_attacks[3] !== 0 && ns_attacks[2] !== 0) ns_attacks[3] = 0; }
             if (chn < 2) uselongblock[chn] = ns_uselongblock; else if (ns_uselongblock === 0) uselongblock[0] = uselongblock[1] = 0;

            energy[chn] = gfc.tot_ener[chn]; // Store previous granule's energy

            // Compute FFTs
            this._compute_ffts(gfp, fftenergy, fftenergy_s, wsamp_L, current_wsamp_l_idx, wsamp_S, current_wsamp_s_idx, gr_out, chn, buffer, bufPos);

            // Compute Masking Thresholds - Short Blocks
            for (let sblock = 0; sblock < 3; sblock++) {
                this._compute_masking_s(gfp, fftenergy_s, eb_s, thr, chn, sblock);
                this._convert_partition2scalefac_s(gfc, eb_s, thr, chn, sblock);
                // Short block pre-echo control
                 for (let sb = 0; sb < Encoder.SBMAX_s; sb++) { let thmm = gfc.thm[chn].s[sb][sblock]; thmm *= NS_PREECHO_ATT0; if (ns_attacks[sblock] >= 2 || ns_attacks[sblock + 1] === 1) { const idx = (sblock !== 0) ? sblock - 1 : 2; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT1 * pcfact); thmm = Math.min(thmm, p); } if (ns_attacks[sblock] === 1) { const idx = (sblock !== 0) ? sblock - 1 : 2; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT2 * pcfact); thmm = Math.min(thmm, p); } else if ((sblock !== 0 && ns_attacks[sblock - 1] === 3) || (sblock === 0 && gfc.nsPsy.lastAttacks[chn] === 3)) { const idx = (sblock !== 2) ? sblock + 1 : 0; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT2 * pcfact); thmm = Math.min(thmm, p); } const sub_en_idx = sblock * 3 + 3; const enn = en_subshort[sub_en_idx] + en_subshort[sub_en_idx + 1] + en_subshort[sub_en_idx + 2]; if (enn > 0 && en_subshort[sub_en_idx + 2] * 6 < enn) { thmm *= 0.5; if (en_subshort[sub_en_idx + 1] * 6 < enn) thmm *= 0.5; } gfc.thm[chn].s[sb][sblock] = thmm; }
            }
            gfc.nsPsy.lastAttacks[chn] = ns_attacks[2]; // Store attack state of last short block

            // Compute Masking Thresholds - Long Blocks
            this._calc_energy(gfc, fftenergy, eb_l, max, avg);
            this._calc_mask_index_l(gfc, max, avg, mask_idx_l);
            let s3_idx = 0;
            for (let b = 0; b < gfc.npart_l; b++) {
                 const first_masker_idx = gfc.s3ind[b][0]; const last_masker_idx = gfc.s3ind[b][1]; let kk = first_masker_idx; let ecb = 0;
                 if (kk <= last_masker_idx) { let masker_energy = eb_l[kk] * tab[mask_idx_l[kk]]; ecb = gfc.s3_ll[s3_idx] * masker_energy; s3_idx++; kk++; while (kk <= last_masker_idx) { masker_energy = eb_l[kk] * tab[mask_idx_l[kk]]; const term = gfc.s3_ll[s3_idx] * masker_energy; ecb = this._mask_add(ecb, term, kk, kk - b, gfc, 0); s3_idx++; kk++; } }
                 ecb *= 0.158489319246111;
                 // Long block pre-echo control
                 if (gfc.blocktype_old[chn & 1] === Encoder.SHORT_TYPE) thr[b] = ecb; else { const min_limit = Math.min(rpelev * gfc.nb_1[chn][b], rpelev2 * gfc.nb_2[chn][b]); const clamped_ecb = Math.min(ecb, min_limit); thr[b] = this._NS_INTERP(clamped_ecb, ecb, pcfact); }
                 gfc.nb_2[chn][b] = gfc.nb_1[chn][b]; gfc.nb_1[chn][b] = ecb;
            }
             for (let b = gfc.npart_l; b <= Encoder.CBANDS; ++b) { eb_l[b] = 0; thr[b] = 0; } if (gfc.npart_l <= Encoder.CBANDS + 1) thr[gfc.npart_l+1]=0;

            this._convert_partition2scalefac_l(gfc, eb_l, thr, chn);

        } // End loop over chn

        // Post-processing
        if (gfp.mode === MPEGMode.STEREO || gfp.mode === MPEGMode.JOINT_STEREO) {
            if (gfp.interChRatio > 0.0) this._calc_interchannel_masking(gfp, gfp.interChRatio);
        }
        if (gfp.mode === MPEGMode.JOINT_STEREO) {
            this._msfix1(gfc);
            const msfix = gfp.msfix;
            if (Math.abs(msfix) > 0.0) {
                const ath_adj_factor = Math.pow(10, gfp.ATHlower * gfc.ATH.adjust); // Calculate linear factor
                this._ns_msfix(gfc, msfix, ath_adj_factor);
            }
        }

        // Determine final block type for previous granule
        this._block_type_set(gfp, uselongblock, blocktype_d, blocktype);

        // Compute PE for previous granule
        for (let chn = 0; chn < numchn; chn++) {
            let ppe_array, ppe_offset, granule_block_type, mr;
            if (chn > 1) { ppe_array = percep_MS_entropy; ppe_offset = -2; granule_block_type = Encoder.NORM_TYPE; if (blocktype_d[0] === Encoder.SHORT_TYPE || blocktype_d[1] === Encoder.SHORT_TYPE) granule_block_type = Encoder.SHORT_TYPE; mr = masking_MS_ratio[gr_out][chn - 2]; }
            else { ppe_array = percep_entropy; ppe_offset = 0; granule_block_type = blocktype_d[chn]; mr = masking_ratio[gr_out][chn]; }
            if (granule_block_type === Encoder.SHORT_TYPE || granule_block_type === Encoder.STOP_TYPE || granule_block_type === Encoder.START_TYPE) ppe_array[ppe_offset + chn] = this._pecalc_s(mr, gfc.masking_lower);
            else ppe_array[ppe_offset + chn] = this._pecalc_l(mr, gfc.masking_lower);
            if (gfp.analysis) gfc.pinfo.pe[gr_out][chn] = ppe_array[ppe_offset + chn];
        }

        return 0; // Success
    }

    /**
     * Performs psychoacoustic analysis using the VBR model.
     * This version optimizes by calculating only the necessary block type (long or short)
     * based on attack detection. It uses VBR-specific masking calculations and adjustments.
     * Returns results delayed by one granule.
     *
     * @public
     * @param {object} gfp - LAME global flags and settings.
     * @param {Array<Float32Array>} buffer - Input PCM buffer [channels][samples]. Contains 1152 samples per channel.
     * @param {number} bufPos - Starting index within the buffer (typically 0).
     * @param {number} gr_out - Granule index (0 or 1) for storing output results.
     * @param {Array<Array<object>>} masking_ratio - Output: Array [2][2] storing masking results (en/thm) for L/R channels of the *previous* granule. Structure same as `L3psycho_anal_ns`.
     * @param {Array<Array<object>>} masking_MS_ratio - Output: Array [2][2] storing masking results (en/thm) for M/S channels of the *previous* granule. Structure same as `L3psycho_anal_ns`.
     * @param {Float32Array} percep_entropy - Output: Array [2] for perceptual entropy (PE) for L/R channels of the *previous* granule.
     * @param {Float32Array} percep_MS_entropy - Output: Array [2] for perceptual entropy (PE) for M/S channels of the *previous* granule.
     * @param {Float32Array} energy - Output: Array [4] for total energy (L, R, M, S) of the *previous* granule.
     * @param {Int32Array} blocktype_d - Output: Array [2] indicating the determined block type (see `Encoder` constants) for L/R channels of the *previous* granule.
     * @returns {number} Status code (0 for success).
     */
    L3psycho_anal_vbr(gfp, buffer, bufPos, gr_out, masking_ratio, masking_MS_ratio, percep_entropy, percep_MS_entropy, energy, blocktype_d) {
        const gfc = gfp.internal_flags;
        const wsamp_L = new_float_n([2, Encoder.BLKSIZE]);
        const wsamp_S = new_float_n([2, 3, Encoder.BLKSIZE_s]);
        const fftenergy = new_float(Encoder.HBLKSIZE);
        const fftenergy_s = new_float_n([3, Encoder.HBLKSIZE_s]);
        const eb = new_float_n([4, Encoder.CBANDS + 1]);
        const thr = new_float_n([4, Encoder.CBANDS + 2]);
        const sub_short_factor = new_float_n([4, 3]);
        const pcfact = 0.6; // Fixed for VBR?
        const ns_attacks = new_array_n([4], () => new_int(4));
        const uselongblock = new_int(2);
        const n_chn_psy = (gfp.mode === MPEGMode.JOINT_STEREO) ? 4 : gfc.channels_out;

        // 1. Attack detection
        this._vbrpsy_attack_detection(gfp, buffer, bufPos, gr_out, masking_ratio, masking_MS_ratio, energy, sub_short_factor, ns_attacks, uselongblock);

        // 2. Compute initial block type choice
        this._vbrpsy_compute_block_type(gfp, uselongblock);

        // 3. Compute Masking - Long Blocks (conditional)
        for (let chn = 0; chn < n_chn_psy; chn++) {
             const ch01 = chn & 0x01;
             this._vbrpsy_compute_fft_l(gfp, buffer, bufPos, chn, gr_out, fftenergy, wsamp_L, ch01);
             this._vbrpsy_compute_loudness_approximation_l(gfp, gr_out, chn, fftenergy);
             if (uselongblock[ch01] !== 0) this._vbrpsy_compute_masking_l(gfc, fftenergy, eb[chn], thr[chn], chn);
             else this._vbrpsy_skip_masking_l(gfc, chn);
        }
        if (uselongblock[0] !== 0 && uselongblock[1] !== 0) {
            if (gfp.mode === MPEGMode.JOINT_STEREO) {
                 const ath_adj_factor = Math.pow(10, gfp.ATHlower * gfc.ATH.adjust); // Linear factor
                 this._vbrpsy_compute_MS_thresholds(eb, thr, gfc.mld_cb_l, gfc.ATH.cb_l, ath_adj_factor, gfp.msfix, gfc.npart_l);
            }
             for (let chn = 0; chn < n_chn_psy; chn++) this._convert_partition2scalefac_l(gfc, eb[chn], thr[chn], chn);
        }


        // 4. Compute Masking - Short Blocks (conditional)
        const compute_short = (uselongblock[0] === 0 || uselongblock[1] === 0);
        if (compute_short) {
            for (let sblock = 0; sblock < 3; sblock++) {
                for (let chn = 0; chn < n_chn_psy; ++chn) {
                    const ch01 = chn & 0x01;
                    if (uselongblock[ch01] === 0) {
                        this._vbrpsy_compute_fft_s(gfp, buffer, bufPos, chn, sblock, fftenergy_s, wsamp_S, ch01);
                        this._vbrpsy_compute_masking_s(gfp, fftenergy_s, eb[chn], thr[chn], chn, sblock);
                    } else { this._vbrpsy_skip_masking_s(gfc, chn, sblock); }
                }
                if (uselongblock[0] === 0 && uselongblock[1] === 0) {
                    if (gfp.mode === MPEGMode.JOINT_STEREO) {
                         const ath_adj_factor = Math.pow(10, gfp.ATHlower * gfc.ATH.adjust) * (Encoder.BLKSIZE_s / Encoder.BLKSIZE); // Scaled linear factor
                         this._vbrpsy_compute_MS_thresholds(eb, thr, gfc.mld_cb_s, gfc.ATH.cb_s, ath_adj_factor, gfp.msfix, gfc.npart_s);
                    }
                     for (let chn = 0; chn < n_chn_psy; ++chn) this._convert_partition2scalefac_s(gfc, eb[chn], thr[chn], chn, sblock);
                } else {
                     for (let chn = 0; chn < n_chn_psy; ++chn) if (uselongblock[chn & 0x01] === 0) this._convert_partition2scalefac_s(gfc, eb[chn], thr[chn], chn, sblock);
                }
            }
            // Short block pre-echo control
            for (let chn = 0; chn < n_chn_psy; chn++) {
                const ch01 = chn & 0x01; if (uselongblock[ch01] === 0) {
                    const new_thmm = new_float(3);
                    for (let sb = 0; sb < Encoder.SBMAX_s; sb++) {
                        for (let sblock = 0; sblock < 3; sblock++) {
                            let thmm = gfc.thm[chn].s[sb][sblock]; thmm *= NS_PREECHO_ATT0;
                             if (ns_attacks[chn][sblock] >= 2 || ns_attacks[chn][sblock + 1] === 1) { const idx = (sblock !== 0) ? sblock - 1 : 2; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT1 * pcfact); thmm = Math.min(thmm, p); }
                             else if (ns_attacks[chn][sblock] === 1) { const idx = (sblock !== 0) ? sblock - 1 : 2; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT2 * pcfact); thmm = Math.min(thmm, p); }
                             else if ((sblock !== 0 && ns_attacks[chn][sblock - 1] === 3) || (sblock === 0 && gfc.nsPsy.lastAttacks[chn] === 3)) { const idx = (sblock !== 2) ? sblock + 1 : 0; const p = this._NS_INTERP(gfc.thm[chn].s[sb][idx], thmm, NS_PREECHO_ATT2 * pcfact); thmm = Math.min(thmm, p); }
                             thmm *= sub_short_factor[chn][sblock]; new_thmm[sblock] = thmm;
                        }
                        for (let sblock = 0; sblock < 3; sblock++) gfc.thm[chn].s[sb][sblock] = new_thmm[sblock];
                    }
                }
            }
        }

        // Store last attack state
        for (let chn = 0; chn < n_chn_psy; chn++) gfc.nsPsy.lastAttacks[chn] = ns_attacks[chn][2];

        // 5. Determine final block type for previous granule
        this._vbrpsy_apply_block_type(gfp, uselongblock, blocktype_d);

        // 6. Compute PE for previous granule
        for (let chn = 0; chn < n_chn_psy; chn++) {
            let ppe_array, ppe_offset, granule_block_type, mr;
            if (chn > 1) { ppe_array = percep_MS_entropy; ppe_offset = -2; granule_block_type = Encoder.NORM_TYPE; if (blocktype_d[0] === Encoder.SHORT_TYPE || blocktype_d[1] === Encoder.SHORT_TYPE) granule_block_type = Encoder.SHORT_TYPE; mr = masking_MS_ratio[gr_out][chn - 2]; }
            else { ppe_array = percep_entropy; ppe_offset = 0; granule_block_type = blocktype_d[chn]; mr = masking_ratio[gr_out][chn]; }
            if (granule_block_type === Encoder.SHORT_TYPE || granule_block_type === Encoder.STOP_TYPE || granule_block_type === Encoder.START_TYPE) ppe_array[ppe_offset + chn] = this._pecalc_s(mr, gfc.masking_lower);
            else ppe_array[ppe_offset + chn] = this._pecalc_l(mr, gfc.masking_lower);
            if (gfp.analysis) gfc.pinfo.pe[gr_out][chn] = ppe_array[ppe_offset + chn];
        }

        return 0; // Success
    }

    /**
     * Initializes the psychoacoustic model constants and workspaces.
     * Must be called before analysis. Calculates partition bands, spreading functions,
     * ATH, MLD, etc., based on sample rate and encoding parameters.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags and settings. Holds encoding parameters like sample rate, mode, quality settings, ATH type etc.
     * @returns {number} Status code (0 for success).
     */
    psymodel_init(gfp) {
        const gfc = gfp.internal_flags;
        let i, j, k;    
        let useOldS3 = true; let bvl_a = 13.0, bvl_b = 24.0;
        let snr_l_a = 0.0, snr_l_b = 0.0; let snr_s_a = -8.25, snr_s_b = -4.5;
        const bval = new_float(Encoder.CBANDS); const bval_width = new_float(Encoder.CBANDS);
        const norm = new_float(Encoder.CBANDS); const sfreq = gfp.out_samplerate;

        switch (gfp.experimentalZ) {
            default: case 0: useOldS3 = true; break;
            case 1: useOldS3 = !(gfp.VBR === VbrMode.vbr_mtrh || gfp.VBR === VbrMode.vbr_mt); break;
            case 2: useOldS3 = false; break;
            case 3: bvl_a = 8; snr_l_a = -1.75; snr_l_b = -0.0125; snr_s_a = -8.25; snr_s_b = -2.25; break;
        }
        gfc.ms_ener_ratio_old = .25; gfc.blocktype_old[0] = gfc.blocktype_old[1] = Encoder.NORM_TYPE;
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < Encoder.CBANDS; ++j) { gfc.nb_1[i][j] = 1e20; gfc.nb_2[i][j] = 1e20; gfc.nb_s1[i][j] = gfc.nb_s2[i][j] = 1.0; }
            if (gfc.en[i] && gfc.thm[i]) { // Check both en and thm
                for (let sb = 0; sb < Encoder.SBMAX_l; sb++) { gfc.en[i].l[sb] = 1e20; gfc.thm[i].l[sb] = 1e20; }
                for (j = 0; j < 3; ++j) for (let sb = 0; sb < Encoder.SBMAX_s; sb++) { gfc.en[i].s[sb][j] = 1e20; gfc.thm[i].s[sb][j] = 1e20; }
            } else {
                console.error(`gfc.en[${i}] or gfc.thm[${i}] is undefined in psymodel_init!`);
                return -1; // Critical error
            }
             gfc.nsPsy.lastAttacks[i] = 0;
             for (j = 0; j < 9; j++) gfc.nsPsy.last_en_subshort[i][j] = 10.;
        }
        gfc.loudness_sq_save[0] = gfc.loudness_sq_save[1] = 0.0;

        // Compute partition mappings, bark values, MLD for LONG blocks
        gfc.npart_l = this._init_numline(gfc.numlines_l, gfc.bo_l, gfc.bm_l, bval, bval_width, gfc.mld_l, gfc.PSY.bo_l_weight, sfreq, Encoder.BLKSIZE, gfc.scalefac_band.l, Encoder.BLKSIZE / (2.0 * 576), Encoder.SBMAX_l);
        assert(gfc.npart_l < Encoder.CBANDS);
        // Compute normalization factors (based on SNR) and reciprocal line counts
        for (i = 0; i < gfc.npart_l; i++) {
            let snr = snr_l_a; if (bval[i] >= bvl_a) snr = snr_l_b * (bval[i] - bvl_a) / (bvl_b - bvl_a) + snr_l_a * (bvl_b - bval[i]) / (bvl_b - bvl_a);
            norm[i] = Math.pow(10.0, snr / 10.0);
            gfc.rnumlines_l[i] = (gfc.numlines_l[i] > 0) ? 1.0 / gfc.numlines_l[i] : 0;
        }
        // Compute spreading function values (flattened array)
        gfc.s3_ll = this._init_s3_values(gfc.npart_l, bval, bval_width, norm, useOldS3);

        // Compute long block ATH and MINVAL per partition
        j = 0;
        for (i = 0; i < gfc.npart_l; i++) {
            let minval = Util.FLOAT_MAX;
            for (let l = 0; l < gfc.numlines_l[i]; l++, j++) {
                const freq = sfreq * j / (1000.0 * Encoder.BLKSIZE);
                let level = this.ATHformula(freq * 1000, gfp) - 20;
                level = Math.pow(10, 0.1 * level) * gfc.numlines_l[i];
                minval = Math.min(minval, level);
            }
            gfc.minval_l[i] = minval;
        }

        // Compute partition mappings, bark values, MLD for SHORT blocks
        gfc.npart_s = this._init_numline(gfc.numlines_s, gfc.bo_s, gfc.bm_s, bval, bval_width, gfc.mld_s, gfc.PSY.bo_s_weight, sfreq, Encoder.BLKSIZE_s, gfc.scalefac_band.s, Encoder.BLKSIZE_s / (2.0 * 192), Encoder.SBMAX_s);
        assert(gfc.npart_s < Encoder.CBANDS);
        // Compute short block normalization factors, ATH, MINVAL per partition
        j = 0;
        for (i = 0; i < gfc.npart_s; i++) {
             let x; let snr = snr_s_a; if (bval[i] >= bvl_a) snr = snr_s_b * (bval[i] - bvl_a) / (bvl_b - bvl_a) + snr_s_a * (bvl_b - bval[i]) / (bvl_b - bvl_a);
             norm[i] = Math.pow(10.0, snr / 10.0);
             x = Float.MAX_VALUE;
             for (k = 0; k < gfc.numlines_s[i]; k++, j++) { const freq = sfreq * j / (1000.0 * Encoder.BLKSIZE_s); let level = this.ATHformula(freq * 1000, gfp) - 20; level = Math.pow(10., 0.1 * level) * gfc.numlines_s[i]; if (x > level) x = level; }
             gfc.ATH.cb_s[i] = x;
             x = (-7.0 + bval[i] * 7.0 / 12.0); if (bval[i] > 12) x *= 1 + Math.log(1 + x) * 3.1; if (bval[i] < 12) x *= 1 + Math.log(1 - x) * 2.3; if (x < -15) x = -15; x -= 8; gfc.minval_s[i] = Math.pow(10.0, x / 10) * gfc.numlines_s[i];
        }
        // Compute short block spreading function values
        gfc.s3_ss = this._init_s3_values(gfc.npart_s, bval, bval_width, norm, useOldS3);

        // Initialize other parameters
        init_mask_add_max_values();
        this.fft.init_fft(gfc);
        gfc.decay = Math.exp(-1.0 * LOG10 / (temporalmask_sustain_sec * sfreq / 192.0));
        let msfix = NS_MSFIX; if ((gfp.exp_nspsytune & 2) != 0) msfix = 1.0; if (Math.abs(gfp.msfix) > 0.0) msfix = gfp.msfix; gfp.msfix = msfix;
        for (let b = 0; b < gfc.npart_l; b++) if (gfc.s3ind[b][1] > gfc.npart_l - 1) gfc.s3ind[b][1] = gfc.npart_l - 1;

        // ATH auto adjustment init
        const frame_duration = (576. * gfc.mode_gr / sfreq);
        gfc.ATH.decay = Math.pow(10., -12. / 10. * frame_duration);
        gfc.ATH.adjust = 0.01; gfc.ATH.adjustLimit = 1.0;

        // Equal loudness weights init
        if (gfp.ATHtype != -1) {
            const freq_inc = gfp.out_samplerate / (Encoder.BLKSIZE); let eql_balance = 0.0; let freq = 0.0;
            for (i = 0; i < Encoder.BLKSIZE / 2; ++i) { freq += freq_inc; gfc.ATH.eql_w[i] = 1. / Math.pow(10, this.ATHformula(freq, gfp) / 10); eql_balance += gfc.ATH.eql_w[i]; }
             if(eql_balance > 0) eql_balance = 1.0 / eql_balance; else eql_balance = 0; // Avoid division by zero
            for (i = Encoder.BLKSIZE / 2; --i >= 0;) gfc.ATH.eql_w[i] *= eql_balance;
        }

        // Partition band MLD calculation (using center frequency)
         j = 0; for (i = 0; i < gfc.npart_l; i++) { const freq = sfreq * (j + gfc.numlines_l[i] / 2.0) / (1.0 * Encoder.BLKSIZE); gfc.mld_cb_l[i] = this._stereo_demask(freq); j += gfc.numlines_l[i]; } for (; i < Encoder.CBANDS; ++i) gfc.mld_cb_l[i] = 1;
         j = 0; for (i = 0; i < gfc.npart_s; i++) { const freq = sfreq * (j + gfc.numlines_s[i] / 2.0) / (1.0 * Encoder.BLKSIZE_s); gfc.mld_cb_s[i] = this._stereo_demask(freq); j += gfc.numlines_s[i]; } for (; i < Encoder.CBANDS; ++i) gfc.mld_cb_s[i] = 1;

        return 0; // Success
    }


    /**
     * Calculates the Absolute Threshold of Hearing (ATH) in dB SPL for a given frequency.
     * Selects the appropriate formula based on `gfp.ATHtype`.
     * Provides options for different ATH curves (e.g., Painter/Spanias, Bouvigne).
     * Includes a special input `f < -0.3` to request the minimum ATH value.
     *
     * @public
     * @param {number} f - Frequency in Hz. Use f < -0.3 to get the minimum ATH.
     * @param {object} gfp - LAME global flags (used to select ATH type and curve).
     * @returns {number} ATH level in dB SPL.
     */
    ATHformula(f, gfp) {
        let ath;
        switch (gfp.ATHtype) {
            case 0: ath = this._ATHformula_GB(f, 9.0); break;
            case 1: ath = this._ATHformula_GB(f, -1.0); break;
            case 2: ath = this._ATHformula_GB(f, 0.0); break;
            case 3: ath = this._ATHformula_GB(f, 1.0) + 6.0; break;
            case 4: ath = this._ATHformula_GB(f, gfp.ATHcurve); break;
            default: ath = this._ATHformula_GB(f, 0.0); break;
        }
        return ath;
    }

}

// Export the class
export { PsyModel };
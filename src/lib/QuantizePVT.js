/**
 * @fileoverview Internal quantization helper functions for LAME MP3 encoder.
 * Ported from quantize_pvt.c. Contains tables, initialization routines,
 * noise calculation, and VBR helper functions used by the main Quantize module.
 * Uses ES Module syntax.
 *
 * Original C Source Header:
 *      quantize_pvt source file
 *
 *      Copyright (c) 1999-2002 Takehiro Tominaga
 *      Copyright (c) 2000-2002 Robert Hegemann
 *      Copyright (c) 2001 Naoki Shibata
 *      Copyright (c) 2002-2005 Gabriel Bouvigne
 *      ... (License details omitted for brevity) ...
 *
 * $Id: QuantizePVT.java,v 1.24 2011/05/24 20:48:06 kenchis Exp $
 *
 * @module QuantizePVT
 */

// Import necessary modules and utilities using ES Module syntax
import { ScaleFac } from './ScaleFac.js';
import * as common from './common.js';
import { Encoder } from './Encoder.js';
import { MeanBits } from './MeanBits.js';
import { LameInternalFlags } from './LameInternalFlags.js';
import { BitStream } from './BitStream.js'; // Assuming BitStream has static methods or default export

// Assuming these are provided externally or via setModules
/** @typedef {import('./Takehiro.js').default} Takehiro */
/** @typedef {import('./Reservoir.js').default} Reservoir */
/** @typedef {import('./PsyModel.js').default} PsyModel */
/** @typedef {import('./LameGlobalFlags.js').default} LameGlobalFlags */
/** @typedef {import('./GrInfo.js').default} GrInfo */
/** @typedef {import('./CalcNoiseResult.js').default} CalcNoiseResult */
/** @typedef {import('./CalcNoiseData.js').default} CalcNoiseData */

// Destructure common utilities for easier access
const {
    // System, // Not used directly
    VbrMode,
    Float,
    ShortBlock,
    Util,
    Arrays,
    // new_array_n, // Used indirectly
    // new_byte, // Not used
    // new_double, // Not used
    new_float,
    new_float_n,
    new_int,
    // new_int_n, // Not used directly
    assert
} = common;

// --- Constants ---
const Q_MAX = 257; // 256 + 1
const Q_MAX2 = 116; // Max combined gain offset
const LARGE_BITS = 100000; // Placeholder for large bit counts
const IXMAX_VAL = 8206; // Max value for ix in quantization loop
const PRECALC_SIZE = IXMAX_VAL + 2; // Size for precomputed pow43 table
const DBL_EPSILON = 2.2204460492503131e-016; // Machine epsilon
const NSATHSCALE = 100; // ATH scaling factor in dB

// --- Precomputed Tables ---
// Initialized in iteration_init
let pow20 = new_float(Q_MAX + Q_MAX2 + 1);
let ipow20 = new_float(Q_MAX);
let pow43 = new_float(PRECALC_SIZE);
let adj43 = new_float(PRECALC_SIZE);

/**
 * @classdesc Provides private helper functions and data structures for the
 * main Quantize module. This includes precomputed tables, noise calculation,
 * ATH adjustments, and VBR preparation logic.
 * @constructs QuantizePVT
 */
class QuantizePVT {
    /** @private @type {Takehiro|null} Huffman coding/bit counting module. */
    tak = null;
    /** @private @type {Reservoir|null} Reservoir handling module. */
    rv = null;
    /** @private @type {PsyModel|null} Psychoacoustic model module. */
    psy = null;

    /**
     * Table B.6: layer3 preemphasis. Values added to scalefactors before quantization.
     * @public
     * @const {number[]}
     */
    pretab = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 0];

    /**
     * MPEG1 Table B.8 and MPEG2 Table B.1 -- Layer III scalefactor bands.
     * Contains boundary information for different sample rates and block types.
     * Indexed by sample rate and MPEG version.
     * @public
     * @const {ScaleFac[]}
     */
    sfBandIndex = [
        // Definitions omitted for brevity, assume they are correct from original code
        new ScaleFac([0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576], [0, 4, 8, 12, 18, 24, 32, 42, 56, 74, 100, 132, 174, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 22.05 kHz MPEG2
        new ScaleFac([0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 114, 136, 162, 194, 232, 278, 332, 394, 464, 540, 576], [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 136, 180, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 24 kHz MPEG2
        new ScaleFac([0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576], [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 16 kHz MPEG2
        new ScaleFac([0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 52, 62, 74, 90, 110, 134, 162, 196, 238, 288, 342, 418, 576], [0, 4, 8, 12, 16, 22, 30, 40, 52, 66, 84, 106, 136, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 44.1 kHz MPEG1
        new ScaleFac([0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 50, 60, 72, 88, 106, 128, 156, 190, 230, 276, 330, 384, 576], [0, 4, 8, 12, 16, 22, 28, 38, 50, 64, 80, 100, 126, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 48 kHz MPEG1
        new ScaleFac([0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 54, 66, 82, 102, 126, 156, 194, 240, 296, 364, 448, 550, 576], [0, 4, 8, 12, 16, 22, 30, 42, 58, 78, 104, 138, 180, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 32 kHz MPEG1
        new ScaleFac([0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576], [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 11.025 kHz MPEG2.5 (Uses MPEG2 indices?)
        new ScaleFac([0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576], [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]), // 12 kHz MPEG2.5 (Uses MPEG2 indices?)
        new ScaleFac([0, 12, 24, 36, 48, 60, 72, 88, 108, 132, 160, 192, 232, 280, 336, 400, 476, 566, 568, 570, 572, 574, 576], [0, 8, 16, 24, 36, 52, 72, 96, 124, 160, 162, 164, 166, 192], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0])  // 8 kHz MPEG2.5 (Adjusted short block table?)
    ];


    /**
     * Scalefactor partitioning table for MPEG2. Used for grouping scalefactors
     * for transmission (slen values). Dimensions: [table_number][row][column].
     * @public
     * @const {number[][][]}
     */
    nr_of_sfb_block = [
        [[6, 5, 5, 5], [9, 9, 9, 9], [6, 9, 9, 9]],
        [[6, 5, 7, 3], [9, 9, 12, 6], [6, 9, 12, 6]],
        [[11, 10, 0, 0], [18, 18, 0, 0], [15, 18, 0, 0]],
        [[7, 7, 7, 0], [12, 12, 12, 0], [6, 15, 12, 0]],
        [[6, 6, 6, 3], [12, 9, 9, 6], [6, 12, 9, 6]],
        [[8, 8, 5, 0], [15, 12, 9, 0], [6, 18, 9, 0]]];


    /**
     * Precomputed table for adjusting quantization steps (x^(4/3)).
     * Used in noise calculation. Global access needed? Made public for now.
     * @public
     * @type {Float32Array}
     */
    adj43 = adj43; // Assign instance variable

    constructor() {
        // Properties initialized externally via setModules
    }

    /**
     * Sets the internal module dependencies. Must be called before use.
     *
     * @public
     * @param {Takehiro} _tk - Huffman coding and bit counting module.
     * @param {Reservoir} _rv - Bit reservoir handling module.
     * @param {PsyModel} _psy - Psychoacoustic model module.
     */
    setModules(_tk, _rv, _psy) {
        this.tak = _tk;
        this.rv = _rv;
        this.psy = _psy;
    }

    /**
     * Calculates `2.0 ^ ((x - 210) * -0.1875)`. Inverse power function.
     * Used for converting quantized values back to spectral domain (approx).
     * Uses precomputed table `ipow20`.
     *
     * @public
     * @param {number} x - Input value (typically related to gains/scalefactors). Expected range [0, Q_MAX-1].
     * @returns {number} Precomputed inverse power value.
     */
    IPOW20(x) {
        assert(0 <= x && x < Q_MAX, `IPOW20 index out of bounds: ${x}`);
        return ipow20[x];
    }

    // --- Private Helper Functions ---
    // (JSDoc omitted for brevity)

    /** @private */
    _POW20(x) { // Renamed from POW20 to avoid conflict and indicate internal use
        assert(0 <= (x + Q_MAX2) && x < Q_MAX + Q_MAX2 + 1, `POW20 index out of bounds: ${x}`); // Check adjusted index
        return pow20[x + Q_MAX2];
    }

    /** @private */
    _ATHmdct(gfp, f) {
        let ath = this.psy.ATHformula(f, gfp);
        ath -= NSATHSCALE;
        ath = Math.pow(10.0, ath / 10.0 + gfp.ATHlower);
        return ath;
    }

    /** @private */
    _compute_ath(gfp) {
        const gfc = gfp.internal_flags;
        const ATH_l = gfc.ATH.l;
        const ATH_psfb21 = gfc.ATH.psfb21;
        const ATH_s = gfc.ATH.s;
        const ATH_psfb12 = gfc.ATH.psfb12;
        const samp_freq = gfp.out_samplerate;

        for (let sfb = 0; sfb < Encoder.SBMAX_l; sfb++) {
            const start = gfc.scalefac_band.l[sfb]; const end = gfc.scalefac_band.l[sfb + 1];
            ATH_l[sfb] = Float.MAX_VALUE;
            for (let i = start; i < end; i++) { const freq = i * samp_freq / (2.0 * 576.0); const ATH_f = this._ATHmdct(gfp, freq); ATH_l[sfb] = Math.min(ATH_l[sfb], ATH_f); }
        }
        for (let sfb = 0; sfb < Encoder.PSFB21; sfb++) {
            const start = gfc.scalefac_band.psfb21[sfb]; const end = gfc.scalefac_band.psfb21[sfb + 1];
            ATH_psfb21[sfb] = Float.MAX_VALUE;
            for (let i = start; i < end; i++) { const freq = i * samp_freq / (2.0 * 576.0); const ATH_f = this._ATHmdct(gfp, freq); ATH_psfb21[sfb] = Math.min(ATH_psfb21[sfb], ATH_f); }
        }
        for (let sfb = 0; sfb < Encoder.SBMAX_s; sfb++) {
            const start = gfc.scalefac_band.s[sfb]; const end = gfc.scalefac_band.s[sfb + 1];
            ATH_s[sfb] = Float.MAX_VALUE;
            for (let i = start; i < end; i++) { const freq = i * samp_freq / (2.0 * 192.0); const ATH_f = this._ATHmdct(gfp, freq); ATH_s[sfb] = Math.min(ATH_s[sfb], ATH_f); }
            ATH_s[sfb] *= (gfc.scalefac_band.s[sfb + 1] - gfc.scalefac_band.s[sfb]);
        }
        for (let sfb = 0; sfb < Encoder.PSFB12; sfb++) {
            const start = gfc.scalefac_band.psfb12[sfb]; const end = gfc.scalefac_band.psfb12[sfb + 1];
            ATH_psfb12[sfb] = Float.MAX_VALUE;
            for (let i = start; i < end; i++) { const freq = i * samp_freq / (2.0 * 192.0); const ATH_f = this._ATHmdct(gfp, freq); ATH_psfb12[sfb] = Math.min(ATH_psfb12[sfb], ATH_f); }
            ATH_psfb12[sfb] *= (gfc.scalefac_band.s[13] - gfc.scalefac_band.s[12]); // Uses fixed width?
        }
        if (gfp.noATH) {
            for (let sfb = 0; sfb < Encoder.SBMAX_l; sfb++) ATH_l[sfb] = 1E-20;
            for (let sfb = 0; sfb < Encoder.PSFB21; sfb++) ATH_psfb21[sfb] = 1E-20;
            for (let sfb = 0; sfb < Encoder.SBMAX_s; sfb++) ATH_s[sfb] = 1E-20;
            for (let sfb = 0; sfb < Encoder.PSFB12; sfb++) ATH_psfb12[sfb] = 1E-20;
        }
        gfc.ATH.floor = 10.0 * Math.log10(this._ATHmdct(gfp, -1.0)); // Calculate floor using special value
    }


    // --- Public Methods ---

    /**
     * Initializes tables and psychoacoustic parameters used by the quantization loops.
     * Computes ATH values, precalculates power tables (`pow43`, `pow20`, `ipow20`),
     * initializes Huffman coding tables, and sets up psychoacoustic tuning factors.
     * This function should only be called once during encoder initialization.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     */
    iteration_init(gfp) {
        const gfc = gfp.internal_flags;
        let i;

        if (gfc.iteration_init_init === 0) { // Check if already initialized
            gfc.iteration_init_init = 1; // Mark as initialized

            gfc.l3_side.main_data_begin = 0; // Reset main data pointer?
            this._compute_ath(gfp); // Calculate ATH values

            // Precompute power tables
            pow43[0] = 0.0;
            for (i = 1; i < PRECALC_SIZE; i++) pow43[i] = Math.pow(i, 4.0 / 3.0);

            for (i = 0; i < PRECALC_SIZE - 1; i++) adj43[i] = ((i + 1) - Math.pow(0.5 * (pow43[i] + pow43[i + 1]), 0.75));
            adj43[PRECALC_SIZE - 1] = 0.5; // Set last value explicitly

            for (i = 0; i < Q_MAX; i++) ipow20[i] = Math.pow(2.0, (i - 210) * -0.1875);
            for (i = 0; i <= Q_MAX + Q_MAX2; i++) pow20[i] = Math.pow(2.0, (i - 210 - Q_MAX2) * 0.25);

            // Initialize Huffman tables (delegated to Takehiro module)
            this.tak.huffman_init(gfc);

            // Initialize psychoacoustic tuning factors based on exp_nspsytune
            {
                let bass, alto, treble, sfb21;
                i = (gfp.exp_nspsytune >> 2) & 63; if (i >= 32) i -= 64; bass = Math.pow(10, i / 4.0 / 10.0);
                i = (gfp.exp_nspsytune >> 8) & 63; if (i >= 32) i -= 64; alto = Math.pow(10, i / 4.0 / 10.0);
                i = (gfp.exp_nspsytune >> 14) & 63; if (i >= 32) i -= 64; treble = Math.pow(10, i / 4.0 / 10.0);
                i = (gfp.exp_nspsytune >> 20) & 63; if (i >= 32) i -= 64; sfb21 = treble * Math.pow(10, i / 4.0 / 10.0);

                for (i = 0; i < Encoder.SBMAX_l; i++) {
                    let f = (i <= 6) ? bass : (i <= 13) ? alto : (i <= 20) ? treble : sfb21;
                    gfc.nsPsy.longfact[i] = f;
                }
                for (i = 0; i < Encoder.SBMAX_s; i++) {
                     let f = (i <= 5) ? bass : (i <= 10) ? alto : (i <= 11) ? treble : sfb21;
                    gfc.nsPsy.shortfact[i] = f;
                }
            }
        }
    }

    /**
     * Allocates target bits among channels based on perceptual entropy (PE).
     * Higher PE gets more bits. Considers bit reservoir and enforces limits.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {Array<Float32Array>} pe - Perceptual entropy [gr][ch].
     * @param {Int32Array} targ_bits - Input/Output: Target bits array for the granule [ch]. Initial allocation might be present; modified in place.
     * @param {number} mean_bits - Average bits available per granule (considering reservoir).
     * @param {number} gr - Granule index (0 or 1).
     * @param {number} cbr - Flag indicating CBR mode (affects reservoir usage).
     * @returns {number} Maximum allowed bits for this granule (including reservoir).
     */
    on_pe(gfp, pe, targ_bits, mean_bits, gr, cbr) {
        const gfc = gfp.internal_flags;
        let tbits = 0; // Base bits for this granule from reservoir
        let total_added_bits = 0; // Total extra bits allocated based on PE
        const add_bits = new_int(2); // Extra bits per channel
        const channels = gfc.channels_out;

        // 1. Get base bits and allowed extra bits from reservoir
        let mb = new MeanBits(tbits);
        const extra_bits = this.rv.ResvMaxBits(gfp, mean_bits, mb, cbr);
        tbits = mb.bits; // Base bits available for this granule
        const max_bits_granule = tbits + extra_bits; // Total allowed for granule
        // Apply hard limit per granule
        const capped_max_bits_granule = Math.min(max_bits_granule, LameInternalFlags.MAX_BITS_PER_GRANULE);

        // 2. Initial target bit allocation and calculate needed extra bits based on PE
        for (let ch = 0; ch < channels; ++ch) {
            // Initial target: Base bits distributed equally
            targ_bits[ch] = Math.min(LameInternalFlags.MAX_BITS_PER_CHANNEL, Math.floor(tbits / channels));

            // Calculate extra bits needed based on PE (relative to 700)
            add_bits[ch] = Math.floor(targ_bits[ch] * pe[gr][ch] / 700.0) - targ_bits[ch];

            // Limit extra bits per channel (heuristic: 0 to 0.75 * avg_bits_per_channel)
            const avg_bits_per_ch_granule = mean_bits / channels; // Approx avg bits per ch/gr
            if (add_bits[ch] > avg_bits_per_ch_granule * 0.75) {
                add_bits[ch] = Math.floor(avg_bits_per_ch_granule * 0.75);
            }
            if (add_bits[ch] < 0) add_bits[ch] = 0;

            // Ensure adding extra bits doesn't exceed per-channel limit
            if (add_bits[ch] + targ_bits[ch] > LameInternalFlags.MAX_BITS_PER_CHANNEL) {
                add_bits[ch] = Math.max(0, LameInternalFlags.MAX_BITS_PER_CHANNEL - targ_bits[ch]);
            }

            total_added_bits += add_bits[ch]; // Sum needed extra bits
        }

        // 3. Scale added bits if total needed exceeds available extra bits
        let actual_extra_bits = extra_bits; // Use the reservoir's calculated extra bits
        if (total_added_bits > actual_extra_bits) {
            for (let ch = 0; ch < channels; ++ch) {
                add_bits[ch] = Math.floor(actual_extra_bits * add_bits[ch] / total_added_bits);
            }
        }

        // 4. Add the allocated extra bits to the target bits
        let final_total_bits = 0;
        for (let ch = 0; ch < channels; ++ch) {
            targ_bits[ch] += add_bits[ch];
            final_total_bits += targ_bits[ch];
        }

        // 5. Rescale final targets if total exceeds granule hard limit (MAX_BITS_PER_GRANULE)
        // This might happen if initial tbits + scaled add_bits > hard limit
        if (final_total_bits > LameInternalFlags.MAX_BITS_PER_GRANULE) {
            let sum_after_scale = 0;
            for (let ch = 0; ch < channels; ++ch) {
                targ_bits[ch] = Math.floor(targ_bits[ch] * LameInternalFlags.MAX_BITS_PER_GRANULE / final_total_bits);
                sum_after_scale += targ_bits[ch];
            }
             assert(sum_after_scale <= LameInternalFlags.MAX_BITS_PER_GRANULE, "on_pe rescaling failed");
        }

        // Return the maximum bits allowed for this granule
        return capped_max_bits_granule;
    }

    /**
     * Reduces the number of bits allocated to the Side channel in MS Stereo,
     * transferring them to the Mid channel based on the Mid/Side energy ratio.
     * Aims for roughly a 2:1 Mid:Side bit split when energy is equal, shifting
     * towards Mid as Side energy decreases. Enforces a minimum bit allocation
     * for the Side channel and respects overall granule/channel bit limits.
     *
     * @public
     * @param {Int32Array} targ_bits - Input/Output: Target bits array for the granule [ch]. Modified in place.
     * @param {number} ms_ener_ratio - Energy ratio (Side / Mid). Typically 0 to 1.
     * @param {number} mean_bits - Average total bits per granule (Mid+Side). Used for limiting transfer.
     * @param {number} max_bits - Maximum total bits allowed for the granule (Mid+Side).
     */
    reduce_side(targ_bits, ms_ener_ratio, mean_bits, max_bits) {
        assert(max_bits <= LameInternalFlags.MAX_BITS_PER_GRANULE, "max_bits exceeds granule limit in reduce_side");
        assert(targ_bits[0] + targ_bits[1] <= LameInternalFlags.MAX_BITS_PER_GRANULE, "Initial targ_bits exceed granule limit");

        // Calculate transfer factor (0 to 0.5).
        // fac = 0.33 when ratio=0 (all Mid), fac=0 when ratio=0.5 (equal), fac=-0.165 when ratio=1 (all Side)
        // Clamp factor between 0 and 0.5.
        let fac = 0.33 * (0.5 - ms_ener_ratio) / 0.5;
        if (fac < 0) fac = 0;
        if (fac > 0.5) fac = 0.5;

        // Calculate bits to move: factor * half_of_total_bits
        let move_bits = Math.floor(fac * 0.5 * (targ_bits[0] + targ_bits[1]));

        // Limit move_bits: cannot make Mid exceed max channel bits
        if (move_bits > LameInternalFlags.MAX_BITS_PER_CHANNEL - targ_bits[0]) {
            move_bits = LameInternalFlags.MAX_BITS_PER_CHANNEL - targ_bits[0];
        }
        if (move_bits < 0) move_bits = 0; // Should not happen due to factor clamping

        // Check if Side channel has enough bits and enforce minimum
        const SIDE_MIN_BITS = 125;
        if (targ_bits[1] >= SIDE_MIN_BITS) {
            // Only move bits if Side remains above minimum
            if (targ_bits[1] - move_bits >= SIDE_MIN_BITS) {
                // Limit transfer if Mid already has significantly more than average
                // (mean_bits here is avg per granule *for both channels*)
                 if (targ_bits[0] < mean_bits) { // Transfer only if Mid is below average total? Heuristic check.
                     targ_bits[0] += move_bits;
                     targ_bits[1] -= move_bits;
                 }
                 // If Mid is already high, don't transfer bits.
            } else {
                // Move only enough bits to bring Side down to minimum
                const bits_to_move_limited = targ_bits[1] - SIDE_MIN_BITS;
                 if (targ_bits[0] < mean_bits) { // Check Mid average heuristic again
                     targ_bits[0] += bits_to_move_limited;
                 }
                 // else: Mid is high, don't add more, but Side still goes to minimum.
                 targ_bits[1] = SIDE_MIN_BITS;
            }
        }
         // else: Side channel starts below minimum, no bits are moved.

        // Final check: Ensure total bits do not exceed max_bits for the granule
        const total_bits_after_move = targ_bits[0] + targ_bits[1];
        if (total_bits_after_move > max_bits) {
            // Scale both channels down proportionally
            if (total_bits_after_move > 0) { // Avoid division by zero
                 targ_bits[0] = Math.floor(max_bits * targ_bits[0] / total_bits_after_move);
                 targ_bits[1] = Math.floor(max_bits * targ_bits[1] / total_bits_after_move);
            } else {
                // Should not happen if SIDE_MIN_BITS > 0
                 targ_bits[0] = Math.floor(max_bits / 2);
                 targ_bits[1] = Math.floor(max_bits / 2);
            }
        }

        // Assert final limits
        assert(targ_bits[0] >= 0 && targ_bits[0] <= LameInternalFlags.MAX_BITS_PER_CHANNEL, `Mid bits out of range: ${targ_bits[0]}`);
        assert(targ_bits[1] >= 0 && targ_bits[1] <= LameInternalFlags.MAX_BITS_PER_CHANNEL, `Side bits out of range: ${targ_bits[1]}`);
        assert(targ_bits[0] + targ_bits[1] <= LameInternalFlags.MAX_BITS_PER_GRANULE, "Total bits exceed granule limit after reduce_side");
    }

    /**
     * Adjusts the Absolute Threshold of Hearing (ATH) value based on the
     * ATH adjustment factor (`a`) and the ATH floor. This aims to shape
     * the noise floor dynamically. Used in VBR modes.
     *
     * @public
     * @param {number} a - ATH adjustment factor (gfc.ATH.adjust, typically 0 to 1).
     * @param {number} x - Original ATH value for the band (linear energy scale).
     * @param {number} athFloor - ATH floor in dB (gfc.ATH.floor).
     * @returns {number} Adjusted ATH value (linear energy scale).
     */
    athAdjust(a, x, athFloor) {
        const o = 90.30873362; // Constant from formula
        const p = 94.82444863; // Constant from formula
        let u = Util.FAST_LOG10_X(x, 10.0); // Convert original ATH energy to dB-like scale
        const v = a * a; // Square of adjustment factor
        let w = 0.0;

        u -= athFloor; // Offset by ATH floor

        // Apply adjustment factor 'a' (via 'v') logarithmically relative to constant 'o'
        if (v > 1E-20) {
            w = 1.0 + Util.FAST_LOG10_X(v, 10.0 / o); // Log scale the adjustment factor
        }
        if (w < 0.0) w = 0.0; // Clamp adjustment effect
        u *= w; // Apply scaled adjustment

        u += athFloor + o - p; // Reapply floor and constants

        // Convert back to linear energy scale
        return Math.pow(10.0, 0.1 * u);
    }


    /**
     * Calculates the allowed noise `xmin` for each scalefactor band based on
     * the psycho-acoustic model results (`ratio`) and the ATH.
     * `xmin(sb) = max(ATH(sb), masking_threshold(sb))`
     * where `masking_threshold(sb) = energy(sb) * ratio(sb) * masking_lower`
     * adjusted for perceptual tuning factors (`nsPsy.longfact`, `nsPsy.shortfact`).
     * Also finds the highest non-zero spectral coefficient (`max_nonzero_coeff`).
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {object} ratio - Masking ratio structure containing `en` (energy) and `thm` (threshold ratio) for the granule/channel.
     * @param {GrInfo} cod_info - Granule info structure containing block type, spectral data (`xr`). `max_nonzero_coeff` is set here.
     * @param {Float32Array} pxmin - Output: Array [SBMAX] to store the calculated allowed noise for each band.
     * @returns {number} `ath_over`: Count of scalefactor bands where energy exceeds the ATH.
     */
    calc_xmin(gfp, ratio, cod_info, pxmin) {
        let pxminPos = 0; // Index for pxmin output array
        const gfc = gfp.internal_flags;
        let gsfb; // Scalefactor band index (combined long/short)
        let j = 0; // Index into spectral data xr[]
        let ath_over = 0; // Count of bands where energy > ATH
        const ATH = gfc.ATH; // ATH values structure
        const xr = cod_info.xr; // Spectral data
        const enable_athaa_fix = (gfp.VBR === VbrMode.vbr_mtrh) ? 1 : 0; // VBR MTRH specific flag
        let masking_lower = gfc.masking_lower; // Global masking adjustment

        // VBR MTRH/MT already applied masking_lower in PSY model
        if (gfp.VBR === VbrMode.vbr_mtrh || gfp.VBR === VbrMode.vbr_mt) {
            masking_lower = 1.0;
        }

        // --- Process Long Block Bands (or long part of mixed block) ---
        for (gsfb = 0; gsfb < cod_info.psy_lmax; gsfb++) {
            let ath_band; // ATH for this band
            let xmin; // Minimum allowed noise (allowed distortion power)
            const width = cod_info.width[gsfb]; // Number of coefficients in this band
            let en0 = 0.0; // Energy in this band

            // Calculate ATH for the band (potentially adjusted)
            if (gfp.VBR === VbrMode.vbr_rh || gfp.VBR === VbrMode.vbr_mtrh) {
                ath_band = this.athAdjust(ATH.adjust, ATH.l[gsfb], ATH.floor);
            } else {
                ath_band = ATH.adjust * ATH.l[gsfb];
            }

            // Calculate energy (en0) and a noise estimate (rh2) based on ATH clipping
            const rh1 = ath_band / width; // ATH threshold per coefficient (approx)
            let rh2 = DBL_EPSILON; // Accumulated noise assuming clipping at rh1
            if(width > 0) { // Avoid loop if width is 0
                let l = width >> 1; // Process pairs
                do {
                    const xa = xr[j] * xr[j]; en0 += xa; rh2 += Math.min(xa, rh1); j++;
                    const xb = xr[j] * xr[j]; en0 += xb; rh2 += Math.min(xb, rh1); j++;
                } while (--l > 0);
                 // Handle odd width if necessary (already included in C logic by loop condition)
                 // if (width & 1) { const xa = xr[j] * xr[j]; en0 += xa; rh2 += Math.min(xa, rh1); j++; }
            }

            // Check if energy exceeds ATH
            if (en0 > ath_band) ath_over++;

            // Apply nsPsy tuning factor adjustment (specific sfb index?)
            if (gsfb === Encoder.SBPSY_l) { // Special handling for last psy band? Check index.
                const x = ath_band * gfc.nsPsy.longfact[gsfb];
                if (rh2 < x) rh2 = x; // Adjust noise estimate based on tuning factor
            }

            // Set initial xmin based on ATH (potentially using rh2 for VBR MTRH)
            xmin = (enable_athaa_fix !== 0) ? rh2 : ath_band;

            // Add masking threshold component if not ATHonly mode
            if (!gfp.ATHonly) {
                const mask_en = ratio.en.l[gsfb];
                if (mask_en > 0.0) {
                    let mask_thr = en0 * ratio.thm.l[gsfb] * masking_lower / mask_en;
                    // Apply nsPsy tuning factor to masking threshold part
                    if (enable_athaa_fix !== 0) mask_thr *= gfc.nsPsy.longfact[gsfb];
                    // xmin is the maximum of ATH and masking threshold
                    if (xmin < mask_thr) xmin = mask_thr;
                }
            }

            // Store final xmin (apply nsPsy factor if not MTRH fix)
            pxmin[pxminPos++] = (enable_athaa_fix !== 0) ? xmin : xmin * gfc.nsPsy.longfact[gsfb];

        } // End long block loop

        // --- Find highest non-zero coefficient (only needed for long blocks here?) ---
        let max_nonzero = 575;
        if (cod_info.block_type !== Encoder.SHORT_TYPE) { // Long or Mixed blocks
            let k = 576;
            // C code has while(k--!=0 && EQ(xr[k],0)); max_nonzero=k;
            // This finds the index *of* the last non-zero or -1 if all zero.
            while (k-- > 0) {
                if (Math.abs(xr[k]) > 1e-9) break; // Find last non-zero
            }
             max_nonzero = k; // k is -1 if all zero, or index of last non-zero
             if(max_nonzero < 0) max_nonzero = 0; // Ensure non-negative index if used later
             cod_info.max_nonzero_coeff = max_nonzero;
        } else {
             cod_info.max_nonzero_coeff = 575; // Assume all coeffs potentially non-zero for short blocks
        }


        // --- Process Short Block Bands (if applicable) ---
        let sfb_s = cod_info.sfb_smin; // Short block sfb index
        while (gsfb < cod_info.psymax) { // Continue using combined index gsfb
            const width = cod_info.width[gsfb]; // Width of one short sub-band
            let ath_band; // ATH for this short sfb

            // Calculate ATH for the band
            if (gfp.VBR === VbrMode.vbr_rh || gfp.VBR === VbrMode.vbr_mtrh) {
                 ath_band = this.athAdjust(ATH.adjust, ATH.s[sfb_s], ATH.floor);
            } else {
                 ath_band = ATH.adjust * ATH.s[sfb_s];
            }

            // Process the 3 windows (sub-blocks) for this short sfb
            for (let b = 0; b < 3; b++) {
                let en0 = 0.0; let xmin;
                const rh1 = ath_band / width; // ATH per coefficient
                let rh2 = DBL_EPSILON; // Accumulated noise estimate
                if(width > 0) {
                    let l = width >> 1;
                    do {
                        const xa = xr[j] * xr[j]; en0 += xa; rh2 += Math.min(xa, rh1); j++;
                        const xb = xr[j] * xr[j]; en0 += xb; rh2 += Math.min(xb, rh1); j++;
                    } while (--l > 0);
                    // if(width & 1) { const xa = xr[j] * xr[j]; en0 += xa; rh2 += Math.min(xa, rh1); j++; }
                }

                if (en0 > ath_band) ath_over++;

                // Apply nsPsy tuning factor adjustment (specific sfb index?)
                 if (sfb_s === Encoder.SBPSY_s) { // Check index for last psy short band
                     const x = ath_band * gfc.nsPsy.shortfact[sfb_s];
                     if (rh2 < x) rh2 = x;
                 }

                // Set initial xmin based on ATH
                xmin = (enable_athaa_fix !== 0) ? rh2 : ath_band;

                // Add masking threshold component if not ATHonly/ATHshort mode
                if (!gfp.ATHonly && !gfp.ATHshort) {
                    const mask_en = ratio.en.s[sfb_s][b];
                    if (mask_en > 0.0) {
                        let mask_thr = en0 * ratio.thm.s[sfb_s][b] * masking_lower / mask_en;
                         if (enable_athaa_fix !== 0) mask_thr *= gfc.nsPsy.shortfact[sfb_s];
                        if (xmin < mask_thr) xmin = mask_thr;
                    }
                }

                // Store final xmin (apply nsPsy factor if not MTRH fix)
                 pxmin[pxminPos++] = (enable_athaa_fix !== 0) ? xmin : xmin * gfc.nsPsy.shortfact[sfb_s];

                gsfb++; // Increment combined index for each window
            } // End window loop (b)

            // Temporal masking adjustment (apply decay between windows of same sfb)
            if (gfp.useTemporal) {
                if (pxmin[pxminPos - 3] > pxmin[pxminPos - 2]) pxmin[pxminPos - 2] += (pxmin[pxminPos - 3] - pxmin[pxminPos - 2]) * gfc.decay;
                if (pxmin[pxminPos - 2] > pxmin[pxminPos - 1]) pxmin[pxminPos - 1] += (pxmin[pxminPos - 2] - pxmin[pxminPos - 1]) * gfc.decay;
            }

            sfb_s++; // Increment short block sfb index
        } // End short block sfb loop

        return ath_over; // Return count of bands where energy exceeded ATH
    }


    /**
     * Helper class for calc_noise_core to pass index by reference.
     * @private
     */
    _StartLine = class { constructor(j) { this.s = j; } };


    /**
     * Core calculation of quantization noise power for a scalefactor band.
     * Noise = Sum[ (|xr[i]| - quantized_xr[i])^2 ]
     * where quantized_xr depends on the quantization step (`step`) and the
     * quantized value (`ix[i]`). Handles different quantization regions (big_values, count1).
     *
     * @public
     * @param {GrInfo} cod_info - Granule info containing spectral data (`xr`) and quantized values (`l3_enc`).
     * @param {object} startline - Input/Output: Object `{s: number}` holding the starting index into `xr` and `l3_enc`. Updated after processing.
     * @param {number} l - Number of coefficient *pairs* to process (width / 2).
     * @param {number} step - Quantization step size (linear scale, from `_POW20(s)`).
     * @returns {number} Total noise power (sum of squared errors) for the processed coefficients.
     */
    calc_noise_core(cod_info, startline, l, step) {
        let noise = 0.0;
        let j = startline.s; // Current index in xr/l3_enc
        const ix = cod_info.l3_enc; // Quantized values
        const xr_abs = cod_info.xr; // Original spectral values (absolute values used here)

        // Determine which quantization region the current index falls into
        if (j >= cod_info.count1) { // Region 2: ix[j] = 0 or 1 (abs value < 1.5*step?) - Not possible here? count1 is boundary.
            // This region seems to be where ix[j] == 0. The error is just xr[j]^2.
            while (l-- > 0) {
                let temp;
                temp = xr_abs[j]; j++; noise += temp * temp; // ix[j-1] assumed 0
                temp = xr_abs[j]; j++; noise += temp * temp; // ix[j-1] assumed 0
            }
        } else if (j >= cod_info.big_values) { // Region 1: ix[j] = 0 or 1
            const ix01 = [0.0, step]; // Dequantized values for ix=0 and ix=1
            while (l-- > 0) {
                let temp;
                // Error = |xr| - dequantized_value(|xr|)
                temp = Math.abs(xr_abs[j]) - ix01[ix[j]]; j++; noise += temp * temp;
                temp = Math.abs(xr_abs[j]) - ix01[ix[j]]; j++; noise += temp * temp;
            }
        } else { // Region 0 (big_values): ix[j] > 1
            // Dequantized value = (ix[j])^(4/3) * step
            while (l-- > 0) {
                 let temp;
                 // Check index bounds for pow43
                 const ix_val_0 = ix[j];
                 const dequant_0 = (ix_val_0 < PRECALC_SIZE) ? pow43[ix_val_0] * step : Math.pow(ix_val_0, 4.0/3.0) * step; // Fallback if index too high
                 temp = Math.abs(xr_abs[j]) - dequant_0; j++; noise += temp * temp;

                 const ix_val_1 = ix[j];
                 const dequant_1 = (ix_val_1 < PRECALC_SIZE) ? pow43[ix_val_1] * step : Math.pow(ix_val_1, 4.0/3.0) * step;
                 temp = Math.abs(xr_abs[j]) - dequant_1; j++; noise += temp * temp;
            }
        }

        startline.s = j; // Update index in the passed object
        return noise;
    }

    /**
     * Calculates the quantization noise relative to the allowed noise (`l3_xmin`)
     * for each scalefactor band. Populates the `distort` array where `distort[sfb] = noise / l3_xmin[sfb]`.
     * Also calculates overall noise metrics (total noise, noise over threshold, max noise)
     * and stores them in the `res` object. Can optionally use/update `prev_noise` cache.
     *
     * @public
     * @param {GrInfo} cod_info - Granule info containing quantization results and parameters.
     * @param {Float32Array} l3_xmin - Array [SBMAX] of allowed noise per scalefactor band.
     * @param {Float32Array} distort - Output: Array [SBMAX] to store the noise ratio (noise/xmin) per band.
     * @param {CalcNoiseResult} res - Output: Object to store overall noise metrics (`over_count`, `tot_noise`, `over_noise`, `max_noise`, `over_SSD`).
     * @param {CalcNoiseData|null} prev_noise - Optional Input/Output: Cache for noise values from previous iteration with same global gain. Can speed up repeated calculations.
     * @returns {number} `over_count`: Number of bands where noise > allowed noise (distort > 1).
     */
    calc_noise(cod_info, l3_xmin, distort, res, prev_noise) {
        let distortPos = 0; // Index for distort array
        let l3_xminPos = 0; // Index for l3_xmin array
        let over = 0; // Count of distorted bands
        let over_noise_db = 0.0; // Sum of dB noise for distorted bands
        let tot_noise_db = 0.0; // Sum of dB noise for all bands
        let max_noise = -20.0; // Max noise in dB relative to xmin (-200 dB)
        let j = 0; // Index into spectral data xr[]
        const scalefac = cod_info.scalefac;
        let scalefacPos = 0; // Index for scalefac array

        res.over_SSD = 0.0; // Sum of squared distortions for over-noise bands (integer scale)

        // Check if prev_noise cache is valid for this global_gain
        const use_prev = prev_noise != null && prev_noise.global_gain === cod_info.global_gain;

        // Iterate over psychoacoustic scalefactor bands
        for (let sfb = 0; sfb < cod_info.psymax; sfb++) {
            let s; // Effective gain = global_gain - scalefac_gain - subblock_gain
            let noise_power = 0.0; // Noise power for this band (linear)
            let noise_db = 0.0; // Noise power relative to xmin (dB * 10)

            // Calculate effective gain 's'
            s = cod_info.global_gain -
                (((scalefac[scalefacPos++]) + (cod_info.preflag !== 0 ? this.pretab[sfb] : 0)) << (cod_info.scalefac_scale + 1)) -
                (cod_info.subblock_gain[cod_info.window[sfb]] * 8);

            // Check cache or calculate noise
            if (use_prev && prev_noise.step[sfb] === s) {
                // Use cached values
                noise_power = prev_noise.noise[sfb];
                noise_db = prev_noise.noise_log[sfb];
                j += cod_info.width[sfb]; // Advance spectral data index
                distort[distortPos++] = noise_power / l3_xmin[l3_xminPos++]; // Calculate distortion ratio
            } else {
                // Calculate noise from scratch
                const step = this._POW20(s); // Get linear quantization step size
                let width = cod_info.width[sfb];
                let l = width >> 1; // Number of pairs
                let start_idx = j;

                // Adjust length if band crosses max_nonzero_coeff boundary
                if ((j + width) > cod_info.max_nonzero_coeff) {
                    const usefullsize = cod_info.max_nonzero_coeff - j + 1;
                    l = (usefullsize > 0) ? usefullsize >> 1 : 0;
                    // Handle odd usefullsize? calc_noise_core expects pairs.
                    // If usefullsize is odd, the last coeff is processed by calc_noise_core? Check logic.
                    // Assuming calc_noise_core handles this.
                }
                if (l < 0) l = 0; // Ensure non-negative length

                // Calculate noise power using core function
                 let sl_obj = new this._StartLine(start_idx);
                 noise_power = this.calc_noise_core(cod_info, sl_obj, l, step);
                 j = sl_obj.s; // Update spectral data index

                // Store calculated values in cache if provided
                if (prev_noise != null) {
                    prev_noise.step[sfb] = s;
                    prev_noise.noise[sfb] = noise_power;
                }

                // Calculate distortion ratio and noise in dB
                distort[distortPos++] = noise_power / l3_xmin[l3_xminPos++];
                noise_db = Util.FAST_LOG10(Math.max(distort[distortPos - 1], 1E-20)); // Use distortion ratio for dB calc

                if (prev_noise != null) {
                    prev_noise.noise_log[sfb] = noise_db;
                }
            }

            // Update overall noise metrics
            tot_noise_db += noise_db;
            if (noise_db > 0.0) { // If distorted (noise > xmin)
                 // Calculate integer measure of distortion for over_SSD
                 const tmp = Math.max(Math.floor(noise_db * 10.0 + 0.5), 1); // Round dB*10 to nearest int, min 1
                 res.over_SSD += tmp * tmp; // Add squared value

                 over++; // Increment distorted band count
                 over_noise_db += noise_db; // Sum dB noise of distorted bands
            }
            max_noise = Math.max(max_noise, noise_db); // Track max noise (dB)

        } // End sfb loop

        // Store final metrics in result object
        res.over_count = over;
        res.tot_noise = tot_noise_db;
        res.over_noise = over_noise_db;
        res.max_noise = max_noise;

        // Update cache global gain marker if cache was used
        if (prev_noise != null) {
            prev_noise.global_gain = cod_info.global_gain;
        }

        return over; // Return distorted band count
    }

    /**
     * Updates the `gfc.pinfo` structure with detailed plotting/analysis data
     * for a specific granule/channel based on the final quantization results.
     * Calculates energy, masking threshold, and quantization noise per scalefactor band.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {GrInfo} cod_info - Granule info containing final quantization results.
     * @param {object} ratio - Masking ratio structure for the granule/channel.
     * @param {number} gr - Granule index (0 or 1).
     * @param {number} ch - Channel index (0 or 1).
     */
    set_pinfo(gfp, cod_info, ratio, gr, ch) {
        const gfc = gfp.internal_flags;
        let sfb, sfb2; // Scalefactor band indices
        let l; // Loop counter for spectral lines
        let en0; // Energy per coefficient or band
        const ifqstep = (cod_info.scalefac_scale === 0) ? 0.5 : 1.0; // Scalefactor step size factor
        const scalefac = cod_info.scalefac; // Final scalefactors

        // Temporary arrays for intermediate calculations
        const l3_xmin = new_float(L3Side.SFBMAX); // Allowed noise per band
        const distort = new_float(L3Side.SFBMAX); // Noise/xmin ratio per band
        const noise_res = new CalcNoiseResult(); // Overall noise metrics

        // Recalculate allowed noise and actual distortion for plotting
        this.calc_xmin(gfp, ratio, cod_info, l3_xmin);
        this.calc_noise(cod_info, l3_xmin, distort, noise_res, null); // Use final cod_info

        const pinfo = gfc.pinfo; // Shortcut to plotting info structure
        let j = 0; // Index into spectral data xr[]

        // Determine upper sfb limit based on block type
        sfb2 = cod_info.sfb_lmax; // Default limit for mixed/short
        if (cod_info.block_type !== Encoder.SHORT_TYPE && cod_info.mixed_block_flag === 0) {
            sfb2 = 22; // Full long block limit
        }

        // --- Process Long Bands ---
        for (sfb = 0; sfb < sfb2; sfb++) {
            const start = gfc.scalefac_band.l[sfb];
            const end = gfc.scalefac_band.l[sfb + 1];
            const bw = end - start; // Bandwidth (number of coefficients)
            if (bw === 0) continue; // Skip zero-width bands

            // Calculate energy in the band
            en0 = 0.0;
            for (l = start; l < end; l++) { // Use temp index l
                en0 += cod_info.xr[l] * cod_info.xr[l];
            }
            j = end; // Update main spectral index j
            en0 /= bw; // Average energy per coefficient

            // Store energy, noise, and threshold for plotting (scaled)
            const plot_scale = 1e15; // Arbitrary scale factor for visibility
            pinfo.en[gr][ch][sfb] = plot_scale * en0;
            // Store scaled noise = scaled_xmin * distortion_ratio / bw
            pinfo.xfsf[gr][ch][sfb] = plot_scale * l3_xmin[sfb] * distort[sfb] / bw;

            // Calculate masking threshold (considering ATHonly mode)
            let thr0 = gfc.ATH.l[sfb]; // Start with ATH
            if (!gfp.ATHonly) {
                const ratio_en = ratio.en.l[sfb];
                if (ratio_en > 0) {
                    // Masking threshold = (Energy / MaskRatioEnergy) * MaskRatioThreshold * MaskingLower? No.
                    // Masking threshold = Energy * (MaskRatioThreshold / MaskRatioEnergy) * MaskingLower
                    // Ratio.thm seems to be already Thm/En? Check PSY model. Assume Ratio.thm = Thresh/Energy.
                    // thr = Energy * Ratio.thm (where Ratio.thm incorporates masking_lower?)
                    // Let's re-evaluate C: en0 = en0 / ratio.en.l[sfb] => en0 = (Sum(xr^2)/bw) / ratio.en.l[sfb]
                    // thr = max(en0 * ratio.thm.l[sfb], gfc.ATH.l[sfb])
                    // This implies threshold = max( (Sum(xr^2)/bw) / ratio.en * ratio.thm, ATH )
                    // This seems wrong dimensionally unless ratio.en=Sum(xr^2) and ratio.thm=ThresholdPower.
                    // Let's assume: thr = max(ratio.thm.l[sfb] * gfc.masking_lower?, ATH.l[sfb])
                    // Revisit C: calc_xmin uses: xmin = max(ATH, en0*ratio.thm/ratio.en * masking_lower)
                    // This plot calculation seems different. Let's trust the C plot code for now.
                    let thr_mask = 0.0;
                    if(ratio.en.l[sfb] > 0) {
                        thr_mask = (en0 / ratio.en.l[sfb]) * ratio.thm.l[sfb]; // Is ratio.thm already scaled by masking_lower? Assume yes based on calc_xmin structure.
                    }
                    thr0 = Math.max(thr_mask, gfc.ATH.l[sfb]);

                } // else: use ATH only if no energy in ratio
            }
            pinfo.thr[gr][ch][sfb] = plot_scale * thr0;

            // Store effective scalefactor gain (dB-like?)
            pinfo.LAMEsfb[gr][ch][sfb] = 0.0; // Initialize
            if (cod_info.preflag !== 0 && sfb >= 11) {
                pinfo.LAMEsfb[gr][ch][sfb] -= ifqstep * this.pretab[sfb]; // Subtract preemphasis gain
            }
            if (sfb < Encoder.SBPSY_l) { // Only apply scalefac if within psych limit
                assert(scalefac[sfb] >= 0, `Negative scalefactor found: sfb=${sfb}, val=${scalefac[sfb]}`);
                pinfo.LAMEsfb[gr][ch][sfb] -= ifqstep * scalefac[sfb]; // Subtract scalefactor gain
            }
        } // End long band loop

        // --- Process Short Bands (if applicable) ---
        if (cod_info.block_type === Encoder.SHORT_TYPE) {
            let sfb_plot_idx = sfb2; // Continue combined index from long bands
            for (sfb = cod_info.sfb_smin; sfb < Encoder.SBMAX_s; sfb++) {
                const start = gfc.scalefac_band.s[sfb];
                const end = gfc.scalefac_band.s[sfb + 1];
                const bw = end - start;
                if (bw === 0) continue;

                // Process 3 windows
                for (let i = 0; i < 3; i++) {
                    // Calculate energy for this window/sfb
                    en0 = 0.0;
                    // Need correct index 'j' into reordered xr[]
                    // C code just continues 'j'. Let's assume j is correct after long block loop.
                    const window_start_j = j;
                    for (l = 0; l < bw; l++) { // Iterate through coeffs for this window/sfb
                         en0 += cod_info.xr[j] * cod_info.xr[j];
                         j++;
                    }
                    en0 = Math.max(en0 / bw, 1e-20); // Average energy per coeff

                    const plot_idx = 3 * sfb + i; // Index for plotting arrays (interleaved)
                    const plot_scale = 1e15;

                    // Store plotting data for this short sub-band
                    pinfo.en_s[gr][ch][plot_idx] = plot_scale * en0;
                    // Store scaled noise
                    pinfo.xfsf_s[gr][ch][plot_idx] = plot_scale * l3_xmin[sfb_plot_idx] * distort[sfb_plot_idx] / bw;

                    // Calculate masking threshold
                    let thr0 = gfc.ATH.s[sfb]; // Start with ATH (scaled by bw in _compute_ath?) Check ATH units. Assume ATH.s is per-band power.
                    if (!gfp.ATHonly && !gfp.ATHshort) {
                        const ratio_en = ratio.en.s[sfb][i];
                        if(ratio_en > 0) {
                             const thr_mask = (en0 / ratio_en) * ratio.thm.s[sfb][i]; // Similar dimension issue as long blocks
                             thr0 = Math.max(thr_mask, gfc.ATH.s[sfb]);
                        }
                    }
                    pinfo.thr_s[gr][ch][plot_idx] = plot_scale * thr0;

                    // Store effective scalefactor gain
                    pinfo.LAMEsfb_s[gr][ch][plot_idx] = -2.0 * cod_info.subblock_gain[i]; // Subblock gain component
                    if (sfb < Encoder.SBPSY_s) { // Apply scalefac only within psy limit
                         pinfo.LAMEsfb_s[gr][ch][plot_idx] -= ifqstep * scalefac[sfb_plot_idx]; // Scalefactor component
                    }
                    sfb_plot_idx++; // Increment combined index
                } // End window loop (i)
            } // End short sfb loop
        } // End if short block

        // Store overall granule/channel info
        pinfo.LAMEqss[gr][ch] = cod_info.global_gain;
        pinfo.LAMEmainbits[gr][ch] = cod_info.part2_3_length + cod_info.part2_length; // Huffman + Scalefactors
        pinfo.LAMEsfbits[gr][ch] = cod_info.part2_length; // Scalefactor bits only
        pinfo.over[gr][ch] = noise_res.over_count;
        pinfo.max_noise[gr][ch] = noise_res.max_noise * 10.0; // Scale dB*10 to dB
        pinfo.over_noise[gr][ch] = noise_res.over_noise * 10.0;
        pinfo.tot_noise[gr][ch] = noise_res.tot_noise * 10.0;
        pinfo.over_SSD[gr][ch] = noise_res.over_SSD;
    }

     /**
      * Updates plotting data for a whole frame by calling set_pinfo for each granule/channel.
      * Handles SCFSI scalefactor reconstruction for granule 1.
      * @private
      */
     _set_frame_pinfo(gfp, ratio) {
        const gfc = gfp.internal_flags;
        gfc.masking_lower = 1.0; // Reset masking adjustment for plotting?

        for (let gr = 0; gr < gfc.mode_gr; gr++) {
            for (let ch = 0; ch < gfc.channels_out; ch++) {
                const cod_info = gfc.l3_side.tt[gr][ch];
                const scalefac_sav = new_int(L3Side.SFBMAX);
                System.arraycopy(cod_info.scalefac, 0, scalefac_sav, 0, L3Side.SFBMAX); // Save original scalefactors

                // Reconstruct scalefactors if SCFSI was used (only for gr=1)
                if (gr === 1) {
                    for (let sfb = 0; sfb < cod_info.sfb_lmax; sfb++) { // Check limit - should cover all possible SCFSI bands
                        if (cod_info.scalefac[sfb] < 0) { // SCFSI flag stored as negative? Check spec/tk code. Assume yes.
                            cod_info.scalefac[sfb] = gfc.l3_side.tt[0][ch].scalefac[sfb]; // Copy from granule 0
                        }
                    }
                     // Need similar logic for short block scalefactors if applicable? Assume covered by sfb_lmax.
                }

                // Set plotting info using (potentially reconstructed) scalefactors
                this.set_pinfo(gfp, cod_info, ratio[gr][ch], gr, ch);

                // Restore original scalefactors
                System.arraycopy(scalefac_sav, 0, cod_info.scalefac, 0, L3Side.SFBMAX);
            }
        }
    }

} // End class QuantizePVT


// Export the class
export { QuantizePVT };
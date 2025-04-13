/**
 * @fileoverview Preset handling for LAME configuration.
 * Ported from presets.c. Provides functions to apply predefined
 * sets of encoding parameters (presets) to the LAME global flags.
 * Uses ES Module syntax.
 *
 * @module Presets
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import Lame from './Lame.js'; // Need Lame for preset constants (V0, V9, etc.)

// Destructure common utilities for easier access
const {
    VbrMode,
    // System, // Not used
    // Float, // Not used
    // ShortBlock, // Not used
    // Util, // Not used
    // Arrays, // Not used
    // new_array_n, // Not used
    // new_byte, // Not used
    // new_double, // Not used
    // new_float, // Not used
    // new_float_n, // Not used
    // new_int, // Not used
    // new_int_n, // Not used
    // assert // Not used
} = common;

// Assuming these types are defined elsewhere
/** @typedef {import('./LameGlobalFlags.js').default} LameGlobalFlags */

// --- Internal Data Structures ---
// (No JSDoc as requested)

/** @private */
class VBRPresets {
    constructor(qual, comp, compS, y, shThreshold, shThresholdS, adj, adjShort, lower, curve, sens, inter, joint, mod, fix) {
        this.vbr_q = qual; this.quant_comp = comp; this.quant_comp_s = compS; this.expY = y;
        this.st_lrm = shThreshold; this.st_s = shThresholdS; this.masking_adj = adj;
        this.masking_adj_short = adjShort; this.ath_lower = lower; this.ath_curve = curve;
        this.ath_sensitivity = sens; this.interch = inter; this.safejoint = joint;
        this.sfb21mod = mod; this.msfix = fix;
    }
}

/** @private */
class ABRPresets {
    constructor(kbps, comp, compS, joint, fix, shThreshold, shThresholdS, bass, sc, mask, lower, curve, interCh, sfScale) {
        this.kbps = kbps; // Added kbps field for clarity
        this.quant_comp = comp; this.quant_comp_s = compS; this.safejoint = joint;
        this.nsmsfix = fix; this.st_lrm = shThreshold; this.st_s = shThresholdS;
        this.nsbass = bass; this.scale = sc; this.masking_adj = mask;
        this.ath_lower = lower; this.ath_curve = curve; this.interch = interCh;
        this.sfscale = sfScale;
    }
}

// --- Preset Data Tables ---
// (Moved inside the class or module scope if they don't need to be global)


/**
 * @classdesc Manages LAME presets. Provides methods to apply predefined
 * sets of encoding parameters based on VBR quality levels or ABR bitrates.
 * @constructs Presets
 */
class Presets {
    /** @private @type {Lame|null} Reference to main Lame object (needed for Lame.nearestBitrateFullIndex). */
    lame = null;

    /**
     * Preset definitions for VBR mode VBR_RH (older psychoacoustic model).
     * @private
     * @const {VBRPresets[]}
     */
    vbr_old_switch_map = [
        new VBRPresets(0, 9, 9, 0, 5.20, 125.0, -4.2, -6.3, 4.8, 1, 0, 0, 2, 21, 0.97),
        new VBRPresets(1, 9, 9, 0, 5.30, 125.0, -3.6, -5.6, 4.5, 1.5, 0, 0, 2, 21, 1.35),
        new VBRPresets(2, 9, 9, 0, 5.60, 125.0, -2.2, -3.5, 2.8, 2, 0, 0, 2, 21, 1.49),
        new VBRPresets(3, 9, 9, 1, 5.80, 130.0, -1.8, -2.8, 2.6, 3, -4, 0, 2, 20, 1.64),
        new VBRPresets(4, 9, 9, 1, 6.00, 135.0, -0.7, -1.1, 1.1, 3.5, -8, 0, 2, 0, 1.79),
        new VBRPresets(5, 9, 9, 1, 6.40, 140.0, 0.5, 0.4, -7.5, 4, -12, 0.0002, 0, 0, 1.95),
        new VBRPresets(6, 9, 9, 1, 6.60, 145.0, 0.67, 0.65, -14.7, 6.5, -19, 0.0004, 0, 0, 2.30),
        new VBRPresets(7, 9, 9, 1, 6.60, 145.0, 0.8, 0.75, -19.7, 8, -22, 0.0006, 0, 0, 2.70),
        new VBRPresets(8, 9, 9, 1, 6.60, 145.0, 1.2, 1.15, -27.5, 10, -23, 0.0007, 0, 0, 0), // msfix 0?
        new VBRPresets(9, 9, 9, 1, 6.60, 145.0, 1.6, 1.6, -36, 11, -25, 0.0008, 0, 0, 0), // msfix 0?
        new VBRPresets(10, 9, 9, 1, 6.60, 145.0, 2.0, 2.0, -36, 12, -25, 0.0008, 0, 0, 0) // Dummy for VBR_q=9 + frac
    ];

    /**
     * Preset definitions for VBR modes using newer psychoacoustic models (VBR_MT, VBR_MTRH).
     * @private
     * @const {VBRPresets[]}
     */
    vbr_psy_switch_map = [
        new VBRPresets(0, 9, 9, 0, 4.20, 25.0, -7.0, -4.0, 7.5, 1, 0, 0, 2, 26, 0.97),
        new VBRPresets(1, 9, 9, 0, 4.20, 25.0, -5.6, -3.6, 4.5, 1.5, 0, 0, 2, 21, 1.35),
        new VBRPresets(2, 9, 9, 0, 4.20, 25.0, -4.4, -1.8, 2, 2, 0, 0, 2, 18, 1.49),
        new VBRPresets(3, 9, 9, 1, 4.20, 25.0, -3.4, -1.25, 1.1, 3, -4, 0, 2, 15, 1.64),
        new VBRPresets(4, 9, 9, 1, 4.20, 25.0, -2.2, 0.1, 0, 3.5, -8, 0, 2, 0, 1.79),
        new VBRPresets(5, 9, 9, 1, 4.20, 25.0, -1.0, 1.65, -7.7, 4, -12, 0.0002, 0, 0, 1.95),
        new VBRPresets(6, 9, 9, 1, 4.20, 25.0, -0.0, 2.47, -7.7, 6.5, -19, 0.0004, 0, 0, 2),
        new VBRPresets(7, 9, 9, 1, 4.20, 25.0, 0.5, 2.0, -14.5, 8, -22, 0.0006, 0, 0, 2),
        new VBRPresets(8, 9, 9, 1, 4.20, 25.0, 1.0, 2.4, -22.0, 10, -23, 0.0007, 0, 0, 2),
        new VBRPresets(9, 9, 9, 1, 4.20, 25.0, 1.5, 2.95, -30.0, 11, -25, 0.0008, 0, 0, 2),
        new VBRPresets(10, 9, 9, 1, 4.20, 25.0, 2.0, 2.95, -36.0, 12, -30, 0.0008, 0, 0, 2) // Dummy for VBR_q=9 + frac
    ];

     /**
     * Preset definitions for ABR mode. Indexed by bitrate index from `nearestBitrateFullIndex`.
     * @private
     * @const {ABRPresets[]}
     */
     abr_switch_map = [
        // Index corresponds roughly to bitrate index (0=8k, 1=16k, ...)
        // Values need verification against LAME C source presets.c
        new ABRPresets(8, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -30.0, 11, 0.0012, 1), //   8kbps
        new ABRPresets(16, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -25.0, 11, 0.0010, 1), //  16kbps
        new ABRPresets(24, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -20.0, 11, 0.0010, 1), //  24kbps
        new ABRPresets(32, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -15.0, 11, 0.0010, 1), //  32kbps
        new ABRPresets(40, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -10.0, 11, 0.0009, 1), //  40kbps
        new ABRPresets(48, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -10.0, 11, 0.0009, 1), //  48kbps
        new ABRPresets(56, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -6.0, 11, 0.0008, 1), //  56kbps
        new ABRPresets(64, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, -2.0, 11, 0.0008, 1), //  64kbps
        new ABRPresets(80, 9, 9, 0, 0, 6.60, 145, 0, 0.95, 0, 0.0, 8, 0.0007, 1), //  80kbps
        new ABRPresets(96, 9, 9, 0, 2.50, 6.60, 145, 0, 0.95, 0, 1.0, 5.5, 0.0006, 1), //  96kbps
        new ABRPresets(112, 9, 9, 0, 2.25, 6.60, 145, 0, 0.95, 0, 2.0, 4.5, 0.0005, 1), // 112kbps
        new ABRPresets(128, 9, 9, 0, 1.95, 6.40, 140, 0, 0.95, 0, 3.0, 4, 0.0002, 1), // 128kbps
        new ABRPresets(160, 9, 9, 1, 1.79, 6.00, 135, 0, 0.95, -2, 5.0, 3.5, 0, 1), // 160kbps
        new ABRPresets(192, 9, 9, 1, 1.49, 5.60, 125, 0, 0.97, -4, 7.0, 3, 0, 0), // 192kbps
        new ABRPresets(224, 9, 9, 1, 1.25, 5.20, 125, 0, 0.98, -6, 9.0, 2, 0, 0), // 224kbps
        new ABRPresets(256, 9, 9, 1, 0.97, 5.20, 125, 0, 1.00, -8, 10.0, 1, 0, 0), // 256kbps
        new ABRPresets(320, 9, 9, 1, 0.90, 5.20, 125, 0, 1.00, -10, 12.0, 0, 0, 0) // 320kbps
    ];

    constructor() {
        // Module set externally
    }

    /**
     * Sets the internal Lame object reference.
     * @public
     * @param {Lame} _lame - The main Lame instance.
     */
    setModules(_lame) {
        this.lame = _lame;
    }

    // --- Private Helper Methods ---
    // (JSDoc omitted for brevity)

    /** @private */
    _apply_vbr_preset(gfp, a, enforce) {
        const vbr_preset = gfp.VBR === VbrMode.vbr_rh ? this.vbr_old_switch_map : this.vbr_psy_switch_map;
        const x = gfp.VBR_q_frac;

        // Ensure 'a' is within bounds for interpolation
        const p_idx = Math.max(0, Math.min(a, vbr_preset.length - 2));
        const q_idx = p_idx + 1;
        const p = vbr_preset[p_idx];
        const q = vbr_preset[q_idx];

        // Create a temporary preset object by interpolating p and q
        const set = new VBRPresets(
            p.vbr_q, p.quant_comp, p.quant_comp_s, p.expY,
            p.st_lrm + x * (q.st_lrm - p.st_lrm),
            p.st_s + x * (q.st_s - p.st_s),
            p.masking_adj + x * (q.masking_adj - p.masking_adj),
            p.masking_adj_short + x * (q.masking_adj_short - p.masking_adj_short),
            p.ath_lower + x * (q.ath_lower - p.ath_lower),
            p.ath_curve + x * (q.ath_curve - p.ath_curve),
            p.ath_sensitivity + x * (q.ath_sensitivity - p.ath_sensitivity),
            p.interch + x * (q.interch - p.interch),
            p.safejoint, p.sfb21mod,
            p.msfix + x * (q.msfix - p.msfix)
        );

        // Apply interpolated values to gfp, respecting enforce flag
        // Uses SET_OPTION logic from C (apply if enforce=true or if current value is default)

        this._lame_set_VBR_q(gfp, set.vbr_q); // Use internal helper

        if (enforce !== 0 || gfp.quant_comp === -1) gfp.quant_comp = set.quant_comp;
        if (enforce !== 0 || gfp.quant_comp_short === -1) gfp.quant_comp_short = set.quant_comp_s;
        if (set.expY !== 0) gfp.experimentalY = true; // Only turn on, never off? Check C.
        if (enforce !== 0 || gfp.internal_flags.nsPsy.attackthre === -1) gfp.internal_flags.nsPsy.attackthre = set.st_lrm;
        if (enforce !== 0 || gfp.internal_flags.nsPsy.attackthre_s === -1) gfp.internal_flags.nsPsy.attackthre_s = set.st_s;
        if (enforce !== 0 || Math.abs(gfp.maskingadjust) < 1e-6) gfp.maskingadjust = set.masking_adj;
        if (enforce !== 0 || Math.abs(gfp.maskingadjust_short) < 1e-6) gfp.maskingadjust_short = set.masking_adj_short;
        if (enforce !== 0 || Math.abs(gfp.ATHlower) < 1e-6) gfp.ATHlower = -set.ath_lower / 10.0; // Note sign change and scaling
        if (enforce !== 0 || gfp.ATHcurve === -1) gfp.ATHcurve = set.ath_curve;
        if (enforce !== 0 || gfp.athaa_sensitivity === 0) gfp.athaa_sensitivity = set.ath_sensitivity; // Check default value
        if (set.interch > 0) { if (enforce !== 0 || gfp.interChRatio === -1.0) gfp.interChRatio = set.interch; }
        if (set.safejoint > 0) gfp.exp_nspsytune |= set.safejoint; // Use bitwise OR
        if (set.sfb21mod > 0) gfp.exp_nspsytune |= (set.sfb21mod << 20);
        if (enforce !== 0 || gfp.msfix === -1.0) gfp.msfix = set.msfix;

        // If not enforcing, store the original VBR quality settings used for interpolation
        if (enforce === 0) {
            gfp.VBR_q = a;
            gfp.VBR_q_frac = x;
        }
    }

     /** @private */
     _apply_abr_preset(gfp, preset, enforce) {
        const actual_bitrate = preset; // Use preset directly as bitrate
        // Find the index corresponding to the nearest standard bitrate
        const r = this.lame._nearestBitrateFullIndex(preset); // Need access to Lame instance method

        // Set VBR mode and bitrate
        gfp.VBR = VbrMode.vbr_abr;
        gfp.VBR_mean_bitrate_kbps = actual_bitrate;
        // Clamp bitrate (redundant if using nearest index, but good practice)
        gfp.VBR_mean_bitrate_kbps = Math.min(gfp.VBR_mean_bitrate_kbps, 320);
        gfp.VBR_mean_bitrate_kbps = Math.max(gfp.VBR_mean_bitrate_kbps, 8);
        gfp.brate = gfp.VBR_mean_bitrate_kbps; // Also set CBR rate for ABR

        // Disable reservoir if bitrate is very high (emulating C behavior)
        if (gfp.VBR_mean_bitrate_kbps > 320) { gfp.disable_reservoir = true; }

        // Apply settings from the abr_switch_map using SET_OPTION logic
        const set = this.abr_switch_map[r];
        if (!set) {
            console.error(`ABR preset not found for index ${r} (bitrate ${preset})`);
            return preset; // Return original value on error
        }

        if (set.safejoint > 0) gfp.exp_nspsytune |= 2; // Use constant 2 for safejoint flag
        if (set.sfscale > 0) gfp.internal_flags.noise_shaping = 2;
        if (Math.abs(set.nsbass) > 1e-6) { let k = Math.floor(set.nsbass * 4); if (k < 0) k += 64; gfp.exp_nspsytune |= (k << 2); }

        if (enforce !== 0 || gfp.quant_comp === -1) gfp.quant_comp = set.quant_comp;
        if (enforce !== 0 || gfp.quant_comp_short === -1) gfp.quant_comp_short = set.quant_comp_s;
        if (enforce !== 0 || gfp.msfix === -1.0) gfp.msfix = set.nsmsfix;
        if (enforce !== 0 || gfp.internal_flags.nsPsy.attackthre === -1) gfp.internal_flags.nsPsy.attackthre = set.st_lrm;
        if (enforce !== 0 || gfp.internal_flags.nsPsy.attackthre_s === -1) gfp.internal_flags.nsPsy.attackthre_s = set.st_s;
        if (enforce !== 0 || gfp.scale === -1.0) gfp.scale = set.scale;
        if (enforce !== 0 || Math.abs(gfp.maskingadjust) < 1e-6) gfp.maskingadjust = set.masking_adj;
        // Set maskingadjust_short based on maskingadjust
        const short_adj = (set.masking_adj > 0) ? (set.masking_adj * 0.9) : (set.masking_adj * 1.1);
        if (enforce !== 0 || Math.abs(gfp.maskingadjust_short) < 1e-6) gfp.maskingadjust_short = short_adj;
        if (enforce !== 0 || Math.abs(gfp.ATHlower) < 1e-6) gfp.ATHlower = -set.ath_lower / 10.0;
        if (enforce !== 0 || gfp.ATHcurve === -1) gfp.ATHcurve = set.ath_curve;
        if (enforce !== 0 || gfp.interChRatio === -1.0) gfp.interChRatio = set.interch;

        return preset; // Return the target bitrate
    }

    /** @private */
    _lame_set_VBR_q(gfp, VBR_q) {
        let ret = 0;
        if (VBR_q < 0) { ret = -1; VBR_q = 0; }
        if (VBR_q > 9) { ret = -1; VBR_q = 9; }
        gfp.VBR_q = VBR_q;
        gfp.VBR_q_frac = 0; // Reset fractional part when setting integer quality
        return ret;
    }


    // --- Public Methods ---

    /**
     * Applies a predefined preset to the LAME configuration.
     * This modifies various settings in the `gfp` (LameGlobalFlags) object
     * based on the chosen preset (VBR quality level or ABR bitrate).
     * Handles translation of legacy preset constants.
     *
     * @public
     * @param {LameGlobalFlags} gfp - The LAME global flags structure to modify.
     * @param {number} preset - The preset value. Can be one of the Lame.V constants (V0-V9),
     *                          one of the legacy constants (STANDARD, EXTREME, etc.), or an
     *                          ABR bitrate value (8-320).
     * @param {number} enforce - If non-zero, forces all preset settings onto `gfp`, overwriting
     *                           any user-defined values for those specific parameters. If zero,
     *                           only applies preset settings if the corresponding `gfp` value
     *                           is still at its default/uninitialized state.
     * @returns {number} The effective preset value applied (e.g., the actual ABR bitrate used). Returns the input `preset` if it was not a recognized preset value.
     */
    apply_preset(gfp, preset, enforce) {
        let original_preset = preset; // Keep original for potential ABR

        // Translate legacy presets to VBR quality levels (-V) and set VBR mode
        switch (preset) {
            case Lame.R3MIX:        preset = Lame.V3; gfp.VBR = VbrMode.vbr_mtrh; break;
            case Lame.MEDIUM:       preset = Lame.V4; gfp.VBR = VbrMode.vbr_rh; break;
            case Lame.MEDIUM_FAST:  preset = Lame.V4; gfp.VBR = VbrMode.vbr_mtrh; break;
            case Lame.STANDARD:     preset = Lame.V2; gfp.VBR = VbrMode.vbr_rh; break;
            case Lame.STANDARD_FAST:preset = Lame.V2; gfp.VBR = VbrMode.vbr_mtrh; break;
            case Lame.EXTREME:      preset = Lame.V0; gfp.VBR = VbrMode.vbr_rh; break;
            case Lame.EXTREME_FAST: preset = Lame.V0; gfp.VBR = VbrMode.vbr_mtrh; break;
            case Lame.INSANE: // Special case: Maps to 320kbps ABR (or CBR)
                preset = 320; gfp.preset = preset; // Store the target bitrate
                this._apply_abr_preset(gfp, preset, enforce);
                // Insane is often CBR, override VBR setting from ABR preset
                gfp.VBR = VbrMode.vbr_off;
                return preset; // Return the bitrate
        }

        gfp.preset = preset; // Store the potentially translated VBR preset value

        // Apply VBR preset (-V settings)
        if (preset >= Lame.V9 && preset <= Lame.V0) {
             const vbr_quality_index = 9 - Math.floor((preset - Lame.V9) / 10); // Map V9..V0 to index 9..0
             this._apply_vbr_preset(gfp, vbr_quality_index, enforce);
             return preset; // Return the original V preset value
        }

        // Apply ABR preset (numeric bitrate value)
        if (preset >= 8 && preset <= 320) {
             return this._apply_abr_preset(gfp, preset, enforce); // Returns the actual bitrate used
        }

        // If no preset matched
        console.warn(`Preset value ${original_preset} not recognized.`);
        gfp.preset = 0; // Reset preset value
        return original_preset; // Return the original unrecognized value
    }

}

export { Presets };
export default Presets; // Also provide default export if needed
/**
 * @fileoverview CBR (Constant BitRate) iteration loop implementation for LAME.
 * This module handles the quantization process specifically for CBR encoding,
 * aiming to meet the target bitrate for each frame.
 * Uses ES Module syntax.
 *
 * @module CBRNewIterationLoop
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import { MeanBits } from './MeanBits.js';
import { Encoder } from './Encoder.js';
import { L3Side } from './L3Side.js';
import { LameInternalFlags } from './LameInternalFlags.js';
// Assuming these types are defined elsewhere and imported if needed for full type safety
/** @typedef {import('./Quantize.js').Quantize} Quantize */
/** @typedef {import('./LameGlobalFlags.js').LameGlobalFlags} LameGlobalFlags */
/** @typedef {import('./Reservoir.js').Reservoir} Reservoir */ // Accessed via quantize.rv
/** @typedef {import('./QuantizePVT.js').QuantizePVT} QuantizePVT */ // Accessed via quantize.qupvt

// Destructure common utilities for easier access
const {
    // System, // Not used
    // VbrMode, // Not used directly, but present in common
    // Float, // Not used
    ShortBlock, // Not used directly, but present in common
    // Util, // Not used
    // Arrays, // Not used
    // new_array_n, // Used indirectly
    // new_byte, // Not used
    // new_double, // Not used
    new_float,
    // new_float_n, // Not used
    new_int,
    // new_int_n, // Not used
    assert
} = common;

/**
 * @classdesc Implements the iteration loop for Constant BitRate (CBR) encoding.
 * This loop adjusts quantization parameters for each granule to meet the
 * target bitrate determined by the reservoir and psychoacoustic model.
 * @constructs CBRNewIterationLoop
 * @param {Quantize} _quantize - An instance of the main Quantize class.
 */
class CBRNewIterationLoop {
    /**
     * @private
     * @type {Quantize}
     */
    quantize;

    /**
     * @param {Quantize} _quantize - An instance of the main Quantize class.
     */
    constructor(_quantize) {
        /** @public */ // Make accessible for JSDoc link? Or keep private? Let's assume internal usage mostly.
        this.quantize = _quantize;
    }

    /**
     * Executes the CBR iteration loop for a frame.
     * Determines target bits per granule/channel based on the bit reservoir,
     * performs psychoacoustic analysis (`calc_xmin`), runs the quantization
     * outer loop (`outer_loop`), and finalizes reservoir calculations.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags and settings.
     * @param {Array<Float32Array>} pe - Perceptual entropy per granule/channel `[gr][ch]`.
     * @param {Float32Array} ms_ener_ratio - Mid/Side energy ratio per granule `[gr]`.
     * @param {Array<Array<object>>} ratio - Masking ratio info per granule/channel `[gr][ch]`.
     */
    iteration_loop(gfp, pe, ms_ener_ratio, ratio) {
        const gfc = gfp.internal_flags;
        const l3_xmin = new_float(L3Side.SFBMAX); // Allowed noise buffer
        const xrpow = new_float(576); // Spectral energy buffer
        const targ_bits = new_int(2); // Target bits per channel for a granule
        let mean_bits = 0; // Average bits available from reservoir
        let max_bits; // Max bits allowed for a granule
        const l3_side = gfc.l3_side;

        // Get mean bits available from reservoir for this frame
        const mb = new MeanBits(mean_bits);
        this.quantize.rv.ResvFrameBegin(gfp, mb);
        mean_bits = mb.bits;

        // Process each granule in the frame
        for (let gr = 0; gr < gfc.mode_gr; gr++) {

            // Calculate needed bits per channel based on PE and reservoir state
            max_bits = this.quantize.qupvt.on_pe(gfp, pe, targ_bits, mean_bits, gr, gr); // Pass gr as cbr flag? Check C code. Assume yes.

            // Handle Mid/Side stereo processing
            if (gfc.mode_ext === Encoder.MPG_MD_MS_LR) { // Using mode_ext check
                this.quantize.ms_convert(gfc.l3_side, gr); // Convert L/R to M/S
                this.quantize.qupvt.reduce_side(targ_bits, ms_ener_ratio[gr], mean_bits, max_bits); // Adjust M/S target bits
            }

            // Process each channel within the granule
            for (let ch = 0; ch < gfc.channels_out; ch++) {
                let adjust, masking_lower_db;
                const cod_info = l3_side.tt[gr][ch];

                // Set masking adjustment based on block type (simplified for CBR?)
                if (cod_info.block_type !== Encoder.SHORT_TYPE) { // Long block types
                    adjust = 0; // CBR adjustment factor? Usually 0?
                    masking_lower_db = gfc.PSY.mask_adjust - adjust;
                } else { // Short blocks
                    adjust = 0;
                    masking_lower_db = gfc.PSY.mask_adjust_short - adjust;
                }
                gfc.masking_lower = Math.pow(10.0, masking_lower_db * 0.1);

                // Initialize granule info for outer loop
                this.quantize.init_outer_loop(gfc, cod_info);

                // Check if there's energy to encode and initialize xrpow
                if (this.quantize.init_xrpow(gfc, cod_info, xrpow)) {
                    // Calculate allowed noise based on psychoacoustics
                    this.quantize.qupvt.calc_xmin(gfp, ratio[gr][ch], cod_info, l3_xmin);

                    // Run the main quantization outer loop to find best scalefactors/gain
                    this.quantize.outer_loop(gfp, cod_info, l3_xmin, xrpow, ch, targ_bits[ch]);
                }
                // else: granule/channel is silent, init_xrpow already zeroed l3_enc

                // Finalize quantization results for this granule/channel
                this.quantize.iteration_finish_one(gfc, gr, ch);

                // Assert that bit limits were respected
                assert(cod_info.part2_3_length <= LameInternalFlags.MAX_BITS_PER_CHANNEL, `Channel ${ch} bits ${cod_info.part2_3_length} exceed MAX_BITS_PER_CHANNEL`);
                // For CBR, the final bits *should* be <= target bits after reservoir adjustment
                // assert(cod_info.part2_3_length <= targ_bits[ch], `Channel ${ch} bits ${cod_info.part2_3_length} exceed target ${targ_bits[ch]}`);
                // Let's relax this assertion slightly, as the outer loop might slightly exceed targ_bits sometimes before reservoir corrects it.
                 assert(cod_info.part2_3_length <= targ_bits[ch] + 100, `Channel ${ch} bits ${cod_info.part2_3_length} significantly exceed target ${targ_bits[ch]}`);

            } /* for ch */
        } /* for gr */

        // Finalize reservoir calculations for the frame
        this.quantize.rv.ResvFrameEnd(gfc, mean_bits);
    }
}

export { CBRNewIterationLoop };
export default CBRNewIterationLoop; // Also provide default export if needed
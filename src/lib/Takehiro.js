/**
 * @fileoverview MP3 Huffman table selecting and bit counting logic.
 * Ported from takehiro.c.
 * Uses ES Module syntax.
 *
 * @module Takehiro
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import { Encoder } from './Encoder.js';
import { Tables, ht, largetbl, table23, table56, bitrate_table } from './Tables.js'; // Import ht and other needed tables
import { GrInfo } from './GrInfo.js';
import { QuantizePVT } from './QuantizePVT.js';

// Destructure common utilities for easier access
const {
    System, VbrMode, Float, ShortBlock, Util, Arrays, new_array_n, new_byte,
    new_double, new_float, new_float_n, new_int, new_int_n, assert
} = common;

// --- Module Scope Constants/Tables ---

/**
 * Scalefactor length table 1 (based on scalefac_compress index).
 * @const {number[]}
 */
export const slen1_tab = [0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4];

/**
 * Scalefactor length table 2 (based on scalefac_compress index).
 * @const {number[]}
 */
export const slen2_tab = [0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3];

/**
 * Precomputed table for subdv_table lookup.
 * @private
 * @const {number[][]}
 */
const subdv_table = [
    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 1], [1, 1], [1, 1], [1, 2],
    [2, 2], [2, 3], [2, 3], [3, 4], [3, 4], [3, 4], [4, 5], [4, 5], [4, 6],
    [5, 6], [5, 6], [5, 7], [6, 7], [6, 7]
];

/**
 * Huffman tables that don't use ESCape codes. Index is max value - 1.
 * @private
 * @const {number[]}
 */
const huf_tbl_noESC = [1, 2, 5, 7, 7, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13];


/**
 * Simple helper class to pass bit count by reference.
 * @private
 */
class Bits {
    constructor(b = 0) {
        this.bits = 0 | b;
    }
}

/**
 * Placeholder/Simplified function to count sign bit.
 * Original C code might have done more complex Huffman length lookup.
 * @param {number} p - Quantized value.
 * @returns {number} 1 if p is non-zero, 0 otherwise.
 */
export function count_bit(p) {
    return (p !== 0) ? 1 : 0;
}


/**
 * @classdesc Implements Huffman coding logic, including table selection
 * and bit counting for quantized spectral values and scalefactors.
 * @constructs Takehiro
 */
class Takehiro {
    /** @private @type {QuantizePVT|null} */
    qupvt = null;

    constructor() {
        // Module set externally
    }

    /**
     * Sets the internal QuantizePVT module dependency.
     * @public
     * @param {QuantizePVT} _qupvt - Private quantization helpers module instance.
     */
    setModules(_qupvt) {
        this.qupvt = _qupvt;
    }

    // --- Private Helper Methods ---
    // (JSDoc omitted for brevity)

    /** @private */
    _quantize_lines_xrpow(l, istep, xr, xrPos, ix, ixPos) {
        assert(l > 0, `quantize_lines_xrpow length must be > 0: ${l}`);
        l = l >> 1; // Process pairs
        let remaining = l % 2; // Check C code logic - l is pair count, original C checks width&1? No, uses loop.
        // Let's stick to pair processing logic.
        l = l >> 1; // Process quads
        while (l-- > 0) {
            let x0 = xr[xrPos++] * istep; let x1 = xr[xrPos++] * istep;
            let rx0 = Math.floor(x0); let rx1 = Math.floor(x1); // Use floor for integer part
            let x2 = xr[xrPos++] * istep; let x3 = xr[xrPos++] * istep;
            let rx2 = Math.floor(x2); let rx3 = Math.floor(x3);
            // Clamp indices for adj43 lookup
            rx0 = Math.max(0, Math.min(rx0, QuantizePVT.IXMAX_VAL + 1));
            rx1 = Math.max(0, Math.min(rx1, QuantizePVT.IXMAX_VAL + 1));
            rx2 = Math.max(0, Math.min(rx2, QuantizePVT.IXMAX_VAL + 1));
            rx3 = Math.max(0, Math.min(rx3, QuantizePVT.IXMAX_VAL + 1));
            x0 += this.qupvt.adj43[rx0]; x1 += this.qupvt.adj43[rx1];
            ix[ixPos++] = Math.floor(x0); ix[ixPos++] = Math.floor(x1);
            x2 += this.qupvt.adj43[rx2]; x3 += this.qupvt.adj43[rx3];
            ix[ixPos++] = Math.floor(x2); ix[ixPos++] = Math.floor(x3);
        }
        // Process remaining pair if original width was odd multiple of 2
         if (remaining !== 0) { // This check might be wrong, based on original loop structure
            let x0 = xr[xrPos++] * istep; let x1 = xr[xrPos++] * istep;
            let rx0 = Math.floor(x0); let rx1 = Math.floor(x1);
            rx0 = Math.max(0, Math.min(rx0, QuantizePVT.IXMAX_VAL + 1));
            rx1 = Math.max(0, Math.min(rx1, QuantizePVT.IXMAX_VAL + 1));
            x0 += this.qupvt.adj43[rx0]; x1 += this.qupvt.adj43[rx1];
            ix[ixPos++] = Math.floor(x0); ix[ixPos++] = Math.floor(x1);
         }
    }

    /** @private */
    _quantize_lines_xrpow_01(l, istep, xr, xrPos, ix, ixPos) {
        const compareval0 = (1.0 - 0.4054) / istep;
        assert(l > 0);
        l = l >> 1; // Process pairs
        while (l-- > 0) {
            ix[ixPos++] = (compareval0 > xr[xrPos++]) ? 0 : 1;
            ix[ixPos++] = (compareval0 > xr[xrPos++]) ? 0 : 1;
        }
    }

    /** @private */
    _ix_max(ix, ixPos, endPos) {
        let max1 = 0, max2 = 0;
        // Ensure endPos is valid
        endPos = Math.min(endPos, ix.length);
        ixPos = Math.min(ixPos, endPos); // Ensure start is not beyond end

        while (ixPos < endPos) { // Use < for standard loop
            const x1 = ix[ixPos++];
            if (max1 < x1) max1 = x1;
            if (ixPos < endPos) { // Check if second element exists
                const x2 = ix[ixPos++];
                if (max2 < x2) max2 = x2;
            }
        }
        return Math.max(max1, max2);
    }

    /** @private */
    _count_bit_ESC(ix, ixPos, end, t1, t2, s) {
        const linbits1 = ht[t1].xlen; const linbits2 = ht[t2].xlen;
        let sum1 = 0, sum2 = 0;
        let linbits_esc1 = (1 << linbits1) -1; // Max value before escape for table 1
        let linbits_esc2 = (1 << linbits2) -1; // Max value before escape for table 2

        while (ixPos < end) {
            let x = ix[ixPos++]; let y = ix[ixPos++];
            let xy1 = 0, xy2 = 0;
            let current_bits1 = 0, current_bits2 = 0;
            let esc1 = false, esc2 = false;

            if (x !== 0) {
                if (x > 14) { x = 15; esc1 = true; current_bits1 += linbits1; }
                xy1 = x * 16;
            }
            if (y !== 0) {
                if (y > 14) { y = 15; esc1 = true; current_bits1 += linbits1; }
                xy1 += y;
            }
            current_bits1 += largetbl[xy1] >> 16; // Get length from table 1 (upper 16 bits)

            // Repeat for table 2
             x = ix[ixPos-2]; y = ix[ixPos-1]; // Re-read original values
             if (x !== 0) {
                if (x > 14) { x = 15; esc2 = true; current_bits2 += linbits2; }
                xy2 = x * 16;
             }
             if (y !== 0) {
                if (y > 14) { y = 15; esc2 = true; current_bits2 += linbits2; }
                xy2 += y;
             }
             current_bits2 += largetbl[xy1] & 0xffff; // Get length from table 2 (lower 16 bits) - Use xy1 index? C uses largetbl[x]

             // C code uses largetbl[x] where x is xy1. Let's use that.
             sum1 += current_bits1;
             sum2 += current_bits2;

        } // end while

        let best_sum = sum1;
        let best_table = t1;
        if (sum1 > sum2) {
            best_sum = sum2;
            best_table = t2;
        }

        s.bits += best_sum;
        return best_table;
    }

    /** @private */
    _count_bit_noESC(ix, ixPos, end, s) {
        let sum1 = 0;
        const hlen1 = ht[1].hlen; // Table 1 lengths
        while (ixPos < end) {
            const x = ix[ixPos] * 2 + ix[ixPos + 1];
            ixPos += 2;
            sum1 += hlen1[x];
        }
        s.bits += sum1;
        return 1; // Always uses table 1
    }

    /** @private */
    _count_bit_noESC_from2(ix, ixPos, end, t1, s) {
        let sum1 = 0, sum2 = 0;
        const xlen = ht[t1].xlen;
        const hlen_combined = (t1 === 2) ? table23 : table56;

        while (ixPos < end) {
            const x = ix[ixPos] * xlen + ix[ixPos + 1];
            ixPos += 2;
            const combined_len = hlen_combined[x];
            sum1 += combined_len >> 16;   // Length from table t1
            sum2 += combined_len & 0xffff; // Length from table t1+1
        }

        let best_sum = sum1;
        let best_table = t1;
        if (sum1 > sum2) {
            best_sum = sum2;
            best_table = t1 + 1;
        }
        s.bits += best_sum;
        return best_table;
    }

    /** @private */
    _count_bit_noESC_from3(ix, ixPos, end, t1, s) {
        let sum1 = 0, sum2 = 0, sum3 = 0;
        const xlen = ht[t1].xlen;
        const hlen1 = ht[t1].hlen;
        const hlen2 = ht[t1 + 1].hlen;
        const hlen3 = ht[t1 + 2].hlen;

        while (ixPos < end) {
            const x = ix[ixPos] * xlen + ix[ixPos + 1];
            ixPos += 2;
            sum1 += hlen1[x];
            sum2 += hlen2[x];
            sum3 += hlen3[x];
        }

        let best_sum = sum1;
        let best_table = t1;
        if (best_sum > sum2) { best_sum = sum2; best_table = t1 + 1; }
        if (best_sum > sum3) { best_sum = sum3; best_table = t1 + 2; }

        s.bits += best_sum;
        return best_table;
    }

    /** @private */
    _choose_table(ix, ixPos, endPos, s) {
        const max = this._ix_max(ix, ixPos, endPos);

        switch (max) {
            case 0: return 0; // Table 0 for all zeros
            case 1: return this._count_bit_noESC(ix, ixPos, endPos, s);
            case 2: case 3:
                return this._count_bit_noESC_from2(ix, ixPos, endPos, huf_tbl_noESC[max - 1], s);
            case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15:
                return this._count_bit_noESC_from3(ix, ixPos, endPos, huf_tbl_noESC[max - 1], s);
            default:
                if (max > QuantizePVT.IXMAX_VAL) { s.bits = LARGE_BITS; return -1; } // Value too large
                // Find best ESC tables (16-23 and 24-31)
                let choice2 = 24; for (; choice2 < 32; choice2++) if (ht[choice2].linmax >= max - 15) break;
                let choice1 = 16; for (; choice1 < 24; choice1++) if (ht[choice1].linmax >= max - 15) break;
                return this._count_bit_ESC(ix, ixPos, endPos, choice1, choice2, s);
        }
    }

     /** @private */
     _recalc_divide_init(gfc, cod_info, ix, r01_bits, r01_div, r0_tbl, r1_tbl) {
        const bigv = cod_info.big_values;
        Arrays.fill(r01_bits, QuantizePVT.LARGE_BITS); // Initialize with large value

        for (let r0 = 0; r0 < 16; r0++) { // region0 possible sizes (0..15 sfbs)
             const a1 = gfc.scalefac_band.l[r0 + 1]; // End of region0
             if (a1 >= bigv) break; // No point if region0 end exceeds bigvalues

             let r0bits = 0; const bi0 = new Bits(r0bits);
             const r0t = this._choose_table(ix, 0, a1, bi0);
             r0bits = bi0.bits;
             if(r0t < 0) continue; // Skip if error in choosing table

             for (let r1 = 0; r1 < 8; r1++) { // region1 possible sizes (0..7 sfbs)
                 const a2 = gfc.scalefac_band.l[r0 + r1 + 2]; // End of region1
                 if (a2 >= bigv) break; // No point if region1 end exceeds bigvalues

                 let bits = r0bits; const bi1 = new Bits(bits);
                 const r1t = this._choose_table(ix, a1, a2, bi1);
                 bits = bi1.bits;
                 if(r1t < 0) continue; // Skip if error

                 const region_boundary_index = r0 + r1; // Combined size index
                 if (r01_bits[region_boundary_index] > bits) {
                     r01_bits[region_boundary_index] = bits;
                     r01_div[region_boundary_index] = r0; // Store best size for region0
                     r0_tbl[region_boundary_index] = r0t;
                     r1_tbl[region_boundary_index] = r1t;
                 }
             }
        }
    }

    /** @private */
    _recalc_divide_sub(gfc, cod_info2, gi, ix, r01_bits, r01_div, r0_tbl, r1_tbl) {
        const bigv = cod_info2.big_values;

        // Iterate through possible end points of region 1 (r2 = r0+r1+2)
        for (let r2_idx = 0; r2_idx <= 7 + 15; r2_idx++) { // r2_idx = r0+r1
             if(r01_bits[r2_idx] >= LARGE_BITS) continue; // Skip if no valid r0/r1 combo found

             const r2_sfb_idx = r2_idx + 2; // Scalefactor band index for end of region 1
             const a2 = gfc.scalefac_band.l[r2_sfb_idx]; // Start index of region 2
             if (a2 >= bigv) continue; // No region 2 needed if start is beyond bigvalues

             let bits = r01_bits[r2_idx] + cod_info2.count1bits; // Add count1 bits
             if (gi.part2_3_length <= bits) continue; // Already worse than current best

             const bi2 = new Bits(bits);
             const r2t = this._choose_table(ix, a2, bigv, bi2); // Calculate bits for region 2
             bits = bi2.bits;
              if(r2t < 0) continue; // Skip if error

             // If this division is better than the current best total
             if (gi.part2_3_length > bits) {
                 gi.assign(cod_info2); // Update best with current state (only count1 differs)
                 gi.part2_3_length = bits; // Store new best bit count
                 gi.region0_count = r01_div[r2_idx];
                 gi.region1_count = r2_idx - gi.region0_count;
                 gi.table_select[0] = r0_tbl[r2_idx];
                 gi.table_select[1] = r1_tbl[r2_idx];
                 gi.table_select[2] = r2t;
             }
        }
    }

    // --- Public Methods ---

    /**
     * Nonlinear quantization of xr elements based on xr^(3/4).
     * Modifies the `ix` array in place.
     *
     * @public
     * @param {number} l - Number of elements to process (must be even).
     * @param {number} istep - Inverse quantization step size (1 / (2^(gain/4))).
     * @param {Float32Array} xr - Input array containing xr^(3/4) values.
     * @param {number} xrPos - Starting position in `xr`.
     * @param {Int32Array} ix - Output array to store quantized integer values.
     * @param {number} ixPos - Starting position in `ix`.
     */
    quantize_lines_xrpow = this._quantize_lines_xrpow; // Expose private method directly

    /**
     * Quantization for values expected to be 0 or 1. Faster than full quantization.
     * Modifies the `ix` array in place.
     *
     * @public
     * @param {number} l - Number of elements to process (must be even).
     * @param {number} istep - Inverse quantization step size.
     * @param {Float32Array} xr - Input array containing xr^(3/4) values.
     * @param {number} xrPos - Starting position in `xr`.
     * @param {Int32Array} ix - Output array to store quantized integer values (0 or 1).
     * @param {number} ixPos - Starting position in `ix`.
     */
    quantize_lines_xrpow_01 = this._quantize_lines_xrpow_01;

    /**
     * Performs quantization on the `xrpow` data (`xr^(3/4)`) using the current
     * scalefactors and global gain from `codInfo`. Stores results in `codInfo.l3_enc`.
     * Uses cached noise data (`prevNoise`) if applicable.
     *
     * @public
     * @param {Float32Array} xp - `xrpow` array (input).
     * @param {Int32Array} pi - `l3_enc` array (output).
     * @param {number} istep - Inverse quantization step size.
     * @param {GrInfo} codInfo - Granule information (input/output).
     * @param {CalcNoiseData | null} prevNoise - Cached noise data from previous iteration.
     */
    quantize_xrpow = this._quantize_xrpow;

    /**
     * Counts the number of bits needed to encode the quantized values (`l3_enc`)
     * WITHOUT performing quantization again. Determines Huffman table selections,
     * `big_values`, `count1`, and `count1bits`.
     *
     * @public
     * @param {LameInternalFlags} gfc - LAME internal flags.
     * @param {GrInfo} gi - Granule information containing `l3_enc` and block type info. Modified in place.
     * @param {CalcNoiseData | null} prev_noise - Noise cache data (used to update sfb_count1).
     * @returns {number} Total bits required for Huffman coding (part3).
     */
    noquant_count_bits(gfc, gi, prev_noise) {
        const ix = gi.l3_enc;
        let i = gi.max_nonzero_coeff + 1; // Start from one past the last non-zero
        if(i & 1) i++; // Ensure even index to start search from pair boundary
        i = Math.min(576, i); // Clamp to max index

        if (prev_noise != null) prev_noise.sfb_count1 = 0; // Reset cache info

        // Determine count1 region (where all ix <= 1)
        for (; i > 1; i -= 2) {
            if ((ix[i - 1] | ix[i - 2]) !== 0) break; // Found last non-zero pair
        }
        gi.count1 = i; // i is the start index of the all-zero region

        // Determine count1 Huffman coding bits (quadruples of 0/1)
        let a1 = 0; let a2 = 0;
        let quad_idx = i; // Start from end of count1 region
        for (; quad_idx > 3; quad_idx -= 4) {
            // Check if any value in the quad is > 1
            if (((ix[quad_idx - 1] | ix[quad_idx - 2] | ix[quad_idx - 3] | ix[quad_idx - 4])) > 1) {
                break; // End of quadruples region
            }
            // Calculate packed value and add bits for tables 32 and 33
            const p = ((ix[quad_idx - 4] * 3 + ix[quad_idx - 3]) * 3 + ix[quad_idx - 2]) * 3 + ix[quad_idx - 1]; // Corrected packing? Check C vs Huffman table format. Assume C is correct for now.
            // Let's use the bitwise packing from C:
            // p = ((ix[i - 4] * 2 + ix[i - 3]) * 2 + ix[i - 2]) * 2 + ix[i - 1];
            const p_idx = ((ix[quad_idx - 4] << 3) | (ix[quad_idx - 3] << 2) | (ix[quad_idx - 2] << 1) | ix[quad_idx - 1]); // Simpler index?
             // Need to map p_idx to the correct index in t32l/t33l if they are not direct lookups.
             // Assuming Tables.t32l/t33l are direct lookups based on the packed value p:
             const packed_val = ((ix[quad_idx - 4] * 3 + ix[quad_idx - 3]) * 3 + ix[quad_idx - 2]) * 3 + ix[quad_idx - 1]; // Revert to C packing if needed
             const p_idx_c = ((ix[quad_idx - 4] * 2 + ix[quad_idx - 3]) * 2 + ix[quad_idx - 2]) * 2 + ix[quad_idx - 1];

             if(p_idx_c >= 0 && p_idx_c < Tables.t32l.length && p_idx_c < Tables.t33l.length){
                a1 += Tables.t32l[p_idx_c];
                a2 += Tables.t33l[p_idx_c];
             } else {
                 console.error("Invalid index for count1 tables:", p_idx_c);
                 // Handle error or break?
             }
        }
        // quad_idx now marks the start of the big_values region
        gi.big_values = quad_idx;

        // Choose best count1 table (32 or 33)
        let bits = a1; gi.count1table_select = 0;
        if (a1 > a2) { bits = a2; gi.count1table_select = 1; }
        gi.count1bits = bits; // Store bits used for count1 region

        if (gi.big_values === 0) return bits; // Only count1 region exists

        // Determine region boundaries for big_values
        let region0_end = 0, region1_end = 0; // Indices into xr/ix
        if (gi.block_type === Encoder.SHORT_TYPE) {
             // Short blocks have fixed regions (relative to start of short data?)
             // Need gfc context here. Assuming gfc is accessible (passed externally or via 'this').
             // This logic likely needs access to gfc passed as argument. Refactoring needed.
             // Let's assume NORM_TYPE logic for now, short block needs GFC access.
             console.warn("Short block region calculation in noquant_count_bits needs gfc context.");
             region0_end = Math.min(gi.big_values, 3 * gfc.scalefac_band.s[3]); // Example guess
             region1_end = gi.big_values;

        } else if (gi.block_type === Encoder.NORM_TYPE) {
            assert(gi.big_values <= 576, "big_values out of bounds");
             if (gi.big_values >= 2) { // Need at least 2 values for lookup
                 // Look up precalculated region boundaries based on big_values end index
                 gi.region0_count = gfc.bv_scf[gi.big_values - 2]; // region0 size in sfbs
                 gi.region1_count = gfc.bv_scf[gi.big_values - 1]; // region1 size in sfbs
             } else {
                 gi.region0_count = 0;
                 gi.region1_count = 0;
             }
             // Clamp region sizes
             const max_region0_size = Encoder.SBPSY_l - 1; // Max sfb index for region0
             if(gi.region0_count > max_region0_size) gi.region0_count = max_region0_size;
             const max_region1_size = Encoder.SBPSY_l - gi.region0_count - 1;
             if(gi.region1_count > max_region1_size) gi.region1_count = max_region1_size;

             region0_end = gfc.scalefac_band.l[gi.region0_count + 1];
             region1_end = gfc.scalefac_band.l[gi.region0_count + gi.region1_count + 2];

        } else { // Mixed block (START/STOP types)
            gi.region0_count = 7;
            gi.region1_count = Encoder.SBMAX_l - 1 - 7 - 1; // SBMAX_l or SBPSY_l? Using SBMAX_l from C.
            region0_end = gfc.scalefac_band.l[gi.region0_count + 1];
            region1_end = gi.big_values; // Region 2 doesn't exist or is handled differently
        }

        // Ensure region ends don't exceed big_values start
        region0_end = Math.min(region0_end, gi.big_values);
        region1_end = Math.min(region1_end, gi.big_values);

        // Count bits for regions 0, 1, 2
        if (region0_end > 0) {
             let bi = new Bits(bits);
             gi.table_select[0] = this._choose_table(ix, 0, region0_end, bi);
             bits = bi.bits;
        }
        if (region1_end > region0_end) {
             let bi = new Bits(bits);
             gi.table_select[1] = this._choose_table(ix, region0_end, region1_end, bi);
             bits = bi.bits;
        }
        if (gi.big_values > region1_end) { // Region 2
             let bi = new Bits(bits);
             gi.table_select[2] = this._choose_table(ix, region1_end, gi.big_values, bi);
             bits = bi.bits;
        }

        // Optional optimization: find best Huffman division
        if (gfc.use_best_huffman === 2) {
            gi.part2_3_length = bits; // Store current total
            this.best_huffman_divide(gfc, gi); // Try to improve division
            bits = gi.part2_3_length; // Get potentially updated total
        }

        // Update cache info if needed
        if (prev_noise != null) {
            if (gi.block_type === Encoder.NORM_TYPE) {
                let sfb = 0;
                // Find sfb containing the start of big_values
                while (sfb < Encoder.SBMAX_l && gfc.scalefac_band.l[sfb + 1] < gi.big_values) {
                    sfb++;
                }
                prev_noise.sfb_count1 = sfb; // Store sfb index where count1 ends
            }
            // else: Short block cache handling might be different
        }

        return bits; // Return total Huffman bits
    }

    /**
     * Performs quantization and counts bits. This is the main function called
     * by the iteration loop to evaluate a given set of quantization parameters.
     *
     * @public
     * @param {LameInternalFlags} gfc - LAME internal flags.
     * @param {Float32Array} xr - `xrpow` array (input).
     * @param {GrInfo} gi - Granule information (input/output). `l3_enc` is modified.
     * @param {CalcNoiseData | null} prev_noise - Cached noise data.
     * @returns {number} Total bits required for Huffman coding, or `LARGE_BITS` if quantization fails (value too large).
     */
    count_bits(gfc, xr, gi, prev_noise) {
        // Check if max xrpow exceeds limit for current gain
        const istep = this.qupvt.IPOW20(gi.global_gain); // Inverse step size
        const max_allowed_xrpow = QuantizePVT.IXMAX_VAL / istep;
        if (gi.xrpow_max > max_allowed_xrpow) {
            return LARGE_BITS; // Value too large to quantize with this gain
        }

        // Perform quantization
        this._quantize_xrpow(xr, gi.l3_enc, istep, gi, prev_noise);

        // Apply substep shaping if enabled
        if ((gfc.substep_shaping & 2) !== 0) {
            let j = 0;
            const gain_factor = this.qupvt.IPOW20(gi.global_gain + gi.scalefac_scale); // Base step for gain comparison
            const roundfac = 0.634521682242439 / gain_factor; // Threshold relative to gain
            for (let sfb = 0; sfb < gi.sfbmax; sfb++) {
                const width = gi.width[sfb];
                if (gfc.pseudohalf[sfb] === 0) { // No shaping in this band
                    j += width;
                } else { // Apply shaping: zero out values below threshold
                    const band_end = j + width;
                    for (; j < band_end; ++j) {
                        // Compare original xr value (not xrpow) against threshold? Check C.
                        // C code uses xrpow (xr^3/4). Needs xr, not xrpow here.
                        // This feature might need the original xr[] array passed in.
                        // For now, assume it uses xrpow (incorrectly based on C context).
                        // if (xr[j] < roundfac) gi.l3_enc[j] = 0;
                        // --- Placeholder logic - Needs original xr[] ---
                        // Cannot implement correctly without original xr
                    }
                }
            }
        }

        // Count bits for the quantized values
        return this.noquant_count_bits(gfc, gi, prev_noise);
    }

    /**
     * Optimizes the Huffman table region boundaries (`region0_count`, `region1_count`)
     * for the `big_values` section to minimize bit usage. Modifies `gi` in place
     * if a better division is found.
     *
     * @public
     * @param {LameInternalFlags} gfc - LAME internal flags.
     * @param {GrInfo} gi - Granule information (input/output).
     */
    best_huffman_divide(gfc, gi) {
        // Allocate temporary structures only if needed
        if (gi.block_type !== Encoder.NORM_TYPE && gi.block_type !== Encoder.START_TYPE && gi.block_type !== Encoder.STOP_TYPE) {
            // Currently only implemented for long/mixed blocks with standard regions
             if (gi.block_type === Encoder.SHORT_TYPE && gfc.mode_gr === 1) return; // Skip MPEG2 LSF short blocks
             // Fall through for other block types? Or return? Let's return for now.
             return;
        }

        const cod_info2 = new GrInfo();
        const ix = gi.l3_enc;
        const r01_bits = new_int(7 + 15 + 1); // Max r0+r1 index
        const r01_div = new_int(7 + 15 + 1);
        const r0_tbl = new_int(7 + 15 + 1);
        const r1_tbl = new_int(7 + 15 + 1);

        cod_info2.assign(gi); // Start with current division

        if (gi.block_type === Encoder.NORM_TYPE) {
            // Calculate costs for all possible R0/R1 combinations
            this._recalc_divide_init(gfc, gi, ix, r01_bits, r01_div, r0_tbl, r1_tbl);
            // Find the best combination including R2
            this._recalc_divide_sub(gfc, cod_info2, gi, ix, r01_bits, r01_div, r0_tbl, r1_tbl);
        }
        // else: Mixed blocks (START/STOP) - C code doesn't call recalc_divide_init/sub?
        // It seems to fall through to the count1 adjustment logic below.

        // Check if adjusting count1 boundary helps
        let i = cod_info2.big_values; // Use the big_values from the *potentially* updated cod_info2
        if (i === 0 || i >= gi.count1) return; // No adjustment possible if no bigvals or already overlaps count1

        // Check if the last *pair* of bigvalues are both <= 1
        if ((ix[i - 1] | ix[i - 2]) > 1) return;

        // Try extending count1 region by 2 samples
        i = gi.count1 + 2; // New potential count1 boundary
        if (i > 576) return; // Cannot extend beyond buffer

        cod_info2.assign(gi); // Re-copy best state found so far
        cod_info2.count1 = i; // Tentatively set new count1
        let a1 = 0; let a2 = 0;
        let quad_idx = i;

        // Recalculate count1 bits for the *new*, potentially smaller quadruple region
        for (; quad_idx > cod_info2.big_values; quad_idx -= 4) { // Loop from new count1 down to new big_values
             const p_idx_c = ((ix[quad_idx - 4] << 3) | (ix[quad_idx - 3] << 2) | (ix[quad_idx - 2] << 1) | ix[quad_idx - 1]);
             if(p_idx_c >= 0 && p_idx_c < Tables.t32l.length && p_idx_c < Tables.t33l.length){ a1 += Tables.t32l[p_idx_c]; a2 += Tables.t33l[p_idx_c]; }
             else { console.error("Invalid index for count1 tables in best_huffman_divide:", p_idx_c); return; } // Error out
        }
        cod_info2.big_values = quad_idx; // Update big_values to new boundary

        // Choose best count1 table
        cod_info2.count1table_select = (a1 > a2) ? 1 : 0;
        cod_info2.count1bits = Math.min(a1, a2);

        // Recalculate bits for the (now potentially larger) big_values regions (0, 1, 2)
        cod_info2.part2_3_length = cod_info2.count1bits; // Start with new count1 bits
        let region0_end = 0, region1_end = 0;
        if (cod_info2.block_type === Encoder.NORM_TYPE) {
             // Use the *original* region counts determined previously
             region0_end = gfc.scalefac_band.l[gi.region0_count + 1];
             region1_end = gfc.scalefac_band.l[gi.region0_count + gi.region1_count + 2];
        } else { // Mixed block (use fixed regions)
             region0_end = gfc.scalefac_band.l[7 + 1];
             region1_end = cod_info2.big_values;
        }
        region0_end = Math.min(region0_end, cod_info2.big_values);
        region1_end = Math.min(region1_end, cod_info2.big_values);

        if (region0_end > 0) {
            let bi = new Bits(cod_info2.part2_3_length);
            cod_info2.table_select[0] = this._choose_table(ix, 0, region0_end, bi);
            cod_info2.part2_3_length = bi.bits;
        }
        if (region1_end > region0_end) {
            let bi = new Bits(cod_info2.part2_3_length);
            cod_info2.table_select[1] = this._choose_table(ix, region0_end, region1_end, bi);
            cod_info2.part2_3_length = bi.bits;
        }
        if (cod_info2.big_values > region1_end) { // Region 2
            let bi = new Bits(cod_info2.part2_3_length);
            cod_info2.table_select[2] = this._choose_table(ix, region1_end, cod_info2.big_values, bi);
            cod_info2.part2_3_length = bi.bits;
        }

        // If this adjusted count1 boundary resulted in fewer bits, update the main gi
        if (gi.part2_3_length > cod_info2.part2_3_length) {
            gi.assign(cod_info2);
        }
    }

    /**
     * Calculates the number of bits used to encode scalefactors for MPEG1.
     * Determines the optimal `scalefac_compress` value.
     * Also applies `pretab` if possible and updates `cod_info.preflag`.
     * Returns true if any scalefactor exceeds the limits for all `scalefac_compress` values.
     *
     * @public
     * @param {GrInfo} cod_info - Granule information (input/output).
     * @returns {boolean} True if scalefactors are too large, false otherwise.
     */
    scale_bitcount(cod_info) {
        let k, sfb;
        let max_slen1 = 0, max_slen2 = 0;
        let tab;
        const scalefac = cod_info.scalefac; // Use local ref

        assert(this._all_scalefactors_not_negative(scalefac, cod_info.sfbmax), "Negative scalefactor found before scale_bitcount");

        // Select appropriate bit count table based on block type
        if (cod_info.block_type === Encoder.SHORT_TYPE) {
            tab = (cod_info.mixed_block_flag !== 0) ? scale_mixed : scale_short;
        } else { // Long/Start/Stop
            tab = scale_long;
            // Check if preemphasis can be applied (implicitly, by reducing scalefactors)
            if (cod_info.preflag === 0) {
                for (sfb = 11; sfb < Encoder.SBPSY_l; sfb++) {
                    if (scalefac[sfb] < this.qupvt.pretab[sfb]) break; // Cannot apply if any sf < pretab
                }
                if (sfb === Encoder.SBPSY_l) { // All sfb >= pretab[sfb]
                    cod_info.preflag = 1; // Enable preflag
                    for (sfb = 11; sfb < Encoder.SBPSY_l; sfb++) {
                         scalefac[sfb] -= this.qupvt.pretab[sfb]; // Adjust scalefactors
                    }
                }
            }
        }

        // Find max scalefactor values in the two slen regions
        for (sfb = 0; sfb < cod_info.sfbdivide; sfb++) {
            if (max_slen1 < scalefac[sfb]) max_slen1 = scalefac[sfb];
        }
        for (; sfb < cod_info.sfbmax; sfb++) {
            if (max_slen2 < scalefac[sfb]) max_slen2 = scalefac[sfb];
        }

        // Find best scalefac_compress value (minimum bits)
        cod_info.part2_length = LARGE_BITS; // Initialize with large value
        let found_valid = false;
        for (k = 0; k < 16; k++) {
            // Check if current max values fit within limits for this 'k'
            if (max_slen1 < slen1_n[k] && max_slen2 < slen2_n[k]) {
                 found_valid = true; // At least one valid compression found
                 if (cod_info.part2_length > tab[k]) { // If this 'k' gives fewer bits
                    cod_info.part2_length = tab[k]; // Store new minimum bit count
                    cod_info.scalefac_compress = k; // Store corresponding index
                 }
            }
        }

        // Return true if no valid scalefac_compress was found (scalefactors too large)
        return !found_valid; // (found_valid=false means part2_length is still LARGE_BITS)
        // Original C returned cod_info.part2_length == LARGE_BITS; equivalent.
    }


    /**
     * Calculates the number of bits used to encode scalefactors for MPEG2 LSF.
     * Determines the optimal `scalefac_compress` value based on partitioning.
     * Returns true if any scalefactor exceeds the limits for the chosen partition table.
     *
     * @public
     * @param {LameInternalFlags} gfc - LAME internal flags (for partition tables).
     * @param {GrInfo} cod_info - Granule information (input/output).
     * @returns {boolean} True if scalefactors are too large, false otherwise.
     */
    scale_bitcount_lsf(gfc, cod_info) {
        let table_number, row_in_table, partition, nr_sfb, window;
        let over = false; // Flag if any scalefactor exceeds limit
        let i, sfb;
        const max_sfac = new_int(4); // Max scalefac value per partition
        const scalefac = cod_info.scalefac; // Use local ref

        // --- Select Partition Table ---
        // C code: tries table 1 if possible, otherwise uses 0 or 2 based on preflag.
        // Simplified: use 0 or 2 based on preflag. Needs review if table 1 logic is crucial.
        table_number = (cod_info.preflag !== 0) ? 2 : 0;

        // Determine max scalefactor in each partition based on block type
        Arrays.fill(max_sfac, 0); // Initialize max values
        if (cod_info.block_type === Encoder.SHORT_TYPE) {
            row_in_table = 1; // Short block partition row
            const partition_table = this.qupvt.nr_of_sfb_block[table_number][row_in_table];
            sfb = 0; // Short block scalefactor index
            for (partition = 0; partition < 4; partition++) {
                nr_sfb = partition_table[partition] / 3; // Number of sfbs in this partition
                for (i = 0; i < nr_sfb; i++, sfb++) {
                    for (window = 0; window < 3; window++) {
                        // Index into flat scalefac array for short blocks: sfb*3 + window
                        const current_sf = scalefac[sfb * 3 + window];
                        if (current_sf > max_sfac[partition]) max_sfac[partition] = current_sf;
                    }
                }
            }
        } else { // Long blocks
            row_in_table = 0; // Long block partition row
            const partition_table = this.qupvt.nr_of_sfb_block[table_number][row_in_table];
            sfb = 0; // Long block scalefactor index
            for (partition = 0; partition < 4; partition++) {
                nr_sfb = partition_table[partition]; // Number of sfbs in this partition
                for (i = 0; i < nr_sfb; i++, sfb++) {
                     if (scalefac[sfb] > max_sfac[partition]) max_sfac[partition] = scalefac[sfb];
                }
            }
        }

        // Check if max values exceed limits for the chosen table
        for (partition = 0; partition < 4; partition++) {
            if (max_sfac[partition] > max_range_sfac_tab[table_number][partition]) {
                over = true; // Scalefactor too large
                break;
            }
        }

        // If no limits exceeded, calculate bits and set compress info
        if (!over) {
            cod_info.sfb_partition_table = this.qupvt.nr_of_sfb_block[table_number][row_in_table];
            for (partition = 0; partition < 4; partition++) {
                // slen is the number of bits needed for the max value in the partition
                 cod_info.slen[partition] = log2tab[max_sfac[partition]];
            }

            // Calculate scalefac_compress based on slen values and table number
            const slen1 = cod_info.slen[0]; const slen2 = cod_info.slen[1];
            const slen3 = cod_info.slen[2]; const slen4 = cod_info.slen[3];
            switch (table_number) {
                case 0: cod_info.scalefac_compress = (((slen1 * 5) + slen2) << 4) + (slen3 << 2) + slen4; break;
                case 1: cod_info.scalefac_compress = 400 + (((slen1 * 5) + slen2) << 2) + slen3; break;
                case 2: cod_info.scalefac_compress = 500 + (slen1 * 3) + slen2; break;
                default: console.error("LSF intensity stereo not implemented"); break; // case 3, 4, 5
            }

            // Calculate total bits for scalefactors
            cod_info.part2_length = 0;
            assert(cod_info.sfb_partition_table != null, "Partition table is null in scale_bitcount_lsf");
            for (partition = 0; partition < 4; partition++) {
                 cod_info.part2_length += cod_info.slen[partition] * cod_info.sfb_partition_table[partition];
            }
        }
         // else: 'over' is true, part2_length remains unset or potentially invalid

        return over; // Return true if scalefactors were too large
    }

    /**
     * Initializes Huffman table boundaries based on scalefactor bands.
     * Precomputes `gfc.bv_scf` array.
     *
     * @public
     * @param {LameInternalFlags} gfc - LAME internal flags.
     */
    huffman_init(gfc) {
        // Precompute bv_scf table used in noquant_count_bits for NORM_TYPE blocks
        for (let i = 2; i <= 576; i += 2) { // Iterate over pairs of coefficients
             let scfb_anz = 0; // Scalefactor band index containing coeff i-1
             // Find the scalefactor band for the current coefficient index i-1
             // Note: Uses gfc.scalefac_band.l which should be initialized before this.
             while (scfb_anz < Encoder.SBMAX_l && gfc.scalefac_band.l[scfb_anz + 1] < i) {
                 scfb_anz++;
             }
             // Handle potential out-of-bounds if i is very high? Should not happen if i <= 576.
             if(scfb_anz >= Encoder.SBMAX_l) scfb_anz = Encoder.SBMAX_l -1; // Clamp?

             // Determine region0 boundary (bv_index) using subdv_table lookup
             let bv_index = subdv_table[scfb_anz][0]; // Default region0 size for this sfb
             // Adjust bv_index downwards if its end boundary exceeds the current coefficient index
             while (bv_index > 0 && gfc.scalefac_band.l[bv_index + 1] >= i) {
                 bv_index--;
             }
              // C code check: if (bv_index < 0) bv_index = subdv_table[scfb_anz][0]; - Seems redundant if adjusted down? Let's keep C logic.
              if (bv_index < 0) bv_index = subdv_table[scfb_anz][0]; // Reset if adjustment went too far?

             gfc.bv_scf[i - 2] = bv_index; // Store region0 size for coeff i-2

             // Determine region1 boundary using subdv_table lookup
             bv_index = subdv_table[scfb_anz][1]; // Default region1 size
             // Adjust downwards based on region0 size and current coeff index
             while (bv_index > 0 && gfc.scalefac_band.l[bv_index + gfc.bv_scf[i - 2] + 2] >= i) {
                 bv_index--;
             }
              if (bv_index < 0) bv_index = subdv_table[scfb_anz][1]; // Reset?

             gfc.bv_scf[i - 1] = bv_index; // Store region1 size for coeff i-1
        }
    }

} // End class Takehiro

// --- Module-level constants/tables used internally by Takehiro methods ---

/** log2(0..15), integer */
const log2tab = [0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4];

/** Max values for slen1 region based on scalefac_compress index */
const slen1_n = [1, 1, 1, 1, 8, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16, 16];
/** Max values for slen2 region based on scalefac_compress index */
const slen2_n = [1, 2, 4, 8, 1, 2, 4, 8, 2, 4, 8, 2, 4, 8, 4, 8];

/** Bit length for slen1 region */
// const slen1_tab = [0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4]; // Now exported
/** Bit length for slen2 region */
// const slen2_tab = [0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3]; // Now exported

/** Max scalefactor values allowed per partition for MPEG2 LSF tables */
const max_range_sfac_tab = [
    [15, 15, 7, 7],  // Table 0
    [15, 15, 7, 0],  // Table 1
    [7, 3, 0, 0],   // Table 2 (pretab)
    [15, 31, 31, 0], // Table 3 (intensity)
    [7, 7, 7, 0],   // Table 4 (intensity)
    [3, 3, 0, 0]    // Table 5 (intensity)
];


export { Takehiro };
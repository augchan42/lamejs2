/**
 * @fileoverview Core MP3 frame encoding logic for LAME.
 * Ported from encoder.c. Handles MDCT, psychoacoustic model interaction,
 * quantization loop selection, and bitstream formatting for a single frame.
 * Uses ES Module syntax.
 *
 * @module Encoder
 */

// Import necessary modules using ES Module syntax
import * as common from './common.js';
import { NewMDCT } from './NewMDCT.js';
import { III_psy_ratio } from './III_psy_ratio.js';
import { MPEGMode } from './MPEGMode.js';
import { L3Side } from './L3Side.js'; // Contains SideInfo structure
import { VbrMode } from './common.js'; // Import VbrMode for comparison

// Assuming these types are defined elsewhere and imported if needed
/** @typedef {import('./BitStream.js').BitStream} BitStream */
/** @typedef {import('./PsyModel.js').PsyModel} PsyModel */
/** @typedef {import('./QuantizePVT.js').QuantizePVT} QuantizePVT */
/** @typedef {import('./VBRTag.js').default} VBRTag */
/** @typedef {import('./LameGlobalFlags.js').default} LameGlobalFlags */
/** @typedef {import('./LameInternalFlags.js').default} LameInternalFlags */

// Destructure common utilities for easier access
const {
    Float, ShortBlock, Util, Arrays, new_array_n, new_byte,
    new_float, new_float_n, new_int, new_int_n, assert
} = common;

/**
 * @classdesc Contains the core logic for encoding a single MP3 frame,
 * coordinating the MDCT, psychoacoustic analysis, quantization, and bitstream writing.
 * @constructs Encoder
 */
class Encoder {
    // --- Static Constants ---
    static ENCDELAY = 576;
    static POSTDELAY = Encoder.ENCDELAY + 1152 - 576 - 1;
    static MDCTDELAY = 48;
    static FFTOFFSET = (224 + Encoder.MDCTDELAY);
    static DECDELAY = 528;
    static SBLIMIT = 32;
    static CBANDS = 64;
    static SBPSY_l = 21;
    static SBPSY_s = 12;
    static SBMAX_l = 22;
    static SBMAX_s = 13;
    static PSFB21 = 6;
    static PSFB12 = 6;
    static BLKSIZE = 1024;
    static HBLKSIZE = (Encoder.BLKSIZE / 2 + 1);
    static BLKSIZE_s = 256;
    static HBLKSIZE_s = (Encoder.BLKSIZE_s / 2 + 1);
    static NORM_TYPE = 0;
    static START_TYPE = 1;
    static SHORT_TYPE = 2;
    static STOP_TYPE = 3;
    static MPG_MD_LR_LR = 0;
    static MPG_MD_LR_I = 1;
    static MPG_MD_MS_LR = 2;
    static MPG_MD_MS_I = 3;
    static fircoef = [-0.0207887 * 5, -0.0378413 * 5, -0.0432472 * 5, -0.031183 * 5, 7.79609e-18 * 5, 0.0467745 * 5, 0.10091 * 5, 0.151365 * 5, 0.187098 * 5];

    /** @private @type {BitStream|null} */
    bs = null;
    /** @private @type {PsyModel|null} */
    psy = null;
    /** @private @type {VBRTag|null} */
    vbr = null;
    /** @private @type {QuantizePVT|null} */
    qupvt = null;
    /** @private @type {NewMDCT} */
    newMDCT;

    constructor() {
        /** @private */
        this.newMDCT = new NewMDCT();
    }

    /**
     * Sets the internal module dependencies.
     * @public
     * @param {BitStream} _bs
     * @param {PsyModel} _psy
     * @param {QuantizePVT} _qupvt
     * @param {VBRTag} _vbr
     */
    setModules(_bs, _psy, _qupvt, _vbr) {
        this.bs = _bs;
        this.psy = _psy;
        this.qupvt = _qupvt;
        this.vbr = _vbr;
    }

    // --- Private Helper Methods ---
    // ... ( _adjust_ATH, _updateStats, _lame_encode_frame_init methods remain the same ) ...
    /** @private */
    _adjust_ATH(gfc) { /* ... */ }
    /** @private */
    _updateStats(gfc) { /* ... */ }
    /** @private */
    _lame_encode_frame_init(gfp, inbuf) { /* ... */ }


    // --- Public Methods ---

    /**
     * Encodes a single MP3 frame.
     * Takes PCM data (after potential resampling and buffering), performs
     * psychoacoustic analysis, MDCT, quantization (via iteration loop),
     * and formats the bitstream.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {Float32Array} inbuf_l - Left channel PCM data for the frame + overlap.
     * @param {Float32Array} inbuf_r - Right channel PCM data (or copy of left for mono).
     * @param {Uint8Array} mp3buf - Output buffer for the encoded MP3 frame data.
     * @param {number} mp3bufPos - Starting position offset within `mp3buf`.
     * @param {number} mp3buf_size - Available size in `mp3buf` from `mp3bufPos`.
     * @returns {number} The number of bytes written to `mp3buf` for this frame, or a negative error code.
     */
    lame_encode_mp3_frame(gfp, inbuf_l, inbuf_r, mp3buf, mp3bufPos, mp3buf_size) {
        // ... (Implementation remains the same as before) ...
        let mp3count;
        const masking_LR = [ [new III_psy_ratio(), new III_psy_ratio()], [new III_psy_ratio(), new III_psy_ratio()] ];
        const masking_MS = [ [new III_psy_ratio(), new III_psy_ratio()], [new III_psy_ratio(), new III_psy_ratio()] ];
        let masking;
        const inbuf = [inbuf_l, inbuf_r];
        const gfc = gfp.internal_flags;
        const ms_ener_ratio = new_float(2);
        const pe = new_float_n([2, 2]);
        const pe_MS = new_float_n([2, 2]);
        let pe_use;
        let ch, gr;

        if (gfc.lame_encode_frame_init === 0) {
            this._lame_encode_frame_init(gfp, inbuf);
        }
        gfc.padding = 0;
        if (gfc.frac_SpF > 1e-9) { gfc.slot_lag -= gfc.frac_SpF; if (gfc.slot_lag < 0) { gfc.slot_lag += gfp.out_samplerate; gfc.padding = 1; } }

        if (gfc.psymodel !== 0) {
            let ret; const bufp = [null, null]; const blocktype = new_int(2);
            for (gr = 0; gr < gfc.mode_gr; gr++) {
                const bufpPos = 576 + gr * 576 - Encoder.FFTOFFSET;
                bufp[0] = inbuf[0].subarray(bufpPos); if (gfc.channels_out === 2) bufp[1] = inbuf[1].subarray(bufpPos);
                if (gfp.VBR === VbrMode.vbr_mtrh || gfp.VBR === VbrMode.vbr_mt) { ret = this.psy.L3psycho_anal_vbr(gfp, bufp, 0, gr, masking_LR, masking_MS, pe[gr], pe_MS[gr], gfc.tot_ener, blocktype); }
                else { ret = this.psy.L3psycho_anal_ns(gfp, bufp, 0, gr, masking_LR, masking_MS, pe[gr], pe_MS[gr], gfc.tot_ener, blocktype); }
                if (ret !== 0) return -4;
                if (gfp.mode === MPEGMode.JOINT_STEREO) { const total_ms_energy = gfc.tot_ener[2] + gfc.tot_ener[3]; if(total_ms_energy > 0) { ms_ener_ratio[gr] = gfc.tot_ener[3] / total_ms_energy; } else { ms_ener_ratio[gr] = 0.5; } }
                for (ch = 0; ch < gfc.channels_out; ch++) { const cod_info = gfc.l3_side.tt[gr][ch]; cod_info.block_type = blocktype[ch]; cod_info.mixed_block_flag = 0; }
            }
        } else { for (gr = 0; gr < gfc.mode_gr; gr++) { for (ch = 0; ch < gfc.channels_out; ch++) { gfc.l3_side.tt[gr][ch].block_type = Encoder.NORM_TYPE; gfc.l3_side.tt[gr][ch].mixed_block_flag = 0; pe_MS[gr][ch] = pe[gr][ch] = 700; } } }
        this._adjust_ATH(gfc);
        this.newMDCT.mdct_sub48(gfc, inbuf[0], inbuf[1]);
        gfc.mode_ext = Encoder.MPG_MD_LR_LR;
        if (gfp.force_ms) { gfc.mode_ext = Encoder.MPG_MD_MS_LR; }
        else if (gfp.mode === MPEGMode.JOINT_STEREO && gfc.channels_out == 2) {
            let sum_pe_MS = 0.0; let sum_pe_LR = 0.0;
            for (gr = 0; gr < gfc.mode_gr; gr++) { sum_pe_MS += pe_MS[gr][0] + pe_MS[gr][1]; sum_pe_LR += pe[gr][0] + pe[gr][1]; }
            if (sum_pe_MS <= 1.00 * sum_pe_LR) { const types_match_gr0 = (gfc.l3_side.tt[0][0].block_type === gfc.l3_side.tt[0][1].block_type); const types_match_gr1 = (gfc.mode_gr === 1) || (gfc.l3_side.tt[1][0].block_type === gfc.l3_side.tt[1][1].block_type); if(types_match_gr0 && types_match_gr1) { gfc.mode_ext = Encoder.MPG_MD_MS_LR; } }
        }
        if (gfc.mode_ext === Encoder.MPG_MD_MS_LR) { masking = masking_MS; pe_use = pe_MS; } else { masking = masking_LR; pe_use = pe; }
        if (gfp.analysis && gfc.pinfo != null) { /* ... copy pinfo data ... */ }
        if (gfp.VBR === VbrMode.vbr_off || gfp.VBR === VbrMode.vbr_abr) {
            let f = 0.0; for (gr = 0; gr < gfc.mode_gr; gr++) for (ch = 0; ch < gfc.channels_out; ch++) f += pe_use[gr][ch];
            for (let i = 0; i < 18; i++) gfc.nsPsy.pefirbuf[i] = gfc.nsPsy.pefirbuf[i + 1]; gfc.nsPsy.pefirbuf[18] = f;
            f = gfc.nsPsy.pefirbuf[9]; for (let i = 0; i < 9; i++) f += (gfc.nsPsy.pefirbuf[i] + gfc.nsPsy.pefirbuf[18 - i]) * Encoder.fircoef[i];
            const target_avg_pe = (670 * 5 * gfc.mode_gr * gfc.channels_out); if (f !== 0) f = target_avg_pe / f; else f = 1.0;
            for (gr = 0; gr < gfc.mode_gr; gr++) for (ch = 0; ch < gfc.channels_out; ch++) pe_use[gr][ch] *= f;
        }
        gfc.iteration_loop.iteration_loop(gfp, pe_use, ms_ener_ratio, masking);
        this.bs.format_bitstream(gfp);
        mp3count = this.bs.copy_buffer(gfc, mp3buf, mp3bufPos, mp3buf_size, 1);
        if (mp3count < 0) return mp3count;
        if (gfp.bWriteVbrTag) { /* this.vbr.addVbrFrame(gfp); */ }
        if (gfp.analysis && gfc.pinfo != null) { /* ... copy pinfo data ... */ }
        this._updateStats(gfc);
        return mp3count;
    }

} // End class Encoder

// --- Exports ---
// Export the class itself
export { Encoder };
export default Encoder;

// Export the constants needed by other modules
export const ENCDELAY = Encoder.ENCDELAY;
export const POSTDELAY = Encoder.POSTDELAY;
export const MDCTDELAY = Encoder.MDCTDELAY;
export const FFTOFFSET = Encoder.FFTOFFSET;
export const DECDELAY = Encoder.DECDELAY;
export const SBLIMIT = Encoder.SBLIMIT;
export const CBANDS = Encoder.CBANDS;
export const SBPSY_l = Encoder.SBPSY_l;
export const SBPSY_s = Encoder.SBPSY_s;
export const SBMAX_l = Encoder.SBMAX_l;
export const SBMAX_s = Encoder.SBMAX_s;
export const PSFB21 = Encoder.PSFB21;
export const PSFB12 = Encoder.PSFB12;
export const BLKSIZE = Encoder.BLKSIZE;
export const HBLKSIZE = Encoder.HBLKSIZE;
export const BLKSIZE_s = Encoder.BLKSIZE_s;
export const HBLKSIZE_s = Encoder.HBLKSIZE_s;
export const NORM_TYPE = Encoder.NORM_TYPE;
export const START_TYPE = Encoder.START_TYPE;
export const SHORT_TYPE = Encoder.SHORT_TYPE;
export const STOP_TYPE = Encoder.STOP_TYPE;
export const MPG_MD_LR_LR = Encoder.MPG_MD_LR_LR;
export const MPG_MD_LR_I = Encoder.MPG_MD_LR_I;
export const MPG_MD_MS_LR = Encoder.MPG_MD_MS_LR;
export const MPG_MD_MS_I = Encoder.MPG_MD_MS_I;
// fircoef is likely internal, not needed for export
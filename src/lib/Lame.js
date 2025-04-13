/**
 * @fileoverview Main LAME encoder class definition and initialization logic.
 * Ported from lame.c. Provides the top-level API for initializing the encoder,
 * setting parameters, and encoding audio frames.
 * Uses ES Module syntax.
 *
 * @module Lame
 */

// Import necessary modules using ES Module syntax
import * as common from './common.js';
import { PsyModel } from './PsyModel.js';
import { LameGlobalFlags } from './LameGlobalFlags.js';
import { LameInternalFlags } from './LameInternalFlags.js';
import { ATH } from './ATH.js';
import { ReplayGain } from './ReplayGain.js';
import { CBRNewIterationLoop } from './CBRNewIterationLoop.js';
// Comment out imports for non-existent classes
// import VBROldIterationLoop from './VBROldIterationLoop.js'; // Assuming this exists
// import ABRIterationLoop from './ABRIterationLoop.js';   // Assuming this exists
import { BitStream } from './BitStream.js';
import { bitrate_table } from './Tables.js';
import { Encoder } from './Encoder.js';
import { MPEGMode } from './MPEGMode.js';

// Assuming these types are defined elsewhere and imported if needed for full type safety
/** @typedef {import('./GainAnalysis.js').default} GainAnalysis */
/** @typedef {import('./Presets.js').default} Presets */
/** @typedef {import('./QuantizePVT.js').QuantizePVT} QuantizePVT */
/** @typedef {import('./Quantize.js').Quantize} Quantize */
/** @typedef {import('./VBRTag.js').default} VBRTag */
/** @typedef {import('./Version.js').default} Version */
/** @typedef {import('./ID3Tag.js').default} ID3Tag */
/** @typedef {import('./MPGLib.js').default} MPGLib */ // Assuming MPGLib module exists
/** @typedef {import('./Takehiro.js').default} Takehiro */ // Assuming Takehiro module exists

// Destructure common utilities for easier access
const {
    VbrMode,
    Float,
    ShortBlock,
    Util,
    Arrays, // Keep for potential fill/sort usage
    new_array_n,
    new_byte,
    // new_double, // Not used
    new_float,
    new_float_n,
    new_int,
    new_int_n,
    new_short_n,
    assert
    // System, // Remove System from destructuring
} = common;

// --- Internal Helper Classes (No JSDoc as requested) ---
/** @private */
class _PSY { // Renamed to match usage
    mask_adjust = 0.0;
    mask_adjust_short = 0.0;
    bo_l_weight = new_float(Encoder.SBMAX_l);
    bo_s_weight = new_float(Encoder.SBMAX_s);
}
/** @private */
class _LowPassHighPass { lowerlimit = 0.0; };
/** @private */
class _BandPass { constructor(bitrate, lPass) { this.lowpass = lPass;} };


/**
 * @classdesc The main LAME MP3 encoder class. Provides methods to initialize
 * the encoder, set parameters, encode PCM data buffers, and finalize the MP3 stream.
 * @constructs Lame
 */
class Lame {
    // --- Static Constants ---
    static LAME_MAXALBUMART = (128 * 1024);
    static V9 = 410; static V8 = 420; static V7 = 430; static V6 = 440; static V5 = 450;
    static V4 = 460; static V3 = 470; static V2 = 480; static V1 = 490; static V0 = 500;
    static R3MIX = 1000; static STANDARD = 1001; static EXTREME = 1002;
    static INSANE = 1003; static STANDARD_FAST = 1004; static EXTREME_FAST = 1005;
    static MEDIUM = 1006; static MEDIUM_FAST = 1007;
    static LAME_MAXMP3BUFFER = (16384 + Lame.LAME_MAXALBUMART);
    static BLACKSIZE = 32; // Maximum filter length + 1

    /** @private @type {GainAnalysis|null} */
    ga = null;
    /** @private @type {BitStream|null} */
    bs = null;
    /** @private @type {Presets|null} */
    p = null;
    /** @private @type {QuantizePVT|null} */
    qupvt = null;
    /** @private @type {Quantize|null} */
    qu = null;
    /** @private @type {PsyModel} */
    psy = new PsyModel();
    /** @private @type {VBRTag|null} */
    vbr = null;
    /** @private @type {Version|null} */
    ver = null;
    /** @private @type {ID3Tag|null} */
    id3 = null;
    /** @private @type {MPGLib|null} */
    mpglib = null;
    /** @public @type {Encoder} Core encoder logic. */
    enc = new Encoder();
    /** @private Internal constant */
    _LAME_ID = 0xFFF88E3B;

    // Internal helper classes (defined outside or passed in)
    /** @private */ _InOut = class InOut { constructor() { this.n_in = 0; this.n_out = 0; } };
    /** @private */ _NumUsed = class NumUsed { constructor() { this.num_used = 0; } };


    constructor() {
        // Modules are set externally via setModules
    }

    /**
     * Sets the internal module dependencies for the Lame instance.
     * Must be called after instantiating Lame and its dependencies.
     *
     * @public
     * @param {GainAnalysis} _ga - Gain Analysis module instance.
     * @param {BitStream} _bs - BitStream module instance.
     * @param {Presets} _p - Presets module instance.
     * @param {QuantizePVT} _qupvt - QuantizePVT module instance.
     * @param {Quantize} _qu - Quantize module instance.
     * @param {VBRTag} _vbr - VBRTag module instance.
     * @param {Version} _ver - Version module instance.
     * @param {ID3Tag} _id3 - ID3Tag module instance.
     * @param {MPGLib} _mpglib - MPGLib module instance.
     */
    setModules(_ga, _bs, _p, _qupvt, _qu, _vbr, _ver, _id3, _mpglib) {
        this.ga = _ga;
        this.bs = _bs;
        this.p = _p;
        this.qupvt = _qupvt;
        this.qu = _qu;
        this.vbr = _vbr;
        this.ver = _ver;
        this.id3 = _id3;
        this.mpglib = _mpglib;
        // Pass necessary modules down to Encoder
        this.enc.setModules(this.bs, this.psy, this.qupvt, this.vbr);
    }

    // --- Internal Helper Methods (Private JSDoc omitted) ---
    /** @private */
    _lame_init_old(gfp) { /* ... (implementation as before) ... */
        let gfc;
        gfp.class_id = this._LAME_ID;
        gfc = gfp.internal_flags = new LameInternalFlags();
        gfp.mode = MPEGMode.NOT_SET; gfp.original = 1; gfp.in_samplerate = 44100;
        gfp.num_channels = 2; gfp.num_samples = -1; gfp.bWriteVbrTag = true;
        gfp.quality = -1; gfp.short_blocks = null; gfc.subblock_gain = -1;
        gfp.lowpassfreq = 0; gfp.highpassfreq = 0; gfp.lowpasswidth = -1; gfp.highpasswidth = -1;
        gfp.VBR = VbrMode.vbr_off; gfp.VBR_q = 4; gfp.ATHcurve = -1;
        gfp.VBR_mean_bitrate_kbps = 128; gfp.VBR_min_bitrate_kbps = 0; gfp.VBR_max_bitrate_kbps = 0;
        gfp.VBR_hard_min = 0; gfc.VBR_min_bitrate = 1; gfc.VBR_max_bitrate = 14;
        gfp.quant_comp = -1; gfp.quant_comp_short = -1; gfp.msfix = -1.0;
        gfc.resample_ratio = 1.0; gfc.OldValue[0] = 180; gfc.OldValue[1] = 180;
        gfc.CurrentStep[0] = 4; gfc.CurrentStep[1] = 4; gfc.masking_lower = 1.0;
        gfc.nsPsy.attackthre = -1; gfc.nsPsy.attackthre_s = -1; gfp.scale = -1.0;
        gfp.athaa_type = -1; gfp.ATHtype = -1; gfp.athaa_loudapprox = -1;
        gfp.athaa_sensitivity = 0.0; gfp.useTemporal = null; gfp.interChRatio = -1.0;
        gfc.mf_samples_to_encode = Encoder.ENCDELAY + Encoder.POSTDELAY;
        gfp.encoder_padding = 0; gfc.mf_size = Encoder.ENCDELAY - Encoder.MDCTDELAY;
        gfp.findReplayGain = false; gfp.decode_on_the_fly = false;
        gfc.decode_on_the_fly = false; gfc.findReplayGain = false; gfc.findPeakSample = false;
        gfc.RadioGain = 0; gfc.AudiophileGain = 0; gfc.noclipGainChange = 0; gfc.noclipScale = -1.0;
        gfp.preset = 0; gfp.write_id3tag_automatic = true;
        return 0;
     }
    /** @private */
    _filter_coef(x) { if (x > 1.0) return 0.0; if (x <= 0.0) return 1.0; return Math.cos(Math.PI / 2.0 * x); }
    /** @private */
    _nearestBitrateFullIndex(bitrate) { /* ... (implementation as before) ... */
        const full_bitrate_table = [8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320];
        let lower_range = 16, lower_range_kbps = 320; let upper_range = 16, upper_range_kbps = 320;
        for (let b = 0; b < 16; b++) { if (full_bitrate_table[b + 1] >= bitrate) { upper_range_kbps = full_bitrate_table[b + 1]; upper_range = b + 1; lower_range_kbps = full_bitrate_table[b]; lower_range = b; break; } }
        return ((upper_range_kbps - bitrate) > (bitrate - lower_range_kbps)) ? lower_range : upper_range;
     }
    /** @private */
    _optimum_samplefreq(lowpassfreq, input_samplefreq) { /* ... (implementation as before) ... */
        let suggested_samplefreq = 44100;
        if (input_samplefreq >= 48000) suggested_samplefreq = 48000;
        else if (input_samplefreq >= 44100) suggested_samplefreq = 44100;
        else if (input_samplefreq >= 32000) suggested_samplefreq = 32000;
        else if (input_samplefreq >= 24000) suggested_samplefreq = 24000;
        else if (input_samplefreq >= 22050) suggested_samplefreq = 22050;
        else if (input_samplefreq >= 16000) suggested_samplefreq = 16000;
        else if (input_samplefreq >= 12000) suggested_samplefreq = 12000;
        else if (input_samplefreq >= 11025) suggested_samplefreq = 11025;
        else if (input_samplefreq >= 8000) suggested_samplefreq = 8000;
        else suggested_samplefreq = 8000;
        if (lowpassfreq === -1) return suggested_samplefreq;
        if (lowpassfreq <= 3970) suggested_samplefreq = 8000;
        else if (lowpassfreq <= 4510) suggested_samplefreq = 11025;
        else if (lowpassfreq <= 5420) suggested_samplefreq = 12000;
        else if (lowpassfreq <= 7230) suggested_samplefreq = 16000;
        else if (lowpassfreq <= 9970) suggested_samplefreq = 22050;
        else if (lowpassfreq <= 11220) suggested_samplefreq = 24000;
        else if (lowpassfreq <= 15250) suggested_samplefreq = 32000;
        else if (lowpassfreq <= 15960) suggested_samplefreq = 44100;
        else suggested_samplefreq = 48000;
        if (input_samplefreq > suggested_samplefreq) {
            if (input_samplefreq > 44100) return 48000; if (input_samplefreq > 32000) return 44100;
            if (input_samplefreq > 24000) return 32000; if (input_samplefreq > 22050) return 24000;
            if (input_samplefreq > 16000) return 22050; if (input_samplefreq > 12000) return 16000;
            if (input_samplefreq > 11025) return 12000; if (input_samplefreq > 8000) return 11025;
            return 8000;
        }
        return suggested_samplefreq;
    }
    /** @private */
    _SmpFrqIndex(sample_freq, gfp) { /* ... (implementation as before) ... */
        switch (sample_freq) {
            case 44100: gfp.version = 1; return 0; case 48000: gfp.version = 1; return 1; case 32000: gfp.version = 1; return 2;
            case 22050: gfp.version = 0; return 0; case 24000: gfp.version = 0; return 1; case 16000: gfp.version = 0; return 2;
            case 11025: gfp.version = 0; return 0; case 12000: gfp.version = 0; return 1; case 8000: gfp.version = 0; return 2;
            default: gfp.version=0; return -1;
        }
     }
    /** @private */
    _FindNearestBitrate(bRate, version, samplerate) { /* ... (implementation as before) ... */
        if (samplerate < 16000) version = 0; else if (samplerate < 32000) version = 0;
        let nearest_bitrate = bitrate_table[version][1]; let min_diff = Math.abs(nearest_bitrate - bRate);
        for (let i = 2; i <= 14; i++) { const current_rate = bitrate_table[version][i]; if (current_rate > 0) { const diff = Math.abs(current_rate - bRate); if (diff < min_diff) { min_diff = diff; nearest_bitrate = current_rate; } } }
        return nearest_bitrate;
    }
    /** @private */
    _BitrateIndex(bRate, version, samplerate) { /* ... (implementation as before) ... */
        if (samplerate < 16000) version = 0; else if (samplerate < 32000) version = 0;
        for (let i = 1; i <= 14; i++) { if (bitrate_table[version][i] === bRate) return i; } return -1;
     }
    /** @private */
    _optimum_bandwidth(lh, bitrate) { /* ... (implementation as before) ... */
         const freq_map = [ new this._BandPass(8, 2000), new this._BandPass(16, 3700), new this._BandPass(24, 3900), new this._BandPass(32, 5500), new this._BandPass(40, 7000), new this._BandPass(48, 7500), new this._BandPass(56, 10000), new this._BandPass(64, 11000), new this._BandPass(80, 13500), new this._BandPass(96, 15100), new this._BandPass(112, 15600), new this._BandPass(128, 17000), new this._BandPass(160, 17500), new this._BandPass(192, 18600), new this._BandPass(224, 19400), new this._BandPass(256, 19700), new this._BandPass(320, 20500) ];
         const table_index = this._nearestBitrateFullIndex(bitrate); lh.lowerlimit = freq_map[table_index].lowpass;
     }
    /** @private */
    _lame_init_params_ppflt(gfp) { /* ... (implementation as before) ... */
        const gfc = gfp.internal_flags;
        let lowpass_band = 32; let highpass_band = -1;
        if (gfc.lowpass1 > 0) { let minband = 999; for (let band = 0; band <= 31; band++) { const freq = band / 31.0; if (freq >= gfc.lowpass2) lowpass_band = Math.min(lowpass_band, band); if (gfc.lowpass1 < freq && freq < gfc.lowpass2) minband = Math.min(minband, band); } gfc.lowpass1 = (minband === 999) ? (lowpass_band - 0.75) / 31.0 : (minband - 0.75) / 31.0; gfc.lowpass2 = lowpass_band / 31.0; }
        if (gfc.highpass2 > 0) { if (gfc.highpass2 < 0.9 * (0.75 / 31.0)) { gfc.highpass1 = 0; gfc.highpass2 = 0; console.warn("Warning: highpass filter disabled. Frequency too small."); } }
        if (gfc.highpass2 > 0) { let maxband = -1; for (let band = 0; band <= 31; band++) { const freq = band / 31.0; if (freq <= gfc.highpass1) highpass_band = Math.max(highpass_band, band); if (gfc.highpass1 < freq && freq < gfc.highpass2) maxband = Math.max(maxband, band); } gfc.highpass1 = highpass_band / 31.0; gfc.highpass2 = (maxband === -1) ? (highpass_band + 0.75) / 31.0 : (maxband + 0.75) / 31.0; }
        for (let band = 0; band < 32; band++) { const freq = band / 31.0; const fc1 = (gfc.highpass2 > gfc.highpass1) ? this._filter_coef((gfc.highpass2 - freq) / (gfc.highpass2 - gfc.highpass1 + 1e-20)) : 1.0; const fc2 = (gfc.lowpass2 > gfc.lowpass1) ? this._filter_coef((freq - gfc.lowpass1) / (gfc.lowpass2 - gfc.lowpass1 + 1e-20)) : 1.0; gfc.amp_filter[band] = fc1 * fc2; }
     }
    /** @private */
    _lame_init_qval(gfp) { /* ... (implementation as before) ... */
        const gfc = gfp.internal_flags; const LAME_DEFAULT_QUALITY = 5;
        if (gfp.quality < 0) gfp.quality = LAME_DEFAULT_QUALITY;
        switch (gfp.quality) {
            default: case 9: gfc.psymodel = 0; gfc.noise_shaping = 0; gfc.noise_shaping_amp = 0; gfc.noise_shaping_stop = 0; gfc.use_best_huffman = 0; gfc.full_outer_loop = 0; break;
            case 8: gfp.quality = 7; // Fallthrough
            case 7: gfc.psymodel = 1; gfc.noise_shaping = 0; gfc.noise_shaping_amp = 0; gfc.noise_shaping_stop = 0; gfc.use_best_huffman = 0; gfc.full_outer_loop = 0; break;
            case 6: case 5: gfc.psymodel = 1; if (gfc.noise_shaping === 0) gfc.noise_shaping = 1; gfc.noise_shaping_amp = 0; gfc.noise_shaping_stop = 0; if (gfc.subblock_gain === -1) gfc.subblock_gain = 1; gfc.use_best_huffman = 0; gfc.full_outer_loop = 0; break;
            case 4: gfc.psymodel = 1; if (gfc.noise_shaping === 0) gfc.noise_shaping = 1; gfc.noise_shaping_amp = 0; gfc.noise_shaping_stop = 0; if (gfc.subblock_gain === -1) gfc.subblock_gain = 1; gfc.use_best_huffman = 1; gfc.full_outer_loop = 0; break;
            case 3: gfc.psymodel = 1; if (gfc.noise_shaping === 0) gfc.noise_shaping = 1; gfc.noise_shaping_amp = 1; gfc.noise_shaping_stop = 1; if (gfc.subblock_gain === -1) gfc.subblock_gain = 1; gfc.use_best_huffman = 1; gfc.full_outer_loop = 0; break;
            case 2: gfc.psymodel = 1; if (gfc.noise_shaping === 0) gfc.noise_shaping = 1; if (gfc.substep_shaping === 0) gfc.substep_shaping = 2; gfc.noise_shaping_amp = 1; gfc.noise_shaping_stop = 1; if (gfc.subblock_gain === -1) gfc.subblock_gain = 1; gfc.use_best_huffman = 1; gfc.full_outer_loop = 0; break;
            case 1: case 0: gfc.psymodel = 1; if (gfc.noise_shaping === 0) gfc.noise_shaping = 1; if (gfc.substep_shaping === 0) gfc.substep_shaping = 2; gfc.noise_shaping_amp = 2; gfc.noise_shaping_stop = 1; if (gfc.subblock_gain === -1) gfc.subblock_gain = 1; gfc.use_best_huffman = 1; gfc.full_outer_loop = 0; break;
        }
     }
    /** @private */
    _lame_init_bitstream(gfp) { /* ... (implementation as before) ... */
        const gfc = gfp.internal_flags; gfp.frameNum = 0;
        if (gfp.write_id3tag_automatic) { /* this.id3.id3tag_write_v2(gfp); */ console.warn("ID3v2 tag writing not implemented."); }
        gfc.bitrate_stereoMode_Hist = new_int_n([16, 5]); gfc.bitrate_blockType_Hist = new_int_n([16, 6]);
        gfc.PeakSample = 0.0;
        if (gfp.bWriteVbrTag) { /* this.vbr.InitVbrTag(gfp); */ console.warn("VBR tag writing not implemented."); }
     }
    /** @private */
    _calcNeeded(gfp) { /* ... (implementation as before) ... */
        let mf_needed = Encoder.BLKSIZE + gfp.framesize - Encoder.FFTOFFSET;
        mf_needed = Math.max(mf_needed, 512 + gfp.framesize - 32);
        assert(LameInternalFlags.MFSIZE >= mf_needed, `MFSIZE ${LameInternalFlags.MFSIZE} too small, need ${mf_needed}`);
        return mf_needed;
     }
    /** @private */
    _lame_encode_frame(gfp, inbuf_l, inbuf_r, mp3buf, mp3bufPos, mp3buf_size) { /* ... (implementation as before) ... */
        const ret = this.enc.lame_encode_mp3_frame(gfp, inbuf_l, inbuf_r, mp3buf, mp3bufPos, mp3buf_size);
        gfp.frameNum++; return ret;
     }
    /** @private */
    _fill_buffer_resample(gfp, outbuf, outbufPos, desired_len, inbuf, in_bufferPos, len, num_used, ch) {
        const gfc = gfp.internal_flags;
        let i, j = 0, k;
        const bpc = Math.min(LameInternalFlags.BPC, gfp.out_samplerate / this._gcd(gfp.out_samplerate, gfp.in_samplerate));
        const intratio = (Math.abs(gfc.resample_ratio - Math.round(gfc.resample_ratio)) < 0.0001) ? 1 : 0;
        let fcn = 1.0 / gfc.resample_ratio;
        if (fcn > 1.0) fcn = 1.0;
        let filter_l = 31;
        if (filter_l % 2 === 0) --filter_l;
        filter_l += intratio;

        if (gfc.fill_buffer_resample_init === 0) {
            gfc.inbuf_old[0] = new_float(LameInternalFlags.INBUF_SIZE);
            gfc.inbuf_old[1] = new_float(LameInternalFlags.INBUF_SIZE);
            for (i = 0; i <= 2 * bpc; ++i) {
                gfc.blackfilt[i] = new_float(Lame.BLACKSIZE);
            }
            gfc.itime[0] = 0;
            gfc.itime[1] = 0;

            for (j = 0; j <= 2 * bpc; j++) {
                let sum = 0.0;
                const offset = (j - bpc) / (2.0 * bpc);
                for (i = 0; i <= filter_l; i++) {
                    sum += gfc.blackfilt[j][i] = this._blackman(i - offset, fcn, filter_l);
                }
                if (sum !== 0) {
                    for (i = 0; i <= filter_l; i++) {
                        gfc.blackfilt[j][i] /= sum;
                    }
                }
            }
            gfc.fill_buffer_resample_init = 1;
        }

        const inbuf_old = gfc.inbuf_old[ch];
        for (k = 0; k < desired_len; k++) {
            const time0 = k * gfc.resample_ratio;
            j = Math.floor(time0 - gfc.itime[ch]);
            if ((filter_l + j - Math.floor(filter_l / 2)) >= len) break;

            const offset = (time0 - gfc.itime[ch] - (j + 0.5 * (filter_l % 2)));
            assert(Math.abs(offset) <= 0.501, `Resample offset error: ${offset}`);
            const joff = Math.floor((offset * 2.0 * bpc) + bpc + 0.5);
            assert(joff >= 0 && joff <= 2 * bpc, `joff out of bounds: ${joff}`);

            let xvalue = 0.0;
            for (i = 0; i <= filter_l; ++i) {
                const j2 = Math.floor(i + j - filter_l / 2);
                const y = (j2 < 0) ? inbuf_old[LameInternalFlags.INBUF_SIZE/2 + j2] : inbuf[in_bufferPos + j2];
                xvalue += y * gfc.blackfilt[joff][i];
            }
            outbuf[outbufPos + k] = xvalue;
        }

        num_used.num_used = Math.min(len, Math.floor(filter_l + j - filter_l / 2));
        if (num_used.num_used < 0) num_used.num_used = 0;

        gfc.itime[ch] += num_used.num_used - k * gfc.resample_ratio;

        if (num_used.num_used >= Lame.BLACKSIZE) {
            for (i = 0; i < Lame.BLACKSIZE; i++) {
                inbuf_old[i] = inbuf[in_bufferPos + num_used.num_used + i - Lame.BLACKSIZE];
            }
        } else {
            const n_shift = Lame.BLACKSIZE - num_used.num_used;
            for (i = 0; i < n_shift; ++i) {
                inbuf_old[i] = inbuf_old[i + num_used.num_used];
            }
            for (j = 0; i < Lame.BLACKSIZE; ++i, ++j) {
                inbuf_old[i] = inbuf[in_bufferPos + j];
            }
            assert(j === num_used.num_used, `Resample buffer fill mismatch: ${j} vs ${num_used.num_used}`);
        }
        return k;
    }
    /** @private */
    _fill_buffer(gfp, mfbuf, in_buffer, in_bufferPos, nsamples, io) {
        const gfc = gfp.internal_flags;
        if (Math.abs(gfc.resample_ratio - 1.0) > 1e-6) {
            // Resampling logic (assuming it handles buffer limits correctly internally)
                for (let ch = 0; ch < gfc.channels_out; ch++) {
                    let numUsed = new this._NumUsed();
                    // Ensure we don't ask resampler to create more samples than fit
                    const space_left_out = LameInternalFlags.MFSIZE - gfc.mf_size;
                    const desired_len = Math.min(gfp.framesize, space_left_out); // Max output needed/possible
                    if (desired_len <= 0) { // Should not happen if mf_size < mf_needed
                    io.n_out = 0;
                    io.n_in = 0;
                    continue; // Skip channel if no space
                    }
                    io.n_out = this._fill_buffer_resample(gfp, mfbuf[ch], gfc.mf_size, desired_len, in_buffer[ch], in_bufferPos, nsamples, numUsed, ch);
                    io.n_in = numUsed.num_used;
                    // If stereo, assume n_in is the same for both, check n_out?
                    // This simplified approach might need review for stereo resampling edge cases.
                }
        } else { // No resampling
            // Determine samples to copy based on input, frame size, AND available buffer space
            const space_left = LameInternalFlags.MFSIZE - gfc.mf_size;
            io.n_out = Math.min(gfp.framesize, nsamples, space_left);
            io.n_in = io.n_out; // Samples consumed = samples output

            // Check if calculation resulted in negative/zero copy, which shouldn't happen if mf_size < mf_needed
            if (io.n_out <= 0) {
                io.n_in = 0; // Consume nothing if outputting nothing
                return; // Nothing to copy
            }

            // Copy the calculated number of samples
            for (let i = 0; i < io.n_out; ++i) {
                    const dest_idx = gfc.mf_size + i;
                    // Double check bounds just in case
                    if (dest_idx >= LameInternalFlags.MFSIZE) {
                        console.error(`Buffer overflow detected in _fill_buffer: index ${dest_idx} >= ${LameInternalFlags.MFSIZE}`);
                        // Adjust n_out/n_in if an error occurs mid-loop? Unlikely with initial check.
                        io.n_out = i; // Record how many were actually copied
                        io.n_in = io.n_out;
                        break;
                    }
                mfbuf[0][dest_idx] = in_buffer[0][in_bufferPos + i];
                if (gfc.channels_out === 2) {
                    mfbuf[1][dest_idx] = in_buffer[1][in_bufferPos + i];
                }
            }
        }
    }
    /** @private */
    _lame_encode_buffer_sample(gfp, buffer_l, buffer_r, nsamples, mp3buf, mp3bufPos, mp3buf_size) { /* ... (implementation as before) ... */
        const gfc = gfp.internal_flags; let mp3size = 0; let ret; let i; let ch; let mp3out; const mfbuf = [gfc.mfbuf[0], gfc.mfbuf[1]]; const in_buffer = [buffer_l, buffer_r];
        if (gfc.Class_ID !== this._LAME_ID) return -3; if (nsamples === 0) return 0;
        mp3out = this.bs.copy_buffer(gfc, mp3buf, mp3bufPos, mp3buf_size, 0); if (mp3out < 0) return mp3out; mp3bufPos += mp3out; mp3size += mp3out;
        const mf_needed = this._calcNeeded(gfp); let in_bufferPos = 0;
        while (nsamples > 0) { const io = new this._InOut(); this._fill_buffer(gfp, mfbuf, in_buffer, in_bufferPos, nsamples, io); const n_in = io.n_in; const n_out = io.n_out;
            if (gfc.findReplayGain && !gfc.decode_on_the_fly) { if (this.ga.AnalyzeSamples(gfc.rgdata, mfbuf[0], gfc.mf_size, mfbuf[1], gfc.mf_size, n_out, gfc.channels_out) === GainAnalysis.GAIN_ANALYSIS_ERROR) return -6; }
            nsamples -= n_in; in_bufferPos += n_in; gfc.mf_size += n_out; assert(gfc.mf_size <= LameInternalFlags.MFSIZE, "mf_size overflow");
            if (gfc.mf_samples_to_encode < 1) gfc.mf_samples_to_encode = Encoder.ENCDELAY + Encoder.POSTDELAY; gfc.mf_samples_to_encode += n_out;
            if (gfc.mf_size >= mf_needed) { const buf_size = (mp3buf_size === 0) ? 0 : mp3buf_size - mp3size; ret = this._lame_encode_frame(gfp, mfbuf[0], mfbuf[1], mp3buf, mp3bufPos, buf_size); if (ret < 0) return ret; mp3bufPos += ret; mp3size += ret;
                gfc.mf_size -= gfp.framesize; gfc.mf_samples_to_encode -= gfp.framesize; for (ch = 0; ch < gfc.channels_out; ch++) {
                    // Use subarray and set for efficient shifting
                    mfbuf[ch].set(mfbuf[ch].subarray(gfp.framesize, gfc.mf_size + gfp.framesize));
                 }
            }
        } assert(nsamples === 0, "Not all input samples processed"); return mp3size;
     }


    // --- Public API Methods ---

    /**
     * Initializes the LAME encoder global flags structure (`LameGlobalFlags`)
     * with default values.
     *
     * @public
     * @returns {LameGlobalFlags | null} A new LameGlobalFlags object with defaults set, or null on failure.
     */
    lame_init() {
        const gfp = new LameGlobalFlags();
        const ret = this._lame_init_old(gfp);
        if (ret !== 0) return null;
        gfp.lame_allocated_gfp = 1;
        return gfp;
    }

    /**
     * Initializes the internal encoder parameters based on the settings in the
     * provided `LameGlobalFlags` structure. Must be called after `lame_init`
     * and setting desired parameters, but before encoding.
     *
     * @public
     * @param {LameGlobalFlags} gfp - The configured LameGlobalFlags structure.
     * @returns {number} 0 on success, negative error code on failure.
     */
    lame_init_params(gfp) {
        if (!gfp || !gfp.internal_flags) { console.error("lame_init_params: lame_init() not called?"); return -1; }
        const gfc = gfp.internal_flags;

        // --- Parameter Validation and Derivation ---
        gfc.Class_ID = this._LAME_ID;
        if (gfc.ATH == null) gfc.ATH = new ATH();
        if (gfc.PSY == null) gfc.PSY = new _PSY(); // Use internal class _PSY
        if (gfc.rgdata == null) gfc.rgdata = new ReplayGain();

        // ...(Rest of lame_init_params implementation remains largely the same)...
        gfc.channels_in = gfp.num_channels;
        if (gfc.channels_in === 1) gfp.mode = MPEGMode.MONO;
        gfc.channels_out = (gfp.mode === MPEGMode.MONO) ? 1 : 2;
        gfc.mode_ext = Encoder.MPG_MD_MS_LR;
        if (gfp.mode === MPEGMode.MONO) gfp.force_ms = false;
        if (gfp.VBR === VbrMode.vbr_off && gfp.VBR_mean_bitrate_kbps !== 128 && gfp.brate === 0) gfp.brate = gfp.VBR_mean_bitrate_kbps;
        if (gfp.VBR === VbrMode.vbr_off && gfp.brate === 0) { if (Math.abs(gfp.compression_ratio) < 1e-6) gfp.compression_ratio = 11.025; }
        if (!(gfp.VBR === VbrMode.vbr_off || gfp.VBR === VbrMode.vbr_mtrh || gfp.VBR === VbrMode.vbr_mt)) gfp.free_format = false;
        if (gfp.VBR === VbrMode.vbr_off && gfp.compression_ratio > 0) {
             if (gfp.out_samplerate === 0) { if (gfp.out_samplerate === 0) gfp.out_samplerate = gfp.in_samplerate; }
             gfp.brate = Math.floor(gfp.out_samplerate * 16 * gfc.channels_out / (1000.0 * gfp.compression_ratio));
             this._SmpFrqIndex(gfp.out_samplerate, gfp);
             if (!gfp.free_format) gfp.brate = this._FindNearestBitrate(gfp.brate, gfp.version, gfp.out_samplerate);
        }
        if (gfp.out_samplerate === 0) {
             if (gfp.lowpassfreq === 0) {
                 let lowpass = 16000.0;
                 if (gfp.VBR === VbrMode.vbr_off || gfp.VBR === VbrMode.vbr_abr) { const lh = new _LowPassHighPass(); this._optimum_bandwidth(lh, gfp.VBR === VbrMode.vbr_off ? gfp.brate : gfp.VBR_mean_bitrate_kbps); lowpass = lh.lowerlimit; }
                 else { const x = [19500, 19000, 18500, 18000, 17500, 16500, 15500, 14500, 12500, 9500, 3950]; if (0 <= gfp.VBR_q && gfp.VBR_q <= 9) { const a = x[gfp.VBR_q], b = x[gfp.VBR_q + 1], m = gfp.VBR_q_frac; lowpass = Util.linear_int(a, b, m); } else lowpass = 19500; }
                 if (gfp.mode === MPEGMode.MONO && (gfp.VBR === VbrMode.vbr_off || gfp.VBR === VbrMode.vbr_abr)) lowpass *= 1.5;
                 gfp.lowpassfreq = Math.floor(lowpass);
             }
             if (2 * gfp.lowpassfreq > gfp.in_samplerate) gfp.lowpassfreq = Math.floor(gfp.in_samplerate / 2);
             gfp.out_samplerate = this._optimum_samplefreq(gfp.lowpassfreq, gfp.in_samplerate);
        }
        gfp.lowpassfreq = Math.min(20500, gfp.lowpassfreq); gfp.lowpassfreq = Math.min(Math.floor(gfp.out_samplerate / 2), gfp.lowpassfreq);
        if (gfp.VBR === VbrMode.vbr_abr) gfp.compression_ratio = gfp.out_samplerate * 16 * gfc.channels_out / (1000.0 * gfp.VBR_mean_bitrate_kbps);
        else if (gfp.VBR !== VbrMode.vbr_off) { const cmp = [5.7, 6.5, 7.3, 8.2, 10.0, 11.9, 13.0, 14.0, 15.0, 16.5]; gfp.compression_ratio = cmp[gfp.VBR_q] || 8.8; }
        else { if (gfp.brate > 0) gfp.compression_ratio = gfp.out_samplerate * 16 * gfc.channels_out / (1000.0 * gfp.brate); }
        gfc.findReplayGain = gfp.findReplayGain && gfp.bWriteVbrTag; gfc.decode_on_the_fly = gfp.decode_on_the_fly && gfp.bWriteVbrTag; gfc.findPeakSample = gfc.decode_on_the_fly;
        if (gfc.findReplayGain) { if (this.ga.InitGainAnalysis(gfc.rgdata, gfp.out_samplerate) === GainAnalysis.GAIN_ANALYSIS_ERROR) return -6; }
        if (gfc.decode_on_the_fly && !gfp.decode_only) { /* Initialize mpglib decoder */ console.warn("Decode on the fly not implemented."); }
        gfc.mode_gr = (gfp.out_samplerate <= 24000 && gfp.version === 0) ? 1 : 2;
        gfp.framesize = (gfp.version === 1 ? 1152 : (576 * gfc.mode_gr));
        gfp.encoder_delay = Encoder.ENCDELAY;
        gfc.resample_ratio = gfp.in_samplerate / gfp.out_samplerate;
        if (gfp.mode === MPEGMode.NOT_SET) gfp.mode = MPEGMode.JOINT_STEREO;
        // this._lame_init_params_ppflt(gfp); // Polyphase filter init
        gfc.samplerate_index = this._SmpFrqIndex(gfp.out_samplerate, gfp); if (gfc.samplerate_index < 0) return -1;
        if (gfp.VBR === VbrMode.vbr_off) { if (gfp.free_format) gfc.bitrate_index = 0; else { gfp.brate = this._FindNearestBitrate(gfp.brate, gfp.version, gfp.out_samplerate); gfc.bitrate_index = this._BitrateIndex(gfp.brate, gfp.version, gfp.out_samplerate); if (gfc.bitrate_index <= 0) return -1; } }
        else { if(gfp.VBR === VbrMode.vbr_abr) { gfc.bitrate_index = this._BitrateIndex(this._FindNearestBitrate(gfp.VBR_mean_bitrate_kbps, gfp.version, gfp.out_samplerate), gfp.version, gfp.out_samplerate); if(gfc.bitrate_index <= 0) gfc.bitrate_index = 1; } else { gfc.bitrate_index = 1; } }
        if (gfp.VBR !== VbrMode.vbr_off) {
            gfc.VBR_min_bitrate = 1; gfc.VBR_max_bitrate = 14;
            if (gfp.out_samplerate < 16000) gfc.VBR_max_bitrate = 8; else if (gfp.out_samplerate < 32000 && gfp.version === 0) gfc.VBR_max_bitrate = 14;
            if (gfp.VBR_min_bitrate_kbps > 0) { gfp.VBR_min_bitrate_kbps = this._FindNearestBitrate(gfp.VBR_min_bitrate_kbps, gfp.version, gfp.out_samplerate); gfc.VBR_min_bitrate = this._BitrateIndex(gfp.VBR_min_bitrate_kbps, gfp.version, gfp.out_samplerate); if (gfc.VBR_min_bitrate < 0) return -1; }
            if (gfp.VBR_max_bitrate_kbps > 0) { gfp.VBR_max_bitrate_kbps = this._FindNearestBitrate(gfp.VBR_max_bitrate_kbps, gfp.version, gfp.out_samplerate); gfc.VBR_max_bitrate = this._BitrateIndex(gfp.VBR_max_bitrate_kbps, gfp.version, gfp.out_samplerate); if (gfc.VBR_max_bitrate < 0) return -1; }
            if (gfc.VBR_min_bitrate > gfc.VBR_max_bitrate) gfc.VBR_min_bitrate = gfc.VBR_max_bitrate;
            gfp.VBR_min_bitrate_kbps = bitrate_table[gfp.version][gfc.VBR_min_bitrate]; gfp.VBR_max_bitrate_kbps = bitrate_table[gfp.version][gfc.VBR_max_bitrate];
            if (gfp.VBR === VbrMode.vbr_abr) { gfp.VBR_mean_bitrate_kbps = Math.min(gfp.VBR_max_bitrate_kbps, gfp.VBR_mean_bitrate_kbps); gfp.VBR_mean_bitrate_kbps = Math.max(gfp.VBR_min_bitrate_kbps, gfp.VBR_mean_bitrate_kbps); }
        }
        let sfb_idx = gfc.samplerate_index + (3 * gfp.version); if (gfp.out_samplerate < 16000) sfb_idx += 6;
        if (sfb_idx >= 0 && sfb_idx < this.qupvt.sfBandIndex.length) {
            const sfb_info = this.qupvt.sfBandIndex[sfb_idx];
            // Use subarray().set() for copying TypedArrays
            gfc.scalefac_band.l.set(sfb_info.l.subarray(0, Encoder.SBMAX_l + 1));
            gfc.scalefac_band.s.set(sfb_info.s.subarray(0, Encoder.SBMAX_s + 1));
            let size_l = Math.floor((gfc.scalefac_band.l[22] - gfc.scalefac_band.l[21]) / Encoder.PSFB21); for (let i = 0; i < Encoder.PSFB21; i++) gfc.scalefac_band.psfb21[i] = gfc.scalefac_band.l[21] + i * size_l; gfc.scalefac_band.psfb21[Encoder.PSFB21] = gfc.scalefac_band.l[22];
            let size_s = Math.floor((gfc.scalefac_band.s[13] - gfc.scalefac_band.s[12]) / Encoder.PSFB12); for (let i = 0; i < Encoder.PSFB12; i++) gfc.scalefac_band.psfb12[i] = gfc.scalefac_band.s[12] + i * size_s; gfc.scalefac_band.psfb12[Encoder.PSFB12] = gfc.scalefac_band.s[13];
        } else { console.error("Invalid samplerate index for scalefactor bands."); return -1; }
        if (gfp.version === 1) gfc.sideinfo_len = (gfc.channels_out === 1) ? 17 + 4 : 32 + 4; else gfc.sideinfo_len = (gfc.channels_out === 1) ? 9 + 4 : 17 + 4; if (gfp.error_protection) gfc.sideinfo_len += 2;
        this.bs.init_bit_stream_w(gfc); this._lame_init_bitstream(gfp);
        // --- Select Iteration Loop ---
        switch (gfp.VBR) {
            case VbrMode.vbr_mt: case VbrMode.vbr_mtrh: console.warn(`Using placeholder CBR loop for VBR mode: ${gfp.VBR}`); gfc.iteration_loop = new CBRNewIterationLoop(this.qu); break;
            case VbrMode.vbr_rh: console.error("VBROldIterationLoop not implemented for VBR mode: vbr_rh"); return -1;
            case VbrMode.vbr_abr: console.error("ABRIterationLoop not implemented for VBR mode: vbr_abr"); return -1;
            case VbrMode.vbr_off: default: gfc.iteration_loop = new CBRNewIterationLoop(this.qu); break;
        }
        this.qupvt.iteration_init(gfp); this.psy.psymodel_init(gfp); this._lame_init_qval(gfp);
        // --- Final Param Checks ---
        if (gfp.scale < 0) gfp.scale = 1.0; if (gfp.scale_left < 0) gfp.scale_left = gfp.scale; if (gfp.scale_right < 0) gfp.scale_right = gfp.scale;
        if (gfp.athaa_type < 0) gfc.ATH.useAdjust = 3; else gfc.ATH.useAdjust = gfp.athaa_type;
        gfc.ATH.aaSensitivityP = Math.pow(10.0, gfp.athaa_sensitivity / -10.0);
        if (gfp.short_blocks == null) gfp.short_blocks = ShortBlock.short_block_allowed;
        if (gfp.short_blocks === ShortBlock.short_block_allowed && (gfp.mode === MPEGMode.JOINT_STEREO || gfp.mode === MPEGMode.STEREO)) gfp.short_blocks = ShortBlock.short_block_coupled;
        if (gfp.quant_comp < 0) gfp.quant_comp = 1; if (gfp.quant_comp_short < 0) gfp.quant_comp_short = 0;
        if (gfp.msfix < -10) gfp.msfix = 0.0;
        gfp.exp_nspsytune |= 1; if (gfc.nsPsy.attackthre < 0) gfc.nsPsy.attackthre = 4.4; if (gfc.nsPsy.attackthre_s < 0) gfc.nsPsy.attackthre_s = 25.0;
        if (gfp.ATHtype < 0) gfp.ATHtype = 4; if (gfp.ATHcurve < -99) gfp.ATHcurve = 0.0;
        if (gfp.athaa_loudapprox < 0) gfp.athaa_loudapprox = 2; if (gfp.interChRatio < 0) gfp.interChRatio = 0.0;
        if (gfp.useTemporal == null) gfp.useTemporal = true;
        gfc.slot_lag = 0; gfc.frac_SpF = 0;
        assert(gfp.scale >= 0, "Scale not initialized");

        // Set sfb21_extra based on VBR mode and out_samplerate
        if (gfp.VBR === VbrMode.vbr_rh || gfp.VBR === VbrMode.vbr_mtrh) {
            if (gfp.experimentalY) {
                gfc.sfb21_extra = false;
            } else {
                gfc.sfb21_extra = (gfp.out_samplerate > 44000);
            }
        } else {
            // CBR/ABR mode
            gfc.sfb21_extra = false;
        }

        return 0; // Success
    }


    /**
     * Encodes the final buffered samples and writes remaining MP3 data,
     * including VBR tag and ID3 tags if configured. Call this after processing
     * all input audio samples.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {Uint8Array} mp3buffer - Output buffer to receive the encoded MP3 data.
     * @param {number} mp3bufPos - Starting position offset within `mp3buffer`.
     * @param {number} mp3buffer_size - Maximum number of bytes available in `mp3buffer` from `mp3bufferPos`.
     * @returns {number} The number of bytes written to `mp3buffer`, or a negative error code.
     */
    lame_encode_flush(gfp, mp3buffer, mp3bufferPos, mp3buffer_size) {
        const gfc = gfp.internal_flags;
        const buffer = [new_float(1152), new_float(1152)];
        let imp3 = 0, mp3count = 0, mp3buffer_size_remaining;

        if (gfc.mf_samples_to_encode < 1) return 0; // Already flushed

        let samples_to_encode = gfc.mf_samples_to_encode - Encoder.POSTDELAY;
        const mf_needed = this._calcNeeded(gfp);
        if (gfp.in_samplerate !== gfp.out_samplerate) samples_to_encode += Math.floor(16.0 * gfp.out_samplerate / gfp.in_samplerate);
        let end_padding = gfp.framesize - (samples_to_encode % gfp.framesize);
        if (end_padding < 576 && gfp.framesize >= 576) end_padding += gfp.framesize;
        else if (end_padding === gfp.framesize) end_padding = 0;
        gfp.encoder_padding = end_padding;
        let frames_left = Math.ceil((samples_to_encode + end_padding) / gfp.framesize);

        while (frames_left > 0 && imp3 >= 0) {
            let bunch = mf_needed - gfc.mf_size;
            if(gfc.resample_ratio > 1e-6) bunch = Math.ceil(bunch * gfc.resample_ratio); else bunch = Math.ceil(bunch);
            bunch = Math.min(bunch, 1152); if (bunch < 1) bunch = 1;
            Arrays.fill(buffer[0], 0, bunch, 0.0); if (gfc.channels_out === 2) Arrays.fill(buffer[1], 0, bunch, 0.0);
            mp3buffer_size_remaining = (mp3buffer_size === 0) ? 0 : mp3buffer_size - mp3count;
            const frame_num_before = gfp.frameNum;
            // Use the public float version here, assumes padding buffer is float
            imp3 = this.lame_encode_buffer_ieee_float(gfp, buffer[0], buffer[1], bunch, mp3buffer, mp3bufferPos, mp3buffer_size_remaining);
            if (imp3 < 0) return imp3;
            mp3bufPos += imp3; mp3count += imp3;
            if (gfp.frameNum > frame_num_before) frames_left--;
        }
        gfc.mf_samples_to_encode = 0;
        if (imp3 < 0) return imp3;

        mp3buffer_size_remaining = (mp3buf_size === 0) ? 0 : mp3buf_size - mp3count;
        this.bs.flush_bitstream(gfp);
        imp3 = this.bs.copy_buffer(gfc, mp3buffer, mp3bufferPos, mp3buffer_size_remaining, 1);
        if (imp3 < 0) return imp3;
        mp3bufPos += imp3; mp3count += imp3;

        if (gfp.write_id3tag_automatic) {
            /* this.id3.id3tag_write_v1(gfp); */ console.warn("ID3v1 tag writing not implemented.");
            mp3buffer_size_remaining = (mp3buf_size === 0) ? 0 : mp3buf_size - mp3count;
            imp3 = this.bs.copy_buffer(gfc, mp3buffer, mp3bufferPos, mp3buffer_size_remaining, 0);
            if (imp3 < 0) return imp3;
            mp3count += imp3;
        }
        return mp3count;
    }

    /**
     * Encodes a buffer of PCM audio samples provided as Float32 arrays.
     * Input samples should be in the range [-1.0, 1.0].
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {Float32Array} buffer_l - Buffer for the left channel (or mono).
     * @param {Float32Array|null} buffer_r - Buffer for the right channel (provide null or same as left if mono).
     * @param {number} nsamples - Number of samples per channel in the input buffers.
     * @param {Uint8Array} mp3buf - Output buffer to receive the encoded MP3 data.
     * @param {number} mp3bufPos - Starting position offset within `mp3buf`.
     * @param {number} mp3buf_size - Maximum number of bytes available in `mp3buf` from `mp3bufPos`.
     * @returns {number} The number of bytes written to `mp3buf`, or a negative error code.
     */
    lame_encode_buffer_ieee_float(gfp, buffer_l, buffer_r, nsamples, mp3buf, mp3bufPos, mp3buf_size) {
        const gfc = gfp.internal_flags;
        if (!gfc || gfc.Class_ID !== this._LAME_ID) return -3;

        let input_r = buffer_r;
        if (gfc.channels_in > 1 && !input_r) { console.error("Right channel buffer required for stereo input."); return -1; }
        if (gfc.channels_in === 1) input_r = buffer_l;

        let working_l = buffer_l; let working_r = input_r; let needs_copy = false;
        if ((Math.abs(gfp.scale - 1.0) > 1e-6 && gfp.scale > 0) || (Math.abs(gfp.scale_left - 1.0) > 1e-6 && gfp.scale_left > 0) || (gfc.channels_out === 2 && Math.abs(gfp.scale_right - 1.0) > 1e-6 && gfp.scale_right > 0) || (gfp.num_channels === 2 && gfc.channels_out === 1)) {
            needs_copy = true; working_l = new_float(nsamples); if (gfc.channels_in > 1) working_r = new_float(nsamples);
        }
        const eff_scale_l = (Math.abs(gfp.scale_left - 1.0) > 1e-6 && gfp.scale_left > 0) ? gfp.scale_left : gfp.scale;
        const eff_scale_r = (gfc.channels_in > 1 && Math.abs(gfp.scale_right - 1.0) > 1e-6 && gfp.scale_right > 0) ? gfp.scale_right : gfp.scale;
        if (gfp.num_channels === 2 && gfc.channels_out === 1) { for (let i = 0; i < nsamples; i++) { working_l[i] = 0.5 * (buffer_l[i] * eff_scale_l + buffer_r[i] * eff_scale_r); } working_r = null; }
        else { const scale_l = (Math.abs(eff_scale_l - 1.0) > 1e-6 && eff_scale_l > 0); const scale_r = (gfc.channels_in > 1 && Math.abs(eff_scale_r - 1.0) > 1e-6 && eff_scale_r > 0); if (needs_copy) { for (let i = 0; i < nsamples; i++) { working_l[i] = scale_l ? buffer_l[i] * eff_scale_l : buffer_l[i]; if (working_r) working_r[i] = scale_r ? buffer_r[i] * eff_scale_r : buffer_r[i]; } } else { working_l = buffer_l; working_r = buffer_r; } }
        if (gfc.channels_out === 1) working_r = working_l;

        return this._lame_encode_buffer_sample(gfp, working_l, working_r, nsamples, mp3buf, mp3bufPos, mp3buf_size);
    }

    /**
     * Encodes a buffer of PCM audio samples provided as Int16 arrays.
     * Input samples should be in the range [-32768, 32767].
     * Converts the input to floating point before encoding.
     *
     * @public
     * @param {LameGlobalFlags} gfp - LAME global flags.
     * @param {Int16Array} buffer_l - Buffer for the left channel (or mono).
     * @param {Int16Array|null} buffer_r - Buffer for the right channel (provide null or same as left if mono).
     * @param {number} nsamples - Number of samples per channel in the input buffers.
     * @param {Uint8Array} mp3buf - Output buffer to receive the encoded MP3 data.
     * @param {number} mp3bufPos - Starting position offset within `mp3buf`.
     * @param {number} mp3buf_size - Maximum number of bytes available in `mp3buf` from `mp3bufPos`.
     * @returns {number} The number of bytes written to `mp3buf`, or a negative error code.
     */
    lame_encode_buffer(gfp, buffer_l, buffer_r, nsamples, mp3buf, mp3bufPos, mp3buf_size) {
        const gfc = gfp.internal_flags;
        if (!gfc || gfc.Class_ID !== this._LAME_ID) return -3;
        if (nsamples === 0) return 0;
        const float_buf_l = new_float(nsamples); let float_buf_r = null;
        if (gfc.channels_in > 1) { if (!buffer_r) { console.error("Right channel Int16 buffer required for stereo input."); return -1; } float_buf_r = new_float(nsamples); }
        const scale_factor = 1.0 / 32768.0;
        for (let i = 0; i < nsamples; i++) { float_buf_l[i] = buffer_l[i] * scale_factor; if (float_buf_r && buffer_r) float_buf_r[i] = buffer_r[i] * scale_factor; }
        return this.lame_encode_buffer_ieee_float(gfp, float_buf_l, float_buf_r, nsamples, mp3buf, mp3bufPos, mp3buf_size);
    }


} // End class Lame

const LAME_MAXMP3BUFFER = Lame.LAME_MAXMP3BUFFER; // Assign static to const
export { Lame, LAME_MAXMP3BUFFER }; // Export both
export default Lame; // Also provide default export if needed
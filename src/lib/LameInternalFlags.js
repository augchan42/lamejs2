/**
 * @fileoverview Internal flags and state variables for the LAME encoder.
 * Ported from internal_flags.h. Contains detailed configuration derived
 * from global flags, psychoacoustic model state, bitstream state, etc.
 * Uses ES Module syntax.
 *
 * @module LameInternalFlags
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import { IIISideInfo } from './IIISideInfo.js';
import { ScaleFac } from './ScaleFac.js';
import { NsPsy } from './NsPsy.js';
import { VBRSeekInfo } from './VBRSeekInfo.js';
import { III_psy_xmin } from './III_psy_xmin.js';
import { Encoder, SBMAX_l, SBMAX_s } from './Encoder.js';
import { L3Side } from './L3Side.js';
// Assuming these types are defined elsewhere and imported if needed for full type safety
/** @typedef {import('./ATH.js').ATH} ATH */
/** @typedef {import('./ReplayGain.js').ReplayGain} ReplayGain */
/** @typedef {import('./IterationLoop.js').IterationLoop} IterationLoop */
/** @typedef {import('./ID3TagSpec.js').ID3TagSpec} ID3TagSpec */
/** @typedef {import('./MPGLib.js').MPGLib} MPGLib */ // Placeholder type
/** @typedef {import('./PlottingData.js').PlottingData} PlottingData */ // Placeholder type

// Destructure common utilities for easier access
const {
    // System, // Not used
    // VbrMode, // Not used directly
    // Float, // Not used
    // ShortBlock, // Not used
    // Util, // Not used
    Arrays, // Keep for potential fill/sort usage
    new_array_n,
    new_byte,
    // new_double, // Use Float64Array directly
    new_float,
    new_float_n,
    new_int,
    new_int_n,
    assert
} = common;

// Helper to create the nested structure for en/thm used in the constructor
/** @private */
function createMaskingInfoInternal() {
    // Helper for internal flags initialization
    return {
        l: new_float(SBMAX_l),
        // Array of [SBMAX_s] Float32Arrays[3]
        s: new_float_n([SBMAX_s, 3])
    };
}


/**
 * @classdesc Holds internal state variables and configuration flags used during encoding.
 * This structure is not meant to be directly manipulated by the application; it's
 * populated and used internally by LAME based on the `LameGlobalFlags`.
 * @constructs LameInternalFlags
 */
class LameInternalFlags {
    // --- Static Constants ---
    /** Size of main filter buffer mfbuf */
    static MFSIZE = (3 * 1152 + Encoder.ENCDELAY - Encoder.MDCTDELAY);
    /** Max size for header buffer array */
    static MAX_HEADER_BUF = 256;
    /** Max bits for a channel in a granule */
    static MAX_BITS_PER_CHANNEL = 4095;
    /** Max bits for a granule (sum of channels) */
    static MAX_BITS_PER_GRANULE = 7680;
    /** Max Coefficients per block for resampling filter calc? */
    static BPC = 320;
    /** Size of one channel's resample history buffer */
    static INBUF_SIZE = 4096; // Size used in original C internal_flags.h for inbuf_old

    // --- Properties ---

    /** Class Identifier (for type checking) @public @type {number} */
    Class_ID = 0;

    /** Initialization flag for lame_encode_frame @public @type {number} */
    lame_encode_frame_init = 0;
    /** Initialization flag for iteration_init @public @type {number} */
    iteration_init_init = 0;
    /** Initialization flag for resampler @public @type {number} */
    fill_buffer_resample_init = 0;

    /** Input buffer for MDCT [2][MFSIZE] @public @type {Array<Float32Array>} */
    mfbuf;

    /** Input buffer history for resampling [2][BLACKSIZE] @public @type {Array<Float32Array>} */
    inbuf_old;

    /** Noise history buffer 1 [4][CBANDS] @public @type {Array<Float32Array>} */
    nb_1;
    /** Noise history buffer 2 [4][CBANDS] @public @type {Array<Float32Array>} */
    nb_2;
    /** Short block noise history buffer 1 [4][CBANDS] @public @type {Array<Float32Array>} */
    nb_s1;
    /** Short block noise history buffer 2 [4][CBANDS] @public @type {Array<Float32Array>} */
    nb_s2;

    /** Loudness calculation state [ch] @public @type {Float32Array} */
    loudness_sq_save;
    /** Loudness squared per granule/channel [gr][ch] @public @type {Array<Float32Array>} */
    loudness_sq;


    // Psychoacoustic model internal arrays (sizes set, filled later)
    /** Number of lines in long partition bands [CBANDS] @public @type {Int32Array} */
    numlines_l;
    /** Number of lines in short partition bands [CBANDS] @public @type {Int32Array} */
    numlines_s;
    /** Boundary partition index for long sfbs [SBMAX_l] @public @type {Int32Array} */
    bo_l;
    /** Boundary partition index for short sfbs [SBMAX_s] @public @type {Int32Array} */
    bo_s;
    /** Mid partition index for long sfbs [SBMAX_l] @public @type {Int32Array} */
    bm_l;
    /** Mid partition index for short sfbs [SBMAX_s] @public @type {Int32Array} */
    bm_s;
    /** MLD factor per long sfb [SBMAX_l] @public @type {Float32Array} */
    mld_l;
    /** MLD factor per short sfb [SBMAX_s] @public @type {Float32Array} */
    mld_s;
    /** MLD factor per long partition band [CBANDS] @public @type {Float32Array} */
    mld_cb_l;
    /** MLD factor per short partition band [CBANDS] @public @type {Float32Array} */
    mld_cb_s;
    /** Reciprocal number of lines in long partition bands [CBANDS] @public @type {Float32Array} */
    rnumlines_l;
    /** Spreading function indices [CBANDS][2] @public @type {Array<Int32Array>} */
    s3ind;
    /** Short block spreading function indices [CBANDS][2] @public @type {Array<Int32Array>} */
    s3ind_s;
    /** Flattened long block spreading function values @public @type {Float32Array | null} */
    s3_ll = null;
    /** Flattened short block spreading function values @public @type {Float32Array | null} */
    s3_ss = null;
    /** Min value table for long blocks [CBANDS] @public @type {Float32Array} */
    minval_l;
    /** Min value table for short blocks [CBANDS] @public @type {Float32Array} */
    minval_s;
    /** State for substep shaping [SBMAX_l]? @public @type {Int32Array} */
    pseudohalf;
    /** Total energy per channel [4] @public @type {Float32Array} */
    tot_ener;
    /** Temporal masking decay factor @public @type {number} */
    decay = 0.0;


    /** Granules per frame (1 or 2) @public @type {number} */
    mode_gr = 0;
    /** Number of input channels @public @type {number} */
    channels_in = 0;
    /** Number of output channels @public @type {number} */
    channels_out = 0;
    /** Resampling ratio (in/out) @public @type {number} */
    resample_ratio = 0.0;

    /** Samples remaining in internal buffers to be encoded @public @type {number} */
    mf_samples_to_encode = 0;
    /** Current number of valid samples in mfbuf @public @type {number} */
    mf_size = 0;
    /** Min VBR bitrate index @public @type {number} */
    VBR_min_bitrate = 0;
    /** Max VBR bitrate index @public @type {number} */
    VBR_max_bitrate = 0;
    /** Current frame bitrate index @public @type {number} */
    bitrate_index = 0;
    /** Output samplerate index @public @type {number} */
    samplerate_index = 0;
    /** Stereo mode extension (e.g., MS/LR) @public @type {number} */
    mode_ext = 0;

    /* lowpass and highpass filter control */
    /** Normalized lower freq bound of lowpass @public @type {number} */
    lowpass1 = 0.0;
    /** Normalized upper freq bound of lowpass @public @type {number} */
    lowpass2 = 0.0;
    /** Normalized lower freq bound of highpass @public @type {number} */
    highpass1 = 0.0;
    /** Normalized upper freq bound of highpass @public @type {number} */
    highpass2 = 0.0;

    /* Noise shaping controls */
    /** @public @type {number} */ noise_shaping = 0;
    /** @public @type {number} */ noise_shaping_amp = 0;
    /** @public @type {number} */ substep_shaping = 0;
    /** @public @type {number} */ noise_shaping_stop = 0;
    /** @public @type {number} */ subblock_gain = 0;

    /** Psychoacoustic model active flag (0=off, 1=on) @public @type {number} */
    psymodel = 0;
    /** Use best Huffman table division @public @type {number} */
    use_best_huffman = 0;
    /** Force full outer loop search @public @type {number} */
    full_outer_loop = 0;

    /** Side info structure for the current frame @public @type {IIISideInfo} */
    l3_side;

    /* Padding state */
    /** Padding flag for current frame @public @type {number} */
    padding = 0;
    /** Fractional samples per frame (for padding calc) @public @type {number} */
    frac_SpF = 0.0;
    /** Slot lag for CBR padding @public @type {number} */
    slot_lag = 0;

    /** ID3 tag specification (details for tag writing) @public @type {ID3TagSpec | null} */
    tag_spec = null;
    /** Music CRC value @public @type {number} */
    nMusicCRC = 0;

    /* Quantization loop state */
    /** Previous global gain values [ch] @public @type {Int32Array} */
    OldValue;
    /** Current gain step size [ch] @public @type {Int32Array} */
    CurrentStep;
    /** Global masking adjustment factor @public @type {number} */
    masking_lower = 0.0;
    /** Scalefactor per coefficient (temporary?) [576] @public @type {Int32Array} */
    bv_scf;

    /** Use sfb21/sfb12 bands beyond standard limits. Will be properly set in lame_init_params based on VBR mode and output sample rate. @public @type {boolean} */
    sfb21_extra;

    /* Resampling state */
    /** Precomputed filter coefficients [joff][filter_l+1] @public @type {Array<Float32Array | null>} */
    blackfilt;
    /** Time offset for resampling [ch] @public @type {Float64Array} */
    itime;

    /** Length of side info in bytes @public @type {number} */
    sideinfo_len = 0;

    /* mdct state */
    /** Polyphase filterbank outputs [ch][gr][sfb][coeff] @public @type {Array<Array<Array<Float32Array>>>} */
    sb_sample;
    /** Polyphase filter amplitude adjustment [32] @public @type {Float32Array} */
    amp_filter;

    /* Header/Ancillary data state (if writing custom headers) */
    /** Header buffer objects @public @type {Array<object>} */
    header;
    /** Current header read pointer @public @type {number} */
    h_ptr = 0;
    /** Current header write pointer @public @type {number} */
    w_ptr = 0;
    /** Flag for ancillary data @public @type {number} */
    ancillary_flag = 0;

    /* Reservoir state */
    /** Current reservoir size (bits) @public @type {number} */
    ResvSize = 0;
    /** Max reservoir size (bits) @public @type {number} */
    ResvMax = 0;

    /** Scalefactor band boundary information @public @type {ScaleFac} */
    scalefac_band;

    /* Masking thresholds and energies (potentially M/S) [4] */
    /** Threshold per sfb @public @type {Array<{l: Float32Array, s: Array<Float32Array>}>} */
    thm;
    /** Energy per sfb @public @type {Array<{l: Float32Array, s: Array<Float32Array>}>} */
    en;

    /* M/S ratio history */
    /** @public @type {number} */ ms_ratio_s_old = 0.0;
    /** @public @type {number} */ ms_ratio_l_old = 0.0;
    /** @public @type {number} */ ms_ener_ratio_old = 0.0;

    /** Previous granule block type [ch] @public @type {Int32Array} */
    blocktype_old;

    /** NS PsyTune specific state @public @type {NsPsy} */
    nsPsy;

    /** VBR tag seek table information @public @type {VBRSeekInfo} */
    VBR_seek_table;

    /** ATH calculation results @public @type {ATH | null} */
    ATH = null;
    /** Psychoacoustic model parameters (internal class) @public @type {_PSY | null} */
    PSY = null; // Refers to internal _PSY class

    /* Gapless encoding state */
    /** Total samples for gapless info @public @type {number} */
    nogap_total = 0;
    /** Current sample count for gapless info @public @type {number} */
    nogap_current = 0;

    /* ReplayGain state */
    /** @public @type {boolean} */ decode_on_the_fly = false; // Default to false
    /** @public @type {boolean} */ findReplayGain = false;
    /** @public @type {boolean} */ findPeakSample = false;
    /** @public @type {number} */ PeakSample = 0.0;
    /** @public @type {number} */ RadioGain = 0;
    /** @public @type {number} */ AudiophileGain = 0;
    /** ReplayGain analysis data @public @type {ReplayGain | null} */
    rgdata = null;
    /** @public @type {number} */ noclipGainChange = 0;
    /** @public @type {number} */ noclipScale = 0.0;

    /* Simple statistics */
    /** Histogram [bitrate_idx][mode+1] @public @type {Array<Int32Array> | null} */
    bitrate_stereoMode_Hist = null;
    /** Histogram [bitrate_idx][blockType+1] @public @type {Array<Int32Array> | null} */
    bitrate_blockType_Hist = null;

    /** Plotting/Analysis info structure @public @type {PlottingData | null} */
    pinfo = null;
    /** mpglib decoder instance (if decode_on_the_fly) @public @type {MPGLib | null} */
    hip = null;

    /* Input buffer state (used by lame_encode_buffer) */
    /** @public @type {number} */ in_buffer_nsamples = 0;
    /** @public @type {Float32Array | null} */ in_buffer_0 = null;
    /** @public @type {Float32Array | null} */ in_buffer_1 = null;

    /** Function pointer to the selected iteration loop @public @type {IterationLoop | null} */
    iteration_loop = null;

    // --- Constructor ---
    constructor() {
        // Initialize arrays
        this.mfbuf = [new_float(LameInternalFlags.MFSIZE), new_float(LameInternalFlags.MFSIZE)];
        // Initialize input buffer history for resampling with INBUF_SIZE to handle negative indices
        this.inbuf_old = [new_float(LameInternalFlags.INBUF_SIZE), new_float(LameInternalFlags.INBUF_SIZE)];
        // Initialize sfb21_extra to false by default - will be updated in lame_init_params based on VBR mode and output sample rate
        this.sfb21_extra = false;

        // Initialize noise history buffers
        this.nb_1 = new_float_n([4, Encoder.CBANDS]);
        this.nb_2 = new_float_n([4, Encoder.CBANDS]);
        this.nb_s1 = new_float_n([4, Encoder.CBANDS]);
        this.nb_s2 = new_float_n([4, Encoder.CBANDS]);
        this.loudness_sq_save = new_float(2);
        this.loudness_sq = new_float_n([2, 2]); // Initialize loudness_sq array
        this.numlines_l = new_int(Encoder.CBANDS);
        this.numlines_s = new_int(Encoder.CBANDS);
        this.bo_l = new_int(Encoder.SBMAX_l);
        this.bo_s = new_int(Encoder.SBMAX_s);
        this.bm_l = new_int(Encoder.SBMAX_l);
        this.bm_s = new_int(Encoder.SBMAX_s);
        this.mld_l = new_float(Encoder.SBMAX_l);
        this.mld_s = new_float(Encoder.SBMAX_s);
        this.mld_cb_l = new_float(Encoder.CBANDS);
        this.mld_cb_s = new_float(Encoder.CBANDS);
        this.rnumlines_l = new_float(Encoder.CBANDS);
        this.s3ind = new_int_n([Encoder.CBANDS, 2]);
        this.s3ind_s = new_int_n([Encoder.CBANDS, 2]);
        this.minval_l = new_float(Encoder.CBANDS);
        this.minval_s = new_float(Encoder.CBANDS);
        this.pseudohalf = new_int(Encoder.SBMAX_l);
        this.tot_ener = new_float(4);
        this.OldValue = new_int(2);
        this.CurrentStep = new_int(2);
        this.bv_scf = new_int(576);
        this.blackfilt = new Array(2 * LameInternalFlags.BPC + 1).fill(null);
        this.itime = new Float64Array(2); // Use Float64Array directly
        this.sb_sample = new_array_n([2, 2, 18, Encoder.SBLIMIT], () => 0.0);
        this.amp_filter = new_float(32);
        this.header = new Array(LameInternalFlags.MAX_HEADER_BUF);
        this.scalefac_band = new ScaleFac();
        this.thm = new Array(4);
        this.en = new Array(4);
        for (let i = 0; i < 4; i++) {
            this.thm[i] = new III_psy_xmin();
            this.en[i] = new III_psy_xmin();
        }
        this.blocktype_old = new_int(2);
        this.nsPsy = new NsPsy();
        this.VBR_seek_table = new VBRSeekInfo();
        this.l3_side = new IIISideInfo();

        // Initialize Header array elements
        for (let i = 0; i < LameInternalFlags.MAX_HEADER_BUF; i++) {
             class Header {
                 constructor() { this.write_timing = 0; this.ptr = 0; this.buf = new_byte(40); }
             }
            this.header[i] = new Header();
        }

        // Ensure histogram arrays are initialized
        if (!this.bitrate_stereoMode_Hist) this.bitrate_stereoMode_Hist = new_int_n([16, 5]);
        if (!this.bitrate_blockType_Hist) this.bitrate_blockType_Hist = new_int_n([16, 6]);

        // Initialize other nullable properties to null explicitly
        this.ATH = null;
        this.PSY = null;
        this.rgdata = null;
        this.iteration_loop = null;
        this.hip = null;
        this.pinfo = null;
        this.tag_spec = null;
        this.in_buffer_0 = null;
        this.in_buffer_1 = null;
        this.s3_ll = null;
        this.s3_ss = null;
    }
}

// Define static constants after class definition if needed
// (Already defined within the class using static keyword)

// At the bottom of LameInternalFlags.js
const MAX_HEADER_BUF = LameInternalFlags.MAX_HEADER_BUF; // Assign static to const
export { LameInternalFlags, MAX_HEADER_BUF }; // Export both
export default LameInternalFlags;
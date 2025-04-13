/**
 * @fileoverview Data structure for ReplayGain analysis state.
 * Ported from gain_analysis.h/c. Stores buffers and running sums for
 * calculating ReplayGain values.
 * Uses ES Module syntax.
 *
 * @module ReplayGain
 */

// Import necessary modules and utilities using ES Module syntax
import * as common from './common.js';
import { GainAnalysis } from './GainAnalysis.js'; // Assuming GainAnalysis exports constants

// Destructure common utilities for easier access
const {
    // System, // Not used
    // VbrMode, // Not used
    // Float, // Not used
    // ShortBlock, // Not used
    // Util, // Not used
    // Arrays, // Not used
    // new_array_n, // Used indirectly
    // new_byte, // Not used
    // new_double, // Not used
    new_float,
    // new_float_n, // Not used
    new_int,
    // new_int_n, // Not used directly
    // assert // Not used
} = common;


/**
 * @classdesc Holds the state variables and buffers required for performing
 * ReplayGain analysis during MP3 encoding or analysis. This includes
 * input buffers with pre-padding for filters, intermediate filter results,
 * windowing information, and running sums for RMS calculations.
 * @constructs ReplayGain
 */
class ReplayGain {
    /**
     * Pre-buffer for left channel input samples, used by filters.
     * Size: `GainAnalysis.MAX_ORDER * 2`
     * @public
     * @type {Float32Array}
     */
    linprebuf;
    /**
     * Current write index within `linprebuf`.
     * @public
     * @type {number}
     */
    linpre = 0;

    /**
     * Buffer for left channel samples after the first filter stage.
     * Size: `GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER`
     * @public
     * @type {Float32Array}
     */
    lstepbuf;
    /**
     * Current write index within `lstepbuf`.
     * @public
     * @type {number}
     */
    lstep = 0;

    /**
     * Buffer for left channel samples after the second filter stage ("out").
     * Size: `GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER`
     * @public
     * @type {Float32Array}
     */
    loutbuf;
    /**
     * Current write index within `loutbuf`.
     * @public
     * @type {number}
     */
    lout = 0;

    /**
     * Pre-buffer for right channel input samples, used by filters.
     * Size: `GainAnalysis.MAX_ORDER * 2`
     * @public
     * @type {Float32Array}
     */
    rinprebuf;
    /**
     * Current write index within `rinprebuf`.
     * @public
     * @type {number}
     */
    rinpre = 0;

    /**
     * Buffer for right channel samples after the first filter stage.
     * Size: `GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER`
     * @public
     * @type {Float32Array}
     */
    rstepbuf;
    /**
     * Current write index within `rstepbuf`.
     * @public
     * @type {number}
     */
    rstep = 0;

    /**
     * Buffer for right channel samples after the second filter stage ("out").
     * Size: `GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER`
     * @public
     * @type {Float32Array}
     */
    routbuf;
    /**
     * Current write index within `routbuf`.
     * @public
     * @type {number}
     */
    rout = 0;

    /**
     * Number of samples required for the RMS window duration (e.g., 50ms).
     * Calculated based on sample rate.
     * @public
     * @type {number}
     */
    sampleWindow = 0;

    /**
     * Total number of samples processed so far.
     * @public
     * @type {number}
     */
    totsamp = 0;

    /**
     * Running sum of squared samples for the left channel within the current window.
     * @public
     * @type {number}
     */
    lsum = 0.0;

    /**
     * Running sum of squared samples for the right channel within the current window.
     * @public
     * @type {number}
     */
    rsum = 0.0;

    /**
     * Index indicating the current frequency weighting step (related to K-weighting).
     * @public
     * @type {number}
     */
    freqindex = 0;

    /**
     * Flag indicating if this is the first block being processed (used for filter initialization).
     * @public
     * @type {number} // Typically 0 or 1
     */
    first = 0;

    /**
     * Histogram A for loudness distribution analysis.
     * Size determined by `GainAnalysis.STEPS_per_dB * GainAnalysis.MAX_dB`.
     * @public
     * @type {Int32Array}
     */
    A;

    /**
     * Histogram B for loudness distribution analysis.
     * Size determined by `GainAnalysis.STEPS_per_dB * GainAnalysis.MAX_dB`.
     * @public
     * @type {Int32Array}
     */
    B;

    constructor() {
        // Initialize arrays
        this.linprebuf = new_float(GainAnalysis.MAX_ORDER * 2);
        this.lstepbuf = new_float(GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER);
        this.loutbuf = new_float(GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER);
        this.rinprebuf = new_float(GainAnalysis.MAX_ORDER * 2);
        this.rstepbuf = new_float(GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER);
        this.routbuf = new_float(GainAnalysis.MAX_SAMPLES_PER_WINDOW + GainAnalysis.MAX_ORDER);

        // Ensure histogram arrays have non-zero size if constants are zero initially
        const histogramSize = Math.max(1, Math.floor(GainAnalysis.STEPS_per_dB * GainAnalysis.MAX_dB));
        this.A = new_int(histogramSize);
        this.B = new_int(histogramSize);

        // Other properties are initialized to defaults (0 or 0.0) by class field initializers
    }
}

export { ReplayGain };
// export default ReplayG
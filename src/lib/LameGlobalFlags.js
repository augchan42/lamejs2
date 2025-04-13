/**
 * @fileoverview Global flags and settings for the LAME encoder.
 * Ported from lame.h. Defines parameters controlling the encoding process.
 * Uses ES Module syntax.
 *
 * @module LameGlobalFlags
 */

import { MPEGMode } from './MPEGMode.js';
// Import common and destructure needed parts
import * as common from './common.js';
const { VbrMode, ShortBlock } = common; // Destructure needed enums/objects

// Assuming LameInternalFlags is defined elsewhere and imported if needed for type info
/** @typedef {import('./LameInternalFlags.js').default} LameInternalFlags */
/** @typedef {import('./VbrMode.js').default} VbrMode */ // Assuming VbrMode enum/object exists
/** @typedef {import('./ShortBlock.js').default} ShortBlock */ // Assuming ShortBlock enum/object exists

/**
 * @classdesc Global flags and settings controlling the LAME encoding process.
 * This structure holds configuration parameters set by the user application,
 * defining input audio properties, encoding mode, quality, VBR settings,
 * psychoacoustic parameters, filtering, etc. It also contains some internal
 * state variables populated by LAME during initialization and encoding.
 * @constructs LameGlobalFlags
 */
class LameGlobalFlags {

    /**
     * Internal flags and state variables used by the encoder.
     * @public
     * @type {LameInternalFlags|null}
     */
    internal_flags = null;

    /**
     * Class identifier (unused in JS version).
     * @public
     * @type {number}
     */
    class_id = 0;

    /* input description */

    /**
     * Total number of PCM samples in the input audio. Set by the application.
     * Default: -1 (unknown)
     * @public
     * @type {number}
     */
    num_samples = 0;

    /**
     * Number of channels in the input audio (e.g., 1 for mono, 2 for stereo).
     * Default: 2
     * @public
     * @type {number}
     */
    num_channels = 0;

    /**
     * Sample rate of the input audio in Hz.
     * Default: 44100
     * @public
     * @type {number}
     */
    in_samplerate = 0;

    /**
     * Sample rate of the output MP3 file in Hz. LAME usually chooses the best rate
     * based on the input rate and MPEG standards unless explicitly set.
     * Default: 0 (LAME chooses)
     * @public
     * @type {number}
     */
    out_samplerate = 0;

    /**
     * Global scaling factor applied to the input PCM data before encoding.
     * Default: 0.95 (historical?) - often set to 1.0.
     * @public
     * @type {number}
     */
    scale = 0.0;

    /**
     * Scaling factor applied specifically to the left channel (channel 0).
     * Overrides `scale` if set.
     * @public
     * @type {number}
     */
    scale_left = 0.0;

    /**
     * Scaling factor applied specifically to the right channel (channel 1).
     * Overrides `scale` if set.
     * @public
     * @type {number}
     */
    scale_right = 0.0;

    /* general control params */

    /**
     * Enable frame analysis for debugging/informational purposes.
     * Default: false
     * @public
     * @type {boolean}
     */
    analysis = false;

    /**
     * Add an Xing VBR header to the MP3 file for accurate duration/seeking.
     * Automatically enabled for VBR modes.
     * Default: false (but often enabled automatically)
     * @public
     * @type {boolean}
     */
    bWriteVbrTag = false;

    /**
     * Use LAME for decoding MP3 to WAV instead of encoding.
     * Default: false
     * @public
     * @type {boolean}
     */
    decode_only = false;

    /**
     * Encoding quality setting. Range: 0 (highest quality, slowest) to 9 (lowest quality, fastest).
     * Affects algorithm choices and speed/quality trade-offs.
     * Default: 5
     * @public
     * @type {number}
     */
    quality = 0;

    /**
     * MPEG encoding mode (e.g., STEREO, JOINT_STEREO, MONO). See `MPEGMode` enum.
     * Default: LAME chooses based on input channels and other settings.
     * @public
     * @type {MPEGMode}
     */
    mode = MPEGMode.STEREO; // Assuming MPEGMode.STEREO is the default intended

    /**
     * Force the use of Mid/Side stereo encoding. Only effective if `mode` is JOINT_STEREO.
     * Default: false
     * @public
     * @type {boolean}
     */
    force_ms = false;

    /**
     * Use the MP3 free format, which allows arbitrary bitrates.
     * Default: false
     * @public
     * @type {boolean}
     */
    free_format = false;

    /**
     * Calculate and store ReplayGain information in the output file.
     * Default: false
     * @public
     * @type {boolean}
     */
    findReplayGain = false;

    /**
     * Decode the MP3 stream during encoding (e.g., for gapless playback calculation).
     * Default: false
     * @public
     * @type {boolean}
     */
    decode_on_the_fly = false;

    /**
     * Automatically write ID3v2 tags (and potentially ID3v1) to the output file.
     * Default: true
     * @public
     * @type {boolean}
     */
    write_id3tag_automatic = false;

    /* Bitrate / Compression settings */

    /**
     * Target bitrate in kbps for CBR (Constant BitRate) mode. Set this *or* `compression_ratio`.
     * Default: 0 (use compression_ratio)
     * @public
     * @type {number}
     */
    brate = 0;

    /**
     * Target compression ratio (input size / output size). Set this *or* `brate`.
     * Default: 11.025 (approximates 128kbps for 44.1kHz stereo)
     * @public
     * @type {number}
     */
    compression_ratio = 0.0;

    /* frame params */

    /**
     * Mark the MP3 frame header's copyright bit.
     * Default: 0 (false)
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    copyright = 0;

    /**
     * Mark the MP3 frame header's original bit.
     * Default: 1 (true)
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    original = 0;

    /**
     * Set the MP3 frame header's 'private extension' bit (rarely used).
     * Default: 0
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    extension = 0;

    /**
     * Indicate that the input PCM data has pre-emphasis applied (rare).
     * Note: LAME's psychoacoustic model does not account for this.
     * Default: 0 (false)
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    emphasis = 0;

    /**
     * Add a CRC checksum to each MP3 frame for error detection. Adds 2 bytes per frame.
     * Default: 0 (false)
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    error_protection = 0;

    /**
     * Enforce strict adherence to ISO MPEG standards where possible.
     * Default: false
     * @public
     * @type {boolean}
     */
    strict_ISO = false;

    /* Bit Reservoir */

    /**
     * Disable the use of the bit reservoir. Primarily for debugging/analysis.
     * Default: false (reservoir enabled)
     * @public
     * @type {boolean}
     */
    disable_reservoir = false;

    /* quantization/noise shaping */

    /**
     * Quantization noise shaping comparison mode (for long blocks). See LAME documentation.
     * Default: 0
     * @public
     * @type {number}
     */
    quant_comp = 0;

    /**
     * Quantization noise shaping comparison mode (for short blocks). See LAME documentation.
     * Default: 0
     * @public
     * @type {number}
     */
    quant_comp_short = 0;

    /**
     * Enable experimental psychoacoustic tuning 'Y'.
     * Default: false
     * @public
     * @type {boolean}
     */
    experimentalY = false;

    /**
     * Select experimental psychoacoustic tuning 'Z' mode (0, 1, 2, 3).
     * Default: 0
     * @public
     * @type {number}
     */
    experimentalZ = 0;

    /**
     * Bitfield for activating various NSPSYTUNE experimental features.
     * Default: 0
     * @public
     * @type {number}
     */
    exp_nspsytune = 0;

    /**
     * LAME preset value used (if any). Informational. See LAME documentation for preset values.
     * Default: 0 (no preset specified)
     * @public
     * @type {number}
     */
    preset = 0;

    /* VBR control */

    /**
     * Variable BitRate mode selection. See `VbrMode` enum.
     * Default: `VbrMode.vbr_off` (CBR)
     * @public
     * @type {VbrMode}
     */
    VBR = null; // Should be initialized, e.g., VbrMode.vbr_off

    /**
     * Fractional part of VBR quality setting (for finer control between VBR_q levels). Range [0, 1).
     * Default: 0.0
     * @public
     * @type {number}
     */
    VBR_q_frac = 0.0;

    /**
     * VBR quality setting. Range [0=highest, 9=lowest]. Lower values use more bits.
     * Default: 4
     * @public
     * @type {number}
     */
    VBR_q = 0;

    /**
     * Target average bitrate in kbps for ABR (Average BitRate) mode.
     * Default: 0
     * @public
     * @type {number}
     */
    VBR_mean_bitrate_kbps = 0;

    /**
     * Minimum allowed bitrate in kbps for VBR/ABR modes.
     * Default: 0 (use LAME default, e.g., 8 or 32)
     * @public
     * @type {number}
     */
    VBR_min_bitrate_kbps = 0;

    /**
     * Maximum allowed bitrate in kbps for VBR/ABR modes.
     * Default: 0 (use LAME default, e.g., 320)
     * @public
     * @type {number}
     */
    VBR_max_bitrate_kbps = 0;

    /**
     * Strictly enforce `VBR_min_bitrate_kbps`. If false (default), the minimum
     * can be violated for periods of analog silence.
     * Default: 0 (false)
     * @public
     * @type {number|boolean} // Typically 0 or 1
     */
    VBR_hard_min = 0;

    /* resampling and filtering */

    /**
     * Lowpass filter cutoff frequency in Hz. 0 means LAME chooses based on bitrate/samplerate. -1 disables filter.
     * Default: 0
     * @public
     * @type {number}
     */
    lowpassfreq = 0;

    /**
     * Highpass filter cutoff frequency in Hz. 0 means LAME chooses. -1 disables filter.
     * Default: 0
     * @public
     * @type {number}
     */
    highpassfreq = 0;

    /**
     * Width of the lowpass filter transition band in Hz.
     * Default: 0 (use LAME default, typically 15% of cutoff freq)
     * @public
     * @type {number}
     */
    lowpasswidth = 0;

    /**
     * Width of the highpass filter transition band in Hz.
     * Default: 0 (use LAME default, typically 15% of cutoff freq)
     * @public
     * @type {number}
     */
    highpasswidth = 0;

    /* Psychoacoustics and advanced settings */

    /**
     * Masking adjustment (dB) for long blocks. Positive values lower thresholds (more masking).
     * Default: 0.0
     * @public
     * @type {number}
     */
    maskingadjust = 0.0;

    /**
     * Masking adjustment (dB) for short blocks.
     * Default: 0.0
     * @public
     * @type {number}
     */
    maskingadjust_short = 0.0;

    /**
     * Use only the Absolute Threshold of Hearing (ATH) for masking, ignore signal masking.
     * Default: false
     * @public
     * @type {boolean}
     */
    ATHonly = false;

    /**
     * Use only ATH for short blocks.
     * Default: false
     * @public
     * @type {boolean}
     */
    ATHshort = false;

    /**
     * Disable the use of the ATH entirely.
     * Default: false
     * @public
     * @type {boolean}
     */
    noATH = false;

    /**
     * Selects the ATH formula to use (0-4). See `ATHformula` documentation.
     * Default: 4 (uses ATHcurve) ? Check LAME defaults. Often defaults based on quality.
     * @public
     * @type {number}
     */
    ATHtype = 0;

    /**
     * Shape adjustment parameter for ATH formula type 4. Range [-100, 10].
     * Default: 0.0
     * @public
     * @type {number}
     */
    ATHcurve = 0.0;

    /**
     * Lower the calculated ATH by this amount in dB. Can compensate for quiet environments.
     * Default: 0.0
     * @public
     * @type {number}
     */
    ATHlower = 0.0;

    /**
     * Selects the ATH auto-adjustment algorithm type.
     * Default: -1 (disabled) ? Check LAME defaults.
     * @public
     * @type {number}
     */
    athaa_type = 0;

    /**
     * Selects the loudness calculation method for ATH auto-adjustment.
     * Default: 0 ? Check LAME defaults.
     * @public
     * @type {number}
     */
    athaa_loudapprox = 0;

    /**
     * Sensitivity parameter for ATH auto-adjustment (dB). Tunes the active region.
     * Default: 0.0
     * @public
     * @type {number}
     */
    athaa_sensitivity = 0.0;

    /**
     * Controls short block switching behavior. See `ShortBlock` enum.
     * Default: `short_block_t.short_block_allowed` ? Check LAME defaults.
     * @public
     * @type {ShortBlock}
     */
    short_blocks = null; // Should be initialized, e.g., ShortBlock.short_block_mixed

    /**
     * Enable the use of temporal masking psychoacoustic effects.
     * Default: true ? Check LAME defaults.
     * @public
     * @type {boolean}
     */
    useTemporal = false;

    /**
     * Inter-channel masking ratio. Reduces masking difference between channels. Range [0, 1].
     * Default: 0.0 (disabled)
     * @public
     * @type {number}
     */
    interChRatio = 0.0;

    /**
     * Naoki Shibata's Mid/Side masking adjustment factor. Range [-10, 10].
     * Positive values increase Side channel masking.
     * Default: 0.0
     * @public
     * @type {number}
     */
    msfix = 0.0;

    /**
     * Enable specific psychoacoustic tuning based on the `--tune` option in command-line LAME.
     * Default: false
     * @public
     * @type {boolean}
     */
    tune = false;

    /**
     * Auxiliary value used by the `--tune` option for specific adjustments.
     * Default: 0.0
     * @public
     * @type {number}
     */
    tune_value_a = 0.0;

    /* Internal variables (read-only for application) */

    /**
     * MPEG version used (0=MPEG-2/2.5, 1=MPEG-1). Set internally by LAME.
     * @public
     * @readonly
     * @type {number}
     */
    version = 0;

    /**
     * Encoder delay in samples (due to MDCT overlap). Set internally.
     * @public
     * @readonly
     * @type {number}
     */
    encoder_delay = 0;

    /**
     * Number of padding samples added at the end of the input stream. Set internally.
     * @public
     * @readonly
     * @type {number}
     */
    encoder_padding = 0;

    /**
     * Number of samples per MP3 frame. Set internally based on version and sample rate.
     * @public
     * @readonly
     * @type {number}
     */
    framesize = 0;

    /**
     * Counter for the number of frames encoded so far. Updated internally.
     * @public
     * @readonly
     * @type {number}
     */
    frameNum = 0;

    /**
     * Internal flag indicating if LAME allocated this structure.
     * @public
     * @readonly
     * @type {number}
     */
    lame_allocated_gfp = 0;

    // Constructor logic can initialize defaults if needed
    constructor() {
        // Set default values matching comments where available
        this.num_samples = -1;
        this.num_channels = 2;
        this.in_samplerate = 44100;
        this.out_samplerate = 0;
        this.scale = 1.0; // Changed from 0 to 1 based on common usage
        this.scale_left = 1.0;
        this.scale_right = 1.0;
        this.analysis = false;
        this.bWriteVbrTag = false; // Will be set automatically for VBR
        this.decode_only = false;
        this.quality = 5;
        this.mode = MPEGMode.STEREO; // Default mode
        this.force_ms = false;
        this.free_format = false;
        this.findReplayGain = false;
        this.decode_on_the_fly = false;
        this.write_id3tag_automatic = true; // Common default
        this.brate = 0;
        this.compression_ratio = 11.025; // ~128kbps default
        this.copyright = 0;
        this.original = 1;
        this.extension = 0;
        this.emphasis = 0;
        this.error_protection = 0;
        this.strict_ISO = false;
        this.disable_reservoir = false;
        this.quant_comp = 0;
        this.quant_comp_short = 0;
        this.experimentalY = false;
        this.experimentalZ = 0;
        this.exp_nspsytune = 0;
        this.preset = 0;
        this.VBR = 0; // Use imported VbrMode
        this.VBR_q_frac = 0.0;
        this.VBR_q = 4; // Common VBR default
        this.VBR_mean_bitrate_kbps = 0;
        this.VBR_min_bitrate_kbps = 0;
        this.VBR_max_bitrate_kbps = 0;
        this.VBR_hard_min = 0;
        this.lowpassfreq = 0;
        this.highpassfreq = 0;
        this.lowpasswidth = 0;
        this.highpasswidth = 0;
        this.maskingadjust = 0.0;
        this.maskingadjust_short = 0.0;
        this.ATHonly = false;
        this.ATHshort = false;
        this.noATH = false;
        this.ATHtype = 4; // Default to type 4 / ATHcurve ? Check LAME behavior
        this.ATHcurve = 0.0;
        this.ATHlower = 0.0;
        this.athaa_type = -1; // Default disabled?
        this.athaa_loudapprox = 0;
        this.athaa_sensitivity = 0.0;
        this.short_blocks = common.ShortBlock.short_block_allowed; // Sensible default? Or mixed?
        this.useTemporal = true; // Often enabled by default?
        this.interChRatio = 0.0;
        this.msfix = 0.0;
        this.tune = false;
        this.tune_value_a = 0.0;

        // Internal vars are initialized by LAME itself
        this.version = 0;
        this.encoder_delay = 0;
        this.encoder_padding = 0;
        this.framesize = 0;
        this.frameNum = 0;
        this.lame_allocated_gfp = 0;
        this.internal_flags = null;
    }
}

export { LameGlobalFlags }; // Use named export
export default LameGlobalFlags; // Also provide default export if needed
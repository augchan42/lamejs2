import { MPEGMode, MPEGModeValues } from './MPEGMode.js';
import type { Mp3Config } from './types.js'; // Assuming this path is correct
import { LameContext } from './LameContext.js';
import * as common from '../lib/common.js'; // Import common for VbrMode
// Import specific types needed from Lame if not relying on context's typing
// import { Lame } from '../lib/Lame.js';
// import { LameGlobalFlags } from '../lib/LameGlobalFlags.js';

// These imports might not be strictly necessary here if only using context
// but keep them if they are used for type hints elsewhere or were part of the original structure
/*
import { BitStream } from '../lib/BitStream.js';
import { GainAnalysis } from '../lib/GainAnalysis.js';
import { Quantize} from '../lib/Quantize.js';
import { VBRTag } from '../lib/VBRTag.js';
import { Version } from '../lib/Version.js';
import { ID3TagSpec } from '../lib/ID3TagSpec.js'; // Assuming this exists if needed
import { Reservoir } from '../lib/Reservoir.js';
import { Takehiro } from '../lib/Takehiro.js';
import { QuantizePVT } from '../lib/QuantizePVT.js';
import { Presets } from '../lib/Presets.js';
*/

export class Mp3Encoder {
    private context: LameContext;
    private gfp: any; // Keep as any, or import LameGlobalFlags type
    private maxSamples: number;
    private mp3buf: Uint8Array;
    private mp3buf_size: number;

    constructor(config: Mp3Config) {
        // Create context
        this.context = new LameContext();

        // Initialize LAME using context
        this.gfp = this.context.lame.lame_init();
        if (!this.gfp) {
            throw new Error('Failed to initialize LAME library.');
        }

        // Set configuration from Mp3Config interface
        this.gfp.num_channels = config.channels;
        this.gfp.in_samplerate = config.sampleRate;
        this.gfp.brate = config.bitRate; // Set bitrate for CBR

        // Map the mode directly
        const modeMap = {
            [MPEGMode.STEREO]: MPEGModeValues.STEREO,
            [MPEGMode.JOINT_STEREO]: MPEGModeValues.JOINT_STEREO,
            [MPEGMode.DUAL_CHANNEL]: MPEGModeValues.DUAL_CHANNEL,
            [MPEGMode.MONO]: MPEGModeValues.MONO,
            [MPEGMode.NOT_SET]: MPEGModeValues.NOT_SET
        };
        // Use STEREO as default if mode is undefined or NOT_SET
        this.gfp.mode = modeMap[config.mode ?? MPEGMode.STEREO];
        // Ensure mode is valid if explicitly set to NOT_SET or based on channels
        if (this.gfp.mode === MPEGModeValues.NOT_SET || this.gfp.num_channels === 1) {
             this.gfp.mode = (config.channels === 1) ? MPEGModeValues.MONO : MPEGModeValues.STEREO;
        }

        this.gfp.quality = config.quality ?? 3; // Use nullish coalescing for default

        // --- Assume CBR by default, no VBR settings from config ---
        this.gfp.VBR = common.VbrMode.vbr_off; // Set to CBR mode
        this.gfp.bWriteVbrTag = false; // VBR tag not needed for CBR

        // Other defaults (could potentially be added to Mp3Config if needed)
        this.gfp.disable_reservoir = true; // Commonly used default? Consider making configurable.
        this.gfp.write_id3tag_automatic = false; // Usually handled separately

        const retcode = this.context.lame.lame_init_params(this.gfp);
        if (retcode !== 0) {
            throw new Error(`Failed to initialize LAME encoder parameters (code ${retcode})`);
        }

        // Set up buffers
        // Allow configuration via maxBuffer, fall back to 1152
        this.maxSamples = config.maxBuffer || 1152;
        this.mp3buf_size = Math.floor(1.25 * this.maxSamples + 7200);
        this.mp3buf = new Uint8Array(this.mp3buf_size);
    }

    /**
     * Encodes a buffer of PCM audio samples.
     * Input can be Int16Array or Float32Array.
     * For stereo, provide equal length arrays for left and right channels.
     * For mono, provide the single channel data as the `left` argument.
     *
     * @public
     * @param {Float32Array | Int16Array} left - Array of PCM samples for the left channel (or mono channel). Float values should be in [-1.0, 1.0].
     * @param {Float32Array | Int16Array | null | undefined} [right] - Array of PCM samples for the right channel. Must be the same length and type as `left`. Required for stereo.
     * @returns {Int8Array} A Int8Array containing the encoded MP3 data for this buffer.
     * @throws {Error} If input arrays have mismatched types or lengths for stereo, or if right channel is missing for stereo.
     */
    public encodeBuffer(left: Float32Array | Int16Array, right?: Float32Array | Int16Array | null): Int8Array {
        let inputRight: Float32Array | Int16Array | null = right ?? null;
        const isFloat = left instanceof Float32Array;

        // Handle mono case and basic validation
        if (this.gfp.num_channels === 1) {
            inputRight = left; // LAME internal processing uses left for both in mono mode
        } else { // Stereo
            if (!inputRight) {
                throw new Error('Right channel buffer required for stereo encoding.');
            }
            if (left.length !== inputRight.length) {
                throw new Error('Left and right channel buffers must have the same length.');
            }
            if (left.constructor !== inputRight.constructor) {
                 throw new Error('Left and right channel buffers must have the same type (Float32Array or Int16Array).');
            }
        }

        // Resize MP3 buffer if needed
        if (left.length > this.maxSamples) {
            console.warn(`Input buffer size (${left.length}) exceeds configured maxSamples (${this.maxSamples}). Resizing MP3 buffer.`);
            this.maxSamples = left.length;
            this.mp3buf_size = Math.floor(1.25 * this.maxSamples + 7200);
            this.mp3buf = new Uint8Array(this.mp3buf_size);
        }

        let encodedSize = 0;

        // Call the appropriate LAME encoding function based on input type
        if (isFloat) {
            const rightFloat = (inputRight instanceof Float32Array) ? inputRight : (this.gfp.num_channels === 1 ? left as Float32Array : null);
            if(this.gfp.num_channels === 2 && !rightFloat){ throw new Error('Right Float32Array buffer needed for stereo float encoding.'); }

            encodedSize = this.context.lame.lame_encode_buffer_ieee_float(
                this.gfp, left as Float32Array, rightFloat, left.length,
                this.mp3buf, 0, this.mp3buf_size
            );
        } else { // Input is Int16Array
            const rightInt = (inputRight instanceof Int16Array) ? inputRight : (this.gfp.num_channels === 1 ? left as Int16Array : null);
            if(this.gfp.num_channels === 2 && !rightInt){ throw new Error('Right Int16Array buffer needed for stereo int encoding.'); }

            encodedSize = this.context.lame.lame_encode_buffer(
                this.gfp, left as Int16Array, rightInt, left.length,
                this.mp3buf, 0, this.mp3buf_size
            );
        }

         if (encodedSize < 0) {
             throw new Error(`LAME encoding failed with error code: ${encodedSize}`);
         }

        // Return the encoded data as Int8Array view
        return new Int8Array(this.mp3buf.buffer, this.mp3buf.byteOffset, encodedSize);
    }

    /**
     * Finalizes the MP3 stream, flushing any remaining data.
     * @public
     * @returns {Int8Array} A Int8Array containing the final MP3 data.
     */
    public flush(): Int8Array {
        if (!this.gfp) {
             console.warn("Encoder already closed or not initialized.");
             return new Int8Array(0);
        }
        const finalSize = this.context.lame.lame_encode_flush(
            this.gfp,
            this.mp3buf,
            0,
            this.mp3buf_size
        );
        if (finalSize < 0) {
             throw new Error(`LAME flush failed with error code: ${finalSize}`);
        }

        return new Int8Array(this.mp3buf.buffer, this.mp3buf.byteOffset, finalSize);
    }

    /**
     * Closes the encoder and releases resources (currently just calls flush).
     * @public
     */
    public close(): void {
        if (this.gfp) {
             try {
                 this.context.lame.lame_encode_flush(this.gfp, this.mp3buf, 0, this.mp3buf_size);
             } catch (e) {
                  console.error("Error during LAME close/flush:", e);
             }
             this.gfp = null; // Mark as closed
        }
    }
}
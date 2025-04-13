/**
 * @fileoverview Main LAME encoder interface and WAV header utilities.
 * Provides the primary `Mp3Encoder` class for encoding PCM data to MP3
 * and a `WavHeader` class/namespace for reading WAV file headers.
 * Uses ES Module syntax.
 *
 * @module Main
 */

// Import necessary modules using ES Module syntax
import * as common from './common.js';
import Lame from './Lame.js';
import Presets from './Presets.js';
import GainAnalysis from './GainAnalysis.js';
import QuantizePVT from './QuantizePVT.js';
import Quantize from './Quantize.js';
import Takehiro from './Takehiro.js';
import Reservoir from './Reservoir.js';
import { MPEGMode } from './MPEGMode.js';
import BitStream from './BitStream.js';
import { Encoder } from './Encoder.js'; // Assuming Encoder is a named export
import Version from './Version.js';
import VBRTag from './VBRTag.js';

// Destructure common utilities for easier access
const {
    // System, // Not used directly
    VbrMode,
    Float,
    ShortBlock,
    Util,
    Arrays,
    // new_array_n, // Used indirectly
    new_byte,
    // new_double, // Not used
    new_float,
    // new_float_n, // Used indirectly
    new_int,
    // new_int_n, // Not used directly
    assert
} = common;


// --- Internal Helper Classes (No JSDoc as requested) ---

class GetAudio {
    parse = null;
    mpg = null;
    setModules(parse2, mpg2) {
        this.parse = parse2;
        this.mpg = mpg2;
    }
}

class Parse {
    ver = null;
    id3 = null;
    pre = null;
    setModules(ver2, id32, pre2) {
        this.ver = ver2;
        this.id3 = id32;
        this.pre = pre2;
    }
}

class MPGLib {
    // Placeholder, implementation not shown/needed for export JSDoc
}

class ID3Tag {
    bits = null;
    ver = null;
    setModules(_bits, _ver) {
        this.bits = _bits;
        this.ver = _ver;
    }
}

// --- Main Encoder Class ---

/**
 * @classdesc Provides an interface to the LAME MP3 encoder.
 * Initializes the encoder with specified channel count, sample rate, and bitrate,
 * and provides methods to encode PCM buffers and finalize the MP3 stream.
 * @constructs Mp3Encoder
 * @param {number} channels - Number of input channels (e.g., 1 for mono, 2 for stereo).
 * @param {number} samplerate - Sample rate of the input audio in Hz (e.g., 44100).
 * @param {number} kbps - Target bitrate in kilobits per second (e.g., 128).
 */
function Mp3Encoder(channels, samplerate, kbps) {
    if (arguments.length !== 3) {
        console.error('WARN: Mp3Encoder(channels, samplerate, kbps) not specified, using defaults.');
        channels = 1;
        samplerate = 44100;
        kbps = 128;
    }

    // Instantiate LAME components
    const lame = new Lame();
    const gaud = new GetAudio(); // Seems unused after setup
    const ga = new GainAnalysis();
    const bs = new BitStream();
    const p = new Presets();
    const qupvt = new QuantizePVT();
    const qu = new Quantize();
    const vbr = new VBRTag();
    const ver = new Version();
    const id3 = new ID3Tag();
    const rv = new Reservoir();
    const tak = new Takehiro();
    const parse = new Parse(); // Seems unused after setup
    const mpg = new MPGLib(); // Seems unused after setup

    // --- Module Wiring ---
    lame.setModules(ga, bs, p, qupvt, qu, vbr, ver, id3, mpg);
    bs.setModules(ga, mpg, ver, vbr); // Pass necessary modules to BitStream
    id3.setModules(bs, ver);
    p.setModules(lame);
    qu.setModules(bs, rv, qupvt, tak);
    qupvt.setModules(tak, rv, lame.enc.psy); // Pass psy model from lame instance
    rv.setModules(bs);
    tak.setModules(qupvt);
    vbr.setModules(lame, bs, ver);
    // gaud.setModules(parse, mpg); // Not used?
    // parse.setModules(ver, id3, p); // Not used?

    // --- LAME Initialization ---
    const gfp = lame.lame_init();

    // Configure LAME flags based on constructor arguments
    gfp.num_channels = channels;
    gfp.in_samplerate = samplerate;
    gfp.brate = kbps;
    gfp.mode = (channels === 1) ? MPEGMode.MONO : MPEGMode.STEREO; // Set mode based on channels
    gfp.quality = 3; // Set a default quality level (0=best, 9=fastest)
    gfp.bWriteVbrTag = false; // Disable VBR tag by default for CBR
    gfp.disable_reservoir = true; // Often recommended for JS streaming? Check implications.
    gfp.write_id3tag_automatic = false; // Don't write ID3 tags by default

    // Initialize LAME parameters
    const retcode = lame.lame_init_params(gfp);
    assert(retcode === 0, `lame_init_params failed with code ${retcode}`);

    // --- Buffer Setup ---
    let maxSamples = 1152; // Initial max samples per encode buffer call (1 MPEG frame)
    // Calculate MP3 buffer size (estimate from LAME examples)
    let mp3buf_size = Math.floor(1.25 * maxSamples + 7200);
    let mp3buf = new_byte(mp3buf_size); // Use common.js utility for Uint8Array

    /**
     * Encodes a buffer of PCM audio samples.
     * For stereo, provide equal length arrays for left and right channels.
     * For mono, provide the single channel data as the `left` argument and null/empty for `right`.
     *
     * @public
     * @param {Float32Array} left - Array of PCM samples for the left channel (or mono channel). Values should be in the range [-1.0, 1.0].
     * @param {Float32Array|null} right - Array of PCM samples for the right channel. Must be the same length as `left`. Ignored if `channels` is 1.
     * @returns {Int8Array} A Int8Array containing the encoded MP3 data for this buffer.
     * @throws {Error} If input array lengths differ for stereo encoding.
     */
    this.encodeBuffer = function (left, right = null) {
        let numSamples = left.length;

        if (gfp.num_channels === 1) {
            // For mono, LAME C API expects the same buffer for left and right,
            // or NULL for right. Here, we can just pass the left buffer again.
            right = left;
        } else {
            // For stereo, ensure right channel is provided and lengths match
            if (!right || left.length !== right.length) {
                throw new Error("Mismatched channel buffer lengths or missing right channel for stereo encoding.");
            }
        }

        // Resize MP3 buffer if input sample count is larger than previous calls
        if (numSamples > maxSamples) {
            maxSamples = numSamples;
            mp3buf_size = Math.floor(1.25 * maxSamples + 7200);
            mp3buf = new_byte(mp3buf_size);
        }

        // Call the core LAME encoding function
        const bytesEncoded = lame.lame_encode_buffer_ieee_float(
            gfp,
            left,
            right,
            numSamples,
            mp3buf,
            0, // mp3buf_offset
            mp3buf_size
        );

        // Return the relevant portion of the MP3 buffer
        // Create a new Int8Array view of the subarray to avoid holding reference to large buffer
        return new Int8Array(mp3buf.buffer, mp3buf.byteOffset, bytesEncoded);
    };

    /**
     * Finalizes the MP3 stream. Call this after processing all audio buffers.
     * Returns any remaining data buffered within LAME (e.g., final frame padding).
     *
     * @public
     * @returns {Int8Array} A Int8Array containing the final MP3 data.
     */
    this.flush = function () {
        const bytesEncoded = lame.lame_encode_flush(
            gfp,
            mp3buf,
            0, // mp3buf_offset
            mp3buf_size
        );

        // Return the relevant portion of the MP3 buffer
        return new Int8Array(mp3buf.buffer, mp3buf.byteOffset, bytesEncoded);
    };
}

// --- WAV Header Class ---

/**
 * @classdesc Represents essential information read from a WAV file header.
 * @constructs WavHeader
 */
function WavHeader() {
    /**
     * Byte offset where the raw audio data begins in the WAV file.
     * @public
     * @type {number}
     */
    this.dataOffset = 0;
    /**
     * Length of the raw audio data chunk in bytes.
     * @public
     * @type {number}
     */
    this.dataLen = 0;
    /**
     * Number of audio channels (e.g., 1 for mono, 2 for stereo).
     * @public
     * @type {number}
     */
    this.channels = 0;
    /**
     * Sample rate in Hz (e.g., 44100).
     * @public
     * @type {number}
     */
    this.sampleRate = 0;
}

/**
 * Helper function to convert a 4-character string (FourCC) to a 32-bit integer.
 * @private
 * @param {string} fourcc - The 4-character string.
 * @returns {number} The integer representation.
 */
function fourccToInt(fourcc) {
    // Ensure string has 4 characters? Or assume valid input.
    return (fourcc.charCodeAt(0) << 24) | (fourcc.charCodeAt(1) << 16) | (fourcc.charCodeAt(2) << 8) | fourcc.charCodeAt(3);
}

/**
 * RIFF chunk identifier (Big-Endian).
 * @public
 * @static
 * @readonly
 * @type {number}
 */
WavHeader.RIFF = fourccToInt("RIFF");

/**
 * WAVE format identifier (Big-Endian).
 * @public
 * @static
 * @readonly
 * @type {number}
 */
WavHeader.WAVE = fourccToInt("WAVE");

/**
 * fmt chunk identifier (Big-Endian).
 * @public
 * @static
 * @readonly
 * @type {number}
 */
WavHeader.fmt_ = fourccToInt("fmt ");

/**
 * data chunk identifier (Big-Endian).
 * @public
 * @static
 * @readonly
 * @type {number}
 */
WavHeader.data = fourccToInt("data");

/**
 * Reads and parses the header of a WAV file from a DataView.
 * Populates and returns a WavHeader object. Returns null if the header is invalid.
 * Supports basic PCM fmt chunks (16 or 18 bytes).
 *
 * @public
 * @static
 * @param {DataView} dataView - A DataView representing the beginning of the WAV file data (at least ~44 bytes).
 * @returns {WavHeader|null} A WavHeader object containing the parsed information, or null if the header is not a valid basic WAV header.
 * @throws {Error} If an unsupported extended fmt chunk is encountered.
 */
WavHeader.readHeader = function (dataView) {
    const w = new WavHeader();

    // Basic checks for RIFF/WAVE header
    if (dataView.byteLength < 44) return null; // Need at least minimum header size
    if (WavHeader.RIFF !== dataView.getUint32(0, false)) return null;
    // const fileLen = dataView.getUint32(4, true); // File length (not strictly needed here)
    if (WavHeader.WAVE !== dataView.getUint32(8, false)) return null;

    // Find 'fmt ' chunk
    let pos = 12;
    while (pos < dataView.byteLength - 8) {
        if (WavHeader.fmt_ === dataView.getUint32(pos, false)) break;
        const chunkLen = dataView.getUint32(pos + 4, true);
        pos += 8 + chunkLen; // Skip chunk identifier, length, and data
        if (chunkLen % 2 !== 0) pos++; // Skip padding byte if chunk length is odd
    }
    if (pos >= dataView.byteLength - 8) return null; // 'fmt ' not found

    const fmtLen = dataView.getUint32(pos + 4, true);
    pos += 8; // Move to start of fmt chunk data

    // Parse basic fmt chunk data
    if (fmtLen === 16 || fmtLen === 18) { // 18 includes cbSize which we ignore
        const audioFormat = dataView.getUint16(pos, true);
        if (audioFormat !== 1) return null; // Only support PCM format 1
        w.channels = dataView.getUint16(pos + 2, true);
        w.sampleRate = dataView.getUint32(pos + 4, true);
        // Skip byteRate, blockAlign, bitsPerSample
    } else {
        // Basic implementation doesn't support extended fmt chunks
        // console.warn(`Extended WAV fmt chunk (length ${fmtLen}) not fully supported.`);
        // Attempt to read basic info anyway if possible
        const audioFormat = dataView.getUint16(pos, true);
        if (audioFormat !== 1) return null;
        w.channels = dataView.getUint16(pos + 2, true);
        w.sampleRate = dataView.getUint32(pos + 4, true);
        // Or throw error as per original code comment?
        // throw new Error('Extended fmt chunk not implemented');
    }
    pos += fmtLen; // Move past fmt chunk data
    if (fmtLen % 2 !== 0) pos++; // Skip padding byte

    // Find 'data' chunk
    while (pos < dataView.byteLength - 8) {
        if (WavHeader.data === dataView.getUint32(pos, false)) break;
        const chunkLen = dataView.getUint32(pos + 4, true);
        pos += 8 + chunkLen; // Skip chunk
         if (chunkLen % 2 !== 0) pos++; // Skip padding
    }
    if (pos >= dataView.byteLength - 8) return null; // 'data' not found

    // Store data chunk info
    w.dataLen = dataView.getUint32(pos + 4, true);
    w.dataOffset = pos + 8; // Offset to the actual audio data

    // Basic validation
    if (w.channels === 0 || w.sampleRate === 0 || w.dataOffset === 0) {
        return null; // Invalid header data
    }

    return w;
};

// --- Exports ---
export { Mp3Encoder, WavHeader };
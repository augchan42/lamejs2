declare const enum MPEGMode {
    STEREO = 0,
    JOINT_STEREO = 1,
    DUAL_CHANNEL = 2,
    MONO = 3,
    NOT_SET = 4
}

interface Mp3Config {
    channels: number;
    sampleRate: number;
    bitRate: number;
    mode?: MPEGMode;
    quality?: number;
    maxBuffer?: number;
}
interface EncodedChunk {
    data: Int8Array;
    sampleRate: number;
    channels: number;
}
interface WavHeader {
    dataOffset: number;
    dataLen: number;
    channels: number;
    sampleRate: number;
}

declare class Mp3Encoder {
    private context;
    private gfp;
    private maxSamples;
    private mp3buf;
    private mp3buf_size;
    constructor(config: Mp3Config);
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
    encodeBuffer(left: Float32Array | Int16Array, right?: Float32Array | Int16Array | null): Int8Array;
    /**
     * Finalizes the MP3 stream, flushing any remaining data.
     * @public
     * @returns {Int8Array} A Int8Array containing the final MP3 data.
     */
    flush(): Int8Array;
    /**
     * Closes the encoder and releases resources (currently just calls flush).
     * @public
     */
    close(): void;
}

export { type EncodedChunk, MPEGMode, type Mp3Config, Mp3Encoder, type WavHeader };

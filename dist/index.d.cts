declare enum MPEGMode {
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
    private lame;
    private gfp;
    private maxSamples;
    private mp3buf;
    private mp3buf_size;
    constructor(config: Mp3Config);
    encodeBuffer(left: Float32Array | Int16Array, right?: Float32Array | Int16Array): Int8Array;
    flush(): Int8Array;
    close(): void;
}

export { type EncodedChunk, MPEGMode, type Mp3Config, Mp3Encoder, type WavHeader };

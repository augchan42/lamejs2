import { MPEGMode, MPEGModeValues } from './MPEGMode.js';
import type { Mp3Config } from './types.js';
import { LameContext } from './LameContext';

const Lame = require('../lib/Lame');
const BitStream = require('../lib/BitStream');
const GainAnalysis = require('../lib/GainAnalysis');
const Quantize = require('../lib/Quantize');
const VBRTag = require('../lib/VBRTag');
const Version = require('../lib/Version');
const ID3TagSpec = require('../lib/ID3TagSpec');
const Reservoir = require('../lib/Reservoir');
const Takehiro = require('../lib/Takehiro');
const QuantizePVT = require('../lib/QuantizePVT');
const Presets = require('../lib/Presets');

export class Mp3Encoder {
  private context: LameContext;
  private gfp: any; // Keep as any since it's internal to lamejs
  private maxSamples: number;
  private mp3buf: Uint8Array;
  private mp3buf_size: number;

  constructor(config: Mp3Config) {
    // Create context
    this.context = new LameContext();
    
    // Initialize LAME using context
    this.gfp = this.context.lame.lame_init();

    // Set configuration
    this.gfp.num_channels = config.channels;
    this.gfp.in_samplerate = config.sampleRate;
    this.gfp.brate = config.bitRate;
    
    // Map the mode directly
    const modeMap = {
      [MPEGMode.STEREO]: MPEGModeValues.STEREO,
      [MPEGMode.JOINT_STEREO]: MPEGModeValues.JOINT_STEREO,
      [MPEGMode.DUAL_CHANNEL]: MPEGModeValues.DUAL_CHANNEL,
      [MPEGMode.MONO]: MPEGModeValues.MONO,
      [MPEGMode.NOT_SET]: MPEGModeValues.NOT_SET
    };
    this.gfp.mode = modeMap[config.mode ?? MPEGMode.STEREO];

    this.gfp.quality = config.quality || 3;
    this.gfp.bWriteVbrTag = false;
    this.gfp.disable_reservoir = true;
    this.gfp.write_id3tag_automatic = false;

    const retcode = this.context.lame.lame_init_params(this.gfp);
    if (retcode !== 0) {
      throw new Error('Failed to initialize LAME encoder');
    }

    // Set up buffers
    this.maxSamples = config.maxBuffer || 1152;
    this.mp3buf_size = Math.floor(1.25 * this.maxSamples + 7200);
    this.mp3buf = new Uint8Array(this.mp3buf_size);
  }

  public encodeBuffer(left: Float32Array | Int16Array, right?: Float32Array | Int16Array): Int8Array {
    // Handle mono case
    if (this.gfp.num_channels === 1) {
      right = left;
    }

    if (!right) {
      throw new Error('Right channel required for stereo');
    }

    if (left.length !== right.length) {
      throw new Error('Left and right channels must be same length');
    }

    // Resize buffer if needed
    if (left.length > this.maxSamples) {
      this.maxSamples = left.length;
      this.mp3buf_size = Math.floor(1.25 * this.maxSamples + 7200);
      this.mp3buf = new Uint8Array(this.mp3buf_size);
    }

    const encodedSize = this.context.lame.lame_encode_buffer(
      this.gfp,
      left,
      right,
      left.length,
      this.mp3buf,
      0,
      this.mp3buf_size
    );

    return new Int8Array(this.mp3buf.subarray(0, encodedSize));
  }

  public flush(): Int8Array {
    const finalSize = this.context.lame.lame_encode_flush(
      this.gfp,
      this.mp3buf,
      0,
      this.mp3buf_size
    );

    return new Int8Array(this.mp3buf.subarray(0, finalSize));
  }

  public close(): void {
    if (this.gfp) {
      this.context.lame.lame_encode_flush(this.gfp, this.mp3buf, 0, this.mp3buf_size);
      this.gfp = null;
    }
  }
} 
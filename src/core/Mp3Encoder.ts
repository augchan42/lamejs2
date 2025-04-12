import { MPEGMode } from './MPEGMode.js';
import type { Mp3Config } from './types.js';
import { createRequire } from 'module';

const require = createRequire(import.meta.url);
const Lame = require('../lib/Lame');
const BitStream = require('../lib/BitStream');
const GainAnalysis = require('../lib/GainAnalysis');
const Quantize = require('../lib/Quantize');
const VBRTag = require('../lib/VBRTag');
const Version = require('../lib/Version');
const ID3Tag = require('../lib/ID3TagSpec');
const Reservoir = require('../lib/Reservoir');
const Takehiro = require('../lib/Takehiro');
const QuantizePVT = require('../lib/QuantizePVT');
const Presets = require('../lib/Presets');

export class Mp3Encoder {
  private lame: typeof Lame;
  private gfp: any; // Keep as any since it's internal to lamejs
  private maxSamples: number;
  private mp3buf: Uint8Array;
  private mp3buf_size: number;

  constructor(config: Mp3Config) {
    // Initialize all required modules
    const lame = new Lame();
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

    // Set up module relationships
    lame.setModules(ga, bs, p, qupvt, qu, vbr, ver, id3);
    bs.setModules(ga, ver, vbr);
    id3.setModules(bs, ver);
    p.setModules(lame);
    qu.setModules(bs, rv, qupvt, tak);
    qupvt.setModules(tak, rv, lame.enc.psy);
    rv.setModules(bs);
    tak.setModules(qupvt);
    vbr.setModules(lame, bs, ver);

    this.lame = lame;
    
    // Initialize LAME
    this.gfp = lame.lame_init();

    // Set configuration
    this.gfp.num_channels = config.channels;
    this.gfp.in_samplerate = config.sampleRate;
    this.gfp.brate = config.bitRate;
    this.gfp.mode = config.mode || MPEGMode.STEREO;
    this.gfp.quality = config.quality || 3;
    this.gfp.bWriteVbrTag = false;
    this.gfp.disable_reservoir = true;
    this.gfp.write_id3tag_automatic = false;

    const retcode = lame.lame_init_params(this.gfp);
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

    const encodedSize = this.lame.lame_encode_buffer(
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
    const finalSize = this.lame.lame_encode_flush(
      this.gfp,
      this.mp3buf,
      0,
      this.mp3buf_size
    );

    return new Int8Array(this.mp3buf.subarray(0, finalSize));
  }

  public close(): void {
    if (this.gfp) {
      this.lame.lame_close(this.gfp);
      this.gfp = null;
    }
  }
} 
"use strict";
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};
var __copyProps = (to, from, except, desc) => {
  if (from && typeof from === "object" || typeof from === "function") {
    for (let key of __getOwnPropNames(from))
      if (!__hasOwnProp.call(to, key) && key !== except)
        __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
  }
  return to;
};
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// src/index.ts
var index_exports = {};
__export(index_exports, {
  MPEGMode: () => MPEGMode,
  Mp3Encoder: () => Mp3Encoder
});
module.exports = __toCommonJS(index_exports);

// src/core/MPEGMode.ts
var MPEGMode = /* @__PURE__ */ ((MPEGMode2) => {
  MPEGMode2[MPEGMode2["STEREO"] = 0] = "STEREO";
  MPEGMode2[MPEGMode2["JOINT_STEREO"] = 1] = "JOINT_STEREO";
  MPEGMode2[MPEGMode2["DUAL_CHANNEL"] = 2] = "DUAL_CHANNEL";
  MPEGMode2[MPEGMode2["MONO"] = 3] = "MONO";
  MPEGMode2[MPEGMode2["NOT_SET"] = 4] = "NOT_SET";
  return MPEGMode2;
})(MPEGMode || {});

// src/core/Mp3Encoder.ts
var import_module = require("module");
var import_meta = {};
var require2 = (0, import_module.createRequire)(import_meta.url);
var Lame = require2("../lib/Lame");
var BitStream = require2("../lib/BitStream");
var GainAnalysis = require2("../lib/GainAnalysis");
var Quantize = require2("../lib/Quantize");
var VBRTag = require2("../lib/VBRTag");
var Version = require2("../lib/Version");
var ID3Tag = require2("../lib/ID3TagSpec");
var Reservoir = require2("../lib/Reservoir");
var Takehiro = require2("../lib/Takehiro");
var QuantizePVT = require2("../lib/QuantizePVT");
var Presets = require2("../lib/Presets");
var Mp3Encoder = class {
  constructor(config) {
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
    this.gfp = lame.lame_init();
    this.gfp.num_channels = config.channels;
    this.gfp.in_samplerate = config.sampleRate;
    this.gfp.brate = config.bitRate;
    this.gfp.mode = config.mode || 0 /* STEREO */;
    this.gfp.quality = config.quality || 3;
    this.gfp.bWriteVbrTag = false;
    this.gfp.disable_reservoir = true;
    this.gfp.write_id3tag_automatic = false;
    const retcode = lame.lame_init_params(this.gfp);
    if (retcode !== 0) {
      throw new Error("Failed to initialize LAME encoder");
    }
    this.maxSamples = config.maxBuffer || 1152;
    this.mp3buf_size = Math.floor(1.25 * this.maxSamples + 7200);
    this.mp3buf = new Uint8Array(this.mp3buf_size);
  }
  encodeBuffer(left, right) {
    if (this.gfp.num_channels === 1) {
      right = left;
    }
    if (!right) {
      throw new Error("Right channel required for stereo");
    }
    if (left.length !== right.length) {
      throw new Error("Left and right channels must be same length");
    }
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
  flush() {
    const finalSize = this.lame.lame_encode_flush(
      this.gfp,
      this.mp3buf,
      0,
      this.mp3buf_size
    );
    return new Int8Array(this.mp3buf.subarray(0, finalSize));
  }
  close() {
    if (this.gfp) {
      this.lame.lame_close(this.gfp);
      this.gfp = null;
    }
  }
};
// Annotate the CommonJS export names for ESM import in node:
0 && (module.exports = {
  MPEGMode,
  Mp3Encoder
});
//# sourceMappingURL=index.cjs.map
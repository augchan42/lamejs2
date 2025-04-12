// First, declare the base module interface
interface LameContext {
  version: Version;
  bitStream: BitStream;
  gainAnalysis: GainAnalysis;
  quantize: Quantize;
  vbrTag: VBRTag;
  id3: ID3Tag;
  reservoir: Reservoir;
  takehiro: Takehiro;
  quantizePVT: QuantizePVT;
  presets: Presets;
  lame: Lame;
}

// Declare each module class
declare module './Lame.js' {
  class Lame {
    public constructor();
    public setModules(ga: any, bs: any, p: any, qupvt: any, qu: any, vbr: any, ver: any, id3: any, mpg: any): void;
  }
  export = Lame;
}

declare module './BitStream.js' {
  class BitStream {
    public constructor(context: LameContext);
    public setModules(ga: any, mpg: any, ver: any, vbr: any): void;
  }
  export = BitStream;
}

declare module './GainAnalysis.js' {
  class GainAnalysis {
    public constructor();
  }
  export = GainAnalysis;
}

declare module './Quantize.js' {
  class Quantize {
    public rv: any;
    public qupvt: any;
    public constructor();
    public setModules(bs: import('./BitStream.js'), rv: import('./Reservoir.js'), qupvt: import('./QuantizePVT.js'), tk: import('./Takehiro.js')): void;
    public ms_convert(l3_side: any, gr: number): void;
    public init_xrpow(gfc: any, cod_info: any, xrpow: any): boolean;
  }
  export = Quantize;
}

declare module './VBRTag.js' {
  class VBRTag {
    public constructor();
    public setModules(lame: any, bs: any, ver: any): void;
  }
  export = VBRTag;
}

declare module './Version.js' {
  class Version {
    public constructor(context: LameContext);
  }
  export = Version;
}

declare module './ID3TagSpec.js' {
  class ID3TagSpec {
    public constructor();
    public setModules(bs: any, ver: any): void;
  }
  export = ID3TagSpec;
}

declare module './Reservoir.js' {
  class Reservoir {
    public constructor();
    public setModules(bs: any): void;
  }
  export = Reservoir;
}

declare module './Takehiro.js' {
  class Takehiro {
    public qupvt: import('./QuantizePVT.js').QuantizePVT;
    public constructor(context: LameContext);
    public setModules(qupvt: import('./QuantizePVT.js').QuantizePVT): void;
    public noquant_count_bits(gfc: any, gi: any, prev_noise: any): number;
    public count_bits(gfc: any, xr: any, gi: any, prev_noise: any): number;
    public best_huffman_divide(gfc: any, gi: any): void;
    public best_scalefac_store(gfc: any, gr: number, ch: number, l3_side: any): void;
    public scale_bitcount(cod_info: any): boolean;
    public scale_bitcount_lsf(gfc: any, cod_info: any): boolean;
    public huffman_init(gfc: any): void;
  }
  export = Takehiro;
}

declare module './QuantizePVT.js' {
  class QuantizePVT {
    public constructor();
    public setModules(tk: any, rv: any, psy: any): void;
  }
  export = QuantizePVT;
}

declare module './Presets.js' {
  class Presets {
    public constructor();
    public setModules(lame: any): void;
  }
  export = Presets;
}

declare module './Encoder.js' {
  class Encoder {
    public constructor(context: LameContext);
    public static NORM_TYPE: number;
    public setModules(...args: any[]): void;
    public lame_encode_mp3_frame(gfc: any): void;
  }
  export = Encoder;
}

declare module './PsyModel.js' {
  class PsyModel {
    public constructor();
  }
  export = PsyModel;
} 
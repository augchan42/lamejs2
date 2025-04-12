declare module '../lib/*' {
  const mod: any;
  export default mod;
}

declare module '../lib/Lame' {
  class Lame {
    setModules(...args: any[]): void;
    lame_init(): any;
    lame_init_params(gfp: any): number;
    lame_encode_buffer(gfp: any, left: any, right: any, len: number, mp3buf: Uint8Array, offset: number, size: number): number;
    lame_encode_flush(gfp: any, mp3buf: Uint8Array, offset: number, size: number): number;
    lame_close(gfp: any): void;
    enc: { psy: any };
  }
  const instance: Lame;
  export = instance;
}

declare module '../lib/BitStream' {
  export default class BitStream {
    setModules(...args: any[]): void;
  }
}

declare module '../lib/GainAnalysis' { export default class GainAnalysis { } }
declare module '../lib/Quantize' {
  export default class Quantize {
    setModules(...args: any[]): void;
  }
}
declare module '../lib/VBRTag' { export default class VBRTag { } }
declare module '../lib/Version' { export default class Version { } }
declare module '../lib/ID3TagSpec' { export default class ID3Tag { } }
declare module '../lib/Reservoir' { export default class Reservoir { } }
declare module '../lib/Takehiro' {
  export default class Takehiro {
    setModules(...args: any[]): void;
  }
}
declare module '../lib/QuantizePVT' { export default class QuantizePVT { } }
declare module '../lib/Presets' { export default class Presets { } } 
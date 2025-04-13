import { Version } from '../lib/Version.js';
import { BitStream } from '../lib/BitStream.js';
import { GainAnalysis } from '../lib/GainAnalysis.js';
import { Quantize} from '../lib/Quantize.js';
import { VBRTag } from '../lib/VBRTag.js';
import { ID3TagSpec } from '../lib/ID3TagSpec.js';
import { Reservoir } from '../lib/Reservoir.js';
import { Takehiro } from '../lib/Takehiro.js';
import { QuantizePVT } from '../lib/QuantizePVT.js';
import { Presets } from '../lib/Presets.js';
import { Lame } from '../lib/Lame.js';
import { PsyModel } from '../lib/PsyModel.js';

export class LameContext {
    // All module instances
    public version: Version;
    public bitStream: BitStream;
    public gainAnalysis: GainAnalysis;
    public quantize: Quantize;
    public vbrTag: VBRTag;
    public id3: ID3TagSpec;
    public reservoir: Reservoir;
    public takehiro: Takehiro;
    public quantizePVT: QuantizePVT;
    public presets: Presets;
    public lame: Lame;
    public psyModel: PsyModel;

    constructor() {
        // Create instances first
        this.id3 = new ID3TagSpec();
        this.reservoir = new Reservoir();
        this.quantizePVT = new QuantizePVT();
        this.takehiro = new Takehiro(this);
        this.quantize = new Quantize();
        this.vbrTag = new VBRTag();
        this.presets = new Presets();
        this.lame = new Lame();
        this.version = new Version(this);
        this.bitStream = new BitStream(this);
        this.gainAnalysis = new GainAnalysis();
        this.psyModel = new PsyModel();

        // Then initialize modules with their dependencies
        this.quantizePVT.setModules(this.takehiro, this.reservoir, this.psyModel);
        this.takehiro.setModules(this.quantizePVT);
        this.quantize.setModules(this.bitStream, this.reservoir, this.quantizePVT, this.takehiro);
        this.presets.setModules(this.lame);
        this.bitStream.setModules(this.gainAnalysis, null, this.version, this.vbrTag);
        this.id3.setModules(this.bitStream, this.version);
        this.reservoir.setModules(this.bitStream);
        this.vbrTag.setModules(this.lame, this.bitStream, this.version);
        this.lame.setModules(
            this.gainAnalysis,
            this.bitStream,
            this.presets,            
            this.quantizePVT,
            this.quantize, 
            this.vbrTag,
            this.version,
            this.id3,
            null
        );
    }
} 
import { MPEGMode } from './MPEGMode.js';

export interface Mp3Config {
  channels: number;
  sampleRate: number;
  bitRate: number;
  mode?: MPEGMode;
  quality?: number; // 0 (best) to 9 (worst)
  maxBuffer?: number;
}

export interface EncodedChunk {
  data: Int8Array;
  sampleRate: number;
  channels: number;
}

export interface WavHeader {
  dataOffset: number;
  dataLen: number;
  channels: number;
  sampleRate: number;
} 
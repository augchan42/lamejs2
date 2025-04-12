import { describe, it, expect } from 'vitest';
import { Mp3Encoder, MPEGMode } from '../src';

describe('Mp3Encoder', () => {
  it('should encode mono audio', () => {
    const encoder = new Mp3Encoder({
      channels: 1,
      sampleRate: 44100,
      bitRate: 128,
      mode: MPEGMode.MONO
    });

    // Create 1 second of 440Hz sine wave
    const samples = new Float32Array(44100);
    for (let i = 0; i < samples.length; i++) {
      samples[i] = Math.sin(2 * Math.PI * 440 * i / 44100);
    }

    const mp3Data = encoder.encodeBuffer(samples);
    const finalData = encoder.flush();

    expect(mp3Data).toBeInstanceOf(Int8Array);
    expect(mp3Data.length).toBeGreaterThan(0);
    expect(finalData.length).toBeGreaterThan(0);
  });

  it('should encode stereo audio', () => {
    const encoder = new Mp3Encoder({
      channels: 2,
      sampleRate: 44100,
      bitRate: 128
    });

    const left = new Float32Array(44100);
    const right = new Float32Array(44100);
    
    const mp3Data = encoder.encodeBuffer(left, right);
    expect(mp3Data).toBeInstanceOf(Int8Array);
  });
}); 
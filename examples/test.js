import { Mp3Encoder, MPEGMode } from 'lamejs2';

// Create a simple mono signal (1 second of 440Hz sine wave)
const sampleRate = 44100;
const samples = new Float32Array(sampleRate);
for (let i = 0; i < samples.length; i++) {
  samples[i] = Math.sin(2 * Math.PI * 440 * i / sampleRate);
}

// Create encoder
const encoder = new Mp3Encoder({
  channels: 1,
  sampleRate: 44100,
  bitRate: 128,
  mode: MPEGMode.MONO
});

// Encode and get MP3 data
const mp3Data = encoder.encodeBuffer(samples);
const finalData = encoder.flush();

// In browser, create blob
const blob = new Blob([mp3Data, finalData], { type: 'audio/mp3' });
console.log('MP3 size:', blob.size, 'bytes'); 
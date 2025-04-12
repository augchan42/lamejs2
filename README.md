# lamejs2

A modern TypeScript rewrite of [lamejs](https://github.com/zhuker/lamejs), the fast JavaScript MP3 encoder. This version maintains the same high performance while adding TypeScript support, better error handling, and a more modern API.

## Features

- 🎯 Full TypeScript support
- 📦 Modern ESM and CommonJS support
- 🚀 Same high performance as lamejs
- 💪 Improved error handling
- 🔧 Modern build tooling
- 📘 Better type definitions

## Installation

```bash
npm install lamejs2
```

## Quick Start

```typescript
import { Mp3Encoder, MPEGMode } from 'lamejs2';

// Create encoder (mono 44.1khz @ 128kbps)
const encoder = new Mp3Encoder(1, 44100, 128);

// Encode some audio samples
const samples = new Int16Array(44100); // your audio data here
const mp3Data = encoder.encodeBuffer(samples);

// Get final MP3 data
const finalData = encoder.flush();

// Create MP3 blob (browser)
const blob = new Blob([mp3Data, finalData], { type: 'audio/mp3' });
```

## Stereo Encoding

```typescript
import { Mp3Encoder, MPEGMode } from 'lamejs2';

// Create stereo encoder
const encoder = new Mp3Encoder(2, 44100, 128);

// Encode stereo data
const left = new Int16Array(44100);  // left channel data
const right = new Int16Array(44100); // right channel data
const mp3Data = encoder.encodeBuffer(left, right);

// Finish encoding
const finalData = encoder.flush();
```

## API Reference

The encoder constructor takes three parameters:

```typescript
new Mp3Encoder(
  channels: number,   // 1 for mono, 2 for stereo
  sampleRate: number, // e.g. 44100
  kbps: number       // e.g. 128
)
```

If no parameters are provided, it defaults to:
- channels: 1 (mono)
- sampleRate: 44100 Hz
- kbps: 128

The encoder automatically configures optimal settings for:
- mode (STEREO/MONO based on channels)
- quality (3)
- VBR tag writing (disabled)
- reservoir (disabled)
- ID3 tag writing (disabled)

## Performance

Like the original lamejs, this encoder is very fast:
- Encodes ~20x faster than realtime
- 132 second sample encodes in ~6.5 seconds
- Works efficiently in both Node.js and browsers

## Credits

This is a TypeScript modernization of [lamejs](https://github.com/zhuker/lamejs) by Alex Zhukov, which was a JavaScript port of LAME MP3 encoder.

## License

LGPL-3.0 (same as original lamejs)
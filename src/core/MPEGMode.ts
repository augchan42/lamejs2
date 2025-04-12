// Define the enum
export const enum MPEGMode {
  STEREO = 0,
  JOINT_STEREO = 1,
  DUAL_CHANNEL = 2,
  MONO = 3,
  NOT_SET = 4
}

// Define the type for mode values
type ModeValue = {
  value: number;
  ordinal: () => number;
};

// Create a mapping from enum to mode values
export const MPEGModeValues: Record<keyof typeof MPEGMode, ModeValue> = {
  STEREO: { value: 0, ordinal: () => 0 },
  JOINT_STEREO: { value: 1, ordinal: () => 1 },
  DUAL_CHANNEL: { value: 2, ordinal: () => 2 },
  MONO: { value: 3, ordinal: () => 3 },
  NOT_SET: { value: 4, ordinal: () => 4 }
} as const; 
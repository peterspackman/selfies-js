/**
 * Custom exception classes for SELFIES library
 */

export class EncoderError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'EncoderError';
    Object.setPrototypeOf(this, EncoderError.prototype);
  }
}

export class DecoderError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'DecoderError';
    Object.setPrototypeOf(this, DecoderError.prototype);
  }
}
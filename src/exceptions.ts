/**
 * Custom exception classes for SELFIES library
 */

export class SMILESParserError extends Error {
  public smiles: string;
  public reason: string;
  public index: number;

  constructor(smiles: string, reason: string = "N/A", index: number = -1) {
    const pointer = index >= 0 ? ' '.repeat(index) + '^' : '';
    const message = `SMILES: ${smiles}${pointer ? '\n        ' + pointer : ''}\nIndex:  ${index}\nReason: ${reason}`;
    
    super(message);
    this.name = 'SMILESParserError';
    this.smiles = smiles;
    this.reason = reason;
    this.index = index;
    Object.setPrototypeOf(this, SMILESParserError.prototype);
  }
}

export class EncoderError extends Error {
  public cause?: Error;

  constructor(message: string, cause?: Error) {
    super(message);
    this.name = 'EncoderError';
    this.cause = cause;
    Object.setPrototypeOf(this, EncoderError.prototype);
  }
}

export class DecoderError extends Error {
  public cause?: Error;

  constructor(message: string, cause?: Error) {
    super(message);
    this.name = 'DecoderError';
    this.cause = cause;
    Object.setPrototypeOf(this, DecoderError.prototype);
  }
}
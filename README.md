# WORLD#

WORLD# is a pure C# porting of [WORLD](https://github.com/mmorise/World)

## Usage

```csharp
var framePeriod = 5.0;
var sampleRate = 44100;
var audio = new double[441000]; // monaural audio data

// f0 estimate
var harvest = new Harvest();
harvest.UseMultiThread = true; // optional
var f0Length = harvest.GetSamplesForHarvest(sampleRate, audio.Length, framePeriod);
var f0 = new double[f0Length];
var timeAxis = new double[f0Length];
harvest.Estimate(audio, sampleRate, timeAxis, f0);

// spectral envelope estimate
var cheapTrick = new CheapTrick(sampleRate);
var fftSize = cheapTrick.GetFFTSizeForCheapTrick(sampleRate);;
cheapTrick.FFTSize = fftSize;
cheapTrick.UseMultiThread = true; // optional
var spectrogram = Enumerable.Range(0, f0Length).Select(_ => new double[fftSize / 2 + 1]).ToArray();
cheapTrick.Estimate(audio, sampleRate, timeAxis, f0, spectrogram);

// aperiodicity estimate
var d4c = new D4C(sampleRate);
d4c.UseMultiThread = true; // optional
var aperiodicity = Enumerable.Range(0, f0Length).Select((i) => new double[fftSize / 2 + 1]).ToArray();
d4c.Estimate(audio, timeAxis, f0, f0Length, fftSize, aperiodicity);

// synthesize
var yLength = (int)((f0.Length - 1) * framePeriod / 1000.0 * sampleRate) + 1;
var synthesis = new Synthesis(sampleRate);
var result = new double[yLength];
synthesis.Synthesize(f0, f0.Length, spectrogram, aperiodicity, fftSize, framePeriod, result);
```

## Caveats

* tools are not ported

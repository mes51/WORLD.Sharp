using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    static class Common
    {
        public static int GetSuitableFFTSize(int sample)
        {
            return (int)(Math.Pow(2.0, (int)(MathUtil.Log2(sample) + 1.0)));
        }

        public static void NuttallWindow(int yLength, double[] y)
        {
            for (int i = 0; i < yLength; i++)
            {
                var tmp = i / (yLength - 1.0);
                y[i] = 0.355768 - 0.487396 * Math.Cos(2.0 * Math.PI * tmp) + 0.144232 * Math.Cos(4.0 * Math.PI * tmp) - 0.012604 * Math.Cos(6.0 * Math.PI * tmp);
            }
        }

        public static void DCCorrection(double[] input, double f0, int fs, int fftSize, double[] output)
        {
            var upperLimit = 2 + (int)(f0 * fftSize / fs);
            var lowFrequencyReplica = new double[upperLimit];
            var lowFrequencyAxis = Enumerable.Range(0, upperLimit).Select((i) => (double)i * fs / fftSize).ToArray();

            var upperLimitReplica = upperLimit - 1;
            MatlabFunctions.Interp1Q(f0 - lowFrequencyAxis[0], -(double)fs / fftSize, input.AsSpan(0, upperLimit + 1), lowFrequencyAxis.AsSpan(0, upperLimitReplica), lowFrequencyReplica);

            for (var i = 0; i < upperLimitReplica; i++)
            {
                output[i] = input[i] + lowFrequencyReplica[i];
            }
        }

        public static void LinearSmoothing(double[] input, double width, int fs, int fftSize, double[] output)
        {
            var boundary = (int)(width * fftSize / fs) + 1;

            // These parameters are set by the other function.
            var mirroringSpectrum = new double[fftSize / 2 + boundary * 2 + 1];
            var mirroringSegment = new double[fftSize / 2 + boundary * 2 + 1];
            var frequencyAxis = new double[fftSize / 2 + 1];
            SetParametersForLinearSmoothing(boundary, fftSize, fs, width, input, mirroringSpectrum, mirroringSegment, frequencyAxis);

            var lowLevels = new double[fftSize / 2 + 1];
            var highLevels = new double[fftSize / 2 + 1];
            var originOfMirroringAxis = -(boundary - 0.5) * fs / fftSize;
            var discreteFrequencyInterval = (double)fs / fftSize;

            MatlabFunctions.Interp1Q(originOfMirroringAxis, discreteFrequencyInterval, mirroringSegment.AsSpan(0, fftSize / 2 + boundary * 2 + 1), frequencyAxis.AsSpan(0, fftSize / 2 + 1), lowLevels);

            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                frequencyAxis[i] += width;
            }

            MatlabFunctions.Interp1Q(originOfMirroringAxis, discreteFrequencyInterval, mirroringSegment.AsSpan(0, fftSize / 2 + boundary * 2 + 1), frequencyAxis.AsSpan(0, fftSize / 2 + 1), highLevels);

            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                output[i] = (highLevels[i] - lowLevels[i]) / width;
            }
        }

        static void SetParametersForLinearSmoothing(int boundary, int fftSize, int fs, double width, double[] powerSpectrum, double[] mirroringSpectrum, double[] mirroringSegment, double[] frequencyAxis)
        {
            for (var i = 0; i < boundary; i++)
            {
                mirroringSpectrum[i] = powerSpectrum[boundary - i];
            }
            for (int i = boundary, limit = fftSize / 2 + boundary; i < limit; i++)
            {
                mirroringSpectrum[i] = powerSpectrum[i - boundary];
            }
            for (int i = fftSize / 2 + boundary, limit = fftSize / 2 + boundary * 2; i <= limit; i++)
            {
                mirroringSpectrum[i] = powerSpectrum[fftSize / 2 - (i - (fftSize / 2 + boundary))];
            }

            mirroringSegment[0] = mirroringSpectrum[0] * fs / fftSize;
            for (int i = 1, limit = fftSize / 2 + boundary * 2 + 1; i < limit; i++)
            {
                mirroringSegment[i] = mirroringSpectrum[i] * fs / fftSize + mirroringSegment[i - 1];
            }

            for (int i = 0, limit = fftSize / 2; i <= limit; ++i)
            {
                frequencyAxis[i] = (double)i / fftSize * fs - width / 2.0;
            }
        }

        public static double GetSafeAperiodicity(double x)
        {
            return Math.Max(0.001, Math.Min(0.999999999999, x));
        }

        internal static void GetMinimumPhaseSpectrum(MinimumPhaseAnalysis minimumPhase)
        {
            var fftSize = minimumPhase.FFTSize;
            for (var i = fftSize / 2 + 1; i < fftSize; i++)
            {
                minimumPhase.LogSpectrum[i] = minimumPhase.LogSpectrum[fftSize - i];
            }

            FFT.Execute(minimumPhase.InverseFFT);
            var cepstrum = minimumPhase.Cepstrum;
            cepstrum[0] = new Complex(cepstrum[0].Real, cepstrum[0].Imaginary);
            for (int i = 1, limit = fftSize / 2; i < limit; i++)
            {
                cepstrum[i] = new Complex(cepstrum[i].Real * 2.0, cepstrum[i].Imaginary * -2.0);
            }
            cepstrum[fftSize / 2] = new Complex(cepstrum[fftSize / 2].Real, cepstrum[fftSize / 2].Imaginary * -1.0);
            for (var i = fftSize / 2 + 1; i < fftSize; i++)
            {
                cepstrum[i] = new Complex();
            }

            FFT.Execute(minimumPhase.ForwardFFT);

            // Since x is complex number, calculation of exp(x) is as following.
            // Note: This FFT library does not keep the aliasing.
            var mininimumPhaseSpectrum = minimumPhase.MinimumPhaseSpectrum;
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                var tmp = Math.Exp(mininimumPhaseSpectrum[i].Real / fftSize);
                var real = tmp * Math.Cos(mininimumPhaseSpectrum[i].Imaginary / fftSize);
                var imaginary = tmp * Math.Sin(mininimumPhaseSpectrum[i].Imaginary / fftSize);
                mininimumPhaseSpectrum[i] = new Complex(real, imaginary);
            }
        }
    }
}

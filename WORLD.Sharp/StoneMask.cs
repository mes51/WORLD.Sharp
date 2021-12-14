using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    public static class StoneMask
    {
        const double FloorF0 = 40.0;

        const double SafeGuardMinimum = 0.000000000001;

        public static void Refine(double[] x, int fs, double[] temporalPositions, double[] f0, double[] refinedF0, bool useMultiThread = false)
        {
            if (useMultiThread)
            {
                Parallel.For(0, f0.Length, i =>
                {
                    refinedF0[i] = GetRefinedF0(x, fs, temporalPositions[i], f0[i]);
                });
            }
            else
            {
                for (var i = 0; i < f0.Length; i++)
                {
                    refinedF0[i] = GetRefinedF0(x, fs, temporalPositions[i], f0[i]);
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetRefinedF0() fixes the F0 estimated by Dio(). This function uses
        // instantaneous frequency.
        //-----------------------------------------------------------------------------
        static double GetRefinedF0(double[] x, int fs, double currentPotision, double initialF0)
        {
            if (initialF0 <= FloorF0 || initialF0 > fs / 12.0)
            {
                return 0.0;
            }

            var halfWindowLength = (int)(1.5 * fs / initialF0 + 1.0);
            var windowLengthInTime = (2.0 * halfWindowLength + 1.0) / fs;
            var baseTime = new double[halfWindowLength * 2 + 1];
            for (var i = 0; i < baseTime.Length; i++)
            {
                baseTime[i] = (-halfWindowLength + i) / (double)fs;
            }

            var fftSize = (int)(Math.Pow(2.0, 2.0 + (int)MathUtil.Log2(halfWindowLength * 2.0 + 1.0)));

            var meanF0 = GetMeanF0(x, fs, currentPotision, initialF0, fftSize, windowLengthInTime, baseTime);

            if (Math.Abs(meanF0 - initialF0) > initialF0 * 0.2)
            {
                meanF0 = initialF0;
            }

            return meanF0;
        }

        //-----------------------------------------------------------------------------
        // GetMeanF0() calculates the instantaneous frequency.
        //-----------------------------------------------------------------------------
        static double GetMeanF0(double[] x, int fs, double currentPotision, double initialF0, int fftSize, double windowLengthInTime, double[] baseTime)
        {
            var forwardRealFFT = ForwardRealFFT.Create(fftSize);
            var mainSpectrum = new Complex[fftSize];
            var diffSpectrum = new Complex[fftSize];

            var indexRaw = new int[baseTime.Length];
            var mainWindow = new double[baseTime.Length];
            var diffWindow = new double[baseTime.Length];

            GetBaseIndex(currentPotision, baseTime, fs, indexRaw);
            GetMainWindow(currentPotision, indexRaw, fs, windowLengthInTime, mainWindow);
            GetDiffWindow(mainWindow, diffWindow);
            GetSpectra(x, fftSize, indexRaw, mainWindow, diffWindow, forwardRealFFT, mainSpectrum, diffSpectrum);

            var powerSpectrum = new double[fftSize / 2 + 1];
            var numeratorI = new double[fftSize / 2 + 1];
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                numeratorI[i] = mainSpectrum[i].Real * diffSpectrum[i].Imaginary - mainSpectrum[i].Imaginary * diffSpectrum[i].Real;
                powerSpectrum[i] = mainSpectrum[i].Real * mainSpectrum[i].Real + mainSpectrum[i].Imaginary * mainSpectrum[i].Imaginary;
            }

            forwardRealFFT.Release();

            return GetTentativeF0(powerSpectrum, numeratorI, fftSize, fs, initialF0);
        }

        //-----------------------------------------------------------------------------
        // GetBaseIndex() calculates the temporal positions for windowing.
        // Since the result includes negative value and the value that exceeds the
        // length of the input signal, it must be modified appropriately.
        //-----------------------------------------------------------------------------
        static void GetBaseIndex(double currentPotision, double[] baseTime, int fs, int[] indexRaw)
        {
            for (var i = 0; i < baseTime.Length; i++)
            {
                indexRaw[i] = MatlabFunctions.MatlabRound((currentPotision + baseTime[i]) * fs);
            }
        }

        //-----------------------------------------------------------------------------
        // GetMainWindow() generates the window function.
        //-----------------------------------------------------------------------------
        static void GetMainWindow(double currentPotision, int[] indexRaw, int fs, double windowLengthInTime, double[] mainWindow)
        {
            for (var i = 0; i < mainWindow.Length; i++)
            {
                var tmp = (indexRaw[i] - 1.0) / fs - currentPotision;
                mainWindow[i] = 0.42 + 0.5 * Math.Cos(2.0 * Math.PI * tmp / windowLengthInTime) + 0.08 * Math.Cos(4.0 * Math.PI * tmp / windowLengthInTime);
            }
        }

        //-----------------------------------------------------------------------------
        // GetDiffWindow() generates the differentiated window.
        // Diff means differential.
        //-----------------------------------------------------------------------------
        static void GetDiffWindow(double[] mainWindow, double[] diffWindow)
        {
            diffWindow[0] = -mainWindow[1] / 2.0;
            for (int i = 1, limit = mainWindow.Length - 1; i < limit; i++)
            {
                diffWindow[i] = -(mainWindow[i + 1] - mainWindow[i - 1]) / 2.0;
            }
            diffWindow[mainWindow.Length - 1] = mainWindow[mainWindow.Length - 2] / 2.0;
        }

        //-----------------------------------------------------------------------------
        // GetSpectra() calculates two spectra of the waveform windowed by windows
        // (main window and diff window).
        //-----------------------------------------------------------------------------
        static void GetSpectra(double[] x, int fftSize, int[] indexRaw, double[] mainWindow, double[] diffWindow, ForwardRealFFT forwardRealFFT, Complex[] mainSpectrum, Complex[] diffSpectrum)
        {
            var index = new int[indexRaw.Length];

            for (var i = 0; i < index.Length; i++)
            {
                index[i] = Math.Max(0, Math.Min(x.Length - 1, indexRaw[i] - 1));
            }
            for (var i = 0; i < index.Length; i++)
            {
                forwardRealFFT.Waveform[i] = x[index[i]] * mainWindow[i];
            }
            forwardRealFFT.Waveform.AsSpan(index.Length, fftSize - index.Length).Clear();

            FFT.Execute(forwardRealFFT.ForwardFFT);
            forwardRealFFT.Spectrum.AsSpan(0, fftSize / 2 + 1).CopyTo(mainSpectrum);

            for (int i = 0; i < index.Length; ++i)
            {
                forwardRealFFT.Waveform[i] = x[index[i]] * diffWindow[i];
            }
            forwardRealFFT.Waveform.AsSpan(index.Length, fftSize - index.Length).Clear();

            FFT.Execute(forwardRealFFT.ForwardFFT);
            forwardRealFFT.Spectrum.AsSpan(0, fftSize / 2 + 1).CopyTo(diffSpectrum);    
        }

        //-----------------------------------------------------------------------------
        // GetTentativeF0() calculates the F0 based on the instantaneous frequency.
        //-----------------------------------------------------------------------------
        static double GetTentativeF0(double[] powerSpectrum, double[] numeratorI, int fftSize, int fs, double initialF0)
        {
            var tentativeF0 = FixF0(powerSpectrum, numeratorI, fftSize, fs, initialF0, 2);

            // If the fixed value is too large, the result will be rejected.
            if (tentativeF0 <= 0.0 || tentativeF0 > initialF0 * 2.0)
            {
                return 0.0;
            }

            return FixF0(powerSpectrum, numeratorI, fftSize, fs, tentativeF0, 6);
        }

        //-----------------------------------------------------------------------------
        // FixF0() fixed the F0 by instantaneous frequency.
        //-----------------------------------------------------------------------------
        static double FixF0(double[] powerSpectrum, double[] numeratorI, int fftSize, int fs, double initialF0, int numberOfHarmonics)
        {
            var amplitudeList = new double[numberOfHarmonics];
            var instantaneousFrequencyList = new double[numberOfHarmonics];
            for (var i = 0; i < numberOfHarmonics; i++)
            {
                var index = Math.Min(MatlabFunctions.MatlabRound(initialF0 * fftSize / fs * (i + 1)), fftSize / 2);
                instantaneousFrequencyList[i] = powerSpectrum[index] == 0.0 ? 0.0 : (double)index * fs / fftSize + numeratorI[index] / powerSpectrum[index] * fs / 2.0 / Math.PI;
                amplitudeList[i] = Math.Sqrt(powerSpectrum[index]);
            }

            var denominator = 0.0;
            var numerator = 0.0;
            for (var i = 0; i < numberOfHarmonics; i++)
            {
                numerator += amplitudeList[i] * instantaneousFrequencyList[i];
                denominator += amplitudeList[i] * (i + 1);
            }

            return numerator / (denominator + SafeGuardMinimum);
        }
    }
}

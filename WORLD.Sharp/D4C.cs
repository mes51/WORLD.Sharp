using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    enum WindowType
    {
        Hanning,
        Blackman
    }

    public class D4C
    {
        const double DefaultThreshold = 0.85;
        const double FloorF0D4C = 47.0;
        const double SafeGuardMinimum = 0.000000000001;
        const double UpperLimit = 15000.0;
        const double FrequencyInterval = 3000.0;

        public double Threshold { get; set; }

        public bool UseMultiThread { get; set; }

        int Fs { get; }

        public D4C(int fs)
        {
            Threshold = DefaultThreshold;
            Fs = fs;
        }

        public void Estimate(double[] x, double[] temporalPositions, double[] f0, int f0Length, int fftSize, double[][] aperiodicity)
        {
            aperiodicity.ForEach((a) => a.Fill(1.0 - SafeGuardMinimum));

            var fftSizeD4C = (int)Math.Pow(2.0, 1.0 + (int)MathUtil.Log2(4.0 * Fs / FloorF0D4C + 1.0));

            var numberOfAperiodicities = (int)(Math.Min(UpperLimit, Fs / 2.0 - FrequencyInterval) / FrequencyInterval);

            // Since the window function is common in D4CGeneralBody(),
            // it is designed here to speed up.
            var windowLength = (int)(FrequencyInterval * fftSizeD4C / Fs) * 2 + 1;
            var window = new double[windowLength];
            Common.NuttallWindow(windowLength, window);

            // D4C Love Train (Aperiodicity of 0 Hz is given by the different algorithm)
            var aperiodicity0 = new double[f0Length];
            var rand = new MVN(Fs);
            D4CLoveTrain(x, f0, f0Length, temporalPositions, rand, aperiodicity0);

            var coarseFrequencyAxis = Enumerable.Range(0, numberOfAperiodicities + 2).Select((i) => i * FrequencyInterval).ToArray();
            coarseFrequencyAxis[numberOfAperiodicities + 1] = Fs / 2.0;

            var frequencyAxis = Enumerable.Range(0, fftSize / 2 + 1).Select((i) => (double)i * Fs / fftSize).ToArray();

            if (UseMultiThread)
            {
                Parallel.For(
                    0,
                    f0.Length,
                    () => new { ForwardRealFFT = ForwardRealFFT.Create(fftSizeD4C), CoarseAperiodicity = new double[numberOfAperiodicities + 2], Rand = new MVN(Fs) },
                    (i, loop, t) =>
                    {
                        t.CoarseAperiodicity[0] = -60.0;
                        t.CoarseAperiodicity[numberOfAperiodicities + 1] = -SafeGuardMinimum;

                        if (f0[i] == 0.0 || aperiodicity0[i] <= Threshold)
                        {
                            return t;
                        }

                        D4CGeneralBody(x, Math.Max(FloorF0D4C, f0[i]), fftSizeD4C, temporalPositions[i], numberOfAperiodicities, window, t.Rand, t.ForwardRealFFT, t.CoarseAperiodicity.AsSpan(1));

                        // Linear interpolation to convert the coarse aperiodicity into its
                        // spectral representation.
                        GetAperiodicity(coarseFrequencyAxis, t.CoarseAperiodicity, numberOfAperiodicities, frequencyAxis, fftSize, aperiodicity[i]);

                        return t;
                    },
                    t =>
                    {
                        t.ForwardRealFFT.Release();
                    }
                );
            }
            else
            {
                var coarseAperiodicity = new double[numberOfAperiodicities + 2];
                coarseAperiodicity[0] = -60.0;
                coarseAperiodicity[numberOfAperiodicities + 1] = -SafeGuardMinimum;
                var forwardRealFFT = ForwardRealFFT.Create(fftSizeD4C);

                for (var i = 0; i < f0Length; i++)
                {
                    if (f0[i] == 0.0 || aperiodicity0[i] <= Threshold)
                    {
                        continue;
                    }

                    D4CGeneralBody(x, Math.Max(FloorF0D4C, f0[i]), fftSizeD4C, temporalPositions[i], numberOfAperiodicities, window, rand, forwardRealFFT, coarseAperiodicity.AsSpan(1));

                    // Linear interpolation to convert the coarse aperiodicity into its
                    // spectral representation.
                    GetAperiodicity(coarseFrequencyAxis, coarseAperiodicity, numberOfAperiodicities, frequencyAxis, fftSize, aperiodicity[i]);
                }

                forwardRealFFT.Release();
            }
        }

        void GetAperiodicity(double[] coarseFrequencyAxis, double[] coarseAperiodicity, int numberOfAperiodicities, double[] frequencyAxis, int fftSize, double[] aperiodicity)
        {
            MatlabFunctions.Interp1(coarseFrequencyAxis.AsSpan(0, numberOfAperiodicities + 2), coarseAperiodicity, frequencyAxis.AsSpan(0, fftSize / 2 + 1), aperiodicity);
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                aperiodicity[i] = Math.Pow(10.0, aperiodicity[i] / 20.0);
            }
        }

        //-----------------------------------------------------------------------------
        // D4CGeneralBody() calculates a spectral envelope at a temporal
        // position. This function is only used in D4C().
        // Caution:
        //   forward_fft is allocated in advance to speed up the processing.
        //-----------------------------------------------------------------------------
        void D4CGeneralBody(double[] x, double currentF0, int fftSize, double currentPosition, int numberOfAperiodicities, double[] window, MVN rand, ForwardRealFFT forwardRealFFT, Span<double> coarseAperiodicity)
        {
            var halfFFTSize = fftSize / 2 + 1;
            var pool = ArrayPool<double>.Shared;

            var staticCentroid = pool.Rent(halfFFTSize);
            var smoothedPowerSpectrum = pool.Rent(halfFFTSize);
            var staticGroupDelay = pool.Rent(halfFFTSize);

            GetStaticCentroid(x, currentF0, fftSize, currentPosition, rand, forwardRealFFT, staticCentroid);
            GetSmoothedPowerSpectrum(x, currentF0, fftSize, currentPosition, rand, forwardRealFFT, smoothedPowerSpectrum);
            GetStaticGroupDelay(staticCentroid, smoothedPowerSpectrum, currentF0, fftSize, staticGroupDelay);

            GetCoarseAperiodicity(staticGroupDelay, fftSize, numberOfAperiodicities, window, forwardRealFFT, coarseAperiodicity);

            // Revision of the result based on the F0
            for (var i = 0; i < numberOfAperiodicities; i++)
            {
                coarseAperiodicity[i] = Math.Min(0.0, coarseAperiodicity[i] + (currentF0 - 100) / 50.0);
            }

            pool.Return(staticCentroid);
            pool.Return(smoothedPowerSpectrum);
            pool.Return(staticGroupDelay);
        }

        //-----------------------------------------------------------------------------
        // GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
        // The upper limit is given based on the sampling frequency.
        //-----------------------------------------------------------------------------
        void GetCoarseAperiodicity(double[] staticGroupDelay, int fftSize, int numberOfAperiodicities, double[] window, ForwardRealFFT forwardRealFFT, Span<double> coarseAperiodicity)
        {
            var boundary = MatlabFunctions.MatlabRound(fftSize * 8.0 / window.Length);
            var halfWindowLength = window.Length / 2;

            Array.Clear(forwardRealFFT.Waveform, 0, fftSize);

            var halfFFTSize = fftSize / 2 + 1;
            var powerSpectrum = ArrayPool<double>.Shared.Rent(halfFFTSize);
            for (var i = 0; i < numberOfAperiodicities; i++)
            {
                var center = (int)(FrequencyInterval * (i + 1) * fftSize / Fs);
                for (int j = 0, limit = halfWindowLength * 2; j <= limit; j++)
                {
                    forwardRealFFT.Waveform[j] = staticGroupDelay[center - halfWindowLength + j] * window[j];
                }
                FFT.Execute(forwardRealFFT.ForwardFFT);

                var spectrum = forwardRealFFT.Spectrum;
                for (var j = 0; j < halfFFTSize; j++)
                {
                    powerSpectrum[j] = spectrum[j].Real * spectrum[j].Real + spectrum[j].Imaginary * spectrum[j].Imaginary;
                }

                coarseAperiodicity[i] = 10.0 * Math.Log10(CalcPowerSpectrumRate(powerSpectrum, fftSize, boundary));
            }

            ArrayPool<double>.Shared.Return(powerSpectrum);
        }

        double CalcPowerSpectrumRate(double[] powerSpectrum, int fftSize, int boundary)
        {
            var halfFFTSize = fftSize / 2 + 1;
            var boundaryPos = fftSize / 2 - boundary - 1;

            var sum = 0.0;
            for (var i = 0; i < halfFFTSize; i++)
            {
                sum += powerSpectrum[i];
            }

            LargeOnlyGrouping.Grouping(powerSpectrum, halfFFTSize, boundaryPos + 1);

            var boundarySum = 0.0;
            for (var i = 0; i <= boundaryPos; i++)
            {
                boundarySum += powerSpectrum[i];
            }

            return boundarySum / sum;
        }

        //-----------------------------------------------------------------------------
        // GetStaticGroupDelay() calculates the temporally static group delay.
        // This is the fundamental parameter in D4C.
        //-----------------------------------------------------------------------------
        void GetStaticGroupDelay(double[] staticCentroid, double[] smoothedPowerSpectrum, double f0, int fftSize, double[] staticGroupDelay)
        {
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                staticGroupDelay[i] = staticCentroid[i] / smoothedPowerSpectrum[i];
            }
            Common.LinearSmoothing(staticGroupDelay, f0 / 2.0, Fs, fftSize, staticGroupDelay);

            var halfFFTSize = fftSize / 2 + 1;
            var smoothedGroupDelay = ArrayPool<double>.Shared.Rent(halfFFTSize);
            Common.LinearSmoothing(staticGroupDelay, f0, Fs, fftSize, smoothedGroupDelay);

            for (var i = 0; i < halfFFTSize; i++)
            {
                staticGroupDelay[i] -= smoothedGroupDelay[i];
            }

            ArrayPool<double>.Shared.Return(smoothedGroupDelay);
        }

        //-----------------------------------------------------------------------------
        // GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
        // The parameters used for smoothing are optimized in davance.
        //-----------------------------------------------------------------------------
        void GetSmoothedPowerSpectrum(double[] x, double currentF0, int fftSize, double currentPosition, MVN rand, ForwardRealFFT forwardRealFFT, double[] smoothedPowerSpectrum)
        {
            Array.Clear(forwardRealFFT.Waveform, 0, fftSize);
            GetWindowedWaveform(x, currentF0, currentPosition, WindowType.Hanning, 4.0, rand, forwardRealFFT.Waveform);

            FFT.Execute(forwardRealFFT.ForwardFFT);
            var spectrum = forwardRealFFT.Spectrum;
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                smoothedPowerSpectrum[i] = spectrum[i].Real * spectrum[i].Real + spectrum[i].Imaginary * spectrum[i].Imaginary;
            }
            Common.DCCorrection(smoothedPowerSpectrum, currentF0, Fs, fftSize, smoothedPowerSpectrum);
            Common.LinearSmoothing(smoothedPowerSpectrum, currentF0, Fs, fftSize, smoothedPowerSpectrum);
        }

        //-----------------------------------------------------------------------------
        // GetStaticCentroid() calculates the temporally static energy centroid.
        // Basic idea was proposed by H. Kawahara.
        //-----------------------------------------------------------------------------
        void GetStaticCentroid(double[] x, double currentF0, int fftSize, double currentPosition, MVN rand, ForwardRealFFT forwardRealFFT, double[] staticCentroid)
        {
            var halfFFTSize = fftSize / 2 + 1;
            var pool = ArrayPool<double>.Shared;

            var centroid1 = pool.Rent(halfFFTSize);
            var centroid2 = pool.Rent(halfFFTSize);

            GetCentroid(x, currentF0, fftSize, currentPosition - 0.25 / currentF0, rand, forwardRealFFT, centroid1);
            GetCentroid(x, currentF0, fftSize, currentPosition + 0.25 / currentF0, rand, forwardRealFFT, centroid2);

            for (var i = 0; i < halfFFTSize; i++)
            {
                staticCentroid[i] = centroid1[i] + centroid2[i];
            }

            Common.DCCorrection(staticCentroid, currentF0, Fs, fftSize, staticCentroid);

            pool.Return(centroid1);
            pool.Return(centroid2);
        }

        void GetCentroid(double[] x, double currentF0, int fftSize, double currentPosition, MVN rand, ForwardRealFFT forwardRealFFT, double[] centroid)
        {
            forwardRealFFT.Waveform.AsSpan().Clear();
            GetWindowedWaveform(x, currentF0, currentPosition, WindowType.Blackman, 4.0, rand, forwardRealFFT.Waveform);

            var windowRange = MatlabFunctions.MatlabRound(2.0 * Fs / currentF0) * 2;
            var waveform = forwardRealFFT.Waveform;
            var power = 0.0;
            for (var i = 0; i < windowRange; i++)
            {
                power += waveform[i] * waveform[i];
            }
            power = Math.Sqrt(power);
            for (var i = 0; i <= windowRange; i++)
            {
                forwardRealFFT.Waveform[i] /= power;
            }

            FFT.Execute(forwardRealFFT.ForwardFFT);
            var halfFFTSize = fftSize / 2 + 1;
            var tmpSpectrum = ArrayPool<Complex>.Shared.Rent(halfFFTSize);
            forwardRealFFT.Spectrum.AsSpan(0, halfFFTSize).CopyTo(tmpSpectrum);

            for (var i = 0; i < fftSize; i++)
            {
                forwardRealFFT.Waveform[i] *= i + 1.0;
            }
            FFT.Execute(forwardRealFFT.ForwardFFT);
            var spectrum = forwardRealFFT.Spectrum;
            for (var i = 0; i < halfFFTSize; i++)
            {
                centroid[i] = spectrum[i].Real * tmpSpectrum[i].Real + spectrum[i].Imaginary * tmpSpectrum[i].Imaginary;
            }

            ArrayPool<Complex>.Shared.Return(tmpSpectrum);
        }

        //-----------------------------------------------------------------------------
        // GetWindowedWaveform() windows the waveform by F0-adaptive window
        // In the variable window_type, 1: hanning, 2: blackman
        //-----------------------------------------------------------------------------
        void GetWindowedWaveform(double[] x, double currentF0, double currentPosition, WindowType windowType, double windowLengthRatio, MVN rand, double[] waveform)
        {
            var halfWindowLength = MatlabFunctions.MatlabRound(windowLengthRatio * Fs / currentF0 / 2.0);
            var windowLength = halfWindowLength * 2 + 1;
            var intPool = ArrayPool<int>.Shared;

            var baseIndex = intPool.Rent(windowLength);
            var safeIndex = intPool.Rent(windowLength);
            var window = ArrayPool<double>.Shared.Rent(windowLength);

            SetParametersForGetWindowedWaveform(halfWindowLength, x.Length, currentPosition, currentF0, windowType, windowLengthRatio, baseIndex, safeIndex, window);

            for (var i = 0; i < windowLength; i++)
            {
                waveform[i] = x[safeIndex[i]] * window[i] + rand.GenerateMVN() * SafeGuardMinimum;
            }

            var tmpWeight1 = 0.0;
            var tmpWeight2 = 0.0;
            for (var i = 0; i < waveform.Length; i++)
            {
                tmpWeight1 += waveform[i];
            }
            for (var i = 0; i < window.Length; i++)
            {
                tmpWeight2 += window[i];
            }
            var weightingCoefficient = tmpWeight1 / tmpWeight2;
            for (int i = 0, limit = halfWindowLength * 2; i <= limit; i++)
            {
                waveform[i] -= window[i] * weightingCoefficient;
            }

            intPool.Return(baseIndex);
            intPool.Return(safeIndex);
            ArrayPool<double>.Shared.Return(window);
        }

        //-----------------------------------------------------------------------------
        // SetParametersForGetWindowedWaveform()
        //-----------------------------------------------------------------------------
        void SetParametersForGetWindowedWaveform(int halfWindowLength, int xLength, double currentPosition, double currentF0, WindowType windowType, double windowLengthRatio, int[] baseIndex, int[] safeIndex, double[] window)
        {
            for (var i = -halfWindowLength; i <= halfWindowLength; i++)
            {
                baseIndex[i + halfWindowLength] = i;
            }
            var origin = MatlabFunctions.MatlabRound(currentPosition * Fs + 0.001);
            for (int i = 0, limit = halfWindowLength * 2; i <= limit; i++)
            {
                safeIndex[i] = Math.Min(xLength - 1, Math.Max(0, origin + baseIndex[i]));
            }

            switch (windowType)
            {
                case WindowType.Hanning:
                    for (int i = 0, limit = halfWindowLength * 2; i <= limit; i++)
                    {
                        var position = (2.0 * baseIndex[i] / windowLengthRatio) / Fs;
                        window[i] = 0.5 * Math.Cos(Math.PI * position * currentF0) + 0.5;
                    }
                    break;
                case WindowType.Blackman:
                    for (int i = 0, limit = halfWindowLength * 2; i <= limit; i++)
                    {
                        var position = (2.0 * baseIndex[i] / windowLengthRatio) / Fs;
                        window[i] = 0.42 + 0.5 * Math.Cos(Math.PI * position * currentF0) + 0.08 * Math.Cos(Math.PI * position * currentF0 * 2);
                    }
                    break;
            }
        }

        //-----------------------------------------------------------------------------
        // D4CLoveTrain() determines the aperiodicity with VUV detection.
        // If a frame was determined as the unvoiced section, aperiodicity is set to
        // very high value as the safeguard.
        // If it was voiced section, the aperiodicity of 0 Hz is set to -60 dB.
        //-----------------------------------------------------------------------------
        void D4CLoveTrain(double[] x, double[] f0, int f0Length, double[] temporalPositions, MVN rand, double[] aperiodicity0)
        {
            var lowestF0 = 40.0;
            var fftSize = (int)Math.Pow(2.0, 1.0 + (int)MathUtil.Log2(3.0 * Fs / lowestF0 + 1));

            // Cumulative powers at 100, 4000, 7900 Hz are used for VUV identification.
            var boundary0 = (int)Math.Ceiling(100.0 * fftSize / Fs);
            var boundary1 = (int)Math.Ceiling(4000.0 * fftSize / Fs);
            var boundary2 = (int)Math.Ceiling(7900.0 * fftSize / Fs);

            if (UseMultiThread)
            {
                Parallel.For(
                    0,
                    f0.Length,
                    () => new { ForwardRealFFT = ForwardRealFFT.Create(fftSize), PowerSpectrum = new double[fftSize] },
                    (i, loop, t) =>
                    {
                        if (f0[i] == 0.0)
                        {
                            aperiodicity0[i] = 0.0;
                        }
                        else
                        {
                            aperiodicity0[i] = D4CLoveTrainSub(x, Math.Max(f0[i], lowestF0), temporalPositions[i], f0Length, fftSize, boundary0, boundary1, boundary2, rand, t.ForwardRealFFT, t.PowerSpectrum);
                            t.PowerSpectrum.AsSpan().Clear();
                        }

                        return t;
                    },
                    t =>
                    {
                        t.ForwardRealFFT.Release();
                    }
                );
            }
            else
            {
                var forwardRealFFT = ForwardRealFFT.Create(fftSize);
                var powerSpectrum = new double[fftSize];
                for (var i = 0; i < f0Length; i++)
                {
                    if (f0[i] == 0.0)
                    {
                        aperiodicity0[i] = 0.0;
                        continue;
                    }

                    aperiodicity0[i] = D4CLoveTrainSub(x, Math.Max(f0[i], lowestF0), temporalPositions[i], f0Length, fftSize, boundary0, boundary1, boundary2, rand, forwardRealFFT, powerSpectrum);
                    powerSpectrum.AsSpan().Clear();
                }

                forwardRealFFT.Release();
            }
        }

        double D4CLoveTrainSub(double[] x, double currentF0, double currentPosition, int f0Length, int fftSize, int boundary0, int boundary1, int boundary2, MVN rand, ForwardRealFFT forwardRealFFT, double[] powerSpectrum)
        {
            var windowLength = MatlabFunctions.MatlabRound(1.5 * Fs / currentF0) * 2 + 1;
            GetWindowedWaveform(x, currentF0, currentPosition, WindowType.Blackman, 3.0, rand, forwardRealFFT.Waveform);

            forwardRealFFT.Waveform.AsSpan(windowLength, fftSize - windowLength).Clear();
            FFT.Execute(forwardRealFFT.ForwardFFT);

            var spectrum = forwardRealFFT.Spectrum;
            for (int i = boundary0 + 1, limit = fftSize / 2 + 1; i < limit; i++)
            {
                powerSpectrum[i] = spectrum[i].Real * spectrum[i].Real + spectrum[i].Imaginary * spectrum[i].Imaginary;
            }
            for (var i = boundary0; i <= boundary2; i++)
            {
                powerSpectrum[i] += powerSpectrum[i - 1];
            }

            return powerSpectrum[boundary1] / powerSpectrum[boundary2];
        }
    }

    static class LargeOnlyGrouping
    {
        const int SwitchThreshold = 64;

        public static void Grouping(double[] values, int end, int limit)
        {
            QuickSort(values, 0, end, limit);
        }

        static void QuickSort(double[] values, int start, int end, int limit)
        {
            if (end - start < SwitchThreshold)
            {
                Array.Sort(values, start, end - start + 1);
                return;
            }

            var pivot = Median(values[start], values[start + (end - start) / 2], values[end]);

            var left = start;
            var right = end;
            while (left <= right)
            {
                for (; left < end && pivot > values[left]; left++) ;
                for (; right > start && pivot <= values[right]; right--) ;

                if (left > right)
                {
                    break;
                }
                Swap(ref values[left], ref values[right]);
                left++;
                right--;
            }

            if (left > limit)
            {
                QuickSort(values, start, left - 1, limit);
            }
            else if (left != limit)
            {
                QuickSort(values, left, end, limit);
            }
        }

        static double Median(double a, double b, double c)
        {
            if (a > b)
            {
                Swap(ref a, ref b);
            }
            if (a > c)
            {
                Swap(ref a, ref c);
            }
            if (b > c)
            {
                Swap(ref b, ref c);
            }
            return b;
        }

        static void Swap(ref double a, ref double b)
        {
            var t = a;
            a = b;
            b = t;
        }
    }
}

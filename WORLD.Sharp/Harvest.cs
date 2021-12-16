using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    public class Harvest
    {
        const double DefaultF0Ceil = 800.0;
        const double DefaultF0Floor = 71.0;
        const double DefaultFramePeriod = 5.0;
        const double SafeGuardMinimum = 0.000000000001;

        public double F0Floor { get; set; }

        public double F0Ceil { get; set; }

        public double FramePeriod { get; set; }

        public bool UseMultiThread { get; set; }

        public Harvest()
        {
            F0Ceil = DefaultF0Ceil;
            F0Floor = DefaultF0Floor;
            FramePeriod = DefaultFramePeriod;
        }

        public void Estimate(double[] x, int fs, double[] temporalPositions, double[] f0)
        {
            const double targetFs = 8000.0;
            const double channelInOctave = 40.0;
            var dimensionRatio = MatlabFunctions.MatlabRound(fs / targetFs);

            if (FramePeriod == 1.0)
            {
                HarvestGeneralBody(x, fs, 1, channelInOctave, dimensionRatio, temporalPositions, f0);
            }
            else
            {
                var basicFramePeriod = 1;
                var basicF0Length = GetSamplesForHarvest(fs, x.Length, basicFramePeriod);
                var basicF0 = new double[basicF0Length];
                var basicTemporalPositions = new double[basicF0Length];
                HarvestGeneralBody(x, fs, basicFramePeriod, channelInOctave, dimensionRatio, basicTemporalPositions, basicF0);

                var f0Length = GetSamplesForHarvest(fs, x.Length, FramePeriod);
                var temporalPeriod = FramePeriod / 1000.0;
                for (var i = 0; i < f0Length; i++)
                {
                    temporalPositions[i] = i * temporalPeriod;
                    f0[i] = basicF0[Math.Min(basicF0Length - 1, MatlabFunctions.MatlabRound(temporalPositions[i] * 1000.0))];
                }
            }
        }

        public int GetSamplesForHarvest(int fs, int length, double framePeriod)
        {
            return (int)(1000.0 * length / fs / framePeriod) + 1;
        }

        //-----------------------------------------------------------------------------
        // HarvestGeneralBody() estimates the F0 contour based on Harvest.
        //-----------------------------------------------------------------------------
        void HarvestGeneralBody(double[] x, int fs, int framePeriod, double channelInOctave, int speed, double[] temporalPositions, double[] f0)
        {
            var adjustedF0Floor = F0Floor * 0.9;
            var adjustedF0Ceil = F0Ceil * 1.1;
            var numberOfChannels = 1 + (int)(MathUtil.Log2(adjustedF0Ceil / adjustedF0Floor) * channelInOctave);
            var boundaryF0List = Enumerable.Range(0, numberOfChannels)
                .Select((i) => adjustedF0Floor * Math.Pow(2.0, (i + 1) / channelInOctave))
                .ToArray();

            // normalization
            var decimationRatio = Math.Max(Math.Min(speed, 12), 1);
            var yLength = (int)Math.Ceiling((double)x.Length / decimationRatio);
            var actualFs = (double)fs / decimationRatio;
            var fftSize = Common.GetSuitableFFTSize(yLength + 5 + 2 * (int)(2.0 * actualFs / boundaryF0List[0]));

            // Calculation of the spectrum used for the f0 estimation
            var y = new double[fftSize];
            var ySpectrum = new Complex[fftSize];
            GetWaveformAndSpectrum(x, yLength, actualFs, fftSize, decimationRatio, y, ySpectrum);

            var f0Length = GetSamplesForHarvest(fs, x.Length, framePeriod);
            for (var i = 0; i < f0Length; i++)
            {
                temporalPositions[i] = i * framePeriod / 1000.0;
            }

            const int OverlapParameter = 7;
            var maxCandidates = MatlabFunctions.MatlabRound(numberOfChannels / 10.0) * OverlapParameter;
            var f0Candidates = Util.Mak2DArray<double>(f0Length, maxCandidates);
            var f0CandidatesScore = Util.Mak2DArray<double>(f0Length, maxCandidates);

            var numberOfCandidates = HarvestGeneralBodySub(boundaryF0List, numberOfChannels, f0Length, actualFs, yLength, temporalPositions, ySpectrum, fftSize, maxCandidates, f0Candidates) * OverlapParameter;

            RefineF0Candidates(y.AsMemory(0, yLength), actualFs, temporalPositions, f0Length, numberOfCandidates, f0Candidates, f0CandidatesScore);
            RemoveUnreliableCandidates(f0Length, numberOfCandidates, f0Candidates, f0CandidatesScore);

            var bestF0Counter = new double[f0Length];
            FixF0Contour(f0Candidates, f0CandidatesScore, f0Length, numberOfCandidates, bestF0Counter);
            SmoothF0Contour(bestF0Counter, f0Length, f0);
        }

        //-----------------------------------------------------------------------------
        // GetWaveformAndSpectrum() calculates the downsampled signal and its spectrum
        //-----------------------------------------------------------------------------
        void GetWaveformAndSpectrum(double[] x, int yLength, double actualFs, int fftSize, int decimationRatio, double[] y, Complex[] ySpectrum)
        {
            Array.Clear(y, 0, fftSize);

            // Processing for the compatibility with MATLAB version
            GetWaveformAndSpectrumSub(x, yLength, actualFs, decimationRatio, y);

            // Removal of the DC component (y = y - mean value of y)
            var meanY = y.Take(yLength).Average();
            for (var i = 0; i < yLength; i++)
            {
                y[i] -= meanY;
            }
            Array.Clear(y, yLength, fftSize - yLength);

            var forwardFFT = FFTPlan.CreatePlanForDftR2C(fftSize, y, ySpectrum, FFTFlag.Estimate);
            FFT.Execute(forwardFFT);
        }

        //-----------------------------------------------------------------------------
        // Since the waveform of beginning and ending after decimate include noise,
        // the input waveform is extended. This is the processing for the
        // compatibility with MATLAB version.
        //-----------------------------------------------------------------------------
        void GetWaveformAndSpectrumSub(double[] x, int yLength, double actualFs, int decimationRatio, double[] y)
        {
            if (decimationRatio == 1)
            {
                for (var i = 0; i < x.Length; i++)
                {
                    y[i] = x[i];
                }
                return;
            }

            var lag = (int)(Math.Ceiling(140.0 / decimationRatio) * decimationRatio);
            var newXLength = x.Length + lag * 2;
            var newY = new double[newXLength];
            var newX = new double[newXLength];
            newX.Fill(x[0], 0, lag);
            x.BlockCopy(0, newX, lag, x.Length);
            newX.Fill(x[x.Length - 1], lag + x.Length, lag);

            MatlabFunctions.Decimate(newX, decimationRatio, newY);
            newY.BlockCopy(lag / decimationRatio, y, 0, yLength);
        }

        //-----------------------------------------------------------------------------
        // HarvestGeneralBodySub() is the subfunction of HarvestGeneralBody()
        //-----------------------------------------------------------------------------
        int HarvestGeneralBodySub(double[] boundaryF0List, int numberOfChannels, int f0Length, double actualFs, int yLength, double[] temporalPositions, Complex[] ySpectrum, int fftSize, int maxCandidates, double[][] f0Candidates)
        {
            var rawF0Candidates = Util.Mak2DArray<double>(numberOfChannels, f0Length);
            GetRawF0Candidates(boundaryF0List, numberOfChannels, actualFs, yLength, temporalPositions, f0Length, ySpectrum, fftSize, rawF0Candidates);

            var numberOfCandidates = DetectOfficialF0Candidates(rawF0Candidates, numberOfChannels, f0Length, maxCandidates, f0Candidates);

            OverlapF0Candidates(f0Length, numberOfCandidates, f0Candidates);

            return numberOfCandidates;
        }

        //-----------------------------------------------------------------------------
        // GetRefinedF0() calculates F0 and its score based on instantaneous frequency.
        //-----------------------------------------------------------------------------
        void GetRefinedF0(Span<double> x, double fs, double currentPosition, double currentF0, ref double refinedF0, ref double refinedScore)
        {
            if (currentF0 <= 0.0)
            {
                refinedF0 = 0.0;
                refinedScore = 0.0;
                return;
            }

            var halfWindowLength = (int)(1.5 * fs / currentF0 + 1.0);
            var windowLengthInTime = (2.0 * halfWindowLength + 1.0) / fs;
            var baseTime = new double[halfWindowLength * 2 + 1];
            for (var i = 0; i < baseTime.Length; i++)
            {
                baseTime[i] = (-halfWindowLength + i) / fs;
            }
            var fftSize = (int)Math.Pow(2.0, 2.0 + (int)MathUtil.Log2(halfWindowLength * 2 + 1));

            GetMeanF0(x, fs, currentPosition, currentF0, fftSize, windowLengthInTime, baseTime, ref refinedF0, ref refinedScore);

            if (refinedF0 < F0Floor || refinedF0 > F0Ceil || refinedScore < 2.5)
            {
                refinedF0 = 0.0;
                refinedScore = 0.0;
            }
        }

        //-----------------------------------------------------------------------------
        // RefineF0() modifies the F0 by instantaneous frequency.
        //-----------------------------------------------------------------------------
        void RefineF0Candidates(Memory<double> x, double fs, double[] temporalPositions, int f0Length, int maxCandidates, double[][] refinedF0Candidates, double[][] f0Scores)
        {
            if (UseMultiThread)
            {
                Parallel.For(0, f0Length, i =>
                {
                    for (var j = 0; j < maxCandidates; j++)
                    {
                        GetRefinedF0(x.Span, fs, temporalPositions[i], refinedF0Candidates[i][j], ref refinedF0Candidates[i][j], ref f0Scores[i][j]);
                    }
                });
            }
            else
            {
                for (var i = 0; i < f0Length; i++)
                {
                    for (var j = 0; j < maxCandidates; j++)
                    {
                        GetRefinedF0(x.Span, fs, temporalPositions[i], refinedF0Candidates[i][j], ref refinedF0Candidates[i][j], ref f0Scores[i][j]);
                    }
                }
            }
        }

        void RemoveUnreliableCandidatesSub(int i, int j, double[][] tmpF0Candidates, int numberOfCandidates, double[][] f0Candidates, double[][] f0Scores)
        {
            const double Threshold = 0.05;
            var referenceF0 = f0Candidates[i][j];
            if (referenceF0 == 0.0)
            {
                return;
            }

            if (tmpF0Candidates.Length <= i + 1 || i < 1)
            {
                return;
            }

            var error1 = 0.0;
            var error2 = 0.0;
            SelectBestF0(referenceF0, tmpF0Candidates[i + 1], numberOfCandidates, 1.0, ref error1);
            SelectBestF0(referenceF0, tmpF0Candidates[i - 1], numberOfCandidates, 1.0, ref error2);
            if (Math.Min(error1, error2) <= Threshold)
            {
                return;
            }

            f0Candidates[i][j] = 0.0;
            f0Scores[i][j] = 0.0;
        }

        //-----------------------------------------------------------------------------
        // RemoveUnreliableCandidates().
        //-----------------------------------------------------------------------------
        void RemoveUnreliableCandidates(int f0Length, int numberOfCandidates, double[][] f0Candidates, double[][] f0Scores)
        {
            var tmpF0Candidates = Util.Mak2DArray<double>(f0Length, numberOfCandidates);
            for (int i = 0, limit = f0Length; i < limit; i++)
            {
                f0Candidates[i].BlockCopy(tmpF0Candidates[i]);
            }

            if (UseMultiThread)
            {
                Parallel.For(0, f0Length, i =>
                {
                    for (var j = 0; j < numberOfCandidates; j++)
                    {
                        RemoveUnreliableCandidatesSub(i, j, tmpF0Candidates, numberOfCandidates, f0Candidates, f0Scores);
                    }
                });
            }
            else
            {
                for (int i = 0, limit = f0Length; i < limit; i++)
                {
                    for (var j = 0; j < numberOfCandidates; j++)
                    {
                        RemoveUnreliableCandidatesSub(i, j, tmpF0Candidates, numberOfCandidates, f0Candidates, f0Scores);
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // FixF0Contour() obtains the likely F0 contour.
        //-----------------------------------------------------------------------------
        void FixF0Contour(double[][] f0Candidates, double[][] f0Scores, int f0Length, int numberOfCandidates, double[] bestF0Counter)
        {
            var tmpF0Counter1 = new double[f0Length];
            var tmpF0Counter2 = new double[f0Length];

            SearchF0Base(f0Candidates, f0Scores, f0Length, numberOfCandidates, tmpF0Counter1);
            FixStep1(tmpF0Counter1, f0Length, 0.008, tmpF0Counter2);
            FixStep2(tmpF0Counter2, f0Length, 6, tmpF0Counter1);
            FixStep3(tmpF0Counter1, f0Length, numberOfCandidates, f0Candidates, 0.18, f0Scores, tmpF0Counter2);
            FixStep4(tmpF0Counter2, f0Length, 9, bestF0Counter);
        }

        //-----------------------------------------------------------------------------
        // SmoothF0Contour() uses the zero-lag Butterworth filter for smoothing.
        //-----------------------------------------------------------------------------
        void SmoothF0Contour(double[] f0, int f0Length, double[] smoothedF0)
        {
            var b = new double[] { 0.0078202080334971724, 0.015640416066994345 };
            var a = new double[] { 1.7347257688092754, -0.76600660094326412 };
            const int Lag = 300;
            var newF0Length = f0Length + Lag * 2;
            var f0Counter = new double[Lag].Concat(f0).Concat(new double[Lag]).ToArray();

            var boundaryList = new int[newF0Length];
            var numberOfBoundaries = GetBoundaryList(f0Counter, newF0Length, boundaryList);
            var multiChannelF0 = Util.Mak2DArray<double>(numberOfBoundaries / 2, newF0Length);
            GetMultiChannelF0(f0Counter, boundaryList, numberOfBoundaries, multiChannelF0);

            for (int i = 0, limit = numberOfBoundaries / 2; i < limit; i++)
            {
                FilteringF0(a, b, multiChannelF0[i], boundaryList[i * 2], boundaryList[i * 2 + 1], f0Counter);
                for (var j = boundaryList[i * 2]; j <= boundaryList[i * 2 + 1]; j++)
                {
                    smoothedF0[j - Lag] = f0Counter[j];
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetRawF0Candidates() calculates f0 candidates in all channels.
        //-----------------------------------------------------------------------------
        void GetRawF0Candidates(double[] boundaryF0List, int numberOfBands, double actualFs, int yLength, double[] temporalPositions, int f0Length, Complex[] ySpectrum, int fftSize, double[][] rawF0Candidates)
        {
            if (UseMultiThread)
            {
                Parallel.For(0, numberOfBands, i =>
                {
                    GetF0CandidateFromRawEvent(boundaryF0List[i], actualFs, yLength, ySpectrum, fftSize, temporalPositions, f0Length, rawF0Candidates[i]);
                });
            }
            else
            {
                for (var i = 0; i < numberOfBands; i++)
                {
                    GetF0CandidateFromRawEvent(boundaryF0List[i], actualFs, yLength, ySpectrum, fftSize, temporalPositions, f0Length, rawF0Candidates[i]);
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetF0CandidateFromRawEvent() f0 candidate contour in 1-ch signal
        //-----------------------------------------------------------------------------
        void GetF0CandidateFromRawEvent(double boundaryF0, double fs, int yLength, Complex[] ySpectrum, int fftSize, double[] temporalPositions, int f0Length, double[] f0Candidates)
        {
            var filteredSignal = new double[fftSize];
            GetFilteredSignal(boundaryF0, fftSize, fs, yLength, ySpectrum, filteredSignal);

            var zeroCrossings = ZeroCrossings.GetFourZeroCrossingIntervals(filteredSignal, yLength, fs);

            GetF0CandidateContour(zeroCrossings, boundaryF0, temporalPositions, f0Length, f0Candidates);

            zeroCrossings.Release();
        }

        //-----------------------------------------------------------------------------
        // GetFilteredSignal() calculates the signal that is the convolution of the
        // input signal and band-pass filter.
        //-----------------------------------------------------------------------------
        void GetFilteredSignal(double boundaryF0, int fftSize, double fs, int yLength, Complex[] ySpectrum, double[] filteredSignal)
        {
            var filterLengthHalf = MatlabFunctions.MatlabRound(fs / boundaryF0 * 2.0);
            var bandPassFilter = new double[fftSize];
            Common.NuttallWindow(filterLengthHalf * 2 + 1, bandPassFilter);
            for (int i = -filterLengthHalf; i <= filterLengthHalf; i++)
            {
                bandPassFilter[i + filterLengthHalf] *= Math.Cos(2.0 * Math.PI * boundaryF0 * i / fs);
            }
            Array.Clear(bandPassFilter, filterLengthHalf * 2 + 1, fftSize - filterLengthHalf * 2 - 1);

            var bandPassFilterSpectrum = new Complex[fftSize];
            var forwardFFT = FFTPlan.CreatePlanForDftR2C(fftSize, bandPassFilter, bandPassFilterSpectrum, FFTFlag.Estimate);
            FFT.Execute(forwardFFT);

            // Convolution
            var tmp = ySpectrum[0].Real * bandPassFilterSpectrum[0].Real - ySpectrum[0].Imaginary * bandPassFilterSpectrum[0].Imaginary;
            bandPassFilterSpectrum[0] = new Complex(tmp, ySpectrum[0].Real * bandPassFilterSpectrum[0].Imaginary + ySpectrum[0].Imaginary * bandPassFilterSpectrum[0].Real);
            for (int i = 1, limit = fftSize / 2; i <= limit; i++)
            {
                tmp = ySpectrum[i].Real * bandPassFilterSpectrum[i].Real - ySpectrum[i].Imaginary * bandPassFilterSpectrum[i].Imaginary;
                bandPassFilterSpectrum[i] = new Complex(tmp, ySpectrum[i].Real * bandPassFilterSpectrum[i].Imaginary + ySpectrum[i].Imaginary * bandPassFilterSpectrum[i].Real);
                bandPassFilterSpectrum[fftSize - i - 1] = bandPassFilterSpectrum[i];
            }

            var inverseFFT = FFTPlan.CreatePlanForDftC2R(fftSize, bandPassFilterSpectrum, filteredSignal, FFTFlag.Estimate);
            FFT.Execute(inverseFFT);

            // Compensation of the delay.
            var indexBias = filterLengthHalf + 1;
            for (int i = 0; i < yLength; ++i)
            {
                filteredSignal[i] = filteredSignal[i + indexBias];
            }
        }

        void GetF0CandidateContour(ZeroCrossings zeroCrossings, double boundaryF0, double[] temporalPositions, int f0Length, double[] f0Candidates)
        {
            if (0 == CheckEvent(zeroCrossings.NumberOfNegatives - 2) *
                CheckEvent(zeroCrossings.NumberOfPositives - 2) *
                CheckEvent(zeroCrossings.NumberOfPeaks - 2) *
                CheckEvent(zeroCrossings.NumberOfDips - 2))
            {
                Array.Clear(f0Candidates, 0, f0Length);
                return;
            }

            var interpolatedF0Set = Util.Mak2DArray<double>(4, f0Length);
            MatlabFunctions.Interp1(
                zeroCrossings.NegativeIntervalLocations.AsSpan(0, zeroCrossings.NumberOfNegatives),
                zeroCrossings.NegativeIntervals,
                temporalPositions.AsSpan(0, f0Length),
                interpolatedF0Set[0]
            );
            MatlabFunctions.Interp1(
                zeroCrossings.PositiveIntervalLocations.AsSpan(0, zeroCrossings.NumberOfPositives),
                zeroCrossings.PositiveIntervals,
                temporalPositions.AsSpan(0, f0Length),
                interpolatedF0Set[1]
            );
            MatlabFunctions.Interp1(
                zeroCrossings.PeakIntervalLocations.AsSpan(0, zeroCrossings.NumberOfPeaks),
                zeroCrossings.PeakIntervals,
                temporalPositions.AsSpan(0, f0Length),
                interpolatedF0Set[2]
            );
            MatlabFunctions.Interp1(
                zeroCrossings.DipIntervalLocations.AsSpan(0, zeroCrossings.NumberOfDips),
                zeroCrossings.DipIntervals,
                temporalPositions.AsSpan(0, f0Length),
                interpolatedF0Set[3]
            );

            GetF0CandidateContourSub(interpolatedF0Set, f0Length, boundaryF0, f0Candidates);
        }

        void GetF0CandidateContourSub(double[][] interpolatedF0Set, int f0Length, double boundaryF0, double[] f0Candidates)
        {
            var upper = Math.Min(boundaryF0 * 1.1, F0Ceil);
            var lower = Math.Max(boundaryF0 * 0.9, F0Floor);

            for (var i = 0; i < f0Length; i++)
            {
                f0Candidates[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] + interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) / 4.0;

                if (f0Candidates[i] > upper || f0Candidates[i] < lower)
                {
                    f0Candidates[i] = 0.0;
                }
            }
        }

        //-----------------------------------------------------------------------------
        // CheckEvent() returns 1, provided that the input value is over 1.
        // This function is for RawEventByDio().
        //-----------------------------------------------------------------------------
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        int CheckEvent(int v)
        {
            return v > 0 ? 1 : 0;
        }

        //-----------------------------------------------------------------------------
        // DetectOfficialF0Candidates() detectes F0 candidates from multi-channel
        // candidates.
        //-----------------------------------------------------------------------------
        int DetectOfficialF0Candidates(double[][] rawF0Candidates, int numberOfChannels, int f0Length, int maxCandidates, double[][] f0Candidates)
        {
            var numberOfCandidates = 0;

            var vuv = new int[numberOfChannels];
            var st = new int[numberOfChannels];
            var ed = new int[numberOfChannels];

            for (var i = 0; i < f0Length; ++i)
            {
                for (var j = 0; j < numberOfChannels; j++)
                {
                    vuv[j] = rawF0Candidates[j][i] > 0 ? 1 : 0;
                }

                vuv[0] = vuv[numberOfChannels - 1] = 0;
                var numberOfVoicedSections = DetectOfficialF0CandidatesSub1(vuv, numberOfChannels, st, ed);
                numberOfCandidates = Math.Max(numberOfCandidates, DetectOfficialF0CandidatesSub2(vuv, rawF0Candidates, i, numberOfVoicedSections, st, ed, maxCandidates, f0Candidates[i]));
            }

            return numberOfCandidates;
        }

        //-----------------------------------------------------------------------------
        // DetectF0CandidatesSub1() calculates VUV areas.
        //-----------------------------------------------------------------------------
        int DetectOfficialF0CandidatesSub1(int[] vuv, int numberOfChannels, int[] st, int[] ed)
        {
            var numberOfVoicedSections = 0;
            for (var i = 1; i < numberOfChannels; i++)
            {
                var tmp = vuv[i] - vuv[i - 1];
                if (tmp == 1)
                {
                    st[numberOfVoicedSections] = i;
                }
                if (tmp == -1)
                {
                    ed[numberOfVoicedSections] = i;
                    numberOfVoicedSections++;
                }
            }

            return numberOfVoicedSections;
        }

        //-----------------------------------------------------------------------------
        // DetectOfficialF0CandidatesSub2() calculates F0 candidates in a frame
        //-----------------------------------------------------------------------------
        int DetectOfficialF0CandidatesSub2(int[] vuv, double[][] rawF0Candidates, int index, int numberOfVoicedSections, int[] st, int[] ed, int maxCandidates, double[] f0List)
        {
            var numberOfCandidates = 0;

            for (var i = 0; i < numberOfVoicedSections; i++)
            {
                var count = ed[i] - st[i];
                if (count < 10)
                {
                    continue;
                }

                var tmpF0 = rawF0Candidates.Skip(st[i]).Take(count).Sum((r) => r[index]);
                tmpF0 /= count;
                f0List[numberOfCandidates] = tmpF0;
                numberOfCandidates++;
            }

            return numberOfCandidates;
        }

        //-----------------------------------------------------------------------------
        // OverlapF0Candidates() spreads the candidates to anteroposterior frames.
        //-----------------------------------------------------------------------------
        void OverlapF0Candidates(int f0Length, int numberOfCandidates, double[][] f0Candidates)
        {
            const int N = 3;
            for (var i = 1; i <= N; i++)
            {
                for (var j = 0; j < numberOfCandidates; j++)
                {
                    for (var k = i; k < f0Length; k++)
                    {
                        f0Candidates[k][j + (numberOfCandidates * i)] = f0Candidates[k - i][j];
                    }
                    for (var k = 0; k < f0Length - i; k++)
                    {
                        f0Candidates[k][j + (numberOfCandidates * (i + N))] = f0Candidates[k + i][j];
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetBoundaryList() detects boundaries between voiced and unvoiced sections.
        //-----------------------------------------------------------------------------
        int GetBoundaryList(double[] f0, int f0Length, int[] boundaryList)
        {
            var numberOfBoundaries = 0;
            var vuv = Enumerable.Range(0, f0Length).Select((x) => f0[x] > 0.0 ? 1 : 0).ToArray();
            vuv[0] = 0;
            vuv[f0Length - 1] = 0;

            for (var i = 1; i < f0Length; i++)
            {
                if (vuv[i] - vuv[i - 1] != 0)
                {
                    boundaryList[numberOfBoundaries] = i - numberOfBoundaries % 2;
                    numberOfBoundaries++;
                }
            }

            return numberOfBoundaries;
        }

        //-----------------------------------------------------------------------------
        // GetMultiChannelF0() separates each voiced section into independent channel.
        //-----------------------------------------------------------------------------
        void GetMultiChannelF0(double[] f0, int[] boundaryList, int numberOfBoundaries, double[][] multiChannelF0)
        {
            for (int i = 0, iLimit = numberOfBoundaries / 2; i < iLimit; i++)
            {
                Array.Clear(multiChannelF0[i], 0, multiChannelF0[i].Length);
                var j = i * 2;
                f0.BlockCopy(boundaryList[j], multiChannelF0[i], boundaryList[j], boundaryList[j + 1] - boundaryList[j] + 1);
            }
        }

        //-----------------------------------------------------------------------------
        // This function uses zero-lag Butterworth filter.
        //-----------------------------------------------------------------------------
        void FilteringF0(double[] a, double[] b, double[] x, int st, int ed, double[] y)
        {
            var w = new double[2];

            x.Fill(x[st], 0, st);
            x.Fill(x[ed], ed + 1);

            var tmpX = new double[x.Length];
            for (var i = 0; i < x.Length; i++)
            {
                var wt = x[i] + a[0] * w[0] + a[1] * w[1];
                tmpX[x.Length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
                w[1] = w[0];
                w[0] = wt;
            }

            w = new double[2];
            for (var i = 0; i < x.Length; i++)
            {
                var wt = tmpX[i] + a[0] * w[0] + a[1] * w[1];
                y[x.Length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
                w[1] = w[0];
                w[0] = wt;
            }
        }

        //-----------------------------------------------------------------------------
        // SearchF0Base() gets the F0 with the highest score.
        //-----------------------------------------------------------------------------
        void SearchF0Base(double[][] f0Candidates, double[][] f0Scores, int f0Length, int numberOfCandidates, double[] baseF0Counter)
        {
            for (var i = 0; i < f0Length; i++)
            {
                var tmpBestScore = 0.0;
                baseF0Counter[i] = 0.0;

                for (var j = 0; j < numberOfCandidates; j++)
                {
                    if (f0Scores[i][j] > tmpBestScore)
                    {
                        baseF0Counter[i] = f0Candidates[i][j];
                        tmpBestScore = f0Scores[i][j];
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // Step 1: Rapid change of F0 contour is replaced by 0.
        //-----------------------------------------------------------------------------
        void FixStep1(double[] f0Base, int f0Length, double allowedRange, double[] f0Step1)
        {
            f0Step1.Fill(0.0);

            for (var i = 2; i < f0Length; i++)
            {
                if (f0Base[i] == 0.0)
                {
                    f0Step1[i] = 0.0;
                    continue;
                }

                var referenceF0 = f0Base[i - 1] * 2 - f0Base[i - 2];
                f0Step1[i] = Math.Abs((f0Base[i] - referenceF0) / referenceF0) > allowedRange && Math.Abs((f0Base[i] - f0Base[i - 1])) / f0Base[i - 1] > allowedRange ? 0.0 : f0Base[i];
            }
        }

        //-----------------------------------------------------------------------------
        // Step 2: Voiced sections with a short period are removed.
        //-----------------------------------------------------------------------------
        void FixStep2(double[] f0Step1, int f0Length, int voiceRangeMinimum, double[] f0Step2)
        {
            f0Step1.BlockCopy(f0Step2);
            var boundaryList = new int[f0Length];
            var numberOfBoundary = GetBoundaryList(f0Step1, f0Length, boundaryList);

            for (var i = 0; i < numberOfBoundary; i++)
            {
                if (boundaryList[i * 2 + 1] - boundaryList[i * 2] >= voiceRangeMinimum)
                {
                    continue;
                }

                Array.Clear(f0Step2, boundaryList[i * 2], boundaryList[i * 2 + 1] - boundaryList[i * 2] + 1);
            }
        }

        //-----------------------------------------------------------------------------
        // Step 3: Voiced sections are extended based on the continuity of F0 contour
        //-----------------------------------------------------------------------------
        void FixStep3(double[] f0Step2, int f0Length, int numberOfCandidates, double[][] f0Candidates, double allowedRange, double[][] f0Scores, double[] f0Step3)
        {
            f0Step2.BlockCopy(f0Step3);
            var boundaryList = new int[f0Length];
            var numberOfBoundaries = GetBoundaryList(f0Step2, f0Length, boundaryList);

            var multiChannelF0 = Util.Mak2DArray<double>(numberOfBoundaries / 2, f0Length);
            GetMultiChannelF0(f0Step2, boundaryList, numberOfBoundaries, multiChannelF0);

            var numberOfChannels = Extend(multiChannelF0, numberOfBoundaries / 2, f0Length, boundaryList, f0Candidates, numberOfCandidates, allowedRange, multiChannelF0, boundaryList);

            if (numberOfChannels != 0)
            {
                MergeF0(multiChannelF0, boundaryList, numberOfChannels, f0Length, f0Candidates, f0Scores, numberOfCandidates, f0Step3);
            }
        }

        //-----------------------------------------------------------------------------
        // Step 4: F0s in short unvoiced section are faked
        //-----------------------------------------------------------------------------
        void FixStep4(double[] f0Step3, int f0Length, int threshold, double[] f0Step4)
        {
            f0Step3.BlockCopy(f0Step4);
            var boundaryList = new int[f0Length];
            var numberOfBoundaries = GetBoundaryList(f0Step3, f0Length, boundaryList);

            for (int i = 0, iLimit = numberOfBoundaries / 2 - 1; i < iLimit; i++)
            {
                var distance = boundaryList[(i + 1) * 2] - boundaryList[i * 2 + 1] - 1;
                if (distance >= threshold)
                {
                    continue;
                }

                var tmp0 = f0Step3[boundaryList[i * 2 + 1]] + 1;
                var tmp1 = f0Step3[boundaryList[(i + 1) * 2]] - 1;
                var coefficient = (tmp1 - tmp0) / (distance + 1.0);
                for (int j = boundaryList[i * 2 + 1] + 1, jLimit = boundaryList[(i + 1) * 2] - 1, count = 1; j <= jLimit; j++, count++)
                {
                    f0Step4[j] = tmp0 + coefficient * count;
                }
            }
        }

        //-----------------------------------------------------------------------------
        // Extend() : The Hand erasing the Space.
        //-----------------------------------------------------------------------------
        int Extend(double[][] multiChannelF0, int numberOfSections, int f0Length, int[] boundaryList, double[][] f0Candidates, int numberOfCandidates, double allowedRange, double[][] extendedF0, int[] shiftedBoundaryList)
        {
            const int Threshold = 100;

            for (var i = 0; i < numberOfSections; i++)
            {
                shiftedBoundaryList[i * 2 + 1] = ExtendF0(multiChannelF0[i], f0Length, boundaryList[i * 2 + 1], Math.Min(f0Length - 2, boundaryList[i * 2 + 1] + Threshold), 1, f0Candidates, numberOfCandidates, allowedRange, extendedF0[i]);
                shiftedBoundaryList[i * 2] = ExtendF0(multiChannelF0[i], f0Length, boundaryList[i * 2], Math.Max(1, boundaryList[i * 2] - Threshold), -1, f0Candidates, numberOfCandidates, allowedRange, extendedF0[i]);
            }

            return ExtendSub(multiChannelF0, shiftedBoundaryList, numberOfSections, extendedF0, shiftedBoundaryList);
        }

        int ExtendSub(double[][] extendedF0, int[] boundaryList, int numberOfSections, double[][] selectedExtendedF0, int[] selectedBoundaryList)
        {
            const double Threshold = 2200.0;

            var count = 0;
            var meanF0 = 0.0;
            for (var i = 0; i < numberOfSections; i++)
            {
                var st = boundaryList[i * 2];
                var ed = boundaryList[i * 2 + 1];
                for (var j = st; j < ed; j++)
                {
                    meanF0 += extendedF0[i][j];
                }
                meanF0 /= ed - st;
                if (Threshold / meanF0 < ed - st)
                {
                    Swap(count, i, selectedExtendedF0, selectedBoundaryList);
                    count++;
                }
            }

            return count;
        }

        //-----------------------------------------------------------------------------
        // Swap the f0 contour and boundary.
        // It is used in ExtendSub() and MergeF0();
        //-----------------------------------------------------------------------------
        void Swap(int index1, int index2, double[][] f0, int[] boundary)
        {
            var tmpF0 = f0[index1];
            f0[index1] = f0[index2];
            f0[index2] = tmpF0;

            var tmpIndex = boundary[index1 * 2];
            boundary[index1 * 2] = boundary[index2 * 2];
            boundary[index2 * 2] = tmpIndex;
            tmpIndex = boundary[index1 * 2 + 1];
            boundary[index1 * 2 + 1] = boundary[index2 * 2 + 1];
            boundary[index2 * 2 + 1] = tmpIndex;
        }

        //-----------------------------------------------------------------------------
        // ExtendF0() : The Hand erasing the Space.
        // The subfunction of Extend().
        //-----------------------------------------------------------------------------
        int ExtendF0(double[] f0, int f0Length, int origin, int lastPoint, int shift, double[][] f0Candidates, int numberOfCandidates, double allowedRange, double[] extendedF0)
        {
            const int Threshold = 4;
            var tmpF0 = extendedF0[origin];
            var shiftedOrigin = origin;

            var distance = Math.Abs(lastPoint - origin);
            var indexList = Enumerable.Range(0, distance + 1).Select((x) => origin + shift * x).ToArray();

            var count = 0;
            var dummy = 0.0;
            for (var i = 0; i <= distance && count != Threshold; i++)
            {
                var target = indexList[i] + shift;
                extendedF0[target] = SelectBestF0(tmpF0, f0Candidates[target], numberOfCandidates, allowedRange, ref dummy);
                if (extendedF0[target] == 0.0)
                {
                    count++;
                }
                else
                {
                    tmpF0 = extendedF0[target];
                    count = 0;
                    shiftedOrigin = target;
                }
            }

            return shiftedOrigin;
        }

        //-----------------------------------------------------------------------------
        // Overlapped F0 contours are merged by the likability score.
        //-----------------------------------------------------------------------------
        void MergeF0(double[][] multiChannelF0, int[] boundaryList, int numberOfChannels, int f0Length, double[][] f0Candidates, double[][] f0Scores, int numberOfCandidates, double[] mergedF0)
        {
            var order = new int[numberOfChannels];
            MakeSortedOrder(boundaryList, numberOfChannels, order);

            multiChannelF0[0].BlockCopy(mergedF0);

            for (var i = 1; i < numberOfChannels; i++)
            {
                if (boundaryList[order[i] * 2] - boundaryList[1] > 0)
                {
                    for (var j = boundaryList[order[i] * 2]; j <= boundaryList[order[i] * 2 + 1]; j++)
                    {
                        mergedF0[j] = multiChannelF0[order[i]][j];
                    }
                    boundaryList[0] = boundaryList[order[i] * 2];
                    boundaryList[1] = boundaryList[order[i] * 2 + 1];
                }
                else
                {
                    boundaryList[1] = MergeF0Sub(mergedF0, f0Length, boundaryList[0], boundaryList[1], multiChannelF0[order[i]], boundaryList[order[i] * 2], boundaryList[order[i] * 2 + 1], f0Candidates, f0Scores, numberOfCandidates, mergedF0);
                }
            }
        }

        //-----------------------------------------------------------------------------
        // Subfunction of MergeF0()
        //-----------------------------------------------------------------------------
        int MergeF0Sub(double[] f01, int f0Length, int st1, int ed1, double[] f02, int st2, int ed2, double[][] f0Candidates, double[][] f0Scores, int numberOfCandidates, double[] mergedF0)
        {
            if (st1 <= st2 && ed1 >= ed2)
            {
                return ed1;
            }

            var score1 = 0.0;
            var score2 = 0.0;
            for (var i = st2; i <= ed1; i++)
            {
                score1 += SearchScore(f01[i], f0Candidates[i], f0Scores[i], numberOfCandidates);
                score2 += SearchScore(f02[i], f0Candidates[i], f0Scores[i], numberOfCandidates);
            }

            if (score1 > score2)
            {
                f02.BlockCopy(ed1, mergedF0, ed1, ed2 - ed1 + 1);
            }
            else
            {
                f02.BlockCopy(st2, mergedF0, st2, ed2 - st2 + 1);
            }

            return ed2;
        }

        double SearchScore(double f0, double[] f0Candidates, double[] f0Scores, int numberOfCandidates)
        {
            var score = 0.0;

            for (var i = 0; i < numberOfCandidates; i++)
            {
                if (f0 == f0Candidates[i] && score < f0Scores[i])
                {
                    score = f0Scores[i];
                }
            }

            return score;
        }

        //-----------------------------------------------------------------------------
        // Indices are sorted.
        //-----------------------------------------------------------------------------
        void MakeSortedOrder(int[] boundaryList, int numberOfSections, int[] order)
        {
            Enumerable.Range(0, numberOfSections).ForEach((x) => order[x] = x);

            for (var i = 0; i < numberOfSections; i++)
            {
                for (var j = i - 1; j >= 0; j--)
                {
                    if (boundaryList[order[j] * 2] > boundaryList[order[i] * 2])
                    {
                        var tmp = order[i];
                        order[i] = order[j];
                        order[j] = tmp;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // SelectBestF0() obtains the nearlest F0 in reference_f0.
        //-----------------------------------------------------------------------------
        double SelectBestF0(double referenceF0, double[] f0Candidates, int numberOfCandidates, double allowedRange, ref double bestError)
        {
            var bestF0 = 0.0;
            bestError = allowedRange;

            for (var i = 0; i < numberOfCandidates; i++)
            {
                var tmp = Math.Abs(referenceF0 - f0Candidates[i]) / referenceF0;
                if (tmp > bestError)
                {
                    continue;
                }

                bestF0 = f0Candidates[i];
                bestError = tmp;
            }


            return bestF0;
        }

        void GetMeanF0(Span<double> x, double fs, double currentPosition, double currentF0, int fftSize, double windowLengthInTime, double[] baseTime, ref double refinedF0, ref double refinedScore)
        {
            var complexPool = ArrayPool<Complex>.Shared;
            var doublePool = ArrayPool<double>.Shared;
            var forwardRealFft = ForwardRealFFT.Create(fftSize);
            var mainSpectrum = complexPool.Rent(fftSize);
            var diffSpectrum = complexPool.Rent(fftSize);
            mainSpectrum.AsSpan().Clear();
            diffSpectrum.AsSpan().Clear();

            var baseIndex = ArrayPool<int>.Shared.Rent(baseTime.Length);
            var mainWindow = doublePool.Rent(baseTime.Length);
            var diffWindow = doublePool.Rent(baseTime.Length);

            GetBaseIndex(currentPosition, baseTime, fs, baseIndex.AsSpan(0, baseTime.Length));
            GetMainWindow(currentPosition, baseIndex, fs, windowLengthInTime, mainWindow.AsSpan(0, baseTime.Length));
            GetDiffWindow(mainWindow.AsSpan(0, baseTime.Length), diffWindow.AsSpan(0, baseTime.Length));

            GetSpectra(x, fftSize, baseIndex.AsSpan(0, baseTime.Length), mainWindow.AsSpan(0, baseTime.Length), diffWindow.AsSpan(0, baseTime.Length), forwardRealFft, mainSpectrum, diffSpectrum);

            var halfFFTSize = fftSize / 2 + 1;
            var powerSpectrum = doublePool.Rent(halfFFTSize);
            var numeratorI = doublePool.Rent(halfFFTSize);
            for (var j = 0; j < halfFFTSize; j++)
            {
                powerSpectrum[j] = mainSpectrum[j].Real * mainSpectrum[j].Real + mainSpectrum[j].Imaginary * mainSpectrum[j].Imaginary;
                numeratorI[j] = mainSpectrum[j].Real * diffSpectrum[j].Imaginary - mainSpectrum[j].Imaginary * diffSpectrum[j].Real;
            }

            var numberOfHarmonics = Math.Min((int)(fs / 2.0 / currentF0), 6);
            FixF0(powerSpectrum, numeratorI, fftSize, fs, currentF0, numberOfHarmonics, ref refinedF0, ref refinedScore);

            forwardRealFft.Release();
            complexPool.Return(mainSpectrum);
            complexPool.Return(diffSpectrum);
            ArrayPool<int>.Shared.Return(baseIndex);
            doublePool.Return(mainWindow);
            doublePool.Return(diffWindow);
            doublePool.Return(powerSpectrum);
            doublePool.Return(numeratorI);
        }

        //-----------------------------------------------------------------------------
        // GetBaseIndex() calculates the temporal positions for windowing.
        //-----------------------------------------------------------------------------
        void GetBaseIndex(double currentPosition, double[] baseTime, double fs, Span<int> baseIndex)
        {
            // First-aid treatment
            var basicIndex = MatlabFunctions.MatlabRound((currentPosition + baseTime[0]) * fs + 0.001);
            for (var i = 0; i < baseTime.Length; i++)
            {
                baseIndex[i] = basicIndex + i;
            }
        }

        //-----------------------------------------------------------------------------
        // GetMainWindow() generates the window function.
        //-----------------------------------------------------------------------------
        void GetMainWindow(double currentPosition, int[] baseIndex, double fs, double windowLengthInTime, Span<double> mainWindow)
        {
            for (var i = 0; i < mainWindow.Length; i++)
            {
                var tmp = (baseIndex[i] - 1) / fs - currentPosition;
                mainWindow[i] = 0.42 + 0.5 * Math.Cos(2.0 * Math.PI * tmp / windowLengthInTime) + 0.08 * Math.Cos(4.0 * Math.PI * tmp / windowLengthInTime);
            }
        }

        //-----------------------------------------------------------------------------
        // GetDiffWindow() generates the differentiated window.
        // Diff means differential.
        //-----------------------------------------------------------------------------
        void GetDiffWindow(ReadOnlySpan<double> mainWindow, Span<double> diffWindow)
        {
            diffWindow[0] = -mainWindow[1] / 2.0;
            for (int i = 1, limit = diffWindow.Length - 1; i < limit; i++)
            {
                diffWindow[i] = -(mainWindow[i + 1] - mainWindow[i - 1]) / 2.0;
            }
            diffWindow[diffWindow.Length - 1] = mainWindow[diffWindow.Length - 2] / 2.0;
        }

        //-----------------------------------------------------------------------------
        // GetSpectra() calculates two spectra of the waveform windowed by windows
        // (main window and diff window).
        //-----------------------------------------------------------------------------
        void GetSpectra(Span<double> x, int fftSize, ReadOnlySpan<int> baseIndex, ReadOnlySpan<double> mainWindow, ReadOnlySpan<double> diffWindow, ForwardRealFFT forwardRealFft, Complex[] mainSpectrum, Complex[] diffSpectrum)
        {
            var xLength = x.Length;
            var safeIndex = ArrayPool<int>.Shared.Rent(baseIndex.Length);
            for (var i = 0; i < baseIndex.Length; i++)
            {
                safeIndex[i] = Math.Max(0, Math.Min(xLength - 1, baseIndex[i] - 1));
            }

            for (var i = 0; i < baseIndex.Length; i++)
            {
                forwardRealFft.Waveform[i] = x[safeIndex[i]] * mainWindow[i];
            }
            forwardRealFft.Waveform.AsSpan(baseIndex.Length, fftSize - baseIndex.Length).Clear();

            FFT.Execute(forwardRealFft.ForwardFFT);
            forwardRealFft.Spectrum.AsSpan(0, fftSize / 2 + 1).CopyTo(mainSpectrum);

            for (var i = 0; i < baseIndex.Length; i++)
            {
                forwardRealFft.Waveform[i] = x[safeIndex[i]] * diffWindow[i];
            }
            forwardRealFft.Waveform.AsSpan(baseIndex.Length, fftSize - baseIndex.Length).Clear()    ;

            FFT.Execute(forwardRealFft.ForwardFFT);
            forwardRealFft.Spectrum.AsSpan(0, fftSize / 2 + 1).CopyTo(diffSpectrum);

            ArrayPool<int>.Shared.Return(safeIndex);
        }

        void FixF0(double[] powerSpectrum, double[] numeratorI, int fftSize, double fs, double currentF0, int numberOfHarmonics, ref double refinedF0, ref double refinedScore)
        {
            var amplitudeList = new double[numberOfHarmonics];
            var instantaneousFrequencyList = new double[numberOfHarmonics];

            for (var i = 0; i < numberOfHarmonics; i++)
            {
                var index = MatlabFunctions.MatlabRound(currentF0 * fftSize / fs * (i + 1));
                if (powerSpectrum[index] != 0.0)
                {
                    instantaneousFrequencyList[i] = index * fs / fftSize + numeratorI[index] / powerSpectrum[index] * fs / 2.0 / Math.PI;
                    amplitudeList[i] = Math.Sqrt(powerSpectrum[index]);
                }
            }

            var denominator = 0.0;
            var numerator = 0.0;
            for (var i = 0; i < numberOfHarmonics; i++)
            {
                numerator += amplitudeList[i] * instantaneousFrequencyList[i];
                denominator += amplitudeList[i] * (i + 1);
                refinedScore += Math.Abs((instantaneousFrequencyList[i] / (i + 1) - currentF0) / currentF0);
            }

            refinedF0 = numerator / (denominator + SafeGuardMinimum);
            refinedScore = 1.0 / (refinedScore / numberOfHarmonics + SafeGuardMinimum);
        }
    }
}

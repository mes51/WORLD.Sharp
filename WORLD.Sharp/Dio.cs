using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    public class Dio
    {
        const double DefaultF0Ceil = 800.0;

        const double DefaultF0Floor = 71.0;

        const double CutOff = 50.0;

        const double SafeGuardMinimum = 0.000000000001;

        const double MaximumValue = 100000.0;

        public double ChannelInOctave { get; set; } = 2.0;

        public double F0Floor { get; set; } = DefaultF0Floor;

        public double F0Ceil { get; set; } = DefaultF0Ceil;

        public double FramePeriod { get; set; } = 5.0;

        public bool UseMultiThread { get; set; }

        /// <summary>
        /// You can use the value from 1 to 12.
        /// Default value 11 is for the fs of 44.1 kHz.
        /// The lower value you use, the better performance you can obtain.
        /// </summary>
        public int Speed { get; set; } = 1;

        /// <summary>
        /// You can give a positive real number as the threshold.
        /// The most strict value is 0, and there is no upper limit.
        /// On the other hand, I think that the value from 0.02 to 0.2 is reasonable.
        /// </summary>
        public double AllowRange { get; set; } = 0.1;

        public void Estimate(double[] x, int fs, double[] temporalPositions, double[] f0)
        {
            DioGeneralBody(x, fs, temporalPositions, f0);
        }

        //-----------------------------------------------------------------------------
        // DioGeneralBody() estimates the F0 based on Distributed Inline-filter
        // Operation.
        //-----------------------------------------------------------------------------
        void DioGeneralBody(double[] x, int fs, double[] temporalPositions, double[] f0)
        {
            var numberOfBands = 1 + (int)(MathUtil.Log2(F0Ceil / F0Floor) * ChannelInOctave);
            var boundaryF0List = new double[numberOfBands];
            for (var i = 0; i < numberOfBands; i++)
            {
                boundaryF0List[i] = F0Floor * Math.Pow(2.0, (i + 1) / ChannelInOctave);
            }

            // normalization
            var decimationRatio = Math.Max(Math.Min(Speed, 12), 1);
            var yLength = (1 + (int)(x.Length / decimationRatio));
            var actualFs = (double)fs / decimationRatio;
            var fftSize = Common.GetSuitableFFTSize(yLength + MatlabFunctions.MatlabRound(actualFs / CutOff) * 2 + 1 + (4 * (int)(1.0 + actualFs / boundaryF0List[0] / 2.0)));

            // Calculation of the spectrum used for the f0 estimation
            var ySpectrum = new Complex[fftSize];
            GetSpectrumForEstimation(x, yLength, actualFs, fftSize, decimationRatio, ySpectrum);

            var f0Length = GetSamplesForDIO(fs, x.Length, FramePeriod);
            var f0Candidates = Util.Mak2DArray<double>(numberOfBands, f0Length);
            var f0Scores = Util.Mak2DArray<double>(numberOfBands, f0Length);

            for (var i = 0; i < f0Length; i++)
            {
                temporalPositions[i] = i * FramePeriod / 1000.0;
            }

            GetF0CandidatesAndScores(boundaryF0List, numberOfBands, actualFs, yLength, temporalPositions, f0Length, ySpectrum, fftSize, f0Candidates, f0Scores);

            // Selection of the best value based on fundamental-ness.
            // This function is related with SortCandidates() in MATLAB.
            var bestF0Contour = new double[f0Length];
            GetBestF0Contour(f0Length, f0Candidates, f0Scores, numberOfBands, bestF0Contour);

            FixF0Contour(numberOfBands, fs, f0Candidates, bestF0Contour, f0Length, f0);
        }

        //-----------------------------------------------------------------------------
        // GetSpectrumForEstimation() calculates the spectrum for estimation.
        // This function carries out downsampling to speed up the estimation process　
        // and calculates the spectrum of the downsampled signal.
        //-----------------------------------------------------------------------------
        void GetSpectrumForEstimation(double[] x, int yLength, double actualFs, int fftSize, int decimationRatio, Complex[] ySpectrum)
        {
            var y = new double[fftSize];

            // Downsampling
            if (decimationRatio != 1)
            {
                MatlabFunctions.Decimate(x, decimationRatio, y);
            }
            else
            {
                x.BlockCopy(y);
            }

            // Removal of the DC component (y = y - mean value of y)
            var meanY = y.Take(yLength).Average();
            for (var i = 0; i < yLength; i++)
            {
                y[i] -= meanY;
            }
            y.AsSpan(yLength).Clear();

            var forwardFFT = FFTPlan.CreatePlanForDftR2C(fftSize, y, ySpectrum, FFTFlag.Estimate);
            FFT.Execute(forwardFFT);

            var cutoffInSample = MatlabFunctions.MatlabRound(actualFs / CutOff);
            DesignLowCutFilter(cutoffInSample * 2 + 1, fftSize, y);

            var filterSpectrum = new Complex[fftSize];
            forwardFFT.COut = filterSpectrum;
            FFT.Execute(forwardFFT);

            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                ySpectrum[i] *= filterSpectrum[i];
            }
        }

        //-----------------------------------------------------------------------------
        // GetF0CandidatesAndScores() calculates all f0 candidates and their scores.
        //-----------------------------------------------------------------------------
        void GetF0CandidatesAndScores(double[] boundaryF0List, int numberOfBands, double actualFs, int yLength, double[] temporalPositions, int f0Length, Complex[] ySpectrum, int fftSize, double[][] rawF0Candidates, double[][] rawF0Scores)
        {
            // Calculation of the acoustics events (zero-crossing)
            if (UseMultiThread)
            {
                Parallel.For(
                    0,
                    numberOfBands,
                    () => new { F0Candidate = new double[f0Length], F0Score = new double[f0Length] },
                    (i, loop, t) =>
                    {
                        var f0Candidate = t.F0Candidate;
                        var f0Score = t.F0Score;

                        GetF0CandidateFromRawEvent(boundaryF0List[i], actualFs, ySpectrum, yLength, fftSize, temporalPositions, f0Length, f0Score, f0Candidate);

                        for (var j = 0; j < f0Length; j++)
                        {
                            rawF0Scores[i][j] = f0Score[j] / (f0Candidate[j] + SafeGuardMinimum); rawF0Candidates[i][j] = f0Candidate[j];
                        }
                        return t;
                    },
                    _ => { }
                );
            }
            else
            {
                var f0Candidate = new double[f0Length];
                var f0Score = new double[f0Length];

                for (var i = 0; i < numberOfBands; i++)
                {
                    GetF0CandidateFromRawEvent(boundaryF0List[i], actualFs, ySpectrum, yLength, fftSize, temporalPositions, f0Length, f0Score, f0Candidate);

                    for (var j = 0; j < f0Length; j++)
                    {
                        rawF0Scores[i][j] = f0Score[j] / (f0Candidate[j] + SafeGuardMinimum); rawF0Candidates[i][j] = f0Candidate[j];
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetBestF0Contour() calculates the best f0 contour based on scores of
        // all candidates. The F0 with highest score is selected.
        //-----------------------------------------------------------------------------
        void GetBestF0Contour(int f0Length, double[][] f0Candidates, double[][] f0Scores, int numberOfBands, double[] bestF0Contour)
        {
            for (var i = 0; i < f0Length; i++)
            {
                var tmp = f0Scores[0][i];
                bestF0Contour[i] = f0Candidates[0][i];
                for (var j = 1; j < numberOfBands; j++)
                {
                    if (tmp > f0Scores[j][i])
                    {
                        tmp = f0Scores[j][i];
                        bestF0Contour[i] = f0Candidates[j][i];
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // FixF0Contour() calculates the definitive f0 contour based on all f0
        // candidates. There are four steps.
        //-----------------------------------------------------------------------------
        void FixF0Contour(int numberOfCandidates, int fs, double[][] f0Candidates, double[] bestF0Contour, int f0Length, double[] fixedF0Contour)
        {
            var voiceRangeMinimum = (int)(0.5 + 1000.0 / FramePeriod / F0Floor) * 2 + 1;

            if (f0Length <= voiceRangeMinimum)
            {
                return;
            }

            var f0Tmp1 = new double[f0Length];
            var f0Tmp2 = new double[f0Length];

            FixStep1(bestF0Contour, f0Length, voiceRangeMinimum, f0Tmp1);
            FixStep2(f0Tmp1, f0Length, voiceRangeMinimum, f0Tmp2);

            int positiveCount;
            int negativeCount;
            var positiveIndex = new int[f0Length];
            var negativeIndex = new int[f0Length];
            GetNumberOfVoicedSections(f0Tmp2, f0Length, positiveIndex, negativeIndex, out positiveCount, out negativeCount);
            FixStep3(f0Tmp2, f0Length, f0Candidates, numberOfCandidates, negativeIndex, negativeCount, f0Tmp1);
            FixStep4(f0Tmp1, f0Length, f0Candidates, numberOfCandidates, positiveIndex, positiveCount, fixedF0Contour);
        }

        //-----------------------------------------------------------------------------
        // DesignLowCutFilter() calculates the coefficients the filter.
        //-----------------------------------------------------------------------------
        void DesignLowCutFilter(int N, int fftSize, double[] lowCutFilter)
        {
            for (var i = 1; i <= N; i++)
            {
                lowCutFilter[i - 1] = 0.5 - 0.5 * Math.Cos(i * 2.0 * Math.PI / (N + 1));
            }
            lowCutFilter.AsSpan(N).Clear();

            var sumOfAmplitude = lowCutFilter.Take(N).Sum();
            for (var i = 0; i < N; i++)
            {
                lowCutFilter[i] = -lowCutFilter[i] / sumOfAmplitude;
            }
            for (int i = 0, limit = (N - 1) / 2; i < limit; i++)
            {
                lowCutFilter[fftSize - limit + i] = lowCutFilter[i];
            }
            for (var i = 0; i < N; i++)
            {
                lowCutFilter[i] = lowCutFilter[i + (N - 1) / 2];
            }
            lowCutFilter[0] += 1.0;
        }

        //-----------------------------------------------------------------------------
        // GetF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
        //-----------------------------------------------------------------------------
        void GetF0CandidateFromRawEvent(double boundaryF0, double fs, Complex[] ySpectrum, int yLength, int fftSize, double[] temporalPositions, int f0Length, double[] f0Score, double[] f0Candidate)
        {
            var filteredSignal = new double[fftSize];
            GetFilteredSignal(MatlabFunctions.MatlabRound(fs / boundaryF0 / 2.0), fftSize, ySpectrum, yLength, filteredSignal);

            var zeroCrossings = ZeroCrossings.GetFourZeroCrossingIntervals(filteredSignal, yLength, fs);
            GetF0CandidateContour(zeroCrossings, boundaryF0, temporalPositions, f0Length, f0Candidate, f0Score);
        }

        //-----------------------------------------------------------------------------
        // FixStep1() is the 1st step of the postprocessing.
        // This function eliminates the unnatural change of f0 based on allowed_range.
        //-----------------------------------------------------------------------------
        void FixStep1(double[] bestF0Contour, int f0Length, int voiceRangeMinimum, double[] f0Step1)
        {
            var f0Base = new double[f0Length];

            // Initialization
            bestF0Contour.BlockCopy(voiceRangeMinimum, f0Base, voiceRangeMinimum, f0Length - voiceRangeMinimum * 2);

            // Processing to prevent the jumping of f0
            f0Step1.AsSpan(0, voiceRangeMinimum).Clear();
            for (int i = voiceRangeMinimum; i < f0Length; i++)
            {
                f0Step1[i] = Math.Abs((f0Base[i] - f0Base[i - 1]) / (SafeGuardMinimum + f0Base[i])) < AllowRange ? f0Base[i] : 0.0;
            }
        }

        //-----------------------------------------------------------------------------
        // FixStep2() is the 2nd step of the postprocessing.
        // This function eliminates the suspected f0 in the anlaut and auslaut.
        //-----------------------------------------------------------------------------
        void FixStep2(double[] f0Step1, int f0Length, int voiceRangeMinimum, double[] f0Step2)
        {
            f0Step1.BlockCopy(0, f0Step2, 0, f0Length);

            var center = (voiceRangeMinimum - 1) / 2;
            for (int i = center, limit = f0Length - center; i < limit; i++)
            {
                for (var j = -center; j <= center; j++)
                {
                    if (f0Step1[i + j] == 0)
                    {
                        f0Step2[i] = 0.0;
                        break;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // FixStep3() is the 3rd step of the postprocessing.
        // This function corrects the f0 candidates from backward to forward.
        //-----------------------------------------------------------------------------
        void FixStep3(double[] f0Step2, int f0Length, double[][] f0Candidates, int numberOfCandidates, int[] negativeIndex, int negativeCount, double[] f0Step3)
        {
            f0Step2.BlockCopy(0, f0Step3, 0, f0Length);

            for (var i = 0; i < negativeCount; i++)
            {
                var limit = i == negativeCount - 1 ? f0Length - 1 : negativeIndex[i + 1];
                for (var j = negativeIndex[i]; j < limit; j++)
                {
                    f0Step3[j + 1] = SelectBestF0(f0Step3[j], f0Step3[j - 1], f0Candidates, numberOfCandidates, j + 1);
                    if (f0Step3[j + 1] == 0)
                    {
                        break;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // FixStep4() is the 4th step of the postprocessing.
        // This function corrects the f0 candidates from forward to backward.
        //-----------------------------------------------------------------------------
        void FixStep4(double[] f0Step3, int f0Length, double[][] f0Candidates, int numberOfCandidates, int[] positiveIndex, int positiveCount, double[] f0Step4)
        {
            f0Step3.BlockCopy(0, f0Step4, 0, f0Length);

            for (var i = positiveCount - 1; i >= 0; i--)
            {
                var limit = i == 0 ? 1 : positiveIndex[i - 1];
                for (var j = positiveIndex[i]; j > limit; j--)
                {
                    f0Step4[j - 1] = SelectBestF0(f0Step4[j], f0Step4[j + 1], f0Candidates, numberOfCandidates, j - 1);
                    if (f0Step4[j - 1] == 0)
                    {
                        break;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetNumberOfVoicedSections() counts the number of voiced sections.
        //-----------------------------------------------------------------------------
        void GetNumberOfVoicedSections(double[] f0, int f0Length, int[] positiveIndex, int[] negativeIndex, out int positiveCount, out int negativeCount)
        {
            positiveCount = 0;
            negativeCount = 0;

            for (var i = 1; i < f0Length; i++)
            {
                if (f0[i] == 0 && f0[i - 1] != 0)
                {
                    negativeIndex[negativeCount++] = i - 1;
                }
                else if (f0[i - 1] == 0 && f0[i] != 0)
                {
                    positiveIndex[positiveCount++] = i;
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetFilteredSignal() calculates the signal that is the convolution of the
        // input signal and low-pass filter.
        // This function is only used in RawEventByDio()
        //-----------------------------------------------------------------------------
        void GetFilteredSignal(int halfAverageLength, int fftSize, Complex[] ySpectrum, int yLength, double[] filteredSignal)
        {
            var lowPassFilter = new double[fftSize];
            // Nuttall window is used as a low-pass filter.
            // Cutoff frequency depends on the window length.
            Common.NuttallWindow(halfAverageLength * 4, lowPassFilter);

            var lowPassFilterSpectrum = new Complex[fftSize];
            var forwardFFT = FFTPlan.CreatePlanForDftR2C(fftSize, lowPassFilter, lowPassFilterSpectrum, FFTFlag.Estimate);
            FFT.Execute(forwardFFT);

            // Convolution
            lowPassFilterSpectrum[0] = ySpectrum[0] * lowPassFilterSpectrum[0];
            for (int i = 1, limit = fftSize / 2; i <= limit; i++)
            {
                lowPassFilterSpectrum[i] = ySpectrum[i] * lowPassFilterSpectrum[i];
                lowPassFilterSpectrum[fftSize - i - 1] = lowPassFilterSpectrum[i];
            }

            var inverseFFT = FFTPlan.CreatePlanForDftC2R(fftSize, lowPassFilterSpectrum, filteredSignal, FFTFlag.Estimate);
            FFT.Execute(inverseFFT);

            // Compensation of the delay.
            var indexBias = halfAverageLength * 2;
            for (var i = 0; i < yLength; i++)
            {
                filteredSignal[i] = filteredSignal[i + indexBias];
            }
        }

        //-----------------------------------------------------------------------------
        // GetF0CandidateContour() calculates the F0 candidates based on the
        // zero-crossings.
        //-----------------------------------------------------------------------------
        void GetF0CandidateContour(ZeroCrossings zeroCrossings, double boundaryF0, double[] temporalPositions, int f0Length, double[] f0Candidate, double[] f0Score)
        {
            if (0 == CheckEvent(zeroCrossings.NumberOfNegatives - 2) *
                CheckEvent(zeroCrossings.NumberOfPositives - 2) *
                CheckEvent(zeroCrossings.NumberOfPeaks - 2) *
                CheckEvent(zeroCrossings.NumberOfDips - 2))
            {
                f0Score.Fill(MaximumValue, 0, f0Length);
                f0Candidate.Fill(0.0, 0, f0Length);
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

            GetF0CandidateContourSub(interpolatedF0Set, f0Length, boundaryF0, f0Candidate, f0Score);
        }

        //-----------------------------------------------------------------------------
        // GetF0CandidateContourSub() calculates the f0 candidates and deviations.
        // This is the sub-function of GetF0Candidates() and assumes the calculation.
        //-----------------------------------------------------------------------------
        void GetF0CandidateContourSub(double[][] interpolatedF0Set, int f0Length, double boundaryF0, double[] f0Candidate, double[] f0Score)
        {
            for (var i = 0; i < f0Length; i++)
            {
                f0Candidate[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] + interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) / 4.0;

                f0Score[i] = Math.Sqrt((
                    (interpolatedF0Set[0][i] - f0Candidate[i]) *
                    (interpolatedF0Set[0][i] - f0Candidate[i]) +
                    (interpolatedF0Set[1][i] - f0Candidate[i]) *
                    (interpolatedF0Set[1][i] - f0Candidate[i]) +
                    (interpolatedF0Set[2][i] - f0Candidate[i]) *
                    (interpolatedF0Set[2][i] - f0Candidate[i]) +
                    (interpolatedF0Set[3][i] - f0Candidate[i]) *
                    (interpolatedF0Set[3][i] - f0Candidate[i])
                ) / 3.0);

                if (f0Candidate[i] > boundaryF0 || f0Candidate[i] < boundaryF0 / 2.0 || f0Candidate[i] > F0Ceil || f0Candidate[i] < F0Floor)
                {
                    f0Candidate[i] = 0.0;
                    f0Score[i] = MaximumValue;
                }
            }
        }

        //-----------------------------------------------------------------------------
        // SelectOneF0() corrects the f0[current_index] based on
        // f0[current_index + sign].
        //-----------------------------------------------------------------------------
        double SelectBestF0(double currentF0, double pastF0, double[][] f0Candidates, int numberOfCandidates, int targetIndex)
        {
            var referenceF0 = (currentF0 * 3.0 - pastF0) / 2.0;
            var minimumError = Math.Abs(referenceF0 - f0Candidates[0][targetIndex]);
            var bestF0 = f0Candidates[0][targetIndex];

            for (var i = 1; i < numberOfCandidates; ++i)
            {
                var currentError = Math.Abs(referenceF0 - f0Candidates[i][targetIndex]);
                if (currentError < minimumError)
                {
                    minimumError = currentError;
                    bestF0 = f0Candidates[i][targetIndex];
                }
            }
            if (Math.Abs(1.0 - bestF0 / referenceF0) > AllowRange)
            {
                return 0.0;
            }

            return bestF0;
        }

        //-----------------------------------------------------------------------------
        // CheckEvent() returns 1, provided that the input value is over 1.
        // This function is for RawEventByDio().
        //-----------------------------------------------------------------------------
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static int CheckEvent(int x)
        {
            return x > 0 ? 1 : 0;
        }

        public static int GetSamplesForDIO(int fs, int xLength, double framePeriod)
        {
            return (int)(1000.0 * xLength / fs / framePeriod) + 1;
        }
    }
}

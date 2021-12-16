using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    class ZeroCrossings
    {
        public int Length { get; }
        public double[] NegativeIntervalLocations { get; }
        public double[] NegativeIntervals { get; }
        public int NumberOfNegatives { get; }
        public double[] PositiveIntervalLocations { get; }
        public double[] PositiveIntervals { get; }
        public int NumberOfPositives { get; }
        public double[] PeakIntervalLocations { get; }
        public double[] PeakIntervals { get; }
        public int NumberOfPeaks { get; }
        public double[] DipIntervalLocations { get; }
        public double[] DipIntervals { get; }
        public int NumberOfDips { get; }

        public ZeroCrossings(
            int length,
            double[] negativeIntervalLocations,
            double[] negativeIntervals,
            int numberOfNegatives,
            double[] positiveIntervalLocations,
            double[] positiveIntervals,
            int numberOfPositives,
            double[] peakIntervalLocations,
            double[] peakIntervals,
            int numberOfPeaks,
            double[] dipIntervalLocations,
            double[] dipIntervals,
            int numberOfDips
            )
        {
            Length = length;
            NegativeIntervalLocations = negativeIntervalLocations;
            NegativeIntervals = negativeIntervals;
            NumberOfNegatives = numberOfNegatives;
            PositiveIntervalLocations = positiveIntervalLocations;
            PositiveIntervals = positiveIntervals;
            NumberOfPositives = numberOfPositives;
            PeakIntervalLocations = peakIntervalLocations;
            PeakIntervals = peakIntervals;
            NumberOfPeaks = numberOfPeaks;
            DipIntervalLocations = dipIntervalLocations;
            DipIntervals = dipIntervals;
            NumberOfDips = numberOfDips;
        }

        public void Release()
        {
            var pool = ArrayPool<double>.Shared;
            pool.Return(NegativeIntervalLocations);
            pool.Return(NegativeIntervals);
            pool.Return(PositiveIntervalLocations);
            pool.Return(PositiveIntervals);
            pool.Return(PeakIntervalLocations);
            pool.Return(PeakIntervals);
            pool.Return(DipIntervalLocations);
            pool.Return(DipIntervals);
        }

        //-----------------------------------------------------------------------------
        // GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
        // (1) Zero-crossing going from negative to positive.
        // (2) Zero-crossing going from positive to negative.
        // (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
        // the differential of waveform.
        //-----------------------------------------------------------------------------
        public static ZeroCrossings GetFourZeroCrossingIntervals(double[] filteredSignal, int yLength, double actualFs)
        {
            var pool = ArrayPool<double>.Shared;
            var negativeIntervalLocations = pool.Rent(yLength);
            var negativeIntervals = pool.Rent(yLength);
            var positiveIntervalLocations = pool.Rent(yLength);
            var positiveIntervals = pool.Rent(yLength);
            var peakIntervals = pool.Rent(yLength);
            var peakIntervalLocations = pool.Rent(yLength);
            var dipIntervalLocations = pool.Rent(yLength);
            var dipIntervals = pool.Rent(yLength);

            var numberOfNegatives = ZeroCrossingEngine(filteredSignal, yLength, actualFs, negativeIntervalLocations, negativeIntervals);

            for (var i = 0; i < yLength; i++)
            {
                filteredSignal[i] = -filteredSignal[i];
            }
            var numberOfPositives = ZeroCrossingEngine(filteredSignal, yLength, actualFs, positiveIntervalLocations, positiveIntervals);

            for (var i = 0; i < yLength; i++)
            {
                filteredSignal[i] = filteredSignal[i] - filteredSignal[i + 1];
            }
            var numberOfPeaks = ZeroCrossingEngine(filteredSignal, yLength - 1, actualFs, peakIntervalLocations, peakIntervals);

            for (var i = 0; i < yLength; i++)
            {
                filteredSignal[i] = -filteredSignal[i];
            }
            var numberOfDips = ZeroCrossingEngine(filteredSignal, yLength - 1, actualFs, dipIntervalLocations, dipIntervals);

            return new ZeroCrossings(
                yLength,
                negativeIntervalLocations,
                negativeIntervals,
                numberOfNegatives,
                positiveIntervalLocations,
                positiveIntervals,
                numberOfPositives,
                peakIntervalLocations,
                peakIntervals,
                numberOfPeaks,
                dipIntervalLocations,
                dipIntervals,
                numberOfDips
            );
        }

        //-----------------------------------------------------------------------------
        // ZeroCrossingEngine() calculates the zero crossing points from positive to
        // negative.
        //-----------------------------------------------------------------------------
        static int ZeroCrossingEngine(double[] filteredSignal, int yLength, double fs, double[] intervalLocations, double[] intervals)
        {
            var pool = ArrayPool<int>.Shared;
            var negativeGoingPoints = pool.Rent(yLength);
            for (int i = 0, limit = yLength - 1; i < limit; i++)
            {
                negativeGoingPoints[i] = 0.0 < filteredSignal[i] && filteredSignal[i + 1] <= 0.0 ? i + 1 : 0;
            }
            negativeGoingPoints[yLength - 1] = 0;

            var edges = pool.Rent(yLength);
            var count = 0;
            for (var i = 0; i < yLength; i++)
            {
                if (negativeGoingPoints[i] > 0)
                {
                    edges[count] = negativeGoingPoints[i];
                    count++;
                }
            }
            if (count < 2)
            {
                return 0;
            }

            var prevFineEdge = edges[0] - filteredSignal[edges[0] - 1] / (filteredSignal[edges[0]] - filteredSignal[edges[0] - 1]);
            for (int i = 0, limit = count - 1; i < limit; ++i)
            {
                var fineEdge = edges[i + 1] - filteredSignal[edges[i + 1] - 1] / (filteredSignal[edges[i + 1]] - filteredSignal[edges[i + 1] - 1]);
                intervals[i] = fs / (fineEdge - prevFineEdge);
                intervalLocations[i] = (prevFineEdge + fineEdge) / 2.0 / fs;

                prevFineEdge = fineEdge;
            }

            pool.Return(negativeGoingPoints);
            pool.Return(edges);

            return count - 1;
        }
    }
}

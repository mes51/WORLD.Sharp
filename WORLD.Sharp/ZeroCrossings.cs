using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    class ZeroCrossings
    {
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

        //-----------------------------------------------------------------------------
        // GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
        // (1) Zero-crossing going from negative to positive.
        // (2) Zero-crossing going from positive to negative.
        // (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
        // the differential of waveform.
        //-----------------------------------------------------------------------------
        public static ZeroCrossings GetFourZeroCrossingIntervals(double[] filteredSignal, int yLength, double actualFs)
        {
            var negativeIntervalLocations = new double[yLength];
            var negativeIntervals = new double[yLength];
            var positiveIntervalLocations = new double[yLength];
            var positiveIntervals = new double[yLength];
            var peakIntervals = new double[yLength];
            var peakIntervalLocations = new double[yLength];
            var dipIntervalLocations = new double[yLength];
            var dipIntervals = new double[yLength];

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
            var negativeGoingPoints = new int[yLength];
            for (int i = 0, limit = yLength - 1; i < limit; i++)
            {
                negativeGoingPoints[i] = 0.0 < filteredSignal[i] && filteredSignal[i + 1] <= 0.0 ? i + 1 : 0;
            }
            negativeGoingPoints[yLength - 1] = 0;

            var edges = new int[yLength];
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

            var fineEdges = new double[count];
            for (var i = 0; i < count; i++)
            {
                fineEdges[i] = edges[i] - filteredSignal[edges[i] - 1] / (filteredSignal[edges[i]] - filteredSignal[edges[i] - 1]);
            }
            for (int i = 0, limit = count - 1; i < limit; ++i)
            {
                intervals[i] = fs / (fineEdges[i + 1] - fineEdges[i]);
                intervalLocations[i] = (fineEdges[i] + fineEdges[i + 1]) / 2.0 / fs;
            }

            return count - 1;
        }
    }
}

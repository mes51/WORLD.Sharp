using System;
using System.Collections.Generic;
using System.Text;

namespace WORLD.Sharp
{
    class MVN
    {
        readonly static int[] LengthList = new int[] { 152, 400 };

        public int BaseFs { get; }

        int CurrentPosition { get; set; } = 152;

        int ChunkLength { get; set; } = 152;

        int Count { get; set; }

        ulong GRandiMvn { get; set; } = 88172645463325252UL;

        //

        int[] Chunk { get; set; } = new int[200];

        int[] Amplitude { get; set; } = new int[200];

        public MVN(int baseFs)
        {
            BaseFs = baseFs;

            LengthList[0] = unchecked((int)(152.0 * baseFs / 384000.0)) * 8;
            LengthList[1] = unchecked((int)(400.0 * baseFs / 384000.0)) * 8;
            CurrentPosition = LengthList[0];
            ChunkLength = LengthList[0];
        }

        public int GenerateMVN()
        {
            if (CurrentPosition == ChunkLength)
            {
                UpdateChunk();
            }

            return Chunk[Count] == CurrentPosition++ ? Amplitude[Count++] : 0;
        }

        int RandI(int imax)
        {
            GRandiMvn = GRandiMvn ^ (GRandiMvn << 7);
            GRandiMvn = GRandiMvn ^ (GRandiMvn >> 9);
            unchecked
            {
                return (int)(GRandiMvn % (ulong)imax);
            };
        }

        void UpdateChunk()
        {
            // Determination of the ChunkLength.
            // td = 4;
            ChunkLength = LengthList[RandI(2)];
            var numberOfPulses = ChunkLength / 4;

            // Amplitude of pulse is 2 or -2.
            Amplitude.Fill(2, 0, numberOfPulses / 2);
            Amplitude.Fill(-2, numberOfPulses / 2, numberOfPulses - numberOfPulses / 2);

            // Randomization of the amplitudes
            for (var i = 0; i < numberOfPulses; i++)
            {
                var tmpIndex = RandI(numberOfPulses);
                var tmp = Amplitude[tmpIndex];
                Amplitude[tmpIndex] = Amplitude[i];
                Amplitude[i] = tmp;
            }

            // Temporal positions with an amplitude other than zero are determined.
            for (var i = 0; i < numberOfPulses; i++)
            {
                Chunk[i] = i * 4 + RandI(4);
            }

            CurrentPosition = 0;
            Count = 0;
        }
    }
}

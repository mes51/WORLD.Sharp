using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    public class Synthesis
    {
        const double DefaultF0 = 500.0;
        const double SafeGuardMinimum = 0.000000000001;
        const double PI2 = 2.0 * Math.PI;

        int Fs { get; }

        private MVN Rand { get; }

        public Synthesis(int fs)
        {
            Fs = fs;
            Rand = new MVN(fs);
        }

        public void Synthesize(double[] f0, int f0Length, double[][] spectrogram, double[][] aperiodicity, int fftSize, double framePeriod, double[] y)
        {
            var minimumPhase = MinimumPhaseAnalysis.Create(fftSize);
            var inverseRealFFT = InverseRealFFT.Create(fftSize);
            var forwardRealFFT = ForwardRealFFT.Create(fftSize);

            var pulseLocations = new double[y.Length];
            var pulseLocationsIndex = new int[y.Length];
            var pulseLocationsTimeShift = new double[y.Length];
            var interpolatedVUV = new double[y.Length];
            var numberOfPulses = GetTimeBase(f0, f0Length, Fs, framePeriod / 1000.0, y.Length, Fs / fftSize + 1.0, pulseLocations, pulseLocationsIndex, pulseLocationsTimeShift, interpolatedVUV);

            var dcRemover = GetDCRemover(fftSize);

            framePeriod /= 1000.0;

            var impulseResponse = new double[fftSize];
            for (var i = 0; i < numberOfPulses; i++)
            {
                var noiseSize = pulseLocationsIndex[Math.Min(numberOfPulses - 1, i + 1)] - pulseLocationsIndex[i];

                GetOneFrameSegment(interpolatedVUV[pulseLocationsIndex[i]], noiseSize, spectrogram, fftSize, aperiodicity, f0Length, framePeriod, pulseLocations[i], pulseLocationsTimeShift[i], Fs, forwardRealFFT, inverseRealFFT, minimumPhase, dcRemover, impulseResponse);

                var offset = pulseLocationsIndex[i] - fftSize / 2 + 1;

                for (int j = Math.Max(0, -offset), limit = Math.Min(fftSize, y.Length - offset); j < limit; j++)
                {
                    var index = j + offset;
                    y[index] += impulseResponse[j];
                }
            }

            inverseRealFFT.Release();
            forwardRealFFT.Release();
        }

        //-----------------------------------------------------------------------------
        // GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
        //-----------------------------------------------------------------------------
        void GetOneFrameSegment(double currentVUV, int noiseSize, double[][] spectrogram, int fftSize, double[][] aperiodicity, int f0Length, double framePeriod, double currentTime, double fractionalTimeShift, int fs, ForwardRealFFT forwardRealFFT, InverseRealFFT inverseRealFFT, MinimumPhaseAnalysis minimumPhase, double[] dcRemover, double[] response)
        {
            var aperiodicResponse = new double[fftSize];
            var periodicResponse = new double[fftSize];

            var spectralEnvelope = new double[fftSize];
            var aperiodicRatio = new double[fftSize];
            GetSpectralEnvelope(currentTime, framePeriod, f0Length, spectrogram, fftSize, spectralEnvelope);
            GetAperiodicRatio(currentTime, framePeriod, f0Length, aperiodicity, fftSize, aperiodicRatio);

            // Synthesis of the periodic response
            GetPeriodicResponse(fftSize, spectralEnvelope, aperiodicRatio, currentVUV, inverseRealFFT, minimumPhase, dcRemover, fractionalTimeShift, fs, periodicResponse);

            // Synthesis of the aperiodic response
            GetAperiodicResponse(noiseSize, fftSize, spectralEnvelope, aperiodicRatio, currentVUV, forwardRealFFT, inverseRealFFT, minimumPhase, aperiodicResponse);

            var sqrtNoiseSize = Math.Sqrt(noiseSize);
            for (var i = 0; i < fftSize; i++)
            {
                response[i] = (periodicResponse[i] * sqrtNoiseSize + aperiodicResponse[i]) / fftSize;
            }
        }

        //-----------------------------------------------------------------------------
        // GetAperiodicResponse() calculates an aperiodic response.
        //-----------------------------------------------------------------------------
        void GetAperiodicResponse(int noiseSize, int fftSize, double[] spectrum, double[] aperiodicRatio, double currentVUV, ForwardRealFFT forwardRealFFT, InverseRealFFT inverseRealFFT, MinimumPhaseAnalysis minimumPhase, double[] aperiodicResponse)
        {
            GetNoiseSpectrum(noiseSize, fftSize, forwardRealFFT);

            var logSpectrum = minimumPhase.LogSpectrum;
            if (currentVUV != 0.0)
            {
                for (int i = 0, limit = minimumPhase.FFTSize / 2; i <= limit; i++)
                {
                    logSpectrum[i] = Math.Log(spectrum[i] * aperiodicRatio[i]) / 2.0;
                }
            }
            else
            {
                for (int i = 0, limit = minimumPhase.FFTSize / 2; i <= limit; i++)
                {
                    logSpectrum[i] = Math.Log(spectrum[i]) / 2.0;
                }
            }
            Common.GetMinimumPhaseSpectrum(minimumPhase);

            var minimumPhaseSpectrum = minimumPhase.MinimumPhaseSpectrum;
            Array.Copy(minimumPhaseSpectrum, 0, inverseRealFFT.Spectrum, 0, fftSize / 2 + 1);
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                var real = minimumPhaseSpectrum[i].Real * forwardRealFFT.Spectrum[i].Real - minimumPhaseSpectrum[i].Imaginary * forwardRealFFT.Spectrum[i].Imaginary;
                var imaginary = minimumPhaseSpectrum[i].Real * forwardRealFFT.Spectrum[i].Imaginary + minimumPhaseSpectrum[i].Imaginary * forwardRealFFT.Spectrum[i].Real;

                inverseRealFFT.Spectrum[i] = new Complex(real, imaginary);
            }
            FFT.Execute(inverseRealFFT.InverseFFT);
            MatlabFunctions.FFTShift(inverseRealFFT.Waveform.AsSpan(0, fftSize), aperiodicResponse);
        }

        void GetNoiseSpectrum(int noiseSize, int fftSize, ForwardRealFFT forwardRealFFT)
        {
            var waveform = forwardRealFFT.Waveform;
            for (var i = 0; i < noiseSize; i++)
            {
                waveform[i] = Rand.GenerateMVN();
            }

            var average = noiseSize > 0 ? waveform.Take(noiseSize).Average() : 1.0;
            for (var i = 0; i < noiseSize; i++)
            {
                waveform[i] -= average;
            }
            Array.Clear(waveform, noiseSize, fftSize - noiseSize);
            FFT.Execute(forwardRealFFT.ForwardFFT);
        }

        //-----------------------------------------------------------------------------
        // GetSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
        // the fractional time shift under 1/fs.
        //-----------------------------------------------------------------------------
        void GetSpectrumWithFractionalTimeShift(int fftSize, double coefficient, InverseRealFFT inverseRealFFT)
        {
            for (var i = fftSize / 2; i > -1; i--)
            {
                var complex = inverseRealFFT.Spectrum[i];

                var re2 = Math.Cos(coefficient * i);
                var im2 = Math.Sqrt(1.0 - re2 * re2);

                inverseRealFFT.Spectrum[i] = new Complex(complex.Real * re2 + complex.Imaginary * im2, complex.Imaginary * re2 - complex.Real * im2);
            }
        }

        //-----------------------------------------------------------------------------
        // GetPeriodicResponse() calculates an aperiodic response.
        //-----------------------------------------------------------------------------
        void GetPeriodicResponse(int fftSize, double[] spectrum, double[] aperiodicRatio, double currentVUV, InverseRealFFT inverseRealFFT, MinimumPhaseAnalysis minimumPhase, double[] dcRemover, double fractionalTimeShift, int fs, double[] periodicResponse)
        {
            if (currentVUV <= 0.5 || aperiodicRatio[0] > 0.999)
            {
                Array.Clear(periodicResponse, 0, fftSize);
                return;
            }

            var logSpectrum = minimumPhase.LogSpectrum;
            for (int i = 0, limit = minimumPhase.FFTSize / 2; i <= limit; i++)
            {
                logSpectrum[i] = Math.Log(spectrum[i] * (1.0 - aperiodicRatio[i]) + SafeGuardMinimum) / 2.0;
            }
            Common.GetMinimumPhaseSpectrum(minimumPhase);
            Array.Copy(minimumPhase.MinimumPhaseSpectrum, 0, inverseRealFFT.Spectrum, 0, fftSize / 2 + 1);

            double coefficient = PI2 * fractionalTimeShift * fs / fftSize;
            GetSpectrumWithFractionalTimeShift(fftSize, coefficient, inverseRealFFT);

            FFT.Execute(inverseRealFFT.InverseFFT);
            MatlabFunctions.FFTShift(inverseRealFFT.Waveform.AsSpan(0, fftSize), periodicResponse);
            RemoveDCComponent(periodicResponse, fftSize, dcRemover, periodicResponse);
        }

        //-----------------------------------------------------------------------------
        // RemoveDCComponent()
        //-----------------------------------------------------------------------------
        void RemoveDCComponent(double[] periodicResponse, int fftSize, double[] dcRemover, double[] newPeriodicResponse)
        {
            var dcComponent = periodicResponse.Take(fftSize).Skip(fftSize / 2).Sum();
            for (int i = 0, limit = fftSize / 2; i < limit; i++)
            {
                newPeriodicResponse[i] = -dcComponent * dcRemover[i];
            }
            for (var i = fftSize / 2; i < fftSize; i++)
            {
                newPeriodicResponse[i] -= dcComponent * dcRemover[i];
            }
        }

        void GetAperiodicRatio(double currentTime, double framePeriod, int f0Length, double[][] aperiodicity, int fftSize, double[] aperiodicSpectrum)
        {
            var currentFrameFloor = Math.Min(f0Length - 1, (int)Math.Floor(currentTime / framePeriod));
            var currentFrameCeil = Math.Min(f0Length - 1, (int)Math.Ceiling(currentTime / framePeriod));
            var interpolation = currentTime / framePeriod - currentFrameFloor;

            if (currentFrameFloor == currentFrameCeil)
            {
                for (int i = 0, limit = fftSize / 2; i <= limit; i++)
                {
                    aperiodicSpectrum[i] = Math.Pow(Common.GetSafeAperiodicity(aperiodicity[currentFrameFloor][i]), 2.0);
                }
            }
            else
            {
                for (int i = 0, limit = fftSize / 2; i <= limit; i++)
                {
                    aperiodicSpectrum[i] = Math.Pow((1.0 - interpolation) * Common.GetSafeAperiodicity(aperiodicity[currentFrameFloor][i]) + interpolation * Common.GetSafeAperiodicity(aperiodicity[currentFrameCeil][i]), 2.0);
                }
            }
        }

        void GetSpectralEnvelope(double currentTime, double framePeriod, int f0Length, double[][] spectrogram, int fftSize, double[] spectralEnvelope)
        {
            int currentFrameFloor = Math.Min(f0Length - 1, (int)Math.Floor(currentTime / framePeriod));
            int currentFrameCeil = Math.Min(f0Length - 1, (int)Math.Ceiling(currentTime / framePeriod));
            double interpolation = currentTime / framePeriod - currentFrameFloor;

            if (currentFrameFloor == currentFrameCeil)
            {
                for (int i = 0, limit = fftSize / 2; i <= limit; i++)
                {
                    spectralEnvelope[i] = Math.Abs(spectrogram[currentFrameFloor][i]);
                }
            }
            else
            {
                for (int i = 0, limit = fftSize / 2; i <= limit; i++)
                {
                    spectralEnvelope[i] = (1.0 - interpolation) * Math.Abs(spectrogram[currentFrameFloor][i]) + interpolation * Math.Abs(spectrogram[currentFrameCeil][i]);
                }
            }
        }

        double[] GetDCRemover(int fftSize)
        {
            var dcRemover = new double[fftSize];
            var dcComponent = 0.0;
            for (int i = 0, limit = fftSize / 2; i < limit; i++)
            {
                dcRemover[i] = 0.5 - 0.5 * Math.Cos(PI2 * (i + 1.0) / (1.0 + fftSize));
                dcRemover[fftSize - i - 1] = dcRemover[i];
                dcComponent += dcRemover[i] * 2.0;
            }
            for (int i = 0, limit = fftSize / 2; i < limit; i++)
            {
                dcRemover[i] /= dcComponent;
                dcRemover[fftSize - i - 1] = dcRemover[i];
            }

            return dcRemover;
        }

        int GetTimeBase(double[] f0, int f0Length, int fs, double framePeriod, int yLength, double lowestF0, double[] pulseLocations, int[] pulseLocationsIndex, double[] pulseLocationsTimeShift, double[] interpolatedVUV)
        {
            var timeAxis = new double[yLength];
            var coarseTimeAxis = new double[f0Length + 1];
            var coarseF0 = new double[f0Length + 1];
            var coarseVUV = new double[f0Length + 1];
            GetTemporalParametersForTimeBase(f0, f0Length, fs, yLength, framePeriod, lowestF0, timeAxis, coarseTimeAxis, coarseF0, coarseVUV);

            var interpolatedF0 = new double[yLength];
            MatlabFunctions.Interp1(coarseTimeAxis, coarseF0, timeAxis, interpolatedF0);
            MatlabFunctions.Interp1(coarseTimeAxis, coarseVUV, timeAxis, interpolatedVUV);
            for (var i = 0; i < yLength; i++)
            {
                interpolatedVUV[i] = interpolatedVUV[i] > 0.5 ? 1.0 : 0.0;
                interpolatedF0[i] = interpolatedVUV[i] == 0.0 ? DefaultF0 : interpolatedF0[i];
            }

            return GetPulseLocationsForTimeBase(interpolatedF0, timeAxis, yLength, fs, pulseLocations, pulseLocationsIndex, pulseLocationsTimeShift);
        }

        int GetPulseLocationsForTimeBase(double[] interpolatedF0, double[] timeAxis, int yLength, int fs, double[] pulseLocations, int[] pulseLocationsIndex, double[] pulseLocationsTimeShift)
        {
            var totalPhase = new double[yLength];
            var warpPhase = new double[yLength];
            var warpPhaseAbs = new double[yLength];
            totalPhase[0] = PI2 * interpolatedF0[0] / fs;
            warpPhase[0] = totalPhase[0] % PI2;
            for (var i = 1; i < yLength; i++)
            {
                totalPhase[i] = totalPhase[i - 1] + PI2 * interpolatedF0[i] / fs;
                warpPhase[i] = totalPhase[i] % PI2;
                warpPhaseAbs[i - 1] = Math.Abs(warpPhase[i] - warpPhase[i - 1]);
            }

            var numberOfPulses = 0;
            for (int i = 0, limit = yLength - 1; i < limit; i++)
            {
                if (warpPhaseAbs[i] > Math.PI)
                {
                    pulseLocations[numberOfPulses] = timeAxis[i];
                    pulseLocationsIndex[numberOfPulses] = i;

                    // calculate the time shift in seconds between exact fractional pulse
                    // position and the integer pulse position (sample i)
                    // as we don't have access to the exact pulse position, we infer it
                    // from the point between sample i and sample i + 1 where the
                    // accummulated phase cross a multiple of 2pi
                    // this point is found by solving y1 + x * (y2 - y1) = 0 for x, where y1
                    // and y2 are the phases corresponding to sample i and i + 1, offset so
                    // they cross zero; x >= 0
                    var y1 = warpPhase[i] - PI2;
                    var y2 = warpPhase[i + 1];
                    var x = -y1 / (y2 - y1);
                    pulseLocationsTimeShift[numberOfPulses] = x / fs;

                    numberOfPulses++;
                }
            }

            return numberOfPulses;
        }

        void GetTemporalParametersForTimeBase(double[] f0, int f0Length, int fs, int yLength, double framePeriod, double lowestF0, double[] timeAxis, double[] coarseTimeAxis, double[] coarseF0, double[] coarseVUV)
        {
            for (var i = 0; i < yLength; i++)
            {
                timeAxis[i] = i / (double)fs;
            }

            // the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
            for (var i = 0; i < f0Length; i++)
            {
                coarseTimeAxis[i] = i * framePeriod;
                coarseF0[i] = f0[i] < lowestF0 ? 0.0 : f0[i];
                coarseVUV[i] = coarseF0[i] == 0.0 ? 0.0 : 1.0;
            }
            coarseTimeAxis[f0Length] = f0Length * framePeriod;
            coarseF0[f0Length] = coarseF0[f0Length - 1] * 2.0 - coarseF0[f0Length - 2];
            coarseVUV[f0Length] = coarseVUV[f0Length - 1] * 2.0 - coarseVUV[f0Length - 2];
        }
    }

    public class SynthesisRealTime
    {
        const double DefaultF0 = 500.0;
        const double PI2 = 2.0 * Math.PI;

        class RingBuffer
        {
            public enum SearchType
            {
                Spectrogram,
                Aperiodicity
            }

            public RingBuffer(int capacity)
            {
                Capacity = capacity;
                F0Length = new int[capacity];
                Spectrogram = new double[capacity][][];
                Aperiodicity = new double[capacity][][];
                InterpolatedVUV = new double[capacity][];
                PulseLocations = new double[capacity][];
                PulseLocationIndex = new int[capacity][];
                NumberOfPulses = new int[capacity];
                F0Origin = new int[capacity];
            }

            public int CurrentPointer { get; set; }

            public int CurrentPointer2 { get; set; }

            public int HeadPointer { get; set; }

            public int HeadIndex => HeadPointer % Capacity;

            public int CurrentIndex => CurrentPointer % Capacity;

            public int CurrentIndex2 => CurrentPointer2 % Capacity;

            public int Capacity { get; }

            public double[][][] Spectrogram { get; }

            public double[][][] Aperiodicity { get; }

            public int[] F0Length { get; }

            public int[] F0Origin { get; }

            public double[][] InterpolatedVUV { get; }

            public int[] NumberOfPulses { get; }

            public double[][] PulseLocations { get; }

            public int[][] PulseLocationIndex { get; }

            public void SearchPointer(int frame, SearchType flag, out double[] front, out double[] next)
            {
                var pointer = CurrentIndex2;
                var f0Length = F0Length[pointer];
                var f0Origin = F0Origin[pointer];
                var index = -1;

                for (var i = 0; i < f0Length; i++)
                {
                    if (f0Origin + i == frame)
                    {
                        index = i;
                        break;
                    }
                }

                var tmpPointer = flag == SearchType.Spectrogram ? Spectrogram : Aperiodicity;

                front = tmpPointer[pointer][index];
                if (index == f0Length - 1)
                {
                    next = tmpPointer[(CurrentPointer2 + 1) % Capacity][0];
                }
                else
                {
                    next = tmpPointer[pointer][index + 1];
                }
            }

            public double GetCurrentVUV(int currentLocation, double framePeriod, int fs)
            {
                var pointer = CurrentPointer % Capacity;
                var startSample = Math.Max(0, (int)Math.Ceiling((F0Origin[pointer] - 1) * framePeriod * fs));

                return InterpolatedVUV[pointer][currentLocation - startSample + 1];
            }

            public int GetNextPulseLocationIndex(int index)
            {
                var pointer = CurrentIndex;
                if (index < NumberOfPulses[pointer] - 1)
                {
                    return PulseLocationIndex[pointer][index + 1];
                }
                else if (CurrentPointer == HeadPointer - 1)
                {
                    return 0;
                }

                for (var i = 1; i < Capacity; i++)
                {
                    var p = (i + CurrentPointer) % Capacity;
                    if (NumberOfPulses[p] != 0)
                    {
                        return PulseLocationIndex[p][0];
                    }
                }

                return 0;
            }

            public void Clear(int start = 0, int end = int.MaxValue)
            {
                if (end == int.MaxValue)
                {
                    end = Capacity;
                }

                for (var i = start; i < end; i++)
                {
                    var p = i % Capacity;
                    NumberOfPulses[p] = 0;
                    PulseLocations[p] = null;
                    PulseLocationIndex[p] = null;
                    InterpolatedVUV[p] = null;
                }
            }

            public bool IsFull()
            {
                return HeadPointer - CurrentPointer2 == Capacity;
            }
        }

        const double SafeGuardMinimum = 0.000000000001;

        public SynthesisRealTime(int sampleRate, double framePeriod, int fftSize, int bufferSize, int ringBufferCapacity)
        {
            SampleRate = sampleRate;
            FramePeriod = framePeriod * 0.001;
            AudioBufferSize = bufferSize;
            AudioBuffer = new double[bufferSize * 2 + fftSize];
            Buffer = new RingBuffer(ringBufferCapacity);
            FFTSize = fftSize;
            Rand = new MVN(sampleRate);

            ImpulseResponse = new double[fftSize];
            DCRemover = GetDCRemover(fftSize / 2);

            RefreshSynthesizer();

            MinimumPhase = MinimumPhaseAnalysis.Create(fftSize);
            InverseRealFFT = InverseRealFFT.Create(fftSize);
            ForwardRealFFT = ForwardRealFFT.Create(fftSize);
        }

        public double[] AudioBuffer { get; }

        public int AudioBufferSize { get; }

        public int SynthesizedSample { get; private set; }

        MVN Rand { get; }

        RingBuffer Buffer { get; }

        int SampleRate { get; }

        double FramePeriod { get; }

        int FFTSize { get; }

        double[] ImpulseResponse { get; }

        double[] DCRemover { get; }

        int CumulativeFrame { get; set; }

        int Handoff { get; set; }

        double HandoffF0 { get; set; }

        double HandoffPhase { get; set; }

        int LastLocation { get; set; }

        int I { get; set; }

        int CurrentFrame { get; set; }

        int CurrentLocation { get; set; }

        MinimumPhaseAnalysis MinimumPhase { get; }

        InverseRealFFT InverseRealFFT { get; }

        ForwardRealFFT ForwardRealFFT { get; }

        public void RefreshSynthesizer()
        {
            Buffer.HeadPointer = 0;
            Buffer.CurrentPointer = 0;
            Buffer.CurrentPointer2 = 0;
            Buffer.Clear();
            HandoffPhase = 0.0;
            HandoffF0 = 0.0;
            Handoff = 0;
            CumulativeFrame = -1;
            LastLocation = 0;
            I = 0;
            CurrentFrame = 0;
            SynthesizedSample = 0;
            AudioBuffer.Fill(0.0);
        }

        public bool IsLocked()
        {
            return Buffer.IsFull() && (SynthesizedSample + AudioBufferSize >= LastLocation);
        }

        public bool AddParameters(double[] f0, int f0Length, double[][] spectrogram, double[][] aperiodicity)
        {
            if (Buffer.IsFull())
            {
                return false;
            }

            var pointer = Buffer.HeadIndex;
            Buffer.F0Length[pointer] = f0Length;
            Buffer.F0Origin[pointer] = CumulativeFrame + 1;
            CumulativeFrame += f0Length;

            Buffer.Spectrogram[pointer] = spectrogram;
            Buffer.Aperiodicity[pointer] = aperiodicity;
            if (CumulativeFrame < 1)
            {
                HandoffF0 = f0[f0Length - 1];
                Buffer.NumberOfPulses[pointer] = 0;
                Buffer.HeadPointer++;
                Handoff = 1;
                return true;
            }

            var startSample = Math.Max(0, (int)Math.Ceiling((CumulativeFrame - f0Length) * FramePeriod * SampleRate));
            var endSample = (int)Math.Ceiling(CumulativeFrame * FramePeriod * SampleRate);
            var numberOfSamples = endSample - startSample;

            Buffer.InterpolatedVUV[pointer] = new double[numberOfSamples + 1];
            Buffer.PulseLocations[pointer] = new double[numberOfSamples];
            Buffer.PulseLocationIndex[pointer] = new int[numberOfSamples];

            GetTimeBase(f0, f0Length, startSample, numberOfSamples);

            HandoffF0 = f0[f0Length - 1];
            Buffer.HeadPointer++;
            Handoff = 1;

            return true;
        }

        public bool Synthesize()
        {
            if (!CheckSynthesizer())
            {
                return false;
            }

            for (int i = 0, limit = AudioBufferSize + FFTSize; i < limit; i++)
            {
                AudioBuffer[i] = AudioBuffer[i + AudioBufferSize];
            }

            var pointer = Buffer.CurrentIndex;
            var currentLocation = Buffer.PulseLocationIndex[pointer][I];
            while (currentLocation < SynthesizedSample + AudioBufferSize)
            {
                var tmp = Buffer.GetNextPulseLocationIndex(I);
                var noiseSize = tmp - currentLocation;

                GetOneFrameSegment(noiseSize, currentLocation);
                var offset = currentLocation - SynthesizedSample - FFTSize / 2 + 1;
                for (var i = Math.Max(0, -offset); i < FFTSize; i++)
                {
                    var index = i + offset;
                    AudioBuffer[index] += ImpulseResponse[i];
                }
                currentLocation = tmp;
                UpdateSynthesizer(currentLocation);
            }

            SynthesizedSample += AudioBufferSize;
            SeekSynthesizer(SynthesizedSample);

            return true;
        }

        void GetNoiseSpectrum(int noiseSize, int fftSize, ForwardRealFFT forwardRealFFT)
        {
            var waveform = forwardRealFFT.Waveform;
            for (var i = 0; i < noiseSize; i++)
            {
                waveform[i] = Rand.GenerateMVN();
            }

            var average = noiseSize > 0 ? waveform.Take(noiseSize).Average() : 1.0;
            for (var i = 0; i < noiseSize; i++)
            {
                waveform[i] -= average;
            }
            Array.Clear(waveform, noiseSize, fftSize - noiseSize);
            FFT.Execute(forwardRealFFT.ForwardFFT);
        }

        //-----------------------------------------------------------------------------
        // GetAperiodicResponse() calculates an aperiodic response.
        //-----------------------------------------------------------------------------
        void GetAperiodicResponse(int noiseSize, int fftSize, double[] spectrum, double[] aperiodicRatio, double currentVUV, double[] aperiodicResponse)
        {
            GetNoiseSpectrum(noiseSize, fftSize, ForwardRealFFT);

            var logSpectrum = MinimumPhase.LogSpectrum;
            if (currentVUV != 0.0)
            {
                for (int i = 0, limit = MinimumPhase.FFTSize / 2; i <= limit; i++)
                {
                    logSpectrum[i] = Math.Log(spectrum[i] * aperiodicRatio[i] + SafeGuardMinimum) / 2.0;
                }
            }
            else
            {
                for (int i = 0, limit = MinimumPhase.FFTSize / 2; i <= limit; i++)
                {
                    logSpectrum[i] = Math.Log(spectrum[i]) / 2.0;
                }
            }
            Common.GetMinimumPhaseSpectrum(MinimumPhase);

            var minimumPhaseSpectrum = MinimumPhase.MinimumPhaseSpectrum;
            for (int i = 0, limit = fftSize / 2; i <= limit; i++)
            {
                var real = minimumPhaseSpectrum[i].Real * ForwardRealFFT.Spectrum[i].Real - minimumPhaseSpectrum[i].Imaginary * ForwardRealFFT.Spectrum[i].Imaginary;
                var imaginary = minimumPhaseSpectrum[i].Real * ForwardRealFFT.Spectrum[i].Imaginary + minimumPhaseSpectrum[i].Imaginary * ForwardRealFFT.Spectrum[i].Real;

                InverseRealFFT.Spectrum[i] = new Complex(real, imaginary);
            }
            FFT.Execute(InverseRealFFT.InverseFFT);
            MatlabFunctions.FFTShift(InverseRealFFT.Waveform.AsSpan(0, fftSize), aperiodicResponse);
        }

        //-----------------------------------------------------------------------------
        // RemoveDCComponent()
        //-----------------------------------------------------------------------------
        void RemoveDCComponent(double[] periodicResponse, int fftSize, double[] dcRemover, double[] newPeriodicResponse)
        {
            var dcComponent = periodicResponse.Take(fftSize).Skip(fftSize / 2).Sum();
            newPeriodicResponse.Fill(0.0, 0, fftSize / 2);
            for (var i = fftSize / 2; i < fftSize; i++)
            {
                newPeriodicResponse[i] -= dcComponent * dcRemover[i - fftSize / 2];
            }
        }

        //-----------------------------------------------------------------------------
        // GetPeriodicResponse() calculates an aperiodic response.
        //-----------------------------------------------------------------------------
        void GetPeriodicResponse(int fftSize, double[] spectrum, double[] aperiodicRatio, double currentVUV, double[] periodicResponse)
        {
            if (currentVUV <= 0.5 || aperiodicRatio[0] > 0.999)
            {
                periodicResponse.Fill(0.0, 0, fftSize);
                return;
            }

            var logSpectrum = MinimumPhase.LogSpectrum;
            for (int i = 0, limit = MinimumPhase.FFTSize / 2; i <= limit; i++)
            {
                logSpectrum[i] = Math.Log(spectrum[i] * (1.0 - aperiodicRatio[i]) + SafeGuardMinimum) / 2.0;
            }
            Common.GetMinimumPhaseSpectrum(MinimumPhase);
            Array.Copy(MinimumPhase.MinimumPhaseSpectrum, 0, InverseRealFFT.Spectrum, 0, fftSize / 2 + 1);

            FFT.Execute(InverseRealFFT.InverseFFT);
            MatlabFunctions.FFTShift(InverseRealFFT.Waveform.AsSpan(0, fftSize), periodicResponse);
            RemoveDCComponent(periodicResponse, fftSize, DCRemover, periodicResponse);
        }

        void GetSpectralEnvelope(double currentLocation, double[] spectralEnvelope)
        {
            var currentFrameFloor = (int)(currentLocation / FramePeriod);
            var currentFrameCeil = (int)Math.Ceiling(currentLocation / FramePeriod);
            var interpolation = currentLocation / FramePeriod - currentFrameFloor;
            double[] front;
            double[] next;
            Buffer.SearchPointer(currentFrameFloor, RingBuffer.SearchType.Spectrogram, out front, out next);

            var halfFFTSize = FFTSize / 2;
            if (currentFrameFloor == currentFrameCeil)
            {
                for (var i = 0; i <= halfFFTSize; i++)
                {
                    spectralEnvelope[i] = Math.Abs(front[i]);
                }
            }
            else
            {
                for (var i = 0; i <= halfFFTSize; i++)
                {
                    spectralEnvelope[i] = (1.0 - interpolation) * Math.Abs(front[i]) + interpolation * Math.Abs(next[i]);
                }
            }
        }

        void GetAperiodicRatio(double currentLocation, double[] aperiodicSpectrum)
        {
            var currentFrameFloor = (int)(currentLocation / FramePeriod);
            var currentFrameCeil = (int)Math.Ceiling(currentLocation / FramePeriod);
            var interpolation = currentLocation / FramePeriod - currentFrameFloor;
            double[] front;
            double[] next;
            Buffer.SearchPointer(currentFrameFloor, RingBuffer.SearchType.Aperiodicity, out front, out next);

            var halfFFTSize = FFTSize / 2;
            if (currentFrameFloor == currentFrameCeil)
            {
                for (var i = 0; i <= halfFFTSize; i++)
                {
                    aperiodicSpectrum[i] = Math.Pow(Common.GetSafeAperiodicity(front[i]), 2.0);
                }
            }
            else
            {
                for (var i = 0; i <= halfFFTSize; i++)
                {
                    aperiodicSpectrum[i] = Math.Pow((1.0 - interpolation) * Common.GetSafeAperiodicity(front[i]) + interpolation * Common.GetSafeAperiodicity(next[i]), 2.0);
                }
            }
        }

        //-----------------------------------------------------------------------------
        // GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
        //-----------------------------------------------------------------------------
        void GetOneFrameSegment(int noiseSize, int currentLocation)
        {
            var aperiodicResponse = new double[FFTSize];
            var periodicResponse = new double[FFTSize];
            var spectralEnvelope = new double[FFTSize];
            var aperiodicRatio = new double[FFTSize];

            var tmpLocation = currentLocation / (double)SampleRate;
            SeekSynthesizer(tmpLocation);
            GetSpectralEnvelope(tmpLocation, spectralEnvelope);
            GetAperiodicRatio(tmpLocation, aperiodicRatio);

            var currentVUV = Buffer.GetCurrentVUV(currentLocation, FramePeriod, SampleRate);

            // Synthesis of the periodic response
            GetPeriodicResponse(FFTSize, spectralEnvelope, aperiodicRatio, currentVUV, periodicResponse);

            // Synthesis of the aperiodic response
            GetAperiodicResponse(noiseSize, FFTSize, spectralEnvelope, aperiodicRatio, currentVUV, aperiodicResponse);

            var sqrtNoiseSize = Math.Sqrt(noiseSize);
            for (var i = 0; i < FFTSize; i++)
            {
                ImpulseResponse[i] = (periodicResponse[i] * sqrtNoiseSize + aperiodicResponse[i]) / FFTSize;
            }
        }

        void GetTemporalParametersForTimeBase(double[] f0, int f0Length, double[] coarseTimeAxis, double[] coarseF0, double[] coarseVUV)
        {
            var cumulativeFrame = Math.Max(0, CumulativeFrame - f0Length);

            coarseF0[0] = HandoffF0;
            coarseTimeAxis[0] = cumulativeFrame * FramePeriod;
            coarseVUV[0] = HandoffF0 == 0.0 ? 0.0 : 1.0;
            for (int i = 0, index = Handoff; i < f0Length; i++, index++)
            {
                coarseTimeAxis[index] = (i + cumulativeFrame + Handoff) * FramePeriod;
                coarseF0[index] = f0[i];
                coarseVUV[index] = f0[i] == 0.0 ? 0.0 : 1.0;
            }
        }

        void GetPulseLocationsForTimeBase(double[] interpolatedF0, double[] timeAxis, int numberOfSamples, double origin)
        {
            var totalPhase = new double[numberOfSamples + Handoff];
            totalPhase[0] = Handoff == 1 ? HandoffPhase : PI2 * interpolatedF0[0] / SampleRate;

            totalPhase[1] = totalPhase[0] + PI2 * interpolatedF0[0] / SampleRate;
            for (int i = 1 + Handoff, limit = numberOfSamples + Handoff; i < limit; i++)
            {
                totalPhase[i] = totalPhase[i - 1] + PI2 * interpolatedF0[i - Handoff] / SampleRate;
            }
            HandoffPhase = totalPhase[numberOfSamples - 1 + Handoff];

            var wrapPhase = totalPhase.Select((t) => t % PI2).ToArray();
            var wrapPhaseABS = wrapPhase.Zip(wrapPhase.Skip(1), (w, nw) => Math.Abs(nw - w)).Append(0.0).ToArray();

            var pointer = Buffer.HeadIndex;
            var numberOfPulses = 0;
            for (int i = 0, limit = numberOfSamples - 1 + Handoff; i < limit; i++)
            {
                if (wrapPhaseABS[i] > Math.PI)
                {
                    Buffer.PulseLocations[pointer][numberOfPulses] = timeAxis[i] - (double)Handoff / SampleRate;
                    Buffer.PulseLocationIndex[pointer][numberOfPulses] = MatlabFunctions.MatlabRound(Buffer.PulseLocations[pointer][numberOfPulses] * SampleRate);
                    numberOfPulses++;
                }
            }
            Buffer.NumberOfPulses[pointer] = numberOfPulses;

            if (numberOfPulses != 0)
            {
                LastLocation = Buffer.PulseLocationIndex[pointer][numberOfPulses - 1];
            }
        }

        void GetTimeBase(double[] f0, int f0Length, int startSample, int numberOfSamples)
        {
            var coarseTimeAxis = new double[f0Length + Handoff];
            var coarseF0 = new double[f0Length + Handoff];
            var coarseVUV = new double[f0Length + Handoff];
            GetTemporalParametersForTimeBase(f0, f0Length, coarseTimeAxis, coarseF0, coarseVUV);

            var interpolatedF0 = new double[numberOfSamples];
            var timeAxis = Enumerable.Range(0, numberOfSamples).Select((i) => (i + startSample) / (double)SampleRate).ToArray();

            var pointer = Buffer.HeadIndex;
            MatlabFunctions.Interp1(coarseTimeAxis, coarseF0, timeAxis, interpolatedF0);
            MatlabFunctions.Interp1(coarseTimeAxis, coarseVUV, timeAxis, Buffer.InterpolatedVUV[pointer]);
            var interpolatedVUV = Buffer.InterpolatedVUV[pointer];
            for (var i = 0; i < numberOfSamples; i++)
            {
                interpolatedVUV[i] = interpolatedVUV[i] > 0.5 ? 1.0 : 0.0;
                interpolatedF0[i] = interpolatedVUV[i] == 0.0 ? DefaultF0 : interpolatedF0[i];
            }

            GetPulseLocationsForTimeBase(interpolatedF0, timeAxis, numberOfSamples, coarseTimeAxis[0]);

            HandoffF0 = interpolatedF0[numberOfSamples - 1];
        }

        bool UpdateSynthesizer(int currentLocation)
        {
            var pointer = Buffer.CurrentIndex;
            if (I < Buffer.NumberOfPulses[pointer] - 1)
            {
                I++;
                return true;
            }
            else if (Buffer.CurrentPointer == Buffer.HeadPointer - 1)
            {
                return false;
            }

            for (var i = 1; i < Buffer.Capacity; i++)
            {
                var p = (i + Buffer.CurrentPointer) % Buffer.Capacity;
                if (Buffer.NumberOfPulses[p] != 0)
                {
                    I = 0;
                    Buffer.CurrentPointer += i;
                    return true;
                }
            }

            return false;
        }

        void SeekSynthesizer(double currentLocation)
        {
            var frameNumber = (int)(currentLocation / FramePeriod);

            var tmpPointer = Buffer.CurrentPointer2;
            for (var i = 0; i < (Buffer.HeadPointer - Buffer.CurrentPointer2); i++)
            {
                var tmp = (tmpPointer + i) % Buffer.Capacity;
                if (Buffer.F0Origin[tmp] <= frameNumber && frameNumber < Buffer.F0Origin[tmp] + Buffer.F0Length[tmp])
                {
                    tmpPointer += i;
                    break;
                }
            }

            Buffer.Clear(Buffer.CurrentPointer2, tmpPointer);
            Buffer.CurrentPointer2 = tmpPointer;
        }

        bool CheckSynthesizer()
        {
            if (SynthesizedSample + AudioBufferSize >= LastLocation)
            {
                return false;
            }

            while (Buffer.NumberOfPulses[Buffer.CurrentIndex] == 0 && Buffer.CurrentPointer != Buffer.HeadPointer)
            {
                Buffer.CurrentPointer++;
            }

            return true;
        }

        double[] GetDCRemover(int fftSize)
        {
            var dcRemover = new double[fftSize];
            var dcComponent = 0.0;
            for (int i = 0, limit = fftSize / 2; i < limit; i++)
            {
                dcRemover[i] = 0.5 - 0.5 * Math.Cos(PI2 * (i + 1.0) / (1.0 + fftSize));
                dcRemover[fftSize - i - 1] = dcRemover[i];
                dcComponent += dcRemover[i] * 2.0;
            }
            for (int i = 0, limit = fftSize / 2; i < limit; i++)
            {
                dcRemover[i] /= dcComponent;
                dcRemover[fftSize - i - 1] = dcRemover[i];
            }

            return dcRemover;
        }

        ~SynthesisRealTime()
        {
            InverseRealFFT.Release();
            ForwardRealFFT.Release();
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace WORLD.Sharp
{
    static class MatlabFunctions
    {
        public static void FilterForDecimate(double[] x, int r, double[] result)
        {
            var a = new double[3];
            var b = new double[2];
            switch (r)
            {
                case 11:  // fs : 44100 (default)
                    a[0] = 2.450743295230728;
                    a[1] = -2.06794904601978;
                    a[2] = 0.59574774438332101;
                    b[0] = 0.0026822508007163792;
                    b[1] = 0.0080467524021491377;
                    break;
                case 12:  // fs : 48000
                    a[0] = 2.4981398605924205;
                    a[1] = -2.1368928194784025;
                    a[2] = 0.62187513816221485;
                    b[0] = 0.0021097275904709001;
                    b[1] = 0.0063291827714127002;
                    break;
                case 10:
                    a[0] = 2.3936475118069387;
                    a[1] = -1.9873904075111861;
                    a[2] = 0.5658879979027055;
                    b[0] = 0.0034818622251927556;
                    b[1] = 0.010445586675578267;
                    break;
                case 9:
                    a[0] = 2.3236003491759578;
                    a[1] = -1.8921545617463598;
                    a[2] = 0.53148928133729068;
                    b[0] = 0.0046331164041389372;
                    b[1] = 0.013899349212416812;
                    break;
                case 8:  // fs : 32000
                    a[0] = 2.2357462340187593;
                    a[1] = -1.7780899984041358;
                    a[2] = 0.49152555365968692;
                    b[0] = 0.0063522763407111993;
                    b[1] = 0.019056829022133598;
                    break;
                case 7:
                    a[0] = 2.1225239019534703;
                    a[1] = -1.6395144861046302;
                    a[2] = 0.44469707800587366;
                    b[0] = 0.0090366882681608418;
                    b[1] = 0.027110064804482525;
                    break;
                case 6:  // fs : 24000 and 22050
                    a[0] = 1.9715352749512141;
                    a[1] = -1.4686795689225347;
                    a[2] = 0.3893908434965701;
                    b[0] = 0.013469181309343825;
                    b[1] = 0.040407543928031475;
                    break;
                case 5:
                    a[0] = 1.7610939654280557;
                    a[1] = -1.2554914843859768;
                    a[2] = 0.3237186507788215;
                    b[0] = 0.021334858522387423;
                    b[1] = 0.06400457556716227;
                    break;
                case 4:  // fs : 16000
                    a[0] = 1.4499664446880227;
                    a[1] = -0.98943497080950582;
                    a[2] = 0.24578252340690215;
                    b[0] = 0.036710750339322612;
                    b[1] = 0.11013225101796784;
                    break;
                case 3:
                    a[0] = 0.95039378983237421;
                    a[1] = -0.67429146741526791;
                    a[2] = 0.15412211621346475;
                    b[0] = 0.071221945171178636;
                    b[1] = 0.21366583551353591;
                    break;
                case 2:  // fs : 8000
                    a[0] = 0.041156734567757189;
                    a[1] = -0.42599112459189636;
                    a[2] = 0.041037215479961225;
                    b[0] = 0.16797464681802227;
                    b[1] = 0.50392394045406674;
                    break;
                default:
                    a[0] = 0.0;
                    a[1] = 0.0;
                    a[2] = 0.0;
                    b[0] = 0.0;
                    b[1] = 0.0;
                    break;
            }

            var w = new double[3];
            var wt = 0.0;
            for (var i = 0; i < x.Length; i++)
            {
                wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
                result[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
                w[2] = w[1];
                w[1] = w[0];
                w[0] = wt;
            }
        }

        public static void FFTShift(Span<double> x, double[] result)
        {
            var halfPoint = x.Length / 2;
            for (var i = 0; i < halfPoint; i++)
            {
                result[i] = x[i + halfPoint];
                result[i + halfPoint] = x[i];
            }
        }

        static void Histc(Span<double> x, Span<double> edges, int[] index)
        {
            var count = 1;
            var i = 0;
            for (; i < edges.Length; i++)
            {
                index[i] = 1;
                if (edges[i] >= x[0])
                {
                    break;
                }
            }
            for (; i < edges.Length; i++)
            {
                if (edges[i] < x[count])
                {
                    index[i] = count;
                }
                else
                {
                    index[i--] = count++;
                }
                if (count == x.Length)
                {
                    break;
                }
            }
            for (count--, i++; i < edges.Length; i++)
            {
                index[i] = count;
            }
        }

        public static void Interp1(Span<double> x, double[] y, Span<double> xi, double[] yi)
        {
            var h = ArrayPoolHolder<double>.Shared.Rent(x.Length - 1);
            for (int i = 0, limit = x.Length - 1; i < limit; i++)
            {
                h[i] = x[i + 1] - x[i];
            }
            var k = ArrayPoolHolder<int>.Shared.Rent(xi.Length);
            Histc(x, xi, k);

            for (var i = 0; i < xi.Length; i++)
            {
                var s = (xi[i] - x[k[i] - 1]) / h[k[i] - 1];
                yi[i] = y[k[i] - 1] + s * (y[k[i]] - y[k[i] - 1]);
            }

            ArrayPoolHolder<double>.Shared.Return(h);
            ArrayPoolHolder<int>.Shared.Return(k);
        }

        public static void Decimate(double[] x, int r, double[] result)
        {
            const int NFact = 9;

            var tmp1 = x.Take(NFact)
                .Reverse()
                .Select((rx) => 2 * x[0] - rx)
                .Concat(x)
                .Concat(x.Skip(x.Length - 1 - NFact).Take(NFact).Select((lx) => 2 * x[x.Length - 1] - lx))
                .ToArray();
            var tmp2 = new double[tmp1.Length];

            FilterForDecimate(tmp1, r, tmp2);
            Array.Reverse(tmp2);
            tmp1 = tmp2;
            tmp2 = new double[tmp1.Length];
            FilterForDecimate(tmp1, r, tmp2);
            Array.Reverse(tmp2);
            tmp1 = tmp2;

            var nbeg = r - r * ((x.Length - 1) / r + 1) + x.Length;
            for (int i = nbeg, count = 0; i < x.Length; i += r, count++)
            {
                result[count] = tmp1[i + NFact - 1];
            }
        }

        public static int MatlabRound(double x)
        {
            return x > 0 ? (int)(x + 0.5) : (int)(x - 0.5);
        }

        static void Diff(Span<double> x, double[] result)
        {
            for (var i = x.Length - 2; i > -1; i--)
            {
                result[i] = x[i + 1] - x[i];
            }
        }

        public static void Interp1Q(double x, double shift, Span<double> y, Span<double> xi, double[] yi)
        {
            var deltaY = ArrayPoolHolder<double>.Shared.Rent(y.Length);
            Diff(y, deltaY);
            deltaY[y.Length - 1] = 0.0;

            for (var i = 0; i < xi.Length; i++)
            {
                var xiBase = (int)((xi[i] - x) / shift);
                var xiFraction = (xi[i] - x) / shift - xiBase;
                yi[i] = y[xiBase] + deltaY[xiBase] * xiFraction;
            }

            ArrayPoolHolder<double>.Shared.Return(deltaY);
        }

        static void FastFFTFilt(double[] x, double[] h, int fftSize, ForwardRealFFT forwardRealFFT, InverseRealFFT inverseRealFFT, double[] y)
        {
            for (var i = 0; i < x.Length; i++)
            {
                forwardRealFFT.Waveform[i] = x[i] / fftSize;
            }
            if (x.Length < fftSize)
            {
                Array.Clear(forwardRealFFT.Waveform, x.Length, fftSize - x.Length);
            }
            FFT.Execute(forwardRealFFT.ForwardFFT);
            var xSpectrum = new Complex[fftSize];
            Array.Copy(forwardRealFFT.Spectrum, xSpectrum, fftSize / 2 + 1);

            for (var i = 0; i < h.Length; i++)
            {
                forwardRealFFT.Waveform[i] = h[i] / fftSize;
            }
            if (h.Length < fftSize)
            {
                Array.Clear(forwardRealFFT.Waveform, h.Length, fftSize - h.Length);
            }
            FFT.Execute(forwardRealFFT.ForwardFFT);

            for (var i = fftSize / 2 + 1; i > -1; i--)
            {
                inverseRealFFT.Spectrum[i] = new Complex(
                    xSpectrum[i].Real * forwardRealFFT.Spectrum[i].Real - xSpectrum[i].Imaginary * forwardRealFFT.Spectrum[i].Imaginary,
                    xSpectrum[i].Real * forwardRealFFT.Spectrum[i].Imaginary + xSpectrum[i].Imaginary * forwardRealFFT.Spectrum[i].Real
                );
            }
            FFT.Execute(inverseRealFFT.InverseFFT);

            inverseRealFFT.Waveform.BlockCopy(0, y, 0, fftSize);
        }

        static double MatlabStd(double[] x)
        {
            var avg = x.Average();
            return Math.Sqrt(x.Sum(xe => Math.Pow(xe - avg, 2.0)) / (x.Length - 1));
        }
    }
}

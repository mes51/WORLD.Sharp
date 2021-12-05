using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace WORLD.Sharp
{
    enum FFTSign
    {
        Forward,
        Backward
    }

    enum FFTFlag
    {
        Forward,
        Backward,
        Estimate
    }

    class FFTPlan
    {
        static Dictionary<int, Tuple<int[], double[]>> DFTCache = new Dictionary<int, Tuple<int[], double[]>>();

        static Dictionary<int, Tuple<int[], double[]>> C2RCache = new Dictionary<int, Tuple<int[], double[]>>();

        static Dictionary<int, Tuple<int[], double[]>> R2CCache = new Dictionary<int, Tuple<int[], double[]>>();

        public int N { get; }
        public FFTSign Sign { get; }
        public FFTFlag Flag { get; }
        public Complex[] CIn { get; }
        public double[] In { get; }
        public Complex[] COut { get; set; }
        public double[] Out { get; }
        public double[] Input { get; }
        public int[] IP { get; }
        public double[] W { get; }

        public FFTPlan(int n, FFTSign sign, FFTFlag flag, Complex[] cin, double[] @in, Complex[] cout, double[] @out, double[] input, int[] ip, double[] w)
        {
            N = n;
            Sign = sign;
            Flag = flag;
            CIn = cin;
            In = @in;
            COut = cout;
            Out = @out;
            Input = input;
            IP = ip;
            W = w;
        }

        public static FFTPlan CreatePlanForDft(int n, Complex[] cin, Complex[] cout, FFTSign sign, FFTFlag flag)
        {
            if (!DFTCache.ContainsKey(n))
            {
                lock(DFTCache)
                {
                    if (!DFTCache.ContainsKey(n))
                    {
                        var tip = new int[n];
                        var tw = new double[n * 5 / 4];
                        MakeWT(n >> 1, tip, tw);
                        DFTCache.Add(n, Tuple.Create(tip, tw));
                    }
                }
            }
            var (ip, w) = DFTCache[n];

            return new FFTPlan(n, sign, flag, cin, null, cout, null, new double[n * 2], ip, w);
        }

        public static FFTPlan CreatePlanForDftC2R(int n, Complex[] cin, double[] @out, FFTFlag flag)
        {
            if (!C2RCache.ContainsKey(n))
            {
                lock (C2RCache)
                {
                    if (!C2RCache.ContainsKey(n))
                    {
                        var tip = new int[n];
                        var tw = new double[n * 5 / 4];
                        var c = new double[tw.Length - (n >> 2)];
                        MakeWT(n >> 2, tip, tw);
                        MakeCT(n >> 2, tip, c);
                        c.BlockCopy(0, tw, n >> 2, c.Length);
                        C2RCache.Add(n, Tuple.Create(tip, tw));
                    }
                }
            }
            var (ip, w) = C2RCache[n];

            return new FFTPlan(n, FFTSign.Backward, flag, cin, null, null, @out, new double[n], ip, w);
        }

        public static FFTPlan CreatePlanForDftR2C(int n, double[] @in, Complex[] cout, FFTFlag flag)
        {
            if (!R2CCache.ContainsKey(n))
            {
                lock (R2CCache)
                {
                    if (!R2CCache.ContainsKey(n))
                    {
                        var tip = new int[n];
                        var tw = new double[n * 5 / 4];
                        var c = new double[tw.Length - (n >> 2)];
                        MakeWT(n >> 2, tip, tw);
                        MakeCT(n >> 2, tip, c);
                        c.BlockCopy(0, tw, n >> 2, c.Length);
                        R2CCache.Add(n, Tuple.Create(tip, tw));
                    }
                }
            }
            var (ip, w) = R2CCache[n];

            return new FFTPlan(n, FFTSign.Forward, flag, null, @in, cout, null, new double[n], ip, w);
        }

        static void MakeWT(int nw, int[] ip, double[] w)
        {
            ip[0] = nw;
            ip[1] = 1;

            var nwh = 0;
            var wn4r = 0.0;
            if (nw > 2)
            {
                nwh = nw >> 1;
                var delta = Math.Atan(1.0) / nwh;
                wn4r = Math.Cos(delta * nwh);
                w[0] = 1.0;
                w[1] = wn4r;
                if (nwh == 4)
                {
                    w[2] = 0.5 / Math.Cos(delta * 2.0);
                    w[3] = 0.5 / Math.Sin(delta * 2.0);
                }
                else if (nwh > 4)
                {
                    MakeIPT(nw, ip);
                    w[2] = 0.5 / Math.Cos(delta * 2.0);
                    w[3] = 0.5 / Math.Cos(delta * 6.0);
                    for (var i = 4; i < nwh; i += 4)
                    {
                        w[i] = Math.Cos(delta * i);
                        w[i + 1] = Math.Sin(delta * i);
                        w[i + 2] = Math.Cos(3.0 * delta * i);
                        w[i + 3] = -Math.Sin(3.0 * delta * i);
                    }
                }
            }
            var nw0 = 0;
            while (nwh > 2)
            {
                var nw1 = nw0 + nwh;
                nwh >>= 1;
                w[nw1] = 1;
                w[nw1 + 1] = wn4r;
                if (nwh == 4)
                {
                    w[nw1 + 2] = w[nw0 + 4];
                    w[nw1 + 3] = w[nw0 + 5];
                }
                else if (nwh > 4)
                {
                    w[nw1 + 2] = 0.5 / w[nw0 + 4];
                    w[nw1 + 3] = 0.5 / w[nw0 + 6];
                    for (var i = 4; i < nwh; i += 4)
                    {
                        w[nw1 + i] = w[nw0 + 2 * i];
                        w[nw1 + i + 1] = w[nw0 + 2 * i + 1];
                        w[nw1 + i + 2] = w[nw0 + 2 * i + 2];
                        w[nw1 + i + 3] = w[nw0 + 2 * i + 3];
                    }
                }
                nw0 = nw1;
            }
        }

        static void MakeCT(int nc, int[] ip, double[] c)
        {
            ip[1] = nc;
            if (nc <= 1)
            {
                return;
            }

            var nch = nc >> 1;
            var delta = Math.Atan(1.0) / nch;
            c[0] = Math.Cos(delta * nch);
            c[nch] = 0.5 * c[0];
            for (var i = 1; i < nch; i++)
            {
                c[i] = 0.5 * Math.Cos(delta * i);
                c[nc - i] = 0.5 * Math.Sin(delta * i);
            }
        }

        static void MakeIPT(int nw, int[] ip)
        {
            ip[2] = 0;
            ip[3] = 16;
            var m = 2;
            for (var i = nw; i > 32; i >>= 2)
            {
                var m2 = m << 1;
                var q = m2 << 3;
                for (var j = m; j < m2; j++)
                {
                    var p = ip[j] << 2;
                    ip[m + j] = p;
                    ip[m2 + j] = p + q;
                }
                m = m2;
            }
        }
    }

    class ForwardRealFFT
    {
        public int FFTSize { get; }
        public double[] Waveform { get; }
        public Complex[] Spectrum { get; }
        public FFTPlan ForwardFFT { get; }

        ForwardRealFFT(int fftSize, double[] waveform, Complex[] spectrum, FFTPlan forwardFFT)
        {
            FFTSize = fftSize;
            Waveform = waveform;
            Spectrum = spectrum;
            ForwardFFT = forwardFFT;
        }

        public static ForwardRealFFT Create(int fftSize)
        {
            var waveform = new double[fftSize];
            var spectrum = new Complex[fftSize];
            return new ForwardRealFFT(fftSize, waveform, spectrum, FFTPlan.CreatePlanForDftR2C(fftSize, waveform, spectrum, FFTFlag.Estimate));
        }
    }

    class InverseRealFFT
    {
        public int FFTSize { get; }
        public double[] Waveform { get; }
        public Complex[] Spectrum { get; }
        public FFTPlan InverseFFT { get; }

        InverseRealFFT(int fftSize, double[] waveform, Complex[] spectrum, FFTPlan inverseFFT)
        {
            FFTSize = fftSize;
            Waveform = waveform;
            Spectrum = spectrum;
            InverseFFT = inverseFFT;
        }

        public static InverseRealFFT Create(int fftSize)
        {
            var waveform = new double[fftSize];
            var spectrum = new Complex[fftSize];
            return new InverseRealFFT(fftSize, waveform, spectrum, FFTPlan.CreatePlanForDftC2R(fftSize, spectrum, waveform, FFTFlag.Estimate));
        }
    }

    class InverseComplexFFT
    {
        public int FFTSize { get; }
        public Complex[] Input { get; }
        public Complex[] Output { get; }
        public FFTPlan InverseFFT { get; }

        InverseComplexFFT(int fftSize, Complex[] input, Complex[] output, FFTPlan inverseFFT)
        {
            FFTSize = fftSize;
            Input = input;
            Output = output;
            InverseFFT = inverseFFT;
        }

        public static InverseComplexFFT Create(int fftSize)
        {
            var input = new Complex[fftSize];
            var output = new Complex[fftSize];
            return new InverseComplexFFT(fftSize, input, output, FFTPlan.CreatePlanForDft(fftSize, input, output, FFTSign.Backward, FFTFlag.Estimate));
        }
    }

    class MinimumPhaseAnalysis
    {
        public int FFTSize { get; }
        public double[] LogSpectrum { get; }
        public Complex[] MinimumPhaseSpectrum { get; }
        public Complex[] Cepstrum { get; }
        public FFTPlan InverseFFT { get; }
        public FFTPlan ForwardFFT { get; }

        MinimumPhaseAnalysis(int fftSize, double[] logSpectrum, Complex[] minimumPhaseSpectrum, Complex[] cepstrum, FFTPlan inverseFFT, FFTPlan forwardFFT)
        {
            FFTSize = fftSize;
            LogSpectrum = logSpectrum;
            MinimumPhaseSpectrum = minimumPhaseSpectrum;
            Cepstrum = cepstrum;
            InverseFFT = inverseFFT;
            ForwardFFT = forwardFFT;
        }

        public static MinimumPhaseAnalysis Create(int fftSize)
        {
            var logSpectrum = new double[fftSize];
            var minimumPhaseSpectrum = new Complex[fftSize];
            var cepstrum = new Complex[fftSize];
            return new MinimumPhaseAnalysis(
                fftSize,
                logSpectrum,
                minimumPhaseSpectrum,
                cepstrum,
                FFTPlan.CreatePlanForDftR2C(fftSize, logSpectrum, cepstrum, FFTFlag.Estimate),
                FFTPlan.CreatePlanForDft(fftSize, cepstrum, minimumPhaseSpectrum, FFTSign.Forward, FFTFlag.Estimate)
            );
        }
    }

    static unsafe class FFT
    {
        public static void Execute(FFTPlan plan)
        {
            switch (plan.Sign)
            {
                case FFTSign.Forward:
                    ForwardFFT(plan);
                    break;
                case FFTSign.Backward:
                    BackwardFFT(plan);
                    break;
            }
        }

        //-----------------------------------------------------------------------
        // The following functions are reffered by
        // http://www.kurims.kyoto-u.ac.jp/~ooura/index.html

        static void ForwardFFT(FFTPlan plan)
        {
            if (plan.CIn == null)
            {
                plan.In.BlockCopy(0, plan.Input, 0, plan.N);

                fixed (double* w = &plan.W[0])
                fixed (double* input = &plan.Input[0])
                {
                    Rdft(plan.N, input, plan.IP, w);
                }

                fixed (Complex* coutPtr = &plan.COut[0])
                {
                    var cout = new Span<double>(coutPtr, plan.N + 2);
                    cout[0] = plan.Input[0];
                    cout[1] = 0.0;
                    plan.Input.AsSpan(2).CopyTo(cout.Slice(2, plan.N - 2));

                    for (var i = 1; i < cout.Length; i += 2)
                    {
                        cout[i] = -cout[i];
                    }

                    cout[plan.N] = plan.Input[1];
                    cout[plan.N + 1] = 0.0;
                }
            }
            else
            {
                fixed (Complex* cinPtr = &plan.CIn[0])
                {
                    new Span<double>(cinPtr, plan.N * 2).CopyTo(plan.Input);
                }

                fixed (double* w = &plan.W[0])
                fixed (double* input = &plan.Input[0])
                {
                    Cdft(plan.N * 2, input, plan.IP, w);
                }

                fixed (Complex* coutPtr = &plan.COut[0])
                {
                    var cout = new Span<double>(coutPtr, plan.N * 2);
                    plan.Input.AsSpan().CopyTo(cout);

                    for (var i = 1; i < cout.Length; i += 2)
                    {
                        cout[i] = -cout[i];
                    }
                }
            }
        }

        static void BackwardFFT(FFTPlan plan)
        {
            if (plan.COut == null)
            {
                plan.Input[0] = plan.CIn[0].Real;
                plan.Input[1] = plan.CIn[plan.N / 2].Real;
                for (int i = 1, ii = 2, limit = plan.N / 2; i < limit; i++, ii += 2)
                {
                    plan.Input[ii] = plan.CIn[i].Real;
                    plan.Input[ii + 1] = -plan.CIn[i].Imaginary;
                }
                fixed (double* w = &plan.W[0])
                fixed (double* input = &plan.Input[0])
                {
                    RdftInv(plan.N, input, plan.IP, w);
                }
                for (int i = 0; i < plan.N; i++)
                {
                    plan.Out[i] = plan.Input[i] * 2.0;
                }
            }
            else
            {
                for (int i = 0, ii = 0; i < plan.N; i++, ii += 2)
                {
                    plan.Input[ii] = plan.CIn[i].Real;
                    plan.Input[ii + 1] = plan.CIn[i].Imaginary;
                }
                fixed (double* w = &plan.W[0])
                fixed (double* input = &plan.Input[0])
                {
                    CdftInv(plan.N * 2, input, plan.IP, w);
                }
                for (int i = 0, ii = 0; i < plan.N; i++, ii += 2)
                {
                    plan.COut[i] = new Complex(plan.Input[ii], -plan.Input[ii + 1]);
                }
            }
        }

        static void Rdft(int n, double* a, int[] ip, double* w)
        {
            var nw = ip[0];
            var nc = ip[1];

            if (n > 4)
            {
                CFTFSub(n, a, ip, nw, w);
                RFTFSub(n, a, nc, &w[nw]);
            }
            else if (n == 4)
            {
                CFTFSub(n, a, ip, nw, w);
            }
            var xi = a[0] - a[1];
            a[0] += a[1];
            a[1] = xi;
        }

        static void RdftInv(int n, double* a, int[] ip, double* w)
        {
            var nw = ip[0];
            var nc = ip[1];

            a[1] = 0.5 * (a[0] - a[1]);
            a[0] -= a[1];
            if (n > 4)
            {
                RFTBSub(n, a, nc, &w[nw]);
                CFTBSub(n, a, ip, nw, w);
            }
            else if (n == 4)
            {
                CFTBSub(n, a, ip, nw, w);
            }
        }

        static void Cdft(int n, double* a, int[] ip, double* w)
        {
            CFTFSub(n, a, ip, ip[0], w);
        }

        static void CdftInv(int n, double* a, int[] ip, double* w)
        {
            CFTBSub(n, a, ip, ip[0], w);
        }

        static void CFTFSub(int n, double* a, int[] ip, int nw, double* w)
        {
            if (n > 8)
            {
                if (n > 32)
                {
                    CftF1st(n, a, &w[nw - (n >> 2)]);
                    if (n > 512)
                    {
                        CFTRec4(n, a, nw, w);
                    }
                    else if (n > 128)
                    {
                        CFTLeaf(n, 1, a, nw, w);
                    }
                    else
                    {
                        CFTX41(n, a, nw, w);
                    }
                    BitRV2(n, ip, a);
                }
                else if (n == 32)
                {
                    CFTF161(a, &w[nw - 8]);
                    BitRV216(a);
                }
                else
                {
                    CFTF081(a, w);
                    BitRV208(a);
                }
            }
            else if (n == 8)
            {
                CFTF040(a);
            }
            else if (n == 4)
            {
                CFTX020(a);
            }
        }

        static void CFTBSub(int n, double* a, int[] ip, int nw, double* w)
        {
            if (n > 8)
            {
                if (n > 32)
                {
                    CftB1st(n, a, &w[nw - (n >> 2)]);
                    if (n > 512)
                    {
                        CFTRec4(n, a, nw, w);
                    }
                    else if (n > 128)
                    {
                        CFTLeaf(n, 1, a, nw, w);
                    }
                    else
                    {
                        CFTX41(n, a, nw, w);
                    }
                    BitRV2Conj(n, ip, a);
                }
                else if (n == 32)
                {
                    CFTF161(a, &w[nw - 8]);
                    BitRV216Neg(a);
                }
                else
                {
                    CFTF081(a, w);
                    BitRV208Neg(a);
                }
            }
            else if (n == 8)
            {
                CFTB040(a);
            }
            else if (n == 4)
            {
                CFTX020(a);
            }
        }

        static void BitRV2(int n, int[] ip, double* a)
        {
            var m = 1;
            var l = n >> 2;
            for (; l > 8; l >>= 2)
            {
                m <<= 1;
            }
            var nh = n >> 1;
            var nm = 4 * m;

            var j1 = 0;
            var k1 = 0;
            var xr = 0.0;
            var xi = 0.0;
            var yr = 0.0;
            var yi = 0.0;
            if (l == 8)
            {
                for (var k = 0; k < m; k++)
                {
                    for (var j = 0; j < k; j++)
                    {
                        j1 = 4 * j + 2 * ip[m + k];
                        k1 = 4 * k + 2 * ip[m + j];
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + 2 * ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 += 2 * nm;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 -= nm;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 -= 2;
                    k1 -= nh;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nh + 2;
                    k1 += nh + 2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 -= nh - nm;
                    k1 += 2 * nm - 2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
            }
            else
            {
                for (var k = 0; k < m; k++)
                {
                    for (var j = 0; j < k; j++)
                    {
                        j1 = 4 * j + ip[m + k];
                        k1 = 4 * k + ip[m + j];
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = a[j1 + 1];
                        yr = a[k1];
                        yi = a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 += nm;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
            }
        }

        static void BitRV2Conj(int n, int[] ip, double* a)
        {
            var m = 1;
            var l = n >> 2;
            for (; l > 8; l >>= 2)
            {
                m <<= 1;
            }
            var nh = n >> 1;
            var nm = 4 * m;

            var j1 = 0;
            var k1 = 0;
            var xr = 0.0;
            var xi = 0.0;
            var yr = 0.0;
            var yi = 0.0;
            if (l == 8)
            {
                for (var k = 0; k < m; k++)
                {
                    for (var j = 0; j < k; j++)
                    {
                        j1 = 4 * j + 2 * ip[m + k];
                        k1 = 4 * k + 2 * ip[m + j];
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + 2 * ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    a[j1 - 1] = -a[j1 - 1];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                    j1 += nm;
                    k1 += 2 * nm;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 -= nm;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 -= 2;
                    k1 -= nh;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += nh + 2;
                    k1 += nh + 2;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 -= nh - nm;
                    k1 += 2 * nm - 2;
                    a[j1 - 1] = -a[j1 - 1];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                }
            }
            else
            {
                for (var k = 0; k < m; k++)
                {
                    for (var j = 0; j < k; j++)
                    {
                        j1 = 4 * j + ip[m + k];
                        k1 = 4 * k + ip[m + j];
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr = a[j1];
                        xi = -a[j1 + 1];
                        yr = a[k1];
                        yi = -a[k1 + 1];
                        a[j1] = yr;
                        a[j1 + 1] = yi;
                        a[k1] = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    a[j1 - 1] = -a[j1 - 1];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                    j1 += nm;
                    k1 += nm;
                    a[j1 - 1] = -a[j1 - 1];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                }
            }

        }

        static void BitRV216(double* a)
        {
            var x1r = a[2];
            var x1i = a[3];
            var x2r = a[4];
            var x2i = a[5];
            var x3r = a[6];
            var x3i = a[7];
            var x4r = a[8];
            var x4i = a[9];
            var x5r = a[10];
            var x5i = a[11];
            var x7r = a[14];
            var x7i = a[15];
            var x8r = a[16];
            var x8i = a[17];
            var x10r = a[20];
            var x10i = a[21];
            var x11r = a[22];
            var x11i = a[23];
            var x12r = a[24];
            var x12i = a[25];
            var x13r = a[26];
            var x13i = a[27];
            var x14r = a[28];
            var x14i = a[29];
            a[2] = x8r;
            a[3] = x8i;
            a[4] = x4r;
            a[5] = x4i;
            a[6] = x12r;
            a[7] = x12i;
            a[8] = x2r;
            a[9] = x2i;
            a[10] = x10r;
            a[11] = x10i;
            a[14] = x14r;
            a[15] = x14i;
            a[16] = x1r;
            a[17] = x1i;
            a[20] = x5r;
            a[21] = x5i;
            a[22] = x13r;
            a[23] = x13i;
            a[24] = x3r;
            a[25] = x3i;
            a[26] = x11r;
            a[27] = x11i;
            a[28] = x7r;
            a[29] = x7i;
        }

        static void BitRV216Neg(double* a)
        {
            var x1r = a[2];
            var x1i = a[3];
            var x2r = a[4];
            var x2i = a[5];
            var x3r = a[6];
            var x3i = a[7];
            var x4r = a[8];
            var x4i = a[9];
            var x5r = a[10];
            var x5i = a[11];
            var x6r = a[12];
            var x6i = a[13];
            var x7r = a[14];
            var x7i = a[15];
            var x8r = a[16];
            var x8i = a[17];
            var x9r = a[18];
            var x9i = a[19];
            var x10r = a[20];
            var x10i = a[21];
            var x11r = a[22];
            var x11i = a[23];
            var x12r = a[24];
            var x12i = a[25];
            var x13r = a[26];
            var x13i = a[27];
            var x14r = a[28];
            var x14i = a[29];
            var x15r = a[30];
            var x15i = a[31];
            a[2] = x15r;
            a[3] = x15i;
            a[4] = x7r;
            a[5] = x7i;
            a[6] = x11r;
            a[7] = x11i;
            a[8] = x3r;
            a[9] = x3i;
            a[10] = x13r;
            a[11] = x13i;
            a[12] = x5r;
            a[13] = x5i;
            a[14] = x9r;
            a[15] = x9i;
            a[16] = x1r;
            a[17] = x1i;
            a[18] = x14r;
            a[19] = x14i;
            a[20] = x6r;
            a[21] = x6i;
            a[22] = x10r;
            a[23] = x10i;
            a[24] = x2r;
            a[25] = x2i;
            a[26] = x12r;
            a[27] = x12i;
            a[28] = x4r;
            a[29] = x4i;
            a[30] = x8r;
            a[31] = x8i;
        }

        static void BitRV208(double* a)
        {
            var x1r = a[2];
            var x1i = a[3];
            var x3r = a[6];
            var x3i = a[7];
            var x4r = a[8];
            var x4i = a[9];
            var x6r = a[12];
            var x6i = a[13];
            a[2] = x4r;
            a[3] = x4i;
            a[6] = x6r;
            a[7] = x6i;
            a[8] = x1r;
            a[9] = x1i;
            a[12] = x3r;
            a[13] = x3i;
        }

        static void BitRV208Neg(double* a)
        {
            var x1r = a[2];
            var x1i = a[3];
            var x2r = a[4];
            var x2i = a[5];
            var x3r = a[6];
            var x3i = a[7];
            var x4r = a[8];
            var x4i = a[9];
            var x5r = a[10];
            var x5i = a[11];
            var x6r = a[12];
            var x6i = a[13];
            var x7r = a[14];
            var x7i = a[15];
            a[2] = x7r;
            a[3] = x7i;
            a[4] = x3r;
            a[5] = x3i;
            a[6] = x5r;
            a[7] = x5i;
            a[8] = x1r;
            a[9] = x1i;
            a[10] = x6r;
            a[11] = x6i;
            a[12] = x2r;
            a[13] = x2i;
            a[14] = x4r;
            a[15] = x4i;
        }

        static void CftF1st(int n, double* a, double* w)
        {
            var mh = n >> 3;
            var m = 2 * mh;
            var j0 = 0;
            var j1 = m;
            var j2 = j1 + m;
            var j3 = j2 + m;
            var x0r = a[0] + a[j2];
            var x0i = a[1] + a[j2 + 1];
            var x1r = a[0] - a[j2];
            var x1i = a[1] - a[j2 + 1];
            var x2r = a[j1] + a[j3];
            var x2i = a[j1 + 1] + a[j3 + 1];
            var x3r = a[j1] - a[j3];
            var x3i = a[j1 + 1] - a[j3 + 1];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            a[j2] = x1r - x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
            var wn4r = w[1];
            var csc1 = w[2];
            var csc3 = w[3];
            var wd1r = 1.0;
            var wd1i = 0.0;
            var wd3r = 1.0;
            var wd3i = 0.0;
            var wk1r = 0.0;
            var wk1i = 0.0;
            var wk3r = 0.0;
            var wk3i = 0.0;
            var k = 0;
            for (var j = 2; j < mh - 2; j += 4)
            {
                k += 4;
                wk1r = csc1 * (wd1r + w[k]);
                wk1i = csc1 * (wd1i + w[k + 1]);
                wk3r = csc3 * (wd3r + w[k + 2]);
                wk3i = csc3 * (wd3i + w[k + 3]);
                wd1r = w[k];
                wd1i = w[k + 1];
                wd3r = w[k + 2];
                wd3i = w[k + 3];
                j1 = j + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j] + a[j2];
                x0i = a[j + 1] + a[j2 + 1];
                x1r = a[j] - a[j2];
                x1i = a[j + 1] - a[j2 + 1];
                var y0r = a[j + 2] + a[j2 + 2];
                var y0i = a[j + 3] + a[j2 + 3];
                var y1r = a[j + 2] - a[j2 + 2];
                var y1i = a[j + 3] - a[j2 + 3];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                var y2r = a[j1 + 2] + a[j3 + 2];
                var y2i = a[j1 + 3] + a[j3 + 3];
                var y3r = a[j1 + 2] - a[j3 + 2];
                var y3i = a[j1 + 3] - a[j3 + 3];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                a[j + 2] = y0r + y2r;
                a[j + 3] = y0i + y2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                a[j1 + 2] = y0r - y2r;
                a[j1 + 3] = y0i - y2i;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j2] = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r = y1r - y3i;
                x0i = y1i + y3r;
                a[j2 + 2] = wd1r * x0r - wd1i * x0i;
                a[j2 + 3] = wd1r * x0i + wd1i * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                x0r = y1r + y3i;
                x0i = y1i - y3r;
                a[j3 + 2] = wd3r * x0r + wd3i * x0i;
                a[j3 + 3] = wd3r * x0i - wd3i * x0r;
                j0 = m - j;
                j1 = j0 + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j0] + a[j2];
                x0i = a[j0 + 1] + a[j2 + 1];
                x1r = a[j0] - a[j2];
                x1i = a[j0 + 1] - a[j2 + 1];
                y0r = a[j0 - 2] + a[j2 - 2];
                y0i = a[j0 - 1] + a[j2 - 1];
                y1r = a[j0 - 2] - a[j2 - 2];
                y1i = a[j0 - 1] - a[j2 - 1];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                y2r = a[j1 - 2] + a[j3 - 2];
                y2i = a[j1 - 1] + a[j3 - 1];
                y3r = a[j1 - 2] - a[j3 - 2];
                y3i = a[j1 - 1] - a[j3 - 1];
                a[j0] = x0r + x2r;
                a[j0 + 1] = x0i + x2i;
                a[j0 - 2] = y0r + y2r;
                a[j0 - 1] = y0i + y2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                a[j1 - 2] = y0r - y2r;
                a[j1 - 1] = y0i - y2i;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j2] = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r = y1r - y3i;
                x0i = y1i + y3r;
                a[j2 - 2] = wd1i * x0r - wd1r * x0i;
                a[j2 - 1] = wd1i * x0i + wd1r * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
                x0r = y1r + y3i;
                x0i = y1i - y3r;
                a[j3 - 2] = wd3i * x0r + wd3r * x0i;
                a[j3 - 1] = wd3i * x0i - wd3r * x0r;
            }
            wk1r = csc1 * (wd1r + wn4r);
            wk1i = csc1 * (wd1i + wn4r);
            wk3r = csc3 * (wd3r - wn4r);
            wk3i = csc3 * (wd3i - wn4r);
            j0 = mh;
            j1 = j0 + m;
            j2 = j1 + m;
            j3 = j2 + m;
            x0r = a[j0 - 2] + a[j2 - 2];
            x0i = a[j0 - 1] + a[j2 - 1];
            x1r = a[j0 - 2] - a[j2 - 2];
            x1i = a[j0 - 1] - a[j2 - 1];
            x2r = a[j1 - 2] + a[j3 - 2];
            x2i = a[j1 - 1] + a[j3 - 1];
            x3r = a[j1 - 2] - a[j3 - 2];
            x3i = a[j1 - 1] - a[j3 - 1];
            a[j0 - 2] = x0r + x2r;
            a[j0 - 1] = x0i + x2i;
            a[j1 - 2] = x0r - x2r;
            a[j1 - 1] = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j2 - 2] = wk1r * x0r - wk1i * x0i;
            a[j2 - 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3 - 2] = wk3r * x0r + wk3i * x0i;
            a[j3 - 1] = wk3r * x0i - wk3i * x0r;
            x0r = a[j0] + a[j2];
            x0i = a[j0 + 1] + a[j2 + 1];
            x1r = a[j0] - a[j2];
            x1i = a[j0 + 1] - a[j2 + 1];
            x2r = a[j1] + a[j3];
            x2i = a[j1 + 1] + a[j3 + 1];
            x3r = a[j1] - a[j3];
            x3i = a[j1 + 1] - a[j3 + 1];
            a[j0] = x0r + x2r;
            a[j0 + 1] = x0i + x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j2] = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
            x0r = a[j0 + 2] + a[j2 + 2];
            x0i = a[j0 + 3] + a[j2 + 3];
            x1r = a[j0 + 2] - a[j2 + 2];
            x1i = a[j0 + 3] - a[j2 + 3];
            x2r = a[j1 + 2] + a[j3 + 2];
            x2i = a[j1 + 3] + a[j3 + 3];
            x3r = a[j1 + 2] - a[j3 + 2];
            x3i = a[j1 + 3] - a[j3 + 3];
            a[j0 + 2] = x0r + x2r;
            a[j0 + 3] = x0i + x2i;
            a[j1 + 2] = x0r - x2r;
            a[j1 + 3] = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j2 + 2] = wk1i * x0r - wk1r * x0i;
            a[j2 + 3] = wk1i * x0i + wk1r * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3 + 2] = wk3i * x0r + wk3r * x0i;
            a[j3 + 3] = wk3i * x0i - wk3r * x0r;
        }

        static void CftB1st(int n, double* a, double* w)
        {
            var mh = n >> 3;
            var m = 2 * mh;
            var j0 = 0;
            var j1 = m;
            var j2 = j1 + m;
            var j3 = j2 + m;
            var x0r = a[0] + a[j2];
            var x0i = -a[1] - a[j2 + 1];
            var x1r = a[0] - a[j2];
            var x1i = -a[1] + a[j2 + 1];
            var x2r = a[j1] + a[j3];
            var x2i = a[j1 + 1] + a[j3 + 1];
            var x3r = a[j1] - a[j3];
            var x3i = a[j1 + 1] - a[j3 + 1];
            a[0] = x0r + x2r;
            a[1] = x0i - x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i + x2i;
            a[j2] = x1r + x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3] = x1r - x3i;
            a[j3 + 1] = x1i - x3r;
            var wn4r = w[1];
            var csc1 = w[2];
            var csc3 = w[3];
            var wd1r = 1.0;
            var wd1i = 0.0;
            var wd3r = 1.0;
            var wd3i = 0.0;
            var wk1r = 0.0;
            var wk1i = 0.0;
            var wk3r = 0.0;
            var wk3i = 0.0;
            var k = 0;
            for (var j = 2; j < mh - 2; j += 4)
            {
                k += 4;
                wk1r = csc1 * (wd1r + w[k]);
                wk1i = csc1 * (wd1i + w[k + 1]);
                wk3r = csc3 * (wd3r + w[k + 2]);
                wk3i = csc3 * (wd3i + w[k + 3]);
                wd1r = w[k];
                wd1i = w[k + 1];
                wd3r = w[k + 2];
                wd3i = w[k + 3];
                j1 = j + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j] + a[j2];
                x0i = -a[j + 1] - a[j2 + 1];
                x1r = a[j] - a[j2];
                x1i = -a[j + 1] + a[j2 + 1];
                var y0r = a[j + 2] + a[j2 + 2];
                var y0i = -a[j + 3] - a[j2 + 3];
                var y1r = a[j + 2] - a[j2 + 2];
                var y1i = -a[j + 3] + a[j2 + 3];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                var y2r = a[j1 + 2] + a[j3 + 2];
                var y2i = a[j1 + 3] + a[j3 + 3];
                var y3r = a[j1 + 2] - a[j3 + 2];
                var y3i = a[j1 + 3] - a[j3 + 3];
                a[j] = x0r + x2r;
                a[j + 1] = x0i - x2i;
                a[j + 2] = y0r + y2r;
                a[j + 3] = y0i - y2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i + x2i;
                a[j1 + 2] = y0r - y2r;
                a[j1 + 3] = y0i + y2i;
                x0r = x1r + x3i;
                x0i = x1i + x3r;
                a[j2] = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r = y1r + y3i;
                x0i = y1i + y3r;
                a[j2 + 2] = wd1r * x0r - wd1i * x0i;
                a[j2 + 3] = wd1r * x0i + wd1i * x0r;
                x0r = x1r - x3i;
                x0i = x1i - x3r;
                a[j3] = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                x0r = y1r - y3i;
                x0i = y1i - y3r;
                a[j3 + 2] = wd3r * x0r + wd3i * x0i;
                a[j3 + 3] = wd3r * x0i - wd3i * x0r;
                j0 = m - j;
                j1 = j0 + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j0] + a[j2];
                x0i = -a[j0 + 1] - a[j2 + 1];
                x1r = a[j0] - a[j2];
                x1i = -a[j0 + 1] + a[j2 + 1];
                y0r = a[j0 - 2] + a[j2 - 2];
                y0i = -a[j0 - 1] - a[j2 - 1];
                y1r = a[j0 - 2] - a[j2 - 2];
                y1i = -a[j0 - 1] + a[j2 - 1];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                y2r = a[j1 - 2] + a[j3 - 2];
                y2i = a[j1 - 1] + a[j3 - 1];
                y3r = a[j1 - 2] - a[j3 - 2];
                y3i = a[j1 - 1] - a[j3 - 1];
                a[j0] = x0r + x2r;
                a[j0 + 1] = x0i - x2i;
                a[j0 - 2] = y0r + y2r;
                a[j0 - 1] = y0i - y2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i + x2i;
                a[j1 - 2] = y0r - y2r;
                a[j1 - 1] = y0i + y2i;
                x0r = x1r + x3i;
                x0i = x1i + x3r;
                a[j2] = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r = y1r + y3i;
                x0i = y1i + y3r;
                a[j2 - 2] = wd1i * x0r - wd1r * x0i;
                a[j2 - 1] = wd1i * x0i + wd1r * x0r;
                x0r = x1r - x3i;
                x0i = x1i - x3r;
                a[j3] = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
                x0r = y1r - y3i;
                x0i = y1i - y3r;
                a[j3 - 2] = wd3i * x0r + wd3r * x0i;
                a[j3 - 1] = wd3i * x0i - wd3r * x0r;
            }
            wk1r = csc1 * (wd1r + wn4r);
            wk1i = csc1 * (wd1i + wn4r);
            wk3r = csc3 * (wd3r - wn4r);
            wk3i = csc3 * (wd3i - wn4r);
            j0 = mh;
            j1 = j0 + m;
            j2 = j1 + m;
            j3 = j2 + m;
            x0r = a[j0 - 2] + a[j2 - 2];
            x0i = -a[j0 - 1] - a[j2 - 1];
            x1r = a[j0 - 2] - a[j2 - 2];
            x1i = -a[j0 - 1] + a[j2 - 1];
            x2r = a[j1 - 2] + a[j3 - 2];
            x2i = a[j1 - 1] + a[j3 - 1];
            x3r = a[j1 - 2] - a[j3 - 2];
            x3i = a[j1 - 1] - a[j3 - 1];
            a[j0 - 2] = x0r + x2r;
            a[j0 - 1] = x0i - x2i;
            a[j1 - 2] = x0r - x2r;
            a[j1 - 1] = x0i + x2i;
            x0r = x1r + x3i;
            x0i = x1i + x3r;
            a[j2 - 2] = wk1r * x0r - wk1i * x0i;
            a[j2 - 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r - x3i;
            x0i = x1i - x3r;
            a[j3 - 2] = wk3r * x0r + wk3i * x0i;
            a[j3 - 1] = wk3r * x0i - wk3i * x0r;
            x0r = a[j0] + a[j2];
            x0i = -a[j0 + 1] - a[j2 + 1];
            x1r = a[j0] - a[j2];
            x1i = -a[j0 + 1] + a[j2 + 1];
            x2r = a[j1] + a[j3];
            x2i = a[j1 + 1] + a[j3 + 1];
            x3r = a[j1] - a[j3];
            x3i = a[j1 + 1] - a[j3 + 1];
            a[j0] = x0r + x2r;
            a[j0 + 1] = x0i - x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i + x2i;
            x0r = x1r + x3i;
            x0i = x1i + x3r;
            a[j2] = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r = x1r - x3i;
            x0i = x1i - x3r;
            a[j3] = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
            x0r = a[j0 + 2] + a[j2 + 2];
            x0i = -a[j0 + 3] - a[j2 + 3];
            x1r = a[j0 + 2] - a[j2 + 2];
            x1i = -a[j0 + 3] + a[j2 + 3];
            x2r = a[j1 + 2] + a[j3 + 2];
            x2i = a[j1 + 3] + a[j3 + 3];
            x3r = a[j1 + 2] - a[j3 + 2];
            x3i = a[j1 + 3] - a[j3 + 3];
            a[j0 + 2] = x0r + x2r;
            a[j0 + 3] = x0i - x2i;
            a[j1 + 2] = x0r - x2r;
            a[j1 + 3] = x0i + x2i;
            x0r = x1r + x3i;
            x0i = x1i + x3r;
            a[j2 + 2] = wk1i * x0r - wk1r * x0i;
            a[j2 + 3] = wk1i * x0i + wk1r * x0r;
            x0r = x1r - x3i;
            x0i = x1i - x3r;
            a[j3 + 2] = wk3i * x0r + wk3r * x0i;
            a[j3 + 3] = wk3i * x0i - wk3r * x0r;
        }

        static void CFTRec4(int n, double* a, int nw, double* w)
        {
            var m = n;
            while (m > 512)
            {
                m >>= 2;
                CFTMdl1(m, &a[n - m], &w[nw - (m >> 1)]);
            }
            CFTLeaf(m, 1, &a[n - m], nw, w);
            var k = 0;
            for (var j = n - m; j > 0; j -= m)
            {
                k++;
                var isplt = CFTTree(m, j, k, a, nw, w);
                CFTLeaf(m, isplt, &a[j - m], nw, w);
            }
        }

        static int CFTTree(int n, int j, int k, double* a, int nw, double* w)
        {
            var isplt = 0;
            if ((k & 3) != 0)
            {
                isplt = k & 1;
                if (isplt != 0)
                {
                    CFTMdl1(n, &a[j - n], &w[nw - (n >> 1)]);
                }
                else
                {
                    CFTMdl2(n, &a[j - n], &w[nw - n]);
                }
            }
            else
            {
                var m = n;
                var i = k;
                for (; (i & 3) == 0; i >>= 2)
                {
                    m <<= 2;
                }
                isplt = i & 1;
                if (isplt != 0)
                {
                    while (m > 128)
                    {
                        CFTMdl1(m, &a[j - m], &w[nw - (m >> 1)]);
                        m >>= 2;
                    }
                }
                else
                {
                    while (m > 128)
                    {
                        CFTMdl2(m, &a[j - m], &w[nw - m]);
                        m >>= 2;
                    }
                }
            }
            return isplt;
        }

        static void CFTLeaf(int n, int isplt, double* a, int nw, double* w)
        {
            if (n == 512)
            {
                CFTMdl1(128, a, &w[nw - 64]);
                CFTF161(a, &w[nw - 8]);
                CFTF162(&a[32], &w[nw - 32]);
                CFTF161(&a[64], &w[nw - 8]);
                CFTF161(&a[96], &w[nw - 8]);
                CFTMdl2(128, &a[128], &w[nw - 128]);
                CFTF161(&a[128], &w[nw - 8]);
                CFTF162(&a[160], &w[nw - 32]);
                CFTF161(&a[192], &w[nw - 8]);
                CFTF162(&a[224], &w[nw - 32]);
                CFTMdl1(128, &a[256], &w[nw - 64]);
                CFTF161(&a[256], &w[nw - 8]);
                CFTF162(&a[288], &w[nw - 32]);
                CFTF161(&a[320], &w[nw - 8]);
                CFTF161(&a[352], &w[nw - 8]);
                if (isplt != 0)
                {
                    CFTMdl1(128, &a[384], &w[nw - 64]);
                    CFTF161(&a[480], &w[nw - 8]);
                }
                else
                {
                    CFTMdl2(128, &a[384], &w[nw - 128]);
                    CFTF162(&a[480], &w[nw - 32]);
                }
                CFTF161(&a[384], &w[nw - 8]);
                CFTF162(&a[416], &w[nw - 32]);
                CFTF161(&a[448], &w[nw - 8]);
            }
            else
            {
                CFTMdl1(64, a, &w[nw - 32]);
                CFTF081(a, &w[nw - 8]);
                CFTF082(&a[16], &w[nw - 8]);
                CFTF081(&a[32], &w[nw - 8]);
                CFTF081(&a[48], &w[nw - 8]);
                CFTMdl2(64, &a[64], &w[nw - 64]);
                CFTF081(&a[64], &w[nw - 8]);
                CFTF082(&a[80], &w[nw - 8]);
                CFTF081(&a[96], &w[nw - 8]);
                CFTF082(&a[112], &w[nw - 8]);
                CFTMdl1(64, &a[128], &w[nw - 32]);
                CFTF081(&a[128], &w[nw - 8]);
                CFTF082(&a[144], &w[nw - 8]);
                CFTF081(&a[160], &w[nw - 8]);
                CFTF081(&a[176], &w[nw - 8]);
                if (isplt != 0)
                {
                    CFTMdl1(64, &a[192], &w[nw - 32]);
                    CFTF081(&a[240], &w[nw - 8]);
                }
                else
                {
                    CFTMdl2(64, &a[192], &w[nw - 64]);
                    CFTF082(&a[240], &w[nw - 8]);
                }
                CFTF081(&a[192], &w[nw - 8]);
                CFTF082(&a[208], &w[nw - 8]);
                CFTF081(&a[224], &w[nw - 8]);
            }
        }

        static void CFTMdl1(int n, double* a, double* w)
        {
            var mh = n >> 3;
            var m = 2 * mh;
            var j0 = 0;
            var j1 = m;
            var j2 = j1 + m;
            var j3 = j2 + m;
            var x0r = a[0] + a[j2];
            var x0i = a[1] + a[j2 + 1];
            var x1r = a[0] - a[j2];
            var x1i = a[1] - a[j2 + 1];
            var x2r = a[j1] + a[j3];
            var x2i = a[j1 + 1] + a[j3 + 1];
            var x3r = a[j1] - a[j3];
            var x3i = a[j1 + 1] - a[j3 + 1];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            a[j2] = x1r - x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
            var wk1r = 0.0;
            var wk1i = 0.0;
            var wn4r = w[1];
            var k = 0;
            for (var j = 2; j < mh; j += 2)
            {
                k += 4;
                wk1r = w[k];
                wk1i = w[k + 1];
                var wk3r = w[k + 2];
                var wk3i = w[k + 3];
                j1 = j + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j] + a[j2];
                x0i = a[j + 1] + a[j2 + 1];
                x1r = a[j] - a[j2];
                x1i = a[j + 1] - a[j2 + 1];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j2] = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                j0 = m - j;
                j1 = j0 + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j0] + a[j2];
                x0i = a[j0 + 1] + a[j2 + 1];
                x1r = a[j0] - a[j2];
                x1i = a[j0 + 1] - a[j2 + 1];
                x2r = a[j1] + a[j3];
                x2i = a[j1 + 1] + a[j3 + 1];
                x3r = a[j1] - a[j3];
                x3i = a[j1 + 1] - a[j3 + 1];
                a[j0] = x0r + x2r;
                a[j0 + 1] = x0i + x2i;
                a[j1] = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j2] = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
            }
            j0 = mh;
            j1 = j0 + m;
            j2 = j1 + m;
            j3 = j2 + m;
            x0r = a[j0] + a[j2];
            x0i = a[j0 + 1] + a[j2 + 1];
            x1r = a[j0] - a[j2];
            x1i = a[j0 + 1] - a[j2 + 1];
            x2r = a[j1] + a[j3];
            x2i = a[j1 + 1] + a[j3 + 1];
            x3r = a[j1] - a[j3];
            x3i = a[j1 + 1] - a[j3 + 1];
            a[j0] = x0r + x2r;
            a[j0 + 1] = x0i + x2i;
            a[j1] = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j2] = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
        }

        static void CFTMdl2(int n, double* a, double* w)
        {
            var mh = n >> 3;
            var m = 2 * mh;
            var wk1r = 0.0;
            var wk1i = 0.0;
            var wn4r = w[1];
            var j0 = 0;
            var j1 = m;
            var j2 = j1 + m;
            var j3 = j2 + m;
            var x0r = a[0] - a[j2 + 1];
            var x0i = a[1] + a[j2];
            var x1r = a[0] + a[j2 + 1];
            var x1i = a[1] - a[j2];
            var x2r = a[j1] - a[j3 + 1];
            var x2i = a[j1 + 1] + a[j3];
            var x3r = a[j1] + a[j3 + 1];
            var x3i = a[j1 + 1] - a[j3];
            var y0r = wn4r * (x2r - x2i);
            var y0i = wn4r * (x2i + x2r);
            var y2r = 0.0;
            var y2i = 0.0;
            a[0] = x0r + y0r;
            a[1] = x0i + y0i;
            a[j1] = x0r - y0r;
            a[j1 + 1] = x0i - y0i;
            y0r = wn4r * (x3r - x3i);
            y0i = wn4r * (x3i + x3r);
            a[j2] = x1r - y0i;
            a[j2 + 1] = x1i + y0r;
            a[j3] = x1r + y0i;
            a[j3 + 1] = x1i - y0r;
            var k = 0;
            var kr = 2 * m;
            for (var j = 2; j < mh; j += 2)
            {
                k += 4;
                wk1r = w[k];
                wk1i = w[k + 1];
                var wk3r = w[k + 2];
                var wk3i = w[k + 3];
                kr -= 4;
                var wd1i = w[kr];
                var wd1r = w[kr + 1];
                var wd3i = w[kr + 2];
                var wd3r = w[kr + 3];
                j1 = j + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j] - a[j2 + 1];
                x0i = a[j + 1] + a[j2];
                x1r = a[j] + a[j2 + 1];
                x1i = a[j + 1] - a[j2];
                x2r = a[j1] - a[j3 + 1];
                x2i = a[j1 + 1] + a[j3];
                x3r = a[j1] + a[j3 + 1];
                x3i = a[j1 + 1] - a[j3];
                y0r = wk1r * x0r - wk1i * x0i;
                y0i = wk1r * x0i + wk1i * x0r;
                y2r = wd1r * x2r - wd1i * x2i;
                y2i = wd1r * x2i + wd1i * x2r;
                a[j] = y0r + y2r;
                a[j + 1] = y0i + y2i;
                a[j1] = y0r - y2r;
                a[j1 + 1] = y0i - y2i;
                y0r = wk3r * x1r + wk3i * x1i;
                y0i = wk3r * x1i - wk3i * x1r;
                y2r = wd3r * x3r + wd3i * x3i;
                y2i = wd3r * x3i - wd3i * x3r;
                a[j2] = y0r + y2r;
                a[j2 + 1] = y0i + y2i;
                a[j3] = y0r - y2r;
                a[j3 + 1] = y0i - y2i;
                j0 = m - j;
                j1 = j0 + m;
                j2 = j1 + m;
                j3 = j2 + m;
                x0r = a[j0] - a[j2 + 1];
                x0i = a[j0 + 1] + a[j2];
                x1r = a[j0] + a[j2 + 1];
                x1i = a[j0 + 1] - a[j2];
                x2r = a[j1] - a[j3 + 1];
                x2i = a[j1 + 1] + a[j3];
                x3r = a[j1] + a[j3 + 1];
                x3i = a[j1 + 1] - a[j3];
                y0r = wd1i * x0r - wd1r * x0i;
                y0i = wd1i * x0i + wd1r * x0r;
                y2r = wk1i * x2r - wk1r * x2i;
                y2i = wk1i * x2i + wk1r * x2r;
                a[j0] = y0r + y2r;
                a[j0 + 1] = y0i + y2i;
                a[j1] = y0r - y2r;
                a[j1 + 1] = y0i - y2i;
                y0r = wd3i * x1r + wd3r * x1i;
                y0i = wd3i * x1i - wd3r * x1r;
                y2r = wk3i * x3r + wk3r * x3i;
                y2i = wk3i * x3i - wk3r * x3r;
                a[j2] = y0r + y2r;
                a[j2 + 1] = y0i + y2i;
                a[j3] = y0r - y2r;
                a[j3 + 1] = y0i - y2i;
            }
            wk1r = w[m];
            wk1i = w[m + 1];
            j0 = mh;
            j1 = j0 + m;
            j2 = j1 + m;
            j3 = j2 + m;
            x0r = a[j0] - a[j2 + 1];
            x0i = a[j0 + 1] + a[j2];
            x1r = a[j0] + a[j2 + 1];
            x1i = a[j0 + 1] - a[j2];
            x2r = a[j1] - a[j3 + 1];
            x2i = a[j1 + 1] + a[j3];
            x3r = a[j1] + a[j3 + 1];
            x3i = a[j1 + 1] - a[j3];
            y0r = wk1r * x0r - wk1i * x0i;
            y0i = wk1r * x0i + wk1i * x0r;
            y2r = wk1i * x2r - wk1r * x2i;
            y2i = wk1i * x2i + wk1r * x2r;
            a[j0] = y0r + y2r;
            a[j0 + 1] = y0i + y2i;
            a[j1] = y0r - y2r;
            a[j1 + 1] = y0i - y2i;
            y0r = wk1i * x1r - wk1r * x1i;
            y0i = wk1i * x1i + wk1r * x1r;
            y2r = wk1r * x3r - wk1i * x3i;
            y2i = wk1r * x3i + wk1i * x3r;
            a[j2] = y0r - y2r;
            a[j2 + 1] = y0i - y2i;
            a[j3] = y0r + y2r;
            a[j3 + 1] = y0i + y2i;
        }

        static void CFTX41(int n, double* a, int nw, double* w)
        {
            if (n == 128)
            {
                CFTF161(a, &w[nw - 8]);
                CFTF162(&a[32], &w[nw - 32]);
                CFTF161(&a[64], &w[nw - 8]);
                CFTF161(&a[96], &w[nw - 8]);
            }
            else
            {
                CFTF081(a, &w[nw - 8]);
                CFTF082(&a[16], &w[nw - 8]);
                CFTF081(&a[32], &w[nw - 8]);
                CFTF081(&a[48], &w[nw - 8]);
            }
        }

        static void CFTF161(double* a, double* w)
        {
            var wn4r = w[1];
            var wk1r = w[2];
            var wk1i = w[3];
            var x0r = a[0] + a[16];
            var x0i = a[1] + a[17];
            var x1r = a[0] - a[16];
            var x1i = a[1] - a[17];
            var x2r = a[8] + a[24];
            var x2i = a[9] + a[25];
            var x3r = a[8] - a[24];
            var x3i = a[9] - a[25];
            var y0r = x0r + x2r;
            var y0i = x0i + x2i;
            var y4r = x0r - x2r;
            var y4i = x0i - x2i;
            var y8r = x1r - x3i;
            var y8i = x1i + x3r;
            var y12r = x1r + x3i;
            var y12i = x1i - x3r;
            x0r = a[2] + a[18];
            x0i = a[3] + a[19];
            x1r = a[2] - a[18];
            x1i = a[3] - a[19];
            x2r = a[10] + a[26];
            x2i = a[11] + a[27];
            x3r = a[10] - a[26];
            x3i = a[11] - a[27];
            var y1r = x0r + x2r;
            var y1i = x0i + x2i;
            var y5r = x0r - x2r;
            var y5i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            var y9r = wk1r * x0r - wk1i * x0i;
            var y9i = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            var y13r = wk1i * x0r - wk1r * x0i;
            var y13i = wk1i * x0i + wk1r * x0r;
            x0r = a[4] + a[20];
            x0i = a[5] + a[21];
            x1r = a[4] - a[20];
            x1i = a[5] - a[21];
            x2r = a[12] + a[28];
            x2i = a[13] + a[29];
            x3r = a[12] - a[28];
            x3i = a[13] - a[29];
            var y2r = x0r + x2r;
            var y2i = x0i + x2i;
            var y6r = x0r - x2r;
            var y6i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            var y10r = wn4r * (x0r - x0i);
            var y10i = wn4r * (x0i + x0r);
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            var y14r = wn4r * (x0r + x0i);
            var y14i = wn4r * (x0i - x0r);
            x0r = a[6] + a[22];
            x0i = a[7] + a[23];
            x1r = a[6] - a[22];
            x1i = a[7] - a[23];
            x2r = a[14] + a[30];
            x2i = a[15] + a[31];
            x3r = a[14] - a[30];
            x3i = a[15] - a[31];
            var y3r = x0r + x2r;
            var y3i = x0i + x2i;
            var y7r = x0r - x2r;
            var y7i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            var y11r = wk1i * x0r - wk1r * x0i;
            var y11i = wk1i * x0i + wk1r * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            var y15r = wk1r * x0r - wk1i * x0i;
            var y15i = wk1r * x0i + wk1i * x0r;
            x0r = y12r - y14r;
            x0i = y12i - y14i;
            x1r = y12r + y14r;
            x1i = y12i + y14i;
            x2r = y13r - y15r;
            x2i = y13i - y15i;
            x3r = y13r + y15r;
            x3i = y13i + y15i;
            a[24] = x0r + x2r;
            a[25] = x0i + x2i;
            a[26] = x0r - x2r;
            a[27] = x0i - x2i;
            a[28] = x1r - x3i;
            a[29] = x1i + x3r;
            a[30] = x1r + x3i;
            a[31] = x1i - x3r;
            x0r = y8r + y10r;
            x0i = y8i + y10i;
            x1r = y8r - y10r;
            x1i = y8i - y10i;
            x2r = y9r + y11r;
            x2i = y9i + y11i;
            x3r = y9r - y11r;
            x3i = y9i - y11i;
            a[16] = x0r + x2r;
            a[17] = x0i + x2i;
            a[18] = x0r - x2r;
            a[19] = x0i - x2i;
            a[20] = x1r - x3i;
            a[21] = x1i + x3r;
            a[22] = x1r + x3i;
            a[23] = x1i - x3r;
            x0r = y5r - y7i;
            x0i = y5i + y7r;
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            x0r = y5r + y7i;
            x0i = y5i - y7r;
            x3r = wn4r * (x0r - x0i);
            x3i = wn4r * (x0i + x0r);
            x0r = y4r - y6i;
            x0i = y4i + y6r;
            x1r = y4r + y6i;
            x1i = y4i - y6r;
            a[8] = x0r + x2r;
            a[9] = x0i + x2i;
            a[10] = x0r - x2r;
            a[11] = x0i - x2i;
            a[12] = x1r - x3i;
            a[13] = x1i + x3r;
            a[14] = x1r + x3i;
            a[15] = x1i - x3r;
            x0r = y0r + y2r;
            x0i = y0i + y2i;
            x1r = y0r - y2r;
            x1i = y0i - y2i;
            x2r = y1r + y3r;
            x2i = y1i + y3i;
            x3r = y1r - y3r;
            x3i = y1i - y3i;
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[2] = x0r - x2r;
            a[3] = x0i - x2i;
            a[4] = x1r - x3i;
            a[5] = x1i + x3r;
            a[6] = x1r + x3i;
            a[7] = x1i - x3r;
        }

        static void CFTF162(double* a, double* w)
        {
            var wn4r = w[1];
            var wk1r = w[4];
            var wk1i = w[5];
            var wk3r = w[6];
            var wk3i = -w[7];
            var wk2r = w[8];
            var wk2i = w[9];
            var x1r = a[0] - a[17];
            var x1i = a[1] + a[16];
            var x0r = a[8] - a[25];
            var x0i = a[9] + a[24];
            var x2r = wn4r * (x0r - x0i);
            var x2i = wn4r * (x0i + x0r);
            var y0r = x1r + x2r;
            var y0i = x1i + x2i;
            var y4r = x1r - x2r;
            var y4i = x1i - x2i;
            x1r = a[0] + a[17];
            x1i = a[1] - a[16];
            x0r = a[8] + a[25];
            x0i = a[9] - a[24];
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            var y8r = x1r - x2i;
            var y8i = x1i + x2r;
            var y12r = x1r + x2i;
            var y12i = x1i - x2r;
            x0r = a[2] - a[19];
            x0i = a[3] + a[18];
            x1r = wk1r * x0r - wk1i * x0i;
            x1i = wk1r * x0i + wk1i * x0r;
            x0r = a[10] - a[27];
            x0i = a[11] + a[26];
            x2r = wk3i * x0r - wk3r * x0i;
            x2i = wk3i * x0i + wk3r * x0r;
            var y1r = x1r + x2r;
            var y1i = x1i + x2i;
            var y5r = x1r - x2r;
            var y5i = x1i - x2i;
            x0r = a[2] + a[19];
            x0i = a[3] - a[18];
            x1r = wk3r * x0r - wk3i * x0i;
            x1i = wk3r * x0i + wk3i * x0r;
            x0r = a[10] + a[27];
            x0i = a[11] - a[26];
            x2r = wk1r * x0r + wk1i * x0i;
            x2i = wk1r * x0i - wk1i * x0r;
            var y9r = x1r - x2r;
            var y9i = x1i - x2i;
            var y13r = x1r + x2r;
            var y13i = x1i + x2i;
            x0r = a[4] - a[21];
            x0i = a[5] + a[20];
            x1r = wk2r * x0r - wk2i * x0i;
            x1i = wk2r * x0i + wk2i * x0r;
            x0r = a[12] - a[29];
            x0i = a[13] + a[28];
            x2r = wk2i * x0r - wk2r * x0i;
            x2i = wk2i * x0i + wk2r * x0r;
            var y2r = x1r + x2r;
            var y2i = x1i + x2i;
            var y6r = x1r - x2r;
            var y6i = x1i - x2i;
            x0r = a[4] + a[21];
            x0i = a[5] - a[20];
            x1r = wk2i * x0r - wk2r * x0i;
            x1i = wk2i * x0i + wk2r * x0r;
            x0r = a[12] + a[29];
            x0i = a[13] - a[28];
            x2r = wk2r * x0r - wk2i * x0i;
            x2i = wk2r * x0i + wk2i * x0r;
            var y10r = x1r - x2r;
            var y10i = x1i - x2i;
            var y14r = x1r + x2r;
            var y14i = x1i + x2i;
            x0r = a[6] - a[23];
            x0i = a[7] + a[22];
            x1r = wk3r * x0r - wk3i * x0i;
            x1i = wk3r * x0i + wk3i * x0r;
            x0r = a[14] - a[31];
            x0i = a[15] + a[30];
            x2r = wk1i * x0r - wk1r * x0i;
            x2i = wk1i * x0i + wk1r * x0r;
            var y3r = x1r + x2r;
            var y3i = x1i + x2i;
            var y7r = x1r - x2r;
            var y7i = x1i - x2i;
            x0r = a[6] + a[23];
            x0i = a[7] - a[22];
            x1r = wk1i * x0r + wk1r * x0i;
            x1i = wk1i * x0i - wk1r * x0r;
            x0r = a[14] + a[31];
            x0i = a[15] - a[30];
            x2r = wk3i * x0r - wk3r * x0i;
            x2i = wk3i * x0i + wk3r * x0r;
            var y11r = x1r + x2r;
            var y11i = x1i + x2i;
            var y15r = x1r - x2r;
            var y15i = x1i - x2i;
            x1r = y0r + y2r;
            x1i = y0i + y2i;
            x2r = y1r + y3r;
            x2i = y1i + y3i;
            a[0] = x1r + x2r;
            a[1] = x1i + x2i;
            a[2] = x1r - x2r;
            a[3] = x1i - x2i;
            x1r = y0r - y2r;
            x1i = y0i - y2i;
            x2r = y1r - y3r;
            x2i = y1i - y3i;
            a[4] = x1r - x2i;
            a[5] = x1i + x2r;
            a[6] = x1r + x2i;
            a[7] = x1i - x2r;
            x1r = y4r - y6i;
            x1i = y4i + y6r;
            x0r = y5r - y7i;
            x0i = y5i + y7r;
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            a[8] = x1r + x2r;
            a[9] = x1i + x2i;
            a[10] = x1r - x2r;
            a[11] = x1i - x2i;
            x1r = y4r + y6i;
            x1i = y4i - y6r;
            x0r = y5r + y7i;
            x0i = y5i - y7r;
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            a[12] = x1r - x2i;
            a[13] = x1i + x2r;
            a[14] = x1r + x2i;
            a[15] = x1i - x2r;
            x1r = y8r + y10r;
            x1i = y8i + y10i;
            x2r = y9r - y11r;
            x2i = y9i - y11i;
            a[16] = x1r + x2r;
            a[17] = x1i + x2i;
            a[18] = x1r - x2r;
            a[19] = x1i - x2i;
            x1r = y8r - y10r;
            x1i = y8i - y10i;
            x2r = y9r + y11r;
            x2i = y9i + y11i;
            a[20] = x1r - x2i;
            a[21] = x1i + x2r;
            a[22] = x1r + x2i;
            a[23] = x1i - x2r;
            x1r = y12r - y14i;
            x1i = y12i + y14r;
            x0r = y13r + y15i;
            x0i = y13i - y15r;
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            a[24] = x1r + x2r;
            a[25] = x1i + x2i;
            a[26] = x1r - x2r;
            a[27] = x1i - x2i;
            x1r = y12r + y14i;
            x1i = y12i - y14r;
            x0r = y13r - y15i;
            x0i = y13i + y15r;
            x2r = wn4r * (x0r - x0i);
            x2i = wn4r * (x0i + x0r);
            a[28] = x1r - x2i;
            a[29] = x1i + x2r;
            a[30] = x1r + x2i;
            a[31] = x1i - x2r;
        }

        static void CFTF081(double* a, double* w)
        {
            var wn4r = w[1];
            var x0r = a[0] + a[8];
            var x0i = a[1] + a[9];
            var x1r = a[0] - a[8];
            var x1i = a[1] - a[9];
            var x2r = a[4] + a[12];
            var x2i = a[5] + a[13];
            var x3r = a[4] - a[12];
            var x3i = a[5] - a[13];
            var y0r = x0r + x2r;
            var y0i = x0i + x2i;
            var y2r = x0r - x2r;
            var y2i = x0i - x2i;
            var y1r = x1r - x3i;
            var y1i = x1i + x3r;
            var y3r = x1r + x3i;
            var y3i = x1i - x3r;
            x0r = a[2] + a[10];
            x0i = a[3] + a[11];
            x1r = a[2] - a[10];
            x1i = a[3] - a[11];
            x2r = a[6] + a[14];
            x2i = a[7] + a[15];
            x3r = a[6] - a[14];
            x3i = a[7] - a[15];
            var y4r = x0r + x2r;
            var y4i = x0i + x2i;
            var y6r = x0r - x2r;
            var y6i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            x2r = x1r + x3i;
            x2i = x1i - x3r;
            var y5r = wn4r * (x0r - x0i);
            var y5i = wn4r * (x0r + x0i);
            var y7r = wn4r * (x2r - x2i);
            var y7i = wn4r * (x2r + x2i);
            a[8] = y1r + y5r;
            a[9] = y1i + y5i;
            a[10] = y1r - y5r;
            a[11] = y1i - y5i;
            a[12] = y3r - y7i;
            a[13] = y3i + y7r;
            a[14] = y3r + y7i;
            a[15] = y3i - y7r;
            a[0] = y0r + y4r;
            a[1] = y0i + y4i;
            a[2] = y0r - y4r;
            a[3] = y0i - y4i;
            a[4] = y2r - y6i;
            a[5] = y2i + y6r;
            a[6] = y2r + y6i;
            a[7] = y2i - y6r;
        }

        static void CFTF082(double* a, double* w)
        {
            var wn4r = w[1];
            var wk1r = w[2];
            var wk1i = w[3];
            var y0r = a[0] - a[9];
            var y0i = a[1] + a[8];
            var y1r = a[0] + a[9];
            var y1i = a[1] - a[8];
            var x0r = a[4] - a[13];
            var x0i = a[5] + a[12];
            var y2r = wn4r * (x0r - x0i);
            var y2i = wn4r * (x0i + x0r);
            x0r = a[4] + a[13];
            x0i = a[5] - a[12];
            var y3r = wn4r * (x0r - x0i);
            var y3i = wn4r * (x0i + x0r);
            x0r = a[2] - a[11];
            x0i = a[3] + a[10];
            var y4r = wk1r * x0r - wk1i * x0i;
            var y4i = wk1r * x0i + wk1i * x0r;
            x0r = a[2] + a[11];
            x0i = a[3] - a[10];
            var y5r = wk1i * x0r - wk1r * x0i;
            var y5i = wk1i * x0i + wk1r * x0r;
            x0r = a[6] - a[15];
            x0i = a[7] + a[14];
            var y6r = wk1i * x0r - wk1r * x0i;
            var y6i = wk1i * x0i + wk1r * x0r;
            x0r = a[6] + a[15];
            x0i = a[7] - a[14];
            var y7r = wk1r * x0r - wk1i * x0i;
            var y7i = wk1r * x0i + wk1i * x0r;
            x0r = y0r + y2r;
            x0i = y0i + y2i;
            var x1r = y4r + y6r;
            var x1i = y4i + y6i;
            a[0] = x0r + x1r;
            a[1] = x0i + x1i;
            a[2] = x0r - x1r;
            a[3] = x0i - x1i;
            x0r = y0r - y2r;
            x0i = y0i - y2i;
            x1r = y4r - y6r;
            x1i = y4i - y6i;
            a[4] = x0r - x1i;
            a[5] = x0i + x1r;
            a[6] = x0r + x1i;
            a[7] = x0i - x1r;
            x0r = y1r - y3i;
            x0i = y1i + y3r;
            x1r = y5r - y7r;
            x1i = y5i - y7i;
            a[8] = x0r + x1r;
            a[9] = x0i + x1i;
            a[10] = x0r - x1r;
            a[11] = x0i - x1i;
            x0r = y1r + y3i;
            x0i = y1i - y3r;
            x1r = y5r + y7r;
            x1i = y5i + y7i;
            a[12] = x0r - x1i;
            a[13] = x0i + x1r;
            a[14] = x0r + x1i;
            a[15] = x0i - x1r;
        }

        static void CFTF040(double* a)
        {
            var x0r = a[0] + a[4];
            var x0i = a[1] + a[5];
            var x1r = a[0] - a[4];
            var x1i = a[1] - a[5];
            var x2r = a[2] + a[6];
            var x2i = a[3] + a[7];
            var x3r = a[2] - a[6];
            var x3i = a[3] - a[7];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[2] = x1r - x3i;
            a[3] = x1i + x3r;
            a[4] = x0r - x2r;
            a[5] = x0i - x2i;
            a[6] = x1r + x3i;
            a[7] = x1i - x3r;
        }

        static void CFTB040(double* a)
        {
            var x0r = a[0] + a[4];
            var x0i = a[1] + a[5];
            var x1r = a[0] - a[4];
            var x1i = a[1] - a[5];
            var x2r = a[2] + a[6];
            var x2i = a[3] + a[7];
            var x3r = a[2] - a[6];
            var x3i = a[3] - a[7];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[2] = x1r + x3i;
            a[3] = x1i - x3r;
            a[4] = x0r - x2r;
            a[5] = x0i - x2i;
            a[6] = x1r - x3i;
            a[7] = x1i + x3r;
        }

        static void CFTX020(double* a)
        {
            var x0r = a[0] - a[2];
            var x0i = a[1] - a[3];
            a[0] += a[2];
            a[1] += a[3];
            a[2] = x0r;
            a[3] = x0i;
        }

        static void RFTFSub(int n, double* a, int nc, double* c)
        {
            var m = n >> 1;
            var ks = 2 * nc / m;
            var kk = 0;
            for (var j = 2; j < m; j += 2)
            {
                var k = n - j;
                kk += ks;
                var wkr = 0.5 - c[nc - kk];
                var wki = c[kk];
                var xr = a[j] - a[k];
                var xi = a[j + 1] + a[k + 1];
                var yr = wkr * xr - wki * xi;
                var yi = wkr * xi + wki * xr;
                a[j] -= yr;
                a[j + 1] -= yi;
                a[k] += yr;
                a[k + 1] -= yi;
            }
        }

        static void RFTBSub(int n, double* a, int nc, double* c)
        {
            var m = n >> 1;
            var ks = 2 * nc / m;
            var kk = 0;

            for (var j = 2; j < m; j += 2)
            {
                var k = n - j;
                kk += ks;
                var wkr = 0.5 - c[nc - kk];
                var wki = c[kk];
                var xr = a[j] - a[k];
                var xi = a[j + 1] + a[k + 1];
                var yr = wkr * xr + wki * xi;
                var yi = wkr * xi - wki * xr;
                a[j] -= yr;
                a[j + 1] -= yi;
                a[k] += yr;
                a[k + 1] -= yi;
            }
        }

        static void DCTSub(int n, double* a, int nc, double* c)
        {
            var m = n >> 1;
            var ks = nc / n;
            var kk = 0;
            for (var j = 1; j < m; j++)
            {
                var k = n - j;
                kk += ks;
                var wkr = c[kk] - c[nc - kk];
                var wki = c[kk] + c[nc - kk];
                var xr = wki * a[j] - wkr * a[k];
                a[j] = wkr * a[j] + wki * a[k];
                a[k] = xr;
            }
            a[m] *= c[0];
        }

        static void DSTSub(int n, double* a, int nc, double* c)
        {
            var m = n >> 1;
            var ks = nc / n;
            var kk = 0;
            for (var j = 1; j < m; j++)
            {
                var k = n - j;
                kk += ks;
                var wkr = c[kk] - c[nc - kk];
                var wki = c[kk] + c[nc - kk];
                var xr = wki * a[k] - wkr * a[j];
                a[k] = wkr * a[k] + wki * a[j];
                a[j] = xr;
            }
            a[m] *= c[0];
        }
    }
}

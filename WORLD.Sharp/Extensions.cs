using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    static class Extensions
    {
        public static void Fill<T>(this T[] array, T value)
        {
            array.Fill<T>(value, 0, array.Length);
        }

        public static void Fill<T>(this T[] array, T value, int begin)
        {
            array.Fill<T>(value, begin, array.Length - begin);
        }

        public static void Fill<T>(this T[] array, T value, int begin, int count)
        {
            for (int i = begin, c = 0; c < count; i++, c++)
            {
                array[i] = value;
            }
        }

        public static void ForEach<T>(this IEnumerable<T> source, Action<T> action)
        {
            foreach (var n in source)
            {
                action(n);
            }
        }

        public static void ForEach<T>(this IEnumerable<T> source, Action<T, int> action)
        {
            var i = 0;
            foreach (var n in source)
            {
                action(n, i);
                i++;
            }
        }

        public static void BlockCopy(this double[] src, double[] dst)
        {
            src.BlockCopy(0, dst, 0, Math.Min(src.Length, dst.Length));
        }

        public static void BlockCopy(this double[] src, double[] dst, int dstOffset)
        {
            src.BlockCopy(0, dst, dstOffset, Math.Min(src.Length, dst.Length - dstOffset));
        }

        public static void BlockCopy(this double[] src, int srcOffset, double[] dst, int dstOffset, int count)
        {
            const int DataSize = sizeof(double);
            Buffer.BlockCopy(src, srcOffset * DataSize, dst, dstOffset * DataSize, count * DataSize);
        }
    }
}

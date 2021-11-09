using System;
using System.Linq;

namespace WORLD.Sharp
{
    static class Util
    {
        public static T[][] Mak2DArray<T>(int row, int col)
        {
            return Enumerable.Range(0, row).Select((x) => new T[col]).ToArray();
        }
    }

    static class MathUtil
    {
        public static double Log2(double d)
        {
            const double Log2 = 0.693147180559945; //Math.Log(2);
            return Math.Log(d) / Log2;
        }
    }
}
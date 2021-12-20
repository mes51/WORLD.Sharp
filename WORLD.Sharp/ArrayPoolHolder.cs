using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace WORLD.Sharp
{
    public class ArrayPoolHolder
    {
        static internal List<Type> ArrayPoolTypes { get; } = new List<Type>();

        static readonly object LockObj = new object();

        public static void Purge()
        {
#if NETCOREAPP3_1_OR_GREATER
#else
            lock (LockObj)
            {
                foreach (var type in ArrayPoolTypes)
                {
                    var arrayPoolType = typeof(ArrayPoolHolder<>).MakeGenericType(type);
                    arrayPoolType.GetProperty(nameof(ArrayPoolHolder<object>.Shared), BindingFlags.Public | BindingFlags.Static).SetValue(null, null);
                    arrayPoolType.GetMethod(nameof(ArrayPool<object>.Create), BindingFlags.Public | BindingFlags.Static).Invoke(null, null);
                }
            }
#endif
        }

        internal static void AddType(Type type)
        {
            if (!ArrayPoolTypes.Contains(type))
            {
                lock (LockObj)
                {
                    if (!ArrayPoolTypes.Contains(type))
                    {
                        ArrayPoolTypes.Add(type);
                    }
                }
            }
        }
    }

    internal class ArrayPoolHolder<T>
    {
        public static ArrayPool<T> Shared { get; internal set; }

        static ArrayPoolHolder()
        {
#if NETCOREAPP3_1_OR_GREATER
            Shared = ArrayPool<T>.Shared;
#else
            Create();
#endif
        }

        public static void Create()
        {
#if NETCOREAPP3_1_OR_GREATER
#else
            Shared = new CustomArrayPool<T>();
            ArrayPoolHolder.AddType(typeof(T));
#endif
        }
    }

    // Re-Implementation ConfigurableArrayPool use lock and allocate array in lock.
    // suppress GC Wait by SpinLock and unnecessary allocation (as possible use returned array).
    // https://github.com/dotnet/runtime/blob/main/src/libraries/System.Private.CoreLib/src/System/Buffers/ConfigurableArrayPool.cs
    internal class CustomArrayPool<T> : ArrayPool<T>
    {
        // 1024 * 1024 * 1024
        const int BucketLength = 27;

        const int MaxBucketsToTry = 2;

        readonly Bucket[] Buckets;

        public CustomArrayPool()
        {
            Buckets = new Bucket[BucketLength];
            for (var i = 0; i < BucketLength; i++)
            {
                Buckets[i] = new Bucket(16 << i);
            }
        }

        public override T[] Rent(int minimumLength)
        {
            if (minimumLength < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(minimumLength));
            }
            else if (minimumLength == 0)
            {
                return Array.Empty<T>();
            }

            var index = SelectBucketIndex(minimumLength);
            if (index < BucketLength)
            {
                for (var c = 0; index < BucketLength && c < MaxBucketsToTry; c++, index++)
                {
                    var result = Buckets[index].Rent();
                    if (result != null)
                    {
                        return result;
                    }
                }

                return new T[Buckets[index].ArrayLength];
            }
            else
            {
                return new T[minimumLength];
            }
        }

        public override void Return(T[] array, bool clearArray = false)
        {
            if (array == null)
            {
                throw new ArgumentNullException(nameof(array));
            }
            else if (array.Length == 0)
            {
                return;
            }

            var index = SelectBucketIndex(array.Length);
            if (index < BucketLength)
            {
                if (clearArray)
                {
                    array.AsSpan().Clear();
                }
                Buckets[index].Return(array);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static int SelectBucketIndex(int bufferSize)
        {
            return FastLog2((uint)bufferSize - 1 | 15) - 3;
        }

        // SEE: https://stackoverflow.com/a/11398748
        static readonly int[] Log2Table = new int[]
        {
             0,  9,  1, 10, 13, 21,  2, 29,
            11, 14, 16, 18, 22, 25,  3, 30,
             8, 12, 20, 28, 15, 17, 24,  7,
            19, 27, 23,  6, 26,  5,  4, 31
        };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static int FastLog2(uint value)
        {
            value |= value >> 1;
            value |= value >> 2;
            value |= value >> 4;
            value |= value >> 8;
            value |= value >> 16;

            return Log2Table[(value * 0x07C4ACDD) >> 27];
        }

        private class Bucket
        {
            // SEE: https://source.dot.net/#System.Private.CoreLib/TlsOverPerCoreLockedStacksArrayPool.cs,2670702e13ff57ca,references
            const int MaxProcessorCount = 64;

            // Thread.GetCurrentProcessorId is not defined in 4.7.2
            static readonly int ThreadCount = Math.Min(Environment.ProcessorCount, MaxProcessorCount);

            internal readonly int ArrayLength;

            readonly LockedStack[] Stacks;

            int Cursor = 0;

            public Bucket(int arrayLength)
            {
                ArrayLength = arrayLength;
                Stacks = new LockedStack[ThreadCount];
                for (var i = 0; i < ThreadCount; i++)
                {
                    Stacks[i] = new LockedStack(arrayLength);
                }
            }

            public T[] Rent()
            {
                var threadIndex = Thread.CurrentThread.ManagedThreadId % ThreadCount;
                return Stacks[threadIndex].Rent();
            }

            public void Return(T[] array)
            {
                if (array.Length != ArrayLength)
                {
                    throw new ArgumentException(nameof(array));
                }

                var threadIndex = Thread.CurrentThread.ManagedThreadId % ThreadCount;
                Stacks[threadIndex].Return(array);
            }
        }

        private class LockedStack
        {
            const int ArrayCount = 16;

            readonly T[][] Arrays;

            readonly int ArrayLength;

            int Cursor;

            public LockedStack(int arrayLength)
            {
                ArrayLength = arrayLength;
                Arrays = new T[ArrayCount][];
            }

            public T[] Rent()
            {
                T[] result = null;

                lock (this)
                {
                    if (Cursor < ArrayCount)
                    {
                        result = Arrays[Cursor];
                        Arrays[Cursor] = null;
                        Cursor++;
                    }

                    if (result == null)
                    {
                        result = new T[ArrayLength];
                    }
                }

                return result;
            }

            public void Return(T[] array)
            {
                lock (this)
                {
                    if (Cursor > 0)
                    {
                        Cursor--;
                        Arrays[Cursor] = array;
                    }
                }
            }
        }
    }
}
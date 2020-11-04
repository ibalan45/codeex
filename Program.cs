using System;
using MathNet.Numerics.LinearAlgebra.Double;

namespace _optimization
{
    class Program
    {
        private static void Main(string[] args)
        {
            var function = new Function(
                (Vector p) =>
                { return (Math.Pow((Math.Pow(p.At(0), 2) + p.At(1) - 11), 2) + Math.Pow((p.At(0) + Math.Pow(p.At(1), 2) - 7), 2)); },
                (Vector p) =>
                { return 4 * p.At(0) * (Math.Pow(p.At(0), 2) + p.At(1) - 11) + 2 * p.At(0) + 2 * (float)Math.Pow(p.At(1), 2) - 14; },
                (Vector p) =>
                { return 2 * Math.Pow(p.At(0), 2) + 4 * p.At(1) * (p.At(0) + (float)Math.Pow(p.At(1), 2) - 7) + 2 * p.At(1) - 22; }
            );
            var gd = new Optimization_1 ( function);
            Console.WriteLine(gd.GradientDescent());
            Console.WriteLine(gd.BFGS());
        }
    }
}
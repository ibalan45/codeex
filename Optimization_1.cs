using System;
using MathNet.Numerics.LinearAlgebra.Double;

namespace _optimization
{
    struct Function
    {
        public Func<Vector, double> function;

        public Func<Vector, double> xDerivative;
        public Func<Vector, double> yDerivative;

        public Function(Func<Vector, double> function,
            Func<Vector, double> xDerivative,
            Func<Vector, double> yDerivative)
        {
            this.function = function;
            this.xDerivative = xDerivative;
            this.yDerivative = yDerivative;
        }
    }

    class Optimization_1
    {
        private Function function;
        private Vector min;

        public Optimization_1(Function function)
        {
            this.function = function;
        }

        //метод наискорейшего спуска
        public Vector GradientDescent()
        {
            var currentPoint = DenseVector.OfArray(new double[] { 20, 20 });
            var antigradient = DenseVector.OfArray(new double[] { 0, 0 });
            do
            {
                //поиск антиградиента в точке
                antigradient = DenseVector.OfArray(new double[] { function.xDerivative(currentPoint), function.yDerivative(currentPoint) }) * -1;
                //поиск оптимального шага(применение одномерной оптимизации в направлении антиградиента)
                Func<float, float> antigradientDirection = (float l) =>
                { return (float)function.function(currentPoint + (l * (antigradient / Math.Sqrt(Math.Pow(antigradient.At(0), 2) + Math.Pow(antigradient.At(1), 2))))); };
                var step = getMinimum(antigradientDirection);
                //прибавляем шаг к текущей точке
                currentPoint += (antigradient / Math.Sqrt(Math.Pow(antigradient.At(0), 2) + Math.Pow(antigradient.At(1), 2))) * step;
            } while (Math.Sqrt(Math.Pow(antigradient.At(0), 2) + Math.Pow(antigradient.At(1), 2)) >= 0.5); //в точке, близкой к екстремуму, длинна вектора антиградиента будет стремиться к нолю
            return currentPoint;
        }

        //метод Бройдена-Флечера-Голдфарба-Шанно
        public Vector BFGS()
        {
            var currentPoint = DenseVector.OfArray(new double[] { 12, -42 });
            var antigradient = DenseVector.OfArray(new double[] { 0, 0 });
            do
            {
                var identityMatrix = DenseMatrix.OfArray(new double[,] {
                {1 , 1 }, {0 , 1} });
                //Задаём начальное приближение к антигессиану ф-ции
                var C = identityMatrix;
                //поиск антиградиента в точке
                var diff = antigradient - DenseVector.OfArray(new double[]
                {function.xDerivative(currentPoint), function.yDerivative(currentPoint) }) * -1;
                antigradient = DenseVector.OfArray(new double[]
                {function.xDerivative(currentPoint), function.yDerivative(currentPoint) }) * -1;
                //умножаем антиградиент на матрицу
                antigradient *= C;
                //поиск оптимального шага(применение одномерной оптимизации в направлении антиградиента)
                Func<float, float> antigradientDirection = (float l) =>
                {
                    return (float)function.function(
                        currentPoint + (l * (antigradient / (Math.Sqrt(Math.Pow(antigradient.At(0), 2) +
                        Math.Pow(antigradient.At(1), 2))))));
                };
                var step = getMinimum(antigradientDirection);
                //прибавляем шаг к текущей точке
                currentPoint += (antigradient / Math.Sqrt(Math.Pow(antigradient.At(0), 2) + Math.Pow(antigradient.At(1), 2))) * step;
                //обновняем матрицу
                var temp = diff * antigradient * (Vector.OuterProduct(antigradient * step, diff));
            } while ((Math.Sqrt(Math.Pow(antigradient.At(0), 2) + Math.Pow(antigradient.At(1), 2))) >= 0.5); //в точке, близкой к екстремуму, длинна вектора антиградиента будет стремиться к нолю
            return currentPoint;
        }


        //одномерная оптимизация
        private float getMinimum(Func<float, float> function)
        {
            float currentPoint = 0;
            float step = 100;
            while (Math.Abs(step) > 0.005)
            {
                float a = function(currentPoint);
                float b = function(currentPoint + step);
                float c = function(currentPoint - step);
                if (a < b && a < c)
                {
                    step = step / 2;
                }
                else if (a < b)
                {
                    step = -step;
                }
                else if (a == b && a == c)
                {
                    return currentPoint;
                }
                currentPoint += step;
            }
            return currentPoint;
            

        }
    }
}

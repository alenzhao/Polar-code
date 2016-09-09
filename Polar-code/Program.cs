using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Polar_code
{
    class Program
    {
        static void Main(string[] args)
        {
            Matrix<double> A = DenseMatrix.OfArray(new double[,]{
                {1, 0 },
                {1, 1 } });


            int N = 16;
            int K = 5;
            int snr = 0;
            PolarCode pc1 = new PolarCode();
            PolarCode pc2 = new PolarCode(N, K, snr);

            Console.WriteLine(pc1.Info());
            Console.WriteLine(pc2.Info());

            Console.WriteLine("fin!");
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Polar_code
{
    /*
    Polar Code with Arikan core
    based on works by Harish Vangala, Emmanuele Viterbo and Yi Hong
    http://www.ecse.monash.edu.au/staff/eviterbo/polarcodes.html
    */
    class PolarCode
    {
        //Arikan core
        static Matrix<double> ArikanCore = DenseMatrix.OfArray(new double[,]{
                {1, 0 },
                {1, 1 } });

        //---data fields
        //codeword Length
        private int N;
        //polarization power
        private int n;
        //infword Length
        private int K;

        //---properties
        //codeword Length
        public int codewordLength
        {
            get { return N; }
            set
            {
                if (Math.Pow(ArikanCore.ColumnCount, Math.Log(value, ArikanCore.ColumnCount)) == value)
                    N = value;
                else
                    Console.WriteLine("ERROR! Codeword length N is not a power of 2!");
            }
        }
        //infword Length
        public int infwordLength
        {
            get { return K; }
            set
            {
                if (value < N)
                    K = value;
                else
                    Console.WriteLine("ERROR! Infword length K must be lower than Codeword length N!");
            }
        }
        //design signal to noise ratio
        public double designSNR { get; set; }
        //frozen bits pozitions
        private List<int> frozenBits;

        //---constructors
        public PolarCode()
            : this(256, 128, 0)
        {
            frozenBits = GetFrozenBits(256, 128, 0);
            n = (int)Math.Log(256, 2);
        }
        public PolarCode(int cword, int iword, double dsnr)
        {
            codewordLength = cword;
            infwordLength = iword;
            designSNR = dsnr;
            frozenBits = GetFrozenBits(cword, iword, dsnr);
            n = (int)Math.Log(cword, 2);
        }


        //---methods
        /*
        PCC-0 Algorythm, based on The Bhattacharyya bounds
        from "Comparative study of polar code constructions for the AWGN channel"
        */
        private List<int> GetFrozenBits(int N, int K, double dSNR)
        {
            double S = Math.Pow(10, dSNR / 10);
            double n = Math.Log(N, 2);
            double[] Z0 = new double[N];
            Z0[0] = Math.Pow(2.7182818284590452354, -S);

            for (int j = 1; j <= n; j++)
            {
                int u = (int)Math.Pow(2, j);
                for (int t = 0; t <= (u / 2 - 1); t++)
                {
                    double T = Z0[t];
                    Z0[t] = 2 * T + Math.Pow(T, 2);
                    Z0[u / 2 + t] = Math.Pow(T, 2);
                }
            }
            List<double> z0 = Z0.ToList();
            List<int> frozenbits = MaxIndexes(z0, N, K);
            return frozenbits;
        }
        //internal method of GetFrozenBits method
        private List<int> MaxIndexes(List<double> Z0, int N, int K)
        {
            List<int> indexes = new List<int>();
            Dictionary<int, double> bitParams = new Dictionary<int, double>();
            for (int i = 0; i < Z0.Count; i++)
            {
                bitParams.Add(i, Z0[i]);
            }
            bitParams = bitParams.OrderByDescending(s => s.Value).ToDictionary(s => s.Key, s => s.Value);
            indexes = bitParams.Keys.ToList();
            List<int> fbits = new List<int>();
            for (int i = 0; i < N - K; i++)
            {
                fbits.Add(indexes[i]);
            }
            fbits.Sort();
            return fbits;
        }

        public void ShowFrozenBits()
        {
            Console.WriteLine(string.Join(",", frozenBits.ToArray()));
        }

        public string Info()
        {
            StringBuilder info = new StringBuilder();

            info.AppendFormat("N = {0}\n", this.codewordLength);
            info.AppendFormat("K = {0}\n", this.infwordLength);
            info.AppendFormat("R = {0}\n", this.codewordLength / this.infwordLength);
            info.AppendFormat("Desing SNR = {0}\n", this.designSNR);
            info.AppendFormat(string.Join(",", frozenBits.ToArray()));
            info.AppendFormat("\n");

            return info.ToString();
        }
        
        /*
        Systematic encoder based on EncoderA algorythm 
        from "Efficient algorythms for systematic polar encoding"
        */
        public List<double> SysEncode(List<double> inf)
        {
            if (inf.Count != this.K)
                Console.WriteLine("ERROR! Incorect size of information word!");
            double[,] X = new double[this.N, this.n + 1];
            /*
            fill first column of X with symbols of inf[] at unfrozen positions
            and with zeros at frozen positions
            */
            int count = 0;
            for (int i = 0; i < this.N; i++)
            {
                if (this.frozenBits.Contains(i))
                {
                    X[i, this.n] = 0;
                }
                else
                {
                    X[i, this.n] = inf[count];
                    count++;
                }
            }
            /*
            encoding procedure
            */
            for(int i = this.N - 1; i >= 0; i--)
            {
                int s; int delta;
                if (this.frozenBits.Contains(i))
                {
                    s = 1; delta = 1;
                }
                else
                {
                    s = this.n + 1; delta = -1;
                }
                int[] b = GetBinaryRepresentation(i);
                for (int j = 0; j < this.n; j++)
                {
                    int t = s + (j + 1) * delta;
                    int l = Math.Min(t, t - delta);
                    int k = (int)Math.Pow(2, this.n - l);
                    if (b[l-1] == 0)
                    {
                        X[i, t] = (X[i, t - delta] + X[i + k, t - delta]) % 2;
                    }
                    else
                    {
                        X[i, t] = X[i, t - delta];
                    }
                }
            }
            List<double> txword = new List<double>();
            for (int i = 0; i < this.N; i++)
            {
                txword.Add(X[i,this.n]);
            }
            return txword;
        }
        //internal method for SysEncode method
        //returns binary representation of size N
        private int[] GetBinaryRepresentation(int i)
        {
            int [] result = new int [this.N];
            int count = 0;
            while (i > 0)
            {
                result[count] = i % 2;
                i = i / 2;
                count++;
            }
            return result;
        }

        public List<double> Decode(List<double> coded)
        {
            if (coded.Count != N)
                Console.WriteLine("ERROR! Incorrect length of coded word!");


            List<double> rxword = new List<double>();
            return rxword;
        }
    }
}

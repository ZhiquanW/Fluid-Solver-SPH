using System;

namespace FluidSolver {
    public class Vector {


        private double _x;
        private double _y;
        private double _z;

        public Vector(double x = 0, double y = 0, double z = 0) {
            _x = x;
            _y = y;
            _z = z;
        }

        public override string ToString() {
            return "(" + _x + "," + _y + "," + _z + ")";
        }

        public double X {
            get => _x;
            set => _x = value;
        }

        public double Y {
            get => _y;
            set => _y = value;
        }

        public double Z {
            get => _z;
            set => _z = value;
        }

        public static double DistanceSquare(Vector v1, Vector v2) {
            return Math.Pow(v1.X - v2.X, 2) + Math.Pow(v1.Y - v2.Y, 2) + Math.Pow(v1.Z - v2.Z, 2);
        }
    }
}
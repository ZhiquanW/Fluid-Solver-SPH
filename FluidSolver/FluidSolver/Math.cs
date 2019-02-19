//using System;
//
//namespace FluidSolver {
//    public class Vector {
//
//
//        private float _x;
//        private float _y;
//        private float _z;
//
//        public Vector(float x = 0, float y = 0, float z = 0) {
//            _x = x;
//            _y = y;
//            _z = z;
//        }
//
//        public override string ToString() {
//            return "(" + _x + "," + _y + "," + _z + ")";
//        }
//
//        public float X {
//            get => _x;
//            set => _x = value;
//        }
//
//        public float Y {
//            get => _y;
//            set => _y = value;
//        }
//
//        public float Z {
//            get => _z;
//            set => _z = value;
//        }
//
//        public static float DistanceSquare(Vector v1, Vector v2) {
//            return Math.Pow(v1.X - v2.X, 2) + Math.Pow(v1.Y - v2.Y, 2) + Math.Pow(v1.Z - v2.Z, 2);
//        }
//    }
//}
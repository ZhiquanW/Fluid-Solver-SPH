using System;
using System.Numerics;

namespace FluidSolver {
    static class Program {
        static void Main(string[] args) {
            Console.WriteLine("Start Generate Particles");
            var fluidSolver = new FluidSolver(new Vector3(100,100,100),0.01f);
            fluidSolver.InitParticles(new Vector3(10, 10, 10), 0.01f, 1);
            fluidSolver.StartImitation(600,0.1f,"test1.fd");
        }
    }
}
using System;

namespace FluidSolver {
    static class Program {
        static void Main(string[] args) {

            Console.WriteLine("Start Generate Particles");
            var tmpRestriction = new Vector(100, 100, 50);
            var fluidSolver = new FluidSolver(new Vector(100,100,100));
            fluidSolver.InitParticles(1000,1,tmpRestriction);
            //fluidSolver.TestInitParticles();
        }
    }
    
}
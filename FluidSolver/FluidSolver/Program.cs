using System;

namespace FluidSolver {
    static class Program {
        static void Main(string[] args) {

            Console.WriteLine("Start Generate Particles");
            var fluidSolver = new FluidSolver(new Vector(100,100,100),30);
            
            //fluidSolver.TestInitParticles();
        }
    }
    
}
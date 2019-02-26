using System;
using System.Numerics;

namespace FluidSolver {
    static class Program {
        static void Main(string[] args) {
            Console.WriteLine("Start Generate Particles");
            var fluidSolver = new FluidSolver(new Vector3(100, 100, 100), 0.01f,
                1000.0f, 1.0f, 1.0f,
                0.075f, 9.8f);
            fluidSolver.InitParticles(new Vector3(10, 10, 10), 0.01f, 0.0004f);
            fluidSolver.TestComputeDensity();
            fluidSolver.StartImitation(600, 0.1f,
                @"/home/zhiquan/Git-Repository/Fluid-Solver-SPH/SPH-Solver-Renderer/Assets/Scripts/test0.fb
            ");
        }
    }
}
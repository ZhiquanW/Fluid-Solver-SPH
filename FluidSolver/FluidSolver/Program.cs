using System;
using System.Collections.Generic;
using System.Numerics;

namespace FluidSolver {
    static class Program {
        static void Main(string[] args) {
//            List<Particle> tmpList = new List<Particle>();
//            tmpList.Add(new Particle(1, 0, new Vector3(0, 0, 0)));
//            tmpList.Add(new Particle(2, 0, new Vector3(0, 0, 0)));
//            tmpList.Add(new Particle(3, 0, new Vector3(0, 0, 0)));
//            tmpList[0].Index = -1;
//            List<Particle> tmpList0 = new List<Particle>();
//            tmpList.ForEach(p => tmpList0.Add(new Particle(p)));
//            Console.WriteLine(tmpList0[0].Index);
//            Console.WriteLine(tmpList[0].Index);
            Console.WriteLine("Start Generate Particles");
            var fluidSolver = new FluidSolver(new Vector3(100, 100, 100), 0.01f,
                1000.0f, 1.0f, 1.0f,
                0.075f, 9.8f,600,0.1f);
           
            fluidSolver.StartImitation(
                @"/home/zhiquan/Git-Repository/Fluid-Solver-SPH/SPH-Solver-Renderer/Assets/Scripts/test0.fb");
            
        }
    }
}
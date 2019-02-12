using System;
using Math;
namespace FluidSolver {
    public class FluidSolver {
        private int particleNum;
        private double particleMass; //set particle
        private double restDensity; // set particle \rho_0
        private double viscosityCoefficient;
        private double coreRadius; // h
        private double gasConstant; // k 
        private double tensionCoefficient; //sigma
        private double gravityAcceleration; //g

        public int ParticleNum { get => particleNum; set => particleNum = value; }
        public double ParticleMass { get => particleMass; set => particleMass = value; }
        public double RestDensity { get => restDensity; set => restDensity = value; }

        public double CoreRadius { get => coreRadius; set => coreRadius = value; }
        public double GasConstant { get => gasConstant; set => gasConstant = value; }
        public double TensionCoefficient { get => tensionCoefficient; set => tensionCoefficient = value; }
    }

    public class Particle {
        private double mass;
        private double density;
        private double restDensity;
        private double pressure;
        private Vector position;
        private Vector velocity;
        private Vector acceleration;
        private Vector pressureForce;
        private Vector viscosityForce;

        public double Mass { get => mass; set => mass = value; }
        public double Density { get => density; set => density = value; }
        public double RestDensity { get => restDensity; set => restDensity = value; }

        public double Pressure { get => pressure; set => pressure = value; }
        internal Vector Position { get => position; set => position = value; }
        internal Vector Velocity { get => velocity; set => velocity = value; }
        internal Vector Acceleration { get => acceleration; set => acceleration = value; }
        internal Vector PressureForce { get => pressureForce; set => pressureForce = value; }
        internal Vector ViscosityForce { get => viscosityForce; set => viscosityForce = value; }

        Particle (){

        }
        Particle(double _mass,double _rest_density){
            
        }
    }

}
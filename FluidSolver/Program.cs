using System;
using Math;
namespace FluidSolver {
    public class FluidSolver {
        private List<Particle> particleList;
        private int particleNum;
        private double particleMass; //set particle
        private double restDensity; // set particle \rho_0
        private double viscosityCoefficient; // mu
        private double coreRadius; // h
        private double tensionCoefficient; //sigma
        private double gravityAcceleration; //g
        public int ParticleNum { get => particleNum; set => particleNum = value; }
        public double ParticleMass { get => particleMass; set => particleMass = value; }
        public double RestDensity { get => restDensity; set => restDensity = value; }
        public double ViscosityCoefficient { get => viscosityCoefficient; set => viscosityCoefficient = value; }
        public double CoreRadius { get => coreRadius; set => coreRadius = value; }
        public double GasConstant { get => gasConstant; set => gasConstant = value; }
        public double TensionCoefficient { get => tensionCoefficient; set => tensionCoefficient = value; }
        public double GravityAcceleration { get => gravityAcceleration; set => gravityAcceleration = value; }

        public FluidSolver(){

        }

        public void initParticles(){

        }
    }

    public class Particle {
        private double mass;
        private Vector position;
        private double density;
        private double restDensity;
        private double gasConstant; // k 
        private double pressure;
        private Vector velocity;
        private Vector acceleration;
        private Vector pressureForce;
        private Vector viscosityForce;

        public double Mass { get => mass; set => mass = value; }
        public double Density { get => density; set => density = value; }
        public double RestDensity { get => restDensity; set => restDensity = value; }
        public double Pressure { get => pressure; set => pressure = value; }
        public Vector Position { get => position; set => position = value; }
        public Vector Velocity { get => velocity; set => velocity = value; }
        public Vector Acceleration { get => acceleration; set => acceleration = value; }
        public Vector PressureForce { get => pressureForce; set => pressureForce = value; }
        public Vector ViscosityForce { get => viscosityForce; set => viscosityForce = value; }
        Particle (){

        }
        Particle(double _mass,double _rest_density,Vector _position){
            mass = _mass;
            restDensity = _rest_density;
            position = _position;
        }

        void calDensity(double _gas_constant,double ){

        }
    }

}
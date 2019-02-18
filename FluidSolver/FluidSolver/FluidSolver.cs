using System;
using System.Collections.Generic;
namespace FluidSolver {
    public class FluidSolver {
        private List<Particle> _particleList;
        private Vector _containerVolume;
        private int _particleNum;
        private double _particleMass; //set particle
        private double _restDensity; // set particle \rho_0
        private double _viscosityCoefficient; // mu
        private double _coreRadius; // h
        private double _gasConstant; // k 
        private double _tensionCoefficient; //sigma
        private double _gravityAcceleration; //g
       
        public FluidSolver(Vector inContainerVolume) {
            this._containerVolume = inContainerVolume;
            this._particleList = new List<Particle>();
            this._particleNum = 0;
        }

        public void InitParticles(int inParticleNum,double inParticleMass,Vector inPosRestriction){
            double posInterval = Math.Pow(inPosRestriction.X * inPosRestriction.Y * inPosRestriction.Z/inParticleNum,1/3f);
            int xNum = (int) (inPosRestriction.X / posInterval);
            int yNum = (int) (inPosRestriction.Y / posInterval);
            int zNum = (int) (inPosRestriction.Z / posInterval);
            this._particleNum += xNum * yNum * zNum;
            for (int i = 0; i < xNum; ++i) {
                for (int j = 0; j < yNum; ++j) {
                    for (int k = 0; k < zNum; ++k) {
                        var tmpPos = new Vector(posInterval*i,posInterval*j,posInterval*k);
                        this._particleList.Add(new Particle(inParticleMass,tmpPos));
                    }
                }
            }
        }

        public void TestInitParticles() {
            Console.WriteLine(this._particleList.Count);
            int counter = 0;
            foreach (var p in this._particleList) {
                Console.WriteLine(counter ++ +" "+p.Position.ToString());
            }
        }
    }

    public class Particle {
        
        private double mass;
        private Vector position;
        private double density;
        private double pressure;
        private Vector velocity;
        private Vector acceleration;
        private Vector pressureForce;
        private Vector viscosityForce;

        
        public Particle (){

        }
        
        public Particle(double inMass,Vector inPosition){
            this.mass = inMass;
            this.position = inPosition;
        }


        void CalDensity(double inGasConstant,double inRestDensity){

        }


        public double Mass {
            get => mass;
            set => mass = value;
        }

        public Vector Position {
            get => position;
            set => position = value;
        }

        public double Density {
            get => density;
            set => density = value;
        }

        public double Pressure {
            get => pressure;
            set => pressure = value;
        }

        public Vector Velocity {
            get => velocity;
            set => velocity = value;
        }

        public Vector Acceleration {
            get => acceleration;
            set => acceleration = value;
        }

        public Vector PressureForce {
            get => pressureForce;
            set => pressureForce = value;
        }

        public Vector ViscosityForce {
            get => viscosityForce;
            set => viscosityForce = value;
        }
    }
}
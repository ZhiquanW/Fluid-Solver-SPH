using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Net.NetworkInformation;
using System.Numerics;

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
       
        public FluidSolver(Vector inContainerVolume,double inCoreRadius) {
            this._containerVolume = inContainerVolume;
            this._particleList = new List<Particle>();
            this._particleNum = 0;
            this._coreRadius = inCoreRadius;
            var tmpRestriction = new Vector(100, 100, 50);
            InitParticles(1000,1,tmpRestriction);
            CalDensity();
            TestCalDensity();
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

        private void CalDensity() {
            double hSquare = this._coreRadius * this._coreRadius;
            double tmpCoefficient = 315 / (64 * Math.PI * Math.Pow(hSquare, 9));
            foreach (var pI in this._particleList) {
                double tmpDensity = 0;
                foreach (var pJ in this._particleList) {
                    double tmpDisSqu = Vector.DistanceSquare(pI.Position, pJ.Position);
                    if (tmpDisSqu < hSquare) {
                        tmpDensity += Math.Pow(hSquare - tmpDisSqu, 3);
                    }
                }

                tmpDensity *= pI.Mass * tmpCoefficient;
                pI.Density = tmpDensity;
            }
        }

        public void TestInitParticles() {
            Console.WriteLine(this._particleList.Count);
            int counter = 0;
            foreach (var p in this._particleList) {
                Console.WriteLine(counter ++ +" "+p.Position.ToString());
            }
        }

        public void TestCalDensity() {
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Density);
            }
        }
    }

    public class Particle {
        
        private double _mass;
        private Vector _position;
        private double _density;
        private double _pressure;
        private Vector _velocity;
        private Vector _acceleration;
        private Vector _pressureForce;
        private Vector _viscosityForce;

        
        public Particle (){

        }
        
        public Particle(double inMass,Vector inPosition){
            this._mass = inMass;
            this._position = inPosition;
            this._density = 0;
        }


        public double Mass {
            get => _mass;
            set => _mass = value;
        }

        public Vector Position {
            get => _position;
            set => _position = value;
        }

        public double Density {
            get => _density;
            set => _density = value;
        }

        public double Pressure {
            get => _pressure;
            set => _pressure = value;
        }

        public Vector Velocity {
            get => _velocity;
            set => _velocity = value;
        }

        public Vector Acceleration {
            get => _acceleration;
            set => _acceleration = value;
        }

        public Vector PressureForce {
            get => _pressureForce;
            set => _pressureForce = value;
        }

        public Vector ViscosityForce {
            get => _viscosityForce;
            set => _viscosityForce = value;
        }
    }
}
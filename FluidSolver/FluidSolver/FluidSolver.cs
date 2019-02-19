using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Net.NetworkInformation;
using System.Numerics;

namespace FluidSolver {
    public class FluidSolver {
        private List<Particle> _particleList;
        private Vector3 _containerVolume;
        private int _particleNum;
        private float _particleMass; //set particle
        private float _restDensity; // set particle \rho_0
        private float _viscosityCoefficient; // mu
        private readonly float _coreRadius; // h
        private float _gasConstant; // k 
        private float _tensionCoefficient; //sigma
        private float _gravityAcceleration; //g
       
        public FluidSolver(Vector3 inContainerVolume,float inCoreRadius) {
            this._containerVolume = inContainerVolume;
            this._particleList = new List<Particle>();
            this._particleNum = 0;
            this._coreRadius = inCoreRadius;
            var tmpRestriction = new Vector3(100, 100, 50);
            InitParticles(1000,1,tmpRestriction);
            TestInitParticles();
//            ComputeDensity();
//            TestCalDensity();
        }

        public void InitParticles(int inParticleNum,float inParticleMass,Vector3 inPosRestriction){
            float posInterval = (float) Math.Pow(inPosRestriction.X * inPosRestriction.Y * inPosRestriction.Z/inParticleNum,1/3f);
            int xNum = (int) (inPosRestriction.X / posInterval);
            int yNum = (int) (inPosRestriction.Y / posInterval);
            int zNum = (int) (inPosRestriction.Z / posInterval);
            this._particleNum += xNum * yNum * zNum;
            int counter = 0;
            for (int i = 0; i < xNum; ++i) {
                for (int j = 0; j < yNum; ++j) {
                    for (int k = 0; k < zNum; ++k) {
                        var tmpPos = new Vector3(posInterval*i,posInterval*j,posInterval*k);
                        this._particleList.Add(new Particle(counter++,inParticleMass,tmpPos));
                    }
                }
            }


        }

        private void ComputeDensity() {
            float hSquare = this._coreRadius * this._coreRadius;
            float tmpCoefficient = (float) (315 / (64 * Math.PI * Math.Pow(hSquare, 9)));
            foreach (var pI in this._particleList) {
                float tmpDensity = 0;
                foreach (var pJ in this._particleList) {
                    float tmpDisSqu = Vector3.DistanceSquared(pI.Position, pJ.Position);
                    if (tmpDisSqu < hSquare) {
                        tmpDensity += (float )Math.Pow(hSquare - tmpDisSqu, 3);
                    }
                }

                tmpDensity *= pI.Mass * tmpCoefficient;
                pI.Density = tmpDensity;
            }
        }

        private void ComputePressure() {
            foreach (var p in this._particleList) {
                p.Pressure = this._gasConstant * (p.Density - _restDensity);
            }
        }

        private void ComputePressureForceAcceleration() {
            float tmpCoefficient =  (float) (45 / (Math.PI * Math.Pow(this._coreRadius, 6)));
            foreach (var pI in this._particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in this._particleList) {
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < this._coreRadius) {
                        tmpVec += (float) ((pI.Pressure + pJ.Pressure) / (2 * pI.Density * pJ.Density) *
                                           Math.Pow(this._coreRadius - tmpDis, 2)) *
                                  Vector3.Normalize(pI.Position - pJ.Position);
                    }
                }
                tmpVec *= pI.Mass * tmpCoefficient * tmpVec;
            }
        }
        public void TestInitParticles() {
            Console.WriteLine(this._particleList.Count);
            int counter = 0;
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Index+" "+p.Position.ToString());
            }
        }

        public void TestCalDensity() {
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Density);
            }
        }
    }

    public class Particle {
        private int _index = 0;
        private float _mass = 0;
        private Vector3 _position = new Vector3();
        private float _density = 0;
        private float _pressure = 0;
        private Vector3 _velocity = new Vector3();
        private Vector3 _totalAcceleration = new Vector3();
        private Vector3 _pressureForceAcceleration = new Vector3();
        private Vector3 _viscosityForceAcceleration = new Vector3();
        
        public Particle (){

        }
        
        public Particle(int inIndex, float inMass,Vector3 inPosition) {
            this.Index = inIndex;
            this.Mass = inMass;
            this._position = inPosition;
        }


        public int Index {
            get => _index;
            set => _index = value;
        }

        public float Mass {
            get => _mass;
            set => _mass = value;
        }

        public Vector3 Position {
            get => _position;
            set => _position = value;
        }

        public float Density {
            get => _density;
            set => _density = value;
        }

        public float Pressure {
            get => _pressure;
            set => _pressure = value;
        }

        public Vector3 Velocity {
            get => _velocity;
            set => _velocity = value;
        }

        public Vector3 TotalAcceleration {
            get => _totalAcceleration;
            set => _totalAcceleration = value;
        }

        public Vector3 PressureForceAcceleration {
            get => _pressureForceAcceleration;
            set => _pressureForceAcceleration = value;
        }

        public Vector3 ViscosityForceAcceleration {
            get => _viscosityForceAcceleration;
            set => _viscosityForceAcceleration = value;
        }
 
    }
}
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Net.NetworkInformation;
using System.Numerics;
using System.Reflection.Metadata;

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

        public FluidSolver(Vector3 inContainerVolume, float inCoreRadius) {
            this._containerVolume = inContainerVolume;
            this._particleList = new List<Particle>();
            this._particleNum = 0;
            this._coreRadius = inCoreRadius;
            var tmpRestriction = new Vector3(100, 100, 50);
            InitParticles(1000, 1, tmpRestriction);
            TestInitParticles();
//            ComputeDensity();
//            TestCalDensity();
        }

        public void InitParticles(int inParticleNum, float inParticleMass, Vector3 inPosRestriction) {
            float posInterval =
                (float) Math.Pow(inPosRestriction.X * inPosRestriction.Y * inPosRestriction.Z / inParticleNum, 1 / 3f);
            int xNum = (int) (inPosRestriction.X / posInterval);
            int yNum = (int) (inPosRestriction.Y / posInterval);
            int zNum = (int) (inPosRestriction.Z / posInterval);
            this._particleNum += xNum * yNum * zNum;
            int counter = 0;
            for (int i = 0; i < xNum; ++i) {
                for (int j = 0; j < yNum; ++j) {
                    for (int k = 0; k < zNum; ++k) {
                        var tmpPos = new Vector3(posInterval * i, posInterval * j, posInterval * k);
                        this._particleList.Add(new Particle(counter++, inParticleMass, tmpPos));
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
                        tmpDensity += (float) Math.Pow(hSquare - tmpDisSqu, 3);
                    }
                }

                pI.Density = pI.Mass * tmpCoefficient * tmpDensity;
            }
        }

        private void ComputePressure() {
            foreach (var p in this._particleList) {
                p.Pressure = this._gasConstant * (p.Density - _restDensity);
            }
        }

        private void ComputePressureForceAcceleration() {
            float tmpCoefficient = (float) (45 / (Math.PI * Math.Pow(this._coreRadius, 6)));
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

                pI.PressureForceAcceleration = pI.Mass * tmpCoefficient * tmpVec;
            }
        }

        private void ComputeViscosityForceAcceleration() {
            float tmpCoefficient =
                (float) (this._viscosityCoefficient * 45 / (Math.PI * Math.Pow(this._coreRadius, 6)));
            foreach (var pI in this._particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in this._particleList) {
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < this._coreRadius) {
                        tmpVec += (pI.Velocity - pJ.Velocity) / (pI.Density * pJ.Density) * this._coreRadius * tmpDis;
                    }
                }

                pI.ViscosityForceAcceleration = pI.Mass * tmpCoefficient * tmpVec;
            }
        }

        private void ComputeGravityAcceleration() {
            foreach (var pI in this._particleList) {
                pI.GravityAcceleration = this._gravityAcceleration * new Vector3(0, 0, -1);
            }
        }

        private void ComputeSurfaceTensionAcceleration() {
            float tmpCoefficient = (float) (-945 / (32 * Math.PI * Math.Pow(this._coreRadius, 9)));
            float tmpCoefficientLaplacian = (float) (945 / (8 * Math.PI * Math.Pow(this._coreRadius, 9)));
            float hSquared = (float) Math.Pow(this._coreRadius, 2);
            foreach (var pI in this._particleList) {
                Vector3 colorFieldGradient = new Vector3();
                float colorFieldLaplacian = 0;
                foreach (var pJ in this._particleList) {
                    float rSquared = Vector3.DistanceSquared(pI.Position, pJ.Position);
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < this._coreRadius) {
                        colorFieldGradient += (float) (1 / pI.Density * Math.Pow(hSquared - rSquared, 2)) *
                                              (pI.Position - pJ.Position);
                        colorFieldLaplacian += 1 / pI.Density * (hSquared - rSquared) *
                                               (rSquared - 3 / 4 * (hSquared - rSquared));
                    }
                }

                pI.SurfaceTensionAcceleration = -1 * _tensionCoefficient * pI.Mass * tmpCoefficientLaplacian /
                                                pI.Density * Vector3.Normalize(colorFieldGradient);
            }
        }

        private void ComputeTotalAcceleration() {
            foreach (var pI in this._particleList) {
                pI.TotalAcceleration = pI.GravityAcceleration + pI.PressureForceAcceleration +
                                       pI.ViscosityForceAcceleration + pI.SurfaceTensionAcceleration;
            }
        }

        public void TestInitParticles() {
            Console.WriteLine(this._particleList.Count);
            int counter = 0;
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Index + " " + p.Position.ToString());
            }
        }

        public void TestCalDensity() {
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Density);
            }
        }
    }

    public class Particle {
        public Particle() {

        }

        public Particle(int inIndex, float inMass, Vector3 inPosition) {
            this.Index = inIndex;
            this.Mass = inMass;
            this.Position = inPosition;
        }


        public int Index { get; set; } = 0;

        public float Mass { get; set; } = 0;

        public Vector3 Position { get; set; } = new Vector3();

        public float Density { get; set; } = 0;

        public float Pressure { get; set; } = 0;

        public Vector3 Velocity { get; set; } = new Vector3();

        public Vector3 TotalAcceleration { get; set; } = new Vector3();

        public Vector3 PressureForceAcceleration { get; set; } = new Vector3();

        public Vector3 ViscosityForceAcceleration { get; set; } = new Vector3();

        public Vector3 GravityAcceleration { get; set; } = new Vector3();

        public Vector3 SurfaceTensionAcceleration { get; set; } = new Vector3();
    }
}
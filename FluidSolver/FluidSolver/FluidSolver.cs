using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Net.NetworkInformation;
using System.Numerics;
using System.Reflection.Metadata;

namespace FluidSolver {
    public class FluidSolver {
        private FluidDatabase _fluidDatabase;
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
            this._fluidDatabase = new FluidDatabase();
        }

        public void StartImitation(float inFrameNum, float inTimeStep, string inFileName) {
            this._fluidDatabase.TimeStep = inTimeStep;
            this._fluidDatabase.TimeDuration = inFrameNum * inTimeStep;
            DateTime startTime = DateTime.Now;
            DateTime endTime = new DateTime();
            for (int i = 0; i < inFrameNum; ++i) {
                Console.SetCursorPosition(0, 1);
                ComputeDensity();
                ComputePressure();
                ComputeGravityAcceleration();
                ComputePressureForceAcceleration();
                ComputeViscosityForceAcceleration();
                ComputeSurfaceTensionAcceleration();
                ComputeTotalAcceleration();
                UpdatePosition(inTimeStep);
                endTime = DateTime.Now;
                float runtime = (float) endTime.Subtract(startTime).TotalMilliseconds / 1000;
                float progressPercentage = (float) (i + 1) / inFrameNum * 100;
                float remainingTime = progressPercentage / runtime * (100 - progressPercentage);
                Console.WriteLine("Progress: " + progressPercentage.ToString("F3") + "% "
                                  + "Runtime: " + runtime.ToString("F3") + "s "
                                  + "EstimatedTime: " + remainingTime.ToString("F3") + "s ");
            }

            this._fluidDatabase.Output(inFileName);
        }

        public void InitParticles(Vector3 inParticleVolume, float inParticleInterval, float inParticleMass) {
            this._particleNum = (int) (inParticleVolume.X * inParticleVolume.Y * inParticleVolume.Z);
            int counter = 0;
            for (int i = 0; i < inParticleVolume.X; ++i) {
                for (int j = 0; j < inParticleVolume.Y; ++j) {
                    for (int k = 0; k < inParticleVolume.Z; ++k) {
                        var tmpPos = new Vector3(inParticleInterval * i, inParticleInterval * j,
                            inParticleInterval * k);
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

        private float IntersectionDetection(Vector3 inVecPos, Vector3 inVecDir, Vector3 inPlanePoint,
            Vector3 inPlaneNormal) {
            float vecNDotPlaneN = Vector3.Dot(inVecDir, inPlaneNormal);
            if (Math.Abs(vecNDotPlaneN) <= TOLERANCE) {
                return 0f;
            } else {
                float tmpValue = Vector3.Dot(inVecDir, inPlanePoint - inVecPos);
                return tmpValue / vecNDotPlaneN;
            }
        }

        private Vector3 Reflection(Vector3 inVec, Vector3 inNormal) {
            return inVec - 2 * Vector3.Dot(inVec, inNormal) * inNormal;
        }
        private void RestrictionParticles() {
            foreach (var pI in this._particleList) {
                bool isOutside = !(0 < pI.Position.X && pI.Position.X < this._containerVolume.X &&
                                   0 < pI.Position.Y && pI.Position.Y < this._containerVolume.Y &&
                                   0 < pI.Position.Z && pI.Position.Z < this._containerVolume.Z);
                if (isOutside) {
                    var vecDir = -Vector3.Normalize(pI.Velocity);
                    float disTopBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(this._containerVolume.X, this._containerVolume.Y, this._containerVolume.Z),
                        new Vector3(0, 1, 0));
                    float disFrontBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(this._containerVolume.X, this._containerVolume.Y, this._containerVolume.Z),
                        new Vector3(1, 0, 0));
                    float disRightBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(this._containerVolume.X, this._containerVolume.Y, this._containerVolume.Z),
                        new Vector3(0, 0, 1));
                    float disBottomBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(0, 0, 0), new Vector3(0, -1, 0));
                    float disBackBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(0, 0, 0), new Vector3(-1, 0, 0));
                    float disLeftBorder = IntersectionDetection(pI.Position, vecDir,
                        new Vector3(0, 0, 0), new Vector3(0, 0, -1));
                    float[] disArr =
                        {disTopBorder, disBottomBorder, disFrontBorder, disBackBorder, disRightBorder, disLeftBorder};
                    float tmpMin = Math.Abs(disTopBorder);
                    float tmpIndex = 0;
                    float counter = 0;
                    foreach (var val in disArr) {
                        ++counter;
                        if (val > 0 && val < tmpMin) {
                            tmpMin = val;
                            tmpIndex = counter;
                        }
                    }
                    Vector3[] normalArr = 

                }
            }
        }

        private void UpdatePosition(float inTimeStep) {
            foreach (var pI in this._particleList) {
                pI.Velocity += pI.TotalAcceleration * inTimeStep;
                pI.Position += pI.Velocity * inTimeStep;
            }

            this._fluidDatabase.AddParticles(this._particleList);
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

    public class FluidDatabase {
        public float TimeStep { get; set; } = 0;

        public float TimeDuration { private get; set; }

        private readonly List<List<Particle>> _fluidMatrix = new List<List<Particle>>();

        public FluidDatabase() {
        }

        public void AddParticles(List<Particle> inParticleList) {
            this._fluidMatrix.Add(inParticleList);
        }

        public void Output(string inFileName) {
            using (System.IO.StreamWriter file =
                new System.IO.StreamWriter(inFileName)) {
                file.WriteLine(TimeDuration.ToString(CultureInfo.InvariantCulture));
                file.WriteLine(TimeStep.ToString(CultureInfo.InvariantCulture));
                int counter = 0;
                foreach (var particleList in this._fluidMatrix) {
                    file.WriteLine(counter.ToString());
                    foreach (var pI in particleList) {
                        file.WriteLine(pI.ToString()); 
                    }
                }
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

        public override string ToString() {
            return Index+ " " +Position.Z+" "+Position.Y+" "+Position.Z;
        }
    }
}
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Net.NetworkInformation;
using System.Numerics;
using System.Reflection.Metadata;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;

namespace FluidSolver {
    public class FluidSolver {
        private FluidDatabase _fluidDatabase;
        private List<Particle> _particleList;
        private Vector3 _containerBox;
        private int _particleNum;
        private readonly float _coreRadius; // h
        private float _particleMass; //set particle
        private float _restDensity; // set particle \rho_0
        private float _viscosityCoefficient; // mu
        private float _gasConstant; // k 
        private float _tensionCoefficient; //sigma
        private float _gravityAcceleration; //g
        private float _frameNum;
        private float _timeStep;
        private float _timeDuration; 
        public FluidSolver(Vector3 inContainerBox, float inCoreRadius,
                           float inRestDensity,float inViscosityCoefficient,
                           float inGasConstant,float inTensionCoefficient,
                           float inGravityAcceleration,float inFrameNum,float inTimeStep) {
            this._particleList = new List<Particle>();
            this._fluidDatabase = new FluidDatabase();
            this._particleNum = 0;
            this._containerBox = inContainerBox; 
            this._coreRadius = inCoreRadius;
            this._restDensity = inRestDensity;
            this._viscosityCoefficient = inViscosityCoefficient;
            this._gasConstant = inGasConstant;
            this._tensionCoefficient = inTensionCoefficient;
            this._gravityAcceleration = inGravityAcceleration;
            this._frameNum = inFrameNum;
            this._timeStep = inTimeStep;
            this._timeDuration = this._frameNum * this._timeStep;
        }

        public void StartImitation(string inFileName) { 
            DisplayParameters();
            this._fluidDatabase.TimeStep = this._timeStep;
            this._fluidDatabase.TimeDuration = this._timeDuration;
            DateTime startTime = DateTime.Now;
            DateTime endTime = new DateTime();
            for (int i = 0; i < this._frameNum; ++i) {
                
                ComputeDensity();
                ComputePressure();
                ComputeGravityAcceleration();
                ComputePressureForceAcceleration();
                ComputeViscosityForceAcceleration();
                ComputeSurfaceTensionAcceleration();
                ComputeTotalAcceleration();
                RestrictionParticles();
                UpdatePosition(this._timeStep);
                
                //Display Computing Progress
                
                endTime = DateTime.Now;
                float runtime = (float) endTime.Subtract(startTime).TotalMilliseconds / 1000;
                float progressPercentage = (float) (i + 1) / this._frameNum * 100;
                float remainingTime = progressPercentage / runtime * (100 - progressPercentage);
                Console.SetCursorPosition(0, Console.CursorTop-1);
                Console.WriteLine("Progress: " + progressPercentage.ToString("F3") + "% "
                                  + "Runtime: " + runtime.ToString("F3") + "s "
                                  + "EstimatedTime: " + remainingTime.ToString("F3") + "s");
            }

            this._fluidDatabase.Output(inFileName);
        }

        public void InitParticles(Vector3 inParticleBox, float inParticleInterval, float inParticleMass) {
            this._particleNum = (int) (inParticleBox.X * inParticleBox.Y * inParticleBox.Z);
            this._particleMass = inParticleMass;
            int counter = 0;
            for (int i = 0; i < inParticleBox.X; ++i) {
                for (int j = 0; j < inParticleBox.Y; ++j) {
                    for (int k = 0; k < inParticleBox.Z; ++k) {
                        var tmpPos = new Vector3(inParticleInterval * i, inParticleInterval * j,
                            inParticleInterval * k);
                        this._particleList.Add(new Particle(counter++, inParticleMass, tmpPos));
                    }
                }
            }
        }
        
        public void TestInitParticles() {
            Console.WriteLine(this._particleList.Count);
            int counter = 0;
            foreach (var p in this._particleList) {
                Console.WriteLine(p.Index + " " + p.Position.ToString());
                Console.WriteLine(p.Velocity);
            }
        }
        
        private void ComputeDensity() {
            float hSquare = this._coreRadius * this._coreRadius;
            float tmpCoefficient = (float) (315 / (64 * Math.PI * Math.Pow(this._coreRadius, 9)));
            List<Particle> preParticleList = this._particleList;
            foreach (var pI in this._particleList) {
                float tmpDensity = 0;
                foreach (var pJ in preParticleList) {
                    float tmpDisSqu = Vector3.DistanceSquared(pI.Position, pJ.Position);
                    if (tmpDisSqu < hSquare) {
                        tmpDensity += (float) Math.Pow(hSquare - tmpDisSqu, 3);
                    }
                }

                pI.Density = pI.Mass * tmpCoefficient * tmpDensity;
//                Console.WriteLine(Vector3.DistanceSquared(pI.Position, pI.Position));
//                Console.WriteLine(pI.Mass + " " + tmpCoefficient + " " + tmpDensity);
            }
        }

        public void TestDensity() {
            foreach (var VARIABLE in _particleList) {
                Console.WriteLine(VARIABLE.Density);

            }
        }

        private void ComputePressure() {
            int a = 1;
            foreach (var p in this._particleList) {
                p.Pressure = this._gasConstant * (p.Density - _restDensity);
//                Console.WriteLine("Density: "+ p.Density+" - Pressure: "+p.Pressure);
            }
            
        }

        private void ComputePressureForceAcceleration() {
            float tmpCoefficient = (float) (45 / (Math.PI * Math.Pow(this._coreRadius, 6)));
            List<Particle> preParticleList = this._particleList;
            foreach (var pI in this._particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in preParticleList) {
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < this._coreRadius) {
                        tmpVec += (float) ((pI.Pressure + pJ.Pressure) / (2 * pI.Density * pJ.Density) *
                                           Math.Pow(this._coreRadius - tmpDis, 2)) *
                                  Vector3.Normalize(pI.Position - pJ.Position);
                    }
                }

                pI.PressureForceAcceleration = pI.Mass * tmpCoefficient * tmpVec;
                Console.WriteLine(pI.PressureForceAcceleration);
            }
        }

        private void ComputeViscosityForceAcceleration() {
            float tmpCoefficient =
                (float) (this._viscosityCoefficient * 45 / (Math.PI * Math.Pow(this._coreRadius, 6)));
            List<Particle> preParticleList = this._particleList;

            foreach (var pI in this._particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in preParticleList) {
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
                pI.GravityAcceleration = this._gravityAcceleration * new Vector3(0, -1, 0);
            }
        }

        private void ComputeSurfaceTensionAcceleration() {
            float tmpCoefficient = (float) (-945 / (32 * Math.PI * Math.Pow(this._coreRadius, 9)));
            float tmpCoefficientLaplacian = (float) (945 / (8 * Math.PI * Math.Pow(this._coreRadius, 9)));
            float hSquared = (float) Math.Pow(this._coreRadius, 2);
            List<Particle> preParticleList = this._particleList;

            foreach (var pI in this._particleList) {
                Vector3 colorFieldGradient = new Vector3();
                float colorFieldLaplacian = 0;
                foreach (var pJ in preParticleList) {
                    float rSquared = Vector3.DistanceSquared(pI.Position, pJ.Position);
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < this._coreRadius) {
                        colorFieldGradient += (float) (1 / pI.Density * Math.Pow(hSquared - rSquared, 2)) *
                                              (pI.Position - pJ.Position);
                        colorFieldLaplacian += 1 / pI.Density * (hSquared - rSquared) *
                                               (rSquared - 3 / 4 * (hSquared - rSquared));
                    }
                }

                pI.SurfaceTensionAcceleration = -1 * _tensionCoefficient * pI.Mass * tmpCoefficientLaplacian*colorFieldLaplacian /
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
            if (Math.Abs(vecNDotPlaneN) <= 0.001f) {
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
            Vector3[] normalArr = {
                new Vector3(0, 1, 0), new Vector3(1, 0, 0),
                new Vector3(0, 0, 1), new Vector3(0, -1, 0),
                new Vector3(-1, 0, 0), new Vector3(0, 0, -1)
            };
            

            foreach (var pI in this._particleList) {
                bool isOutside = !(0 < pI.Position.X && pI.Position.X < this._containerBox.X &&
                                   0 < pI.Position.Y && pI.Position.Y < this._containerBox.Y &&
                                   0 < pI.Position.Z && pI.Position.Z < this._containerBox.Z);
                if (isOutside) {
                    var vecDir = -Vector3.Normalize(pI.Velocity);
                    Vector3[] anchorPoints = {this._containerBox, new Vector3(0, 0, 0)};
                    float[] disArr= new float[6];
                    for (int i = 0; i < 6; ++i) {
                        if (i < 3) {
                            disArr[i] = IntersectionDetection(pI.Position, vecDir, anchorPoints[0], normalArr[i]);
                        } else {
                            disArr[i] = IntersectionDetection(pI.Position, vecDir, anchorPoints[1], normalArr[i]);

                        }
                    }

                    float tmpMin = Math.Abs(disArr[0]);
                    int tmpIndex = 0;
                    int counter = 0;
                    foreach (var val in disArr) {
                        ++counter;
                        if (val > 0 && val < tmpMin) {
                            tmpMin = val;
                            tmpIndex = counter;
                        }
                    }

                    var tmpReflection = Reflection(vecDir, normalArr[tmpIndex]);
                    pI.Velocity = tmpReflection * pI.Velocity.Length();
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

        public void DisplayParameters() {
//        private int _particleNum;
//        private readonly float _coreRadius; // h
//        private float _particleMass; //set particle
//        private float _restDensity; // set particle \rho_0
//        private float _viscosityCoefficient; // mu
//        private float _gasConstant; // k 
//        private float _tensionCoefficient; //sigma
//        private float _gravityAcceleration; //g

            Console.WriteLine("===================="+"Parameters Setting:"+"====================");
            Console.WriteLine("Particle Number: "+this._particleNum);
            Console.WriteLine("Particle Mass(m): "+this._particleMass);
            Console.WriteLine("Frame Number: "+this._frameNum);
            Console.WriteLine("Time Step: "+this._timeStep);
            Console.WriteLine("Time Duration: "+this._timeDuration);
            Console.WriteLine("Container Box: "+this._containerBox);
            Console.WriteLine("Core Radius(h): "+this._coreRadius);
            Console.WriteLine("Rest Density(rho_0): "+this._restDensity);
            Console.WriteLine("Viscosity Coefficient(mu): "+this._viscosityCoefficient);
            Console.WriteLine("Gas Constant: "+ this._gasConstant);
            Console.WriteLine("Tension Coefficient(sigma): "+this._tensionCoefficient);
            Console.WriteLine("Gravity Acceleration Coefficient(g): "+this._gravityAcceleration);
            Console.WriteLine("====================" + "====================" + "====================\n");
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
                    counter++;
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
            this.Density = 0;
            this.Velocity = new Vector3();
            this.TotalAcceleration = new Vector3();
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
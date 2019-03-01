using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Numerics;

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
            _particleList = new List<Particle>();

            _particleNum = 0;
            _containerBox = inContainerBox; 
            _coreRadius = inCoreRadius;
            _restDensity = inRestDensity;
            _viscosityCoefficient = inViscosityCoefficient;
            _gasConstant = inGasConstant;
            _tensionCoefficient = inTensionCoefficient;
            _gravityAcceleration = inGravityAcceleration;
            _frameNum = inFrameNum;
            _timeStep = inTimeStep;
            _timeDuration = _frameNum * _timeStep;
            _fluidDatabase = new FluidDatabase();

        }

        public void StartImitation(string inFileName) {
            DisplayParameters();
            _fluidDatabase.TimeStep = _timeStep;
            _fluidDatabase.TimeDuration = _timeDuration;
            _fluidDatabase.FileName = inFileName;

            _fluidDatabase.OutputTitle(_timeDuration,_timeDuration);
            DateTime startTime = DateTime.Now;
            InitParticles(new Vector3(10, 10, 10), 0.001f, 0.0004f);
            for (int i = 0; i < _frameNum; ++i) {
                
                ComputeDensity();
                ComputePressure();
                ComputeGravityAcceleration();
                ComputePressureForceAcceleration();
                ComputeViscosityForceAcceleration();
                ComputeSurfaceTensionAcceleration();
                ComputeTotalAcceleration();
                RestrictionParticles();
                UpdatePosition(i+1,_timeStep);
                
                //Display Computing Progress
                
                var endTime = DateTime.Now;
                float runtime = (float) endTime.Subtract(startTime).TotalMilliseconds / 1000;
                float progressPercentage = (i + 1) / _frameNum * 100;
                float remainingTime = progressPercentage / runtime * (100 - progressPercentage);
                Console.SetCursorPosition(0, Console.CursorTop-1);
                Console.WriteLine("Progress: " + progressPercentage.ToString("F3") + "% "
                                  + "Runtime: " + runtime.ToString("F3") + "s "
                                  + "EstimatedTime: " + remainingTime.ToString("F3") + "s");
            }

//            _fluidDatabase.Output(inFileName);
        }

        public void InitParticles(Vector3 inParticleBox, float inParticleInterval, float inParticleMass) {
           
            _particleNum = (int) (inParticleBox.X * inParticleBox.Y * inParticleBox.Z);
            _particleMass = inParticleMass;
            int counter = 0;
            for (int i = 0; i < inParticleBox.X; ++i) {
                for (int j = 0; j < inParticleBox.Y; ++j) {
                    for (int k = 0; k < inParticleBox.Z; ++k) {
                        var tmpPos = new Vector3(inParticleInterval * i, inParticleInterval * j,
                            inParticleInterval * k);
                        _particleList.Add(new Particle(counter++, inParticleMass, tmpPos));
                    }
                }
            }
            _fluidDatabase.AddParticles(0,_particleList);
        }
        
        public void TestInitParticles() {
            Console.WriteLine(_particleList.Count);
            foreach (var p in _particleList) {
                Console.WriteLine(p.Index + " " + p.Position);
                Console.WriteLine(p.Velocity);
            }
        }
        
        private void ComputeDensity() {
            float hSquare = _coreRadius * _coreRadius;
            float tmpCoefficient = (float) (315 / (64 * Math.PI * Math.Pow(_coreRadius, 9)));
            List<Particle> preParticleList = new List<Particle>(_particleList);
            foreach (var pI in _particleList) {
                float tmpDensity = 0;
                foreach (var pJ in preParticleList) {
                    if (pI.Index == pJ.Index) {
                        continue;
                    }
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
            foreach (var p in _particleList) {
                p.Pressure = _gasConstant * (p.Density - _restDensity);
//                Console.WriteLine("Density: "+ p.Density+" - Pressure: "+p.Pressure);
            }
            
        }


        private void ComputePressureForceAcceleration() {
            float tmpCoefficient = (float) (45 / (Math.PI * Math.Pow(_coreRadius, 6)));
            List<Particle> preParticleList = new List<Particle>(_particleList);
            foreach (var pI in _particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in preParticleList) {
                    if (pI.Index == pJ.Index) {
                        continue;
                    }
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < _coreRadius) {
                        tmpVec += (float) ((pI.Pressure + pJ.Pressure) / (2 * pI.Density * pJ.Density) *
                                           Math.Pow(_coreRadius - tmpDis, 2)) *
                                  Vector3.Normalize(pI.Position - pJ.Position);
                    }
                }

                pI.PressureForceAcceleration = pI.Mass * tmpCoefficient * tmpVec;
//                Console.WriteLine(pI.PressureForceAcceleration);
            }
        }

        private void ComputeViscosityForceAcceleration() {
            float tmpCoefficient =
                (float) (_viscosityCoefficient * 45 / (Math.PI * Math.Pow(_coreRadius, 6)));
            List<Particle> preParticleList = new List<Particle>(_particleList);

            foreach (var pI in _particleList) {
                var tmpVec = new Vector3();
                foreach (var pJ in preParticleList) {
                    if (pI.Index == pJ.Index) {
                        continue;
                    }
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < _coreRadius) {
                        tmpVec += (pI.Velocity - pJ.Velocity) / (pI.Density * pJ.Density) * _coreRadius * tmpDis;
                    }
                }

                pI.ViscosityForceAcceleration = pI.Mass * tmpCoefficient * tmpVec;
//                Console.WriteLine(pI.ViscosityForceAcceleration);
            }
        }

        private void ComputeGravityAcceleration() {
            foreach (var pI in _particleList) {
                pI.GravityAcceleration = _gravityAcceleration * new Vector3(0, -1, 0);
            }
        }

        private void ComputeSurfaceTensionAcceleration() {
            float tmpCoefficient = (float) (-945 / (32 * Math.PI * Math.Pow(_coreRadius, 9)));
            float tmpCoefficientLaplacian = (float) (945 / (8 * Math.PI * Math.Pow(_coreRadius, 9)));
            float hSquared = (float) Math.Pow(_coreRadius, 2);

            foreach (var pI in _particleList) {
                Vector3 colorFieldGradient = new Vector3();
                float colorFieldLaplacian = 0;
                
                foreach (var pJ in _particleList) {
                    if (pI.Index == pJ.Index) {
                        continue;
                    }
                    float rSquared = Vector3.DistanceSquared(pI.Position, pJ.Position);
                    float tmpDis = Vector3.Distance(pI.Position, pJ.Position);
                    if (tmpDis < _coreRadius) {
                        colorFieldGradient += (float) (1 / pI.Density * Math.Pow(hSquared - rSquared, 2)) *
                                              (pI.Position - pJ.Position);
                        colorFieldLaplacian += 1 / pI.Density * (hSquared - rSquared) *
                                               (rSquared - 3f / 4f * (hSquared - rSquared));
                    }
                }

                pI.SurfaceTensionAcceleration = -1 * _tensionCoefficient * pI.Mass * tmpCoefficientLaplacian*colorFieldLaplacian /
                                                pI.Density * Vector3.Normalize(colorFieldGradient);
//                Console.WriteLine(pI.SurfaceTensionAcceleration);
            }
        }

        private void ComputeTotalAcceleration() {
            foreach (var pI in _particleList) {
                pI.TotalAcceleration = pI.GravityAcceleration + pI.PressureForceAcceleration +
                                       pI.ViscosityForceAcceleration + pI.SurfaceTensionAcceleration;
//                Console.WriteLine(pI.TotalAcceleration);
            }
        }

        private float IntersectionDetection(Vector3 inVecPos, Vector3 inVecDir, Vector3 inPlanePoint,
            Vector3 inPlaneNormal) {
            float vecNDotPlaneN = Vector3.Dot(inVecDir, inPlaneNormal);
            if (Math.Abs(vecNDotPlaneN) <= 0.001f) {
                return 0f;
            }

            float tmpValue = Vector3.Dot(inVecDir, inPlanePoint - inVecPos);
            return tmpValue / vecNDotPlaneN;
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
            

            foreach (var pI in _particleList) {
                if (pI.Velocity == new Vector3()) {
                    continue;
                }
                bool isOutside = !(0 < pI.Position.X && pI.Position.X < _containerBox.X &&
                                   0 < pI.Position.Y && pI.Position.Y < _containerBox.Y &&
                                   0 < pI.Position.Z && pI.Position.Z < _containerBox.Z);
                if (isOutside) {
                    
                    var vecDir = -Vector3.Normalize(pI.Velocity);
                   
                    Vector3[] anchorPoints = {_containerBox, new Vector3(0, 0, 0)};
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
                        if (val > 0 && val < tmpMin) {
                            tmpMin = val;
                            tmpIndex = counter;
                        }
                        ++counter;
                    }

                    var tmpReflection = Reflection(vecDir, normalArr[tmpIndex]);
//                    Console.WriteLine("Ref");
//                    Console.WriteLine(pI.Velocity);
                    pI.Velocity = tmpReflection * pI.Velocity.Length();
//                    Console.WriteLine("After");
//                    Console.WriteLine(pI.Velocity);
                }
            }
        }

        private void UpdatePosition(int inFrameIndex,float inTimeStep) {
            
            foreach (var pI in _particleList) {
//                Console.WriteLine("Pre: "+pI.Index);
//                Console.WriteLine(pI.Velocity);
//                Console.WriteLine(pI.Position);
                pI.Velocity += pI.TotalAcceleration * inTimeStep;
                pI.Position += pI.Velocity * inTimeStep;
//                Console.WriteLine("latest");
//                Console.WriteLine(pI.Position);
//                Console.WriteLine(pI.Velocity);
            }
            
            _fluidDatabase.AddParticles(inFrameIndex,_particleList);
        }

        private void DisplayParameters() {
//        private int _particleNum;
//        private readonly float _coreRadius; // h
//        private float _particleMass; //set particle
//        private float _restDensity; // set particle \rho_0
//        private float _viscosityCoefficient; // mu
//        private float _gasConstant; // k 
//        private float _tensionCoefficient; //sigma
//        private float _gravityAcceleration; //g

            Console.WriteLine("===================="+"Parameters Setting:"+"====================");
            Console.WriteLine("Particle Number: "+_particleNum);
            Console.WriteLine("Particle Mass(m): "+_particleMass);
            Console.WriteLine("Frame Number: "+_frameNum);
            Console.WriteLine("Time Step: "+_timeStep);
            Console.WriteLine("Time Duration: "+_timeDuration);
            Console.WriteLine("Container Box: "+_containerBox);
            Console.WriteLine("Core Radius(h): "+_coreRadius);
            Console.WriteLine("Rest Density(rho_0): "+_restDensity);
            Console.WriteLine("Viscosity Coefficient(mu): "+_viscosityCoefficient);
            Console.WriteLine("Gas Constant: "+ _gasConstant);
            Console.WriteLine("Tension Coefficient(sigma): "+_tensionCoefficient);
            Console.WriteLine("Gravity Acceleration Coefficient(g): "+_gravityAcceleration);
            Console.WriteLine("====================" + "====================" + "====================\n");
        }

    }

    public class FluidDatabase {
        public float TimeStep { private get; set; }

        public float TimeDuration { private get; set; }
        public string FileName { get; set; }

        public FluidDatabase() {
            
        }
        public void OutputTitle(float inTimeDuration,float inTimeStep) {
            Console.WriteLine(FileName);
            using (StreamWriter file =
                new StreamWriter(FileName)) {
                file.WriteLine(TimeDuration.ToString(CultureInfo.InvariantCulture));
                file.WriteLine(TimeStep.ToString(CultureInfo.InvariantCulture));
            }
        }
        public void AddParticles(int inFrameIndex,List<Particle> inParticleList) {
            using (System.IO.StreamWriter file =
                new System.IO.StreamWriter(FileName, true)) {
                file.WriteLine(inFrameIndex.ToString());
                foreach (var pI in inParticleList) {
                    file.WriteLine(pI.ToString());
                }
            }
        }


    }

    public class Particle {
        public Particle() {

        }

        public Particle(int inIndex, float inMass, Vector3 inPosition) {
            Index = inIndex;
            Mass = inMass;
            Position = inPosition;
            Density = 0;
            Velocity = new Vector3();
            TotalAcceleration = new Vector3();
        }

        public int Index { get; set; }

        public float Mass { get; }

        public Vector3 Position { get; set; }

        public float Density { get; set; }

        public float Pressure { get; set; }

        public Vector3 Velocity { get; set; }

        public Vector3 TotalAcceleration { get; set; }

        public Vector3 PressureForceAcceleration { get; set; }

        public Vector3 ViscosityForceAcceleration { get; set; }

        public Vector3 GravityAcceleration { get; set; }

        public Vector3 SurfaceTensionAcceleration { get; set; }

        public override string ToString() {
            return Index+ " " +Position.Z+" "+Position.Y+" "+Position.Z;
        }
    }
}
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Serialization;

public class FluidManager : MonoBehaviour {
    public ParticleController particlePerfab;
    public int particleNum;
    public int fps;
    public int frameNum;
    public int currentFrame;
    public double frameInterval;
    public double animationDuration;
    public List<List<Vector3>> particleDatabase;
    public List<ParticleController> particleList;

    private void Awake() {
        Debug.Log("123");
        LoadData();
    }

    // Start is called before the first frame update
    void Start() {
        for (int i = 0; i < particleNum; ++i) {
            var tmpObject = GameObject.Instantiate(particlePerfab, this.transform, true);
            tmpObject.GetComponent<ParticleController>().index = i;
            tmpObject.transform.position = particleDatabase[0][i];
            particleList.Add(particlePerfab);
        }
    }

    // Update is called once per frame
    void FixedUpdate() {
        ++currentFrame;
        if (currentFrame >= frameNum) {
            return;
        }
        for (int i = 0; i < particleNum;++ i) {
            particleList[i].transform.position = particleDatabase[currentFrame][i];
        }
    }

    private void LoadData() {
        var lines = System.IO.File.ReadAllLines(
            @"/home/zhiquan/Git-Repository/Fluid-Solver-SPH/SPH-Fluid-Solver/cmake-build-debug/2019-03-22-17:34:03");
        string[] parameters = lines[1].Split(' ');
        particleNum = int.Parse(parameters[0]);
        frameNum = int.Parse(parameters[1]);
        frameInterval = Double.Parse(parameters[2]);
        animationDuration = Double.Parse(parameters[3]);
        fps = (int) (1 / frameInterval);


        particleDatabase = new List<List<Vector3>>();
        for (int i = 0; i < frameNum; ++i) {
            List<Vector3> tmp_list = new List<Vector3>();
            for (int j = 0; j < particleNum; ++j) {
                tmp_list.Add(new Vector3());
            }

            particleDatabase.Add(tmp_list);
        }

        for (int i = 0; i < frameNum; i++) {
            int startIndex = i * (particleNum + 1) + 3;
            int endIndex = startIndex + particleNum;
            for (int j = startIndex; j < endIndex; ++j) {
                var tmp_parameters = lines[j].Replace("(", "").Replace(")", "").Replace(",", " ").Split();
                var tmp_index = int.Parse(tmp_parameters[0]);
                var tmp_pos = new Vector3(float.Parse(tmp_parameters[1]), float.Parse(tmp_parameters[2]),
                    float.Parse(tmp_parameters[3]));
                particleDatabase[i][tmp_index] = tmp_pos;
            }
        }

    }
}
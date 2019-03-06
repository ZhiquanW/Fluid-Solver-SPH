//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_RESTRICTIONBOX_H
#define SPH_FLUID_SOLVER_RESTRICTIONBOX_H

#include "Vector3.h"
class RestrictionBox {
    Vector3 vertexs[2];
    RestrictionBox(const Vector3 &vertex0, const Vector3 &vertex1){
        vertexs[0] = vertex0;
        vertexs[1] = vertex1;
    }
};


#endif //SPH_FLUID_SOLVER_RESTRICTIONBOX_H

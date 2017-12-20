#pragma once

#include <vector>

#include <glmSupport.h>

#include "Skeleton.h"
#include "BezierSpline.h"

class BSplineSkinner {
	std::vector<glm::vec3> curvePoints;

public:
	BSplineSkinner(Joint *joint);

	void generateCurves(Joint *joint);
	void generateCurve(Joint *joint, int link);


};

std::vector<bezier<glm::vec4>> getCurveSet(Joint *joint, int link, float radius, float theta);
glm::vec3 generatePoint(Joint *joint, int link, float s, float radius, float theta);
glm::vec3 blendPair(Skeleton *skeleton, int pivot, int linkA, int linkB, float u, float theta);
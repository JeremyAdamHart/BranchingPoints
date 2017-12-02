#pragma once

#include <vector>

#include <glmSupport.h>

#include "Skeleton.h"

class BSplineSkinner {
	std::vector<glm::vec3> curvePoints;

public:
	BSplineSkinner(Joint *joint);

	void generateCurves(Joint *joint);
	void generateCurve(Joint *joint, int link);


};
#pragma once

#include <glmSupport.h>
#include<vector>
#include"BezierSpline.h"

struct Link;
struct Joint;

struct Link {
	Link(Joint *a, Joint *b);

	Joint *a, *b;
	glm::vec3 dir() const;
};

struct Joint {
	std::vector<Link> links;
	glm::vec3 pos;

	Joint(glm::vec3 pos);

	void addLink(Joint *b);
};

struct BasisPair {
	glm::vec3 bx, by1, by2;
};

float getAngle(glm::vec3 v1, glm::vec3 v2, glm::vec3 n);
BasisPair getFrame(glm::vec3 v1, glm::vec3 v2);
glm::vec3 getCirclePoint(glm::vec3 center, glm::vec3 bx, glm::vec3 by, float radius, float theta);

class Skeleton {
public: 
	Joint *joint;
	std::vector<BasisPair> bases;
	float radius;

	Skeleton(Joint *joint, float radius=0.2f);

	std::vector<bezier<glm::vec4>> getCurveSet(int link, float theta);
	glm::vec3 getEndpoint(int link, float theta);
	float convertAngle(int linkA, int linkB, float theta);
	glm::vec3 getDir(int link, float theta);
	glm::vec3 getDir(int link);
};

//Performs line intersection by looking for closest point to deal with floating point inaccuracy
template<class V> V lineIntersection(V p1, V l1, V p2, V l2) {
	float s = (dot(l1, p2) - dot(l1, p1) - dot(l1, l2)*(dot(l2, p2) - dot(l2, p1))) /
		(1 - dot(l1, l2)*dot(l1, l2));
	return p1 + l1*s;
}

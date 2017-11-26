#pragma once

#include<glmSupport.h>
#include <vector>

template <class V> class bezier {
public:
	std::vector<V> control;

	bezier();
	bezier(V *points, int numPoints);
	bezier(const vector<V>& points);

	V getPoint(float u);
	V getQuadPoint(float u);		//Faster implementation for quadratic curves
	std::vector<V> getQuadPoints(int divisions, float uStart = 0.f, float uEnd = 1.f);
	std::vector<V> getPoints(int divisions, float uStart = 0.f, float uEnd = 1.f);
};
#pragma once

#include <glmSupport.h>
#include<vector>

struct Link;
struct Joint;

struct Link {
	Link(Joint *a, Joint *b);

	Joint *a, *b;
	glm::vec3 dir();
};

struct Joint {
	std::vector<Link> links;
	glm::vec3 pos;

	Joint(glm::vec3 pos);

	void addLink(Joint *b);
};

class Skeleton {

};
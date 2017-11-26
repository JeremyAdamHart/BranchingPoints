#include "Skeleton.h"

using namespace glm;

Link::Link(Joint *a, Joint *b) :a(a), b(b) {}

vec3 Link::dir() {
	return b->pos - a->pos;
}

Joint::Joint(glm::vec3 pos) :pos(pos) {}

void Joint::addLink(Joint *b) {
	links.push_back(Link(this, b));
}
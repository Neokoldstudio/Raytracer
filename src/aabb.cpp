#include "aabb.h" 

//Basically, on veut construire un volume à l'entour d'un obj.
//1-Trouver 8 coins associés au AABB
//2-Construire le AABB
//3-Intersection ---- https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection

//CE que jai en tete:
// -faut prendre la shape de l'obj
// -savoir ya cmb de pts ............ NVM

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {
       return true;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {

	return std::vector<double3>{};
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {

	return AABB{double3{-DBL_MAX,-DBL_MAX,-DBL_MAX},double3{DBL_MAX,DBL_MAX,DBL_MAX}};
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};
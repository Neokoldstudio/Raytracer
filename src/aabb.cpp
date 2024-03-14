#include "aabb.h" 

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max) {
    for (int i = 0; i < 3; ++i) {
        double invD = 1.0 / ray.direction[i];
        double t0 = (min[i] - ray.origin[i]) * invD;
        double t1 = (max[i] - ray.origin[i]) * invD;
        if (invD < 0.0)
            std::swap(t0, t1);
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max <= t_min)
            return false;
    }
    return true;
}
// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
    // à l'aide des points min et max, on peut facilement retrouver tout les coins de la boite englobante
    // en combinant leurs composantes.
	std::vector<double3> corners;
    corners.push_back(aabb.min);
    corners.push_back(double3{aabb.min.x, aabb.min.y, aabb.max.z});
    corners.push_back(double3{aabb.min.x, aabb.max.y, aabb.min.z});
    corners.push_back(double3{aabb.min.x, aabb.max.y, aabb.max.z});
    corners.push_back(double3{aabb.max.x, aabb.min.y, aabb.min.z});
    corners.push_back(double3{aabb.max.x, aabb.min.y, aabb.max.z});
    corners.push_back(double3{aabb.max.x, aabb.max.y, aabb.min.z});
    corners.push_back(aabb.max);
    return corners;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
    double3 min_point{DBL_MAX, DBL_MAX, DBL_MAX};
    double3 max_point{-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (const auto& point : points) {
        min_point = min(min_point, point);
        max_point = max(max_point, point);
    }
    return AABB{min_point, max_point};
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};
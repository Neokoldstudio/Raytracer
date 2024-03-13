#include "object.h"

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray, 
							 double t_min, double t_max, 
							 Intersection *hit) 
{

    auto a = dot(ray.direction, ray.direction);
    auto b = 2.0 * dot(ray.origin, ray.direction);
    auto c = dot(ray.origin, ray.origin) - this->radius * this->radius;
    auto discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return false;

    auto sqrt_discriminant = sqrt(discriminant);
    auto root = (-b - sqrt_discriminant) / (2.0 * a);
    if (root <= t_min || root >= t_max) {
        root = (-b + sqrt_discriminant) / (2.0 * a);
        if (root < t_min || root > t_max) {
            return false;
        }
    }

    hit->depth = root;
    hit->position = ray.origin + root*ray.direction;
    hit->normal = normalize(hit->position);
    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray ray, 
							double t_min, double t_max, 
							Intersection *hit)
{
    double t = -ray.origin.z / ray.direction.z;

    if (t < t_min || t > t_max) {
        return false;
    }

    double3 point = { ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, 0.0 };

    if (std::abs(point.x) > half_size || std::abs(point.y) > half_size) {
        return false;
    }

    double3 normal = { 0.0, 0.0, 1.0 };

    hit->depth = t;
    hit->position = point;
    hit->normal = normalize(normal);
    hit->uv = {(point.x/half_size), (point.y/half_size)};

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
	return Object::compute_aabb();
	//return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray, 
							   double t_min, double t_max, 
							   Intersection *hit)
{
    double a = ray.direction.x * ray.direction.x + ray.direction.z * ray.direction.z;
    double b = 2 * (ray.origin.x * ray.direction.x + ray.origin.z * ray.direction.z);
    double c = ray.origin.x * ray.origin.x + ray.origin.z * ray.origin.z - radius * radius;

    // Résolution de l'équation quadratique pour trouver les points d'intersection
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return false;
    }

    auto sqrt_discriminant = sqrt(discriminant);
    auto root = (-b - sqrt_discriminant) / (2.0 * a);
    auto root2 = (-b + sqrt_discriminant) / (2.0 * a);

    if (root <= t_min || root >= t_max) {
        if (root2 < t_min || root2 > t_max) {
            return false;
        }
    }

    double3 point = ray.origin + root*ray.direction;
    // Vérification de l'appartenance du point au cylindre (si le premier point est trop "loin", check le deuxième)
    // évite des problèmes qui pourrait s'apparenter à du backface culling (ce n'en est pas en réalité, l'effet est trompeur)

    double3 normal = {point.x, 0.0, point.z};

    if (point.y < -half_height || point.y > half_height) {
        root = root2;
        point = ray.origin + root * ray.direction;
        if (point.y < -half_height || point.y > half_height) {
            return false;
        }
    }
    else
        normal = -normal;

    hit->depth = root;
    hit->position = point;
    hit->normal = normalize(normal);

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
	return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,  
						   double t_min, double t_max, 
						   Intersection* hit)
{
	bool intersected = false;

    // Parcourir chaque triangle du maillage et tester s'il y a une intersection
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        if (intersect_triangle(ray, t_min, t_max, triangles[i], hit))
        {
            // Mettre à jour t_max avec la profondeur de l'intersection actuelle
            t_max = hit->depth;
            intersected = true;
        }
    }

    return intersected;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray  ray, 
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

	// Calcul du vecteur de l'arête BA et de l'arête AC
	double3 AB = p1 - p0;
	double3 AC = p2 - p0;

	// Calcul du produit vectoriel de l'arrête AC et du vecteur de direction du rayon
	double3 h = cross(ray.direction, AC);
    // On projette le vecteur h sur l'arrête AB
	double a = dot(AB, h);

	// Si a est proche de zéro, le rayon est parallèle au triangle
	if (a > -EPSILON && a < EPSILON)
		return false;

	double f = 1.0 / a;
	double3 s = ray.origin - p0;
	double u = f * dot(s, h);

	if (u < 0.0 || u > 1.0)
		return false;

	double3 q = cross(s, AB);
	double v = f * dot(ray.direction, q);

	if (v < 0.0 || u + v > 1.0)
		return false;

	// Calcul de t pour trouver l'intersection du rayon et du triangle
	double t = f * dot(AC, q);

	// Vérifier si l'intersection est dans la plage spécifiée
	if (t < t_min || t > t_max)
		return false;

	// Si l'intersection est plus profonde que la valeur actuelle de hit->depth, retourner false -> on ne veut pas override une intersection plus proche
	if (hit->depth < t)
		return false;

	// Mettre à jour les informations de l'intersection
	hit->depth = t;
	hit->normal = normalize(cross(AB, AC));
	hit->position = ray.origin + t * ray.direction;

	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	return Object::compute_aabb();
}
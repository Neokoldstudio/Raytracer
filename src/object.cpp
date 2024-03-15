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
    if(length(ray.origin)<=radius + EPSILON) return false;

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

    // Calculate spherical coordinates
    double phi = atan2(hit->normal.z, hit->normal.x) - (PI / 2.0);;
    double theta = acos(hit->normal.y);

    // Map spherical coordinates to UV coordinates
    double u = 1-(phi + PI) / (2.0 * PI); // Map phi to range [0, 1]
    double v = (theta / PI); // Map theta to range [0, 1]

    // Shift and scale UV coordinates for desired pole positions
    u = u + 0.5;
    v = v;
    hit->uv = { u, v };
    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
    double3 center = mul(transform, double4{0.0, 0.0, 0.0,1.0}).xyz();
    double4 min_point = {center.x - radius, center.y - radius, center.z - radius, 1.0};
    double4 max_point = {center.x + radius, center.y + radius, center.z + radius, 1.0};
    return AABB{min_point.xyz(), max_point.xyz()};
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
    // si la composante Z est presque nulle, le rayon peut être considéré parallère au plan
    if (std::abs(ray.direction.z) <= EPSILON) {
        return false;
    }

    double t = -ray.origin.z / ray.direction.z;

    if (t < t_min || t > t_max) {
        return false;
    }

    double3 point = { ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, 0.0 };

    if (std::abs(point.x) > half_size || std::abs(point.y) > half_size) {
        return false;
    }

    double3 normal = { 0.0, 0.0, 1.0 };

    if(dot(ray.direction, normal)<0)
        normal = -normal;

    double u = (point.x + half_size) / (2.0 * half_size);
    double v = (point.y - half_size) / (2.0 * -half_size);

    hit->depth = t;
    hit->position = point;
    hit->normal = -normal;
    hit->uv = {u,v};

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
    double3 p0{-half_size, -half_size, -EPSILON};  // Bottom left corner
    double3 p1{-half_size, half_size, -EPSILON};   // Top left corner
    double3 p2{half_size, half_size, EPSILON};   // Top right corner
    double3 p3{half_size, -half_size, EPSILON};  // Bottom right corner
    std::vector<double3> corners = {p0, p1, p2, p3};
    std::vector<double3> transformed_corners;
    for (const auto& corner : corners) {
        transformed_corners.push_back(mul(transform, {corner,1.0}).xyz());
    }
    return construct_aabb(transformed_corners);
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
    // If 'a' is close to zero, the ray is nearly parallel to the cylinder's lateral surface
    if (std::abs(a) < EPSILON) {
        return false; // Ray is parallel to the cylinder
    }

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

    double3 point = ray.origin + root * ray.direction;
    double3 normal = {point.x, 0.0, point.z};

    if (point.y < -half_height || point.y > half_height) {
        // le point est à l'interieur, on utilise l'autre solution
        root = root2;
        point = ray.origin + root * ray.direction;
        normal = -normal;
        if (point.y < -half_height || point.y > half_height) {
            return false;
        }
    }

    hit->depth = root;
    hit->position = point;
    hit->normal = normalize(normal);

    double u = atan2(hit->position.z, -hit->position.x) / (2.0 * PI) + 0.5;
    double v = (hit->position.y - half_height) / (2.0 * -half_height);

    hit->uv = { u, v };

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
    // Get the transformation matrix
    double4x4 trans = Cylinder::transform;

    // Calculate the world-space coordinates of the base center (0, 0, 0) of the cylinder
    //centre global
    double3 center = mul(trans, double4{0.0, 0.0, 0.0, 1.0}).xyz();

    //What if on fait juste comme la sphere, mais on rajoute le half_height
    double4 min_point_try = {center.x - Cylinder::radius - Cylinder::half_height, center.y - radius - Cylinder::half_height, center.z - radius - Cylinder::half_height, 1.0};
    double4 max_point_try = {center.x + Cylinder::radius + Cylinder::half_height, center.y + radius + Cylinder::half_height, center.z + radius + Cylinder::half_height, 1.0};


    /*
    // Calculate the directions of the cylinder's axes after rotation
    double3 xAxis = normalize(mul(trans, double4{1.0, 0.0, 0.0, 0.0}).xyz());
    double3 yAxis = normalize(mul(trans, double4{0.0, 1.0, 0.0, 0.0}).xyz());
    double3 zAxis = normalize(mul(trans, double4{0.0, 0.0, 1.0, 0.0}).xyz());


    // Calculate the radius along each axis after rotation
    double3 radiusVector = {Cylinder::radius,Cylinder::radius, Cylinder::half_height}; // Half_height is along z-axis
    double3 rotatedRadius = {
        dot(radiusVector, xAxis),
        dot(radiusVector, yAxis),
        dot(radiusVector, zAxis)
    };
    */

    //ESSAIES:
    //1-Si fonctionne pas, mettre w min à -1 --MARCHE PAS
    //2-Si fonctionne pas non plus, mettre double3 ---Faut changer le return aussi
    //double3 min_point = {center.x - rotatedRadius.x, center.y - rotatedRadius.y, center.z - rotatedRadius.z/*, 1.0*/};
    //double3 max_point = {center.x + rotatedRadius.x, center.y + rotatedRadius.y, center.z + rotatedRadius.z/*,1.0*/};
    //return AABB{min_point/*.xyz()*/, max_point/*.xyz()*/};

    return AABB{min_point_try.xyz(), max_point_try.xyz()};

    //NORMAL:
    // Calculate the minimum and maximum points of the AABB
    //double4 min_point = {center.x - rotatedRadius.x, center.y - rotatedRadius.y, center.z - rotatedRadius.z, 1.0};
    //double4 max_point = {center.x + rotatedRadius.x, center.y + rotatedRadius.y, center.z + rotatedRadius.z,1.0};


    // Return the AABB
    //return AABB{min_point.xyz(), max_point.xyz()};

    //2-essayer de changer pour match ca : return construct_aabb(transformed_corners);
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

	 // Interpolation
    double2 const &ti0 = tex_coords[tri[0].ti];
    double2 const &ti1 = tex_coords[tri[1].ti];
    double2 const &ti2 = tex_coords[tri[2].ti];

    double2 interpolated_uv = ti0 + u * (ti1 - ti0) + v * (ti2 - ti0);

    // Interpolate normals
    double3 const &ni0 = normals[tri[1].ni];
    double3 const &ni1 = normals[tri[1].ni];
    double3 const &ni2 = normals[tri[1].ni];

    double3 interpolated_normal = normalize(ni0 + u * (ni1 - ni0) + v * (ni2 - ni0));

    // Update intersection information
    hit->depth = t;
    hit->position = ray.origin + t * ray.direction;
    hit->normal = interpolated_normal;
    hit->uv = interpolated_uv;

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
    AABB meshAABB;
    if (!positions.empty()) {
        std::vector<double3> transformed_points;// on génère une liste de tout les points transformés dans l'espace global
        for (const auto& point : positions) {
            transformed_points.push_back(mul(transform,{point,1.0}).xyz());
        }
        meshAABB = construct_aabb(transformed_points);
    } else {
        // cas ou le mesh est vide
        meshAABB = AABB();
    }
    return meshAABB;
}
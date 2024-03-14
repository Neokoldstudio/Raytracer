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

    //ICI, 2 options:
    // (1) Soit on calc pour de vrai tous les pts et c'est pas tant précis et ca fait des problèmes et c'est long
    // (2) Ou tout simplement, on part du point milieu qui est 0,0,0 et on additionne le rayon

    ////On fait un conteneur pts dans lequel on va mettre les pts
    std::vector<double3> pts_Sphere;

    /*
    //(1)
    //On fait un conteneur pts dans lequel on va mettre les pts
    //std::vector<double3> pts_Sphere;
    //On a le rayon, alors on peut calculer tous les pts de la sphère


    //Alors, on va calculer les pts avec une formule
    const double theta_max = 2*PI; // Nombre d'étapes pour theta, initialement 20
    const double phi_max = PI; // Nombre d'étapes pour phi, initialement 40

    for (int i = 0; i <= theta_max*10; ++i) {
        double theta = PI * i / theta_max;
        for (int j = 0; j <= phi_max*10; ++j) {
            double phi = 2 * PI * j / phi_max;
            double x = radius * sin(theta) * cos(phi);
            double y = radius * sin(theta) * sin(phi);
            double z = radius * cos(theta);

            pts_Sphere.emplace_back(x, y, z); //remplacer par push_back???
        }
    }
    //ENCORE des problèmes avec, mais on est presque là
    //Faire des mods pour érduire les redondances !!!
    */

    //(2)
    //On peut juste faire ca et ca fonctionne déjà, JSP trop
    double3 pt_milieu_Sphere = {0, 0, 0};
    //ptmilieu + radius? -> jsp si change qlq chose
    /*
    double3 pt_externeX_plus = {pt_milieu_Sphere.x + radius, 0, 0};
    double3 pt_externeX_moins = {-(pt_milieu_Sphere.x + radius), 0, 0};
    double3 pt_externeY_plus = {0, pt_milieu_Sphere.y + radius, 0};
    double3 pt_externeY_moins = {0, -(pt_milieu_Sphere.y + radius), 0};
    double3 pt_externeZ_plus = {0, 0, pt_milieu_Sphere.z + radius};
    double3 pt_externeZ_moins = {0, 0, -(pt_milieu_Sphere.z + radius)};
    */

    double3 pt_Spheremin = {-(pt_milieu_Sphere.x + Sphere::radius), -(pt_milieu_Sphere.y + Sphere::radius), -(pt_milieu_Sphere.z + Sphere::radius)};
    double3 pt_Spheremax = {pt_milieu_Sphere.x + Sphere::radius, pt_milieu_Sphere.y + Sphere::radius, pt_milieu_Sphere.z + Sphere::radius};

    //pts.push_back(pt_milieu_Sphere);
    /*
    pts_Sphere.push_back(pt_externeX_plus);
    pts_Sphere.push_back(pt_externeX_moins);
    pts_Sphere.push_back(pt_externeY_plus);
    pts_Sphere.push_back(pt_externeY_moins);
    pts_Sphere.push_back(pt_externeZ_plus);
    pts_Sphere.push_back(pt_externeZ_moins);
    */
    pts_Sphere.push_back(pt_Spheremin);
    pts_Sphere.push_back(pt_Spheremax);

    //Nouv tentative
    //On créé un objet spécifique sphere locale qui va calculer les coor min et max locales
    //On appel les fonctions de aabb par la suite
    AABB sphere_locale = construct_aabb(pts_Sphere);
    retrieve_corners(sphere_locale);

    //Sphere::transform[0].w;

    //Transformer espace objet vers espace global, on multiplie les pts par la matrice de transformation
    //double minX_Shpere = (Sphere::transform[0].x * sphere_locale.min.x) + (Sphere::transform[0].y * sphere_locale.min.y) + (Sphere::transform[0].z * sphere_locale.min.z) + (Sphere::transform[0].w * 1);
    //double minY_Shpere = (Sphere::transform[1].x * sphere_locale.min.x) + (Sphere::transform[1].y * sphere_locale.min.y) + (Sphere::transform[1].z * sphere_locale.min.z) + (Sphere::transform[1].w * 1);
    //double minZ_Sphere = (Sphere::transform[2].x * sphere_locale.min.x) + (Sphere::transform[2].y * sphere_locale.min.y) + (Sphere::transform[2].z * sphere_locale.min.z) + (Sphere::transform[2].w * 1);
    //double3 pt_Shpere_min_globale = { minX_Shpere, minY_Shpere, minZ_Sphere};
    for (int i = 0; i < 4; ++i){
        std::cout<<"transform x : " << Sphere::transform[i].x <<"    " << "transform y : " <<Sphere::transform[i].y <<"    " <<"transform z : " << Sphere::transform[i].z <<"    " <<"transform w : " << Sphere::transform[i].w<< "\n";
    }

    //1-On fait un nouvel objet AABB local (dans object.cpp) en appelant la fonction construct_aabb (dans aabb.cpp)
    //  la fonction nous retourne un AABB
    //2-Avec l'info retournée, on calcul les coins locaux en faisant appel dans la sphere(dans object.cpp) à la fonction retrieve_corners (dans aabb.cpp)
    //  la fonction nous retourne une liste de coins
    //3-Puisqu'on a les points locaux des formes (dans object.cpp), on veut mtn faire la transformation pour l'espace global
    //  DANS object.h sont définies les fonctions suivantes : double4x4 transform;   // Transformation de l'espace de l'objet à l'espace global (local --> global).
    //                                                        double4x4 i_transform; // Transformation de l'espace de global à l'espace de l'objet (global --> local).
    //
    //                                                        double3x3 n_transform; // Transformation de l'espace de l'objet à l'espace global pour les normales (local --> global).
    //LIVE j'essaie de mult les points trouvés par la matrice transform
    //  mais jsp si c'est comme ca qu'on doit faire

    //compute_aabb() renvoit un aabb, JAI l'impression qu'on va devoir le changer pour renvoyer le AABB de la sphère globale
    //AUSSI, ya pleins de choses définies dans object.h qui semblent important, mais jps trop comment m'en servir

    //Faudrait transformer le pt_milieu PEUT pas fonctionner pcq c'est tout à 0, à moins que le w transforme comme il faut
    //double3 point_milieu_Sphere_Global;
    //double minX_Shpere = (Sphere::transform[0].x * sphere_locale.min.x) + (Sphere::transform[0].y * sphere_locale.min.y) + (Sphere::transform[0].z * sphere_locale.min.z) + (Sphere::transform[0].w * 1);
    //double minY_Shpere = (Sphere::transform[1].x * sphere_locale.min.x) + (Sphere::transform[1].y * sphere_locale.min.y) + (Sphere::transform[1].z * sphere_locale.min.z) + (Sphere::transform[1].w * 1);
    //double minZ_Sphere = (Sphere::transform[2].x * sphere_locale.min.x) + (Sphere::transform[2].y * sphere_locale.min.y) + (Sphere::transform[2].z * sphere_locale.min.z) + (Sphere::transform[2].w * 1);



    //avant, je faisais : construct_aabb(pts_Sphere); mais je crois pas que ce soit bon

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

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
    //half-size -> demi largeur
    //C'EST quoi la hauteur????
    //C'EST UN RECTANGLE, TOUTES les faces sont des rectangles?
    //ICI, QUAD est fait en carré !!!!!!!!

    std::vector<double3> pts_Quad;

    double3 pt_milieu_Quad = {0, 0, 0};
    //ptmilieu + half_size? -> jsp si change qlq chose
    //         /|__/|
    //         |/__|/
    /*
    double3 pt_milieuFaceDroit = {pt_milieu_Quad.x + half_size, 0, 0};
    double3 pt_milieuFaceGauche = {-(pt_milieu_Quad.x + half_size), 0, 0};
    double3 pt_milieu_FaceDevant = {0, pt_milieu_Quad.y + half_size, 0};
    double3 pt_milieu_FaceDerriere = {0, -(pt_milieu_Quad.y + half_size), 0};
    double3 pt_milieu_FaceHaut = {0, 0, pt_milieu_Quad.z + half_size};
    double3 pt_milieu_FaceBas = {0,0, -(pt_milieu_Quad.z + half_size)};
     */

    double3 pt_Quadmin = {-(pt_milieu_Quad.x + half_size), -(pt_milieu_Quad.y + half_size), -(pt_milieu_Quad.z + half_size)};
    double3 pt_Quadmax = {pt_milieu_Quad.x + half_size, pt_milieu_Quad.y + half_size, pt_milieu_Quad.z + half_size};

    /*
    pts_Quad.push_back(pt_milieuFaceDroit);
    pts_Quad.push_back(pt_milieuFaceGauche);
    pts_Quad.push_back(pt_milieu_FaceDevant);
    pts_Quad.push_back(pt_milieu_FaceDerriere);
    pts_Quad.push_back(pt_milieu_FaceHaut);
    pts_Quad.push_back(pt_milieu_FaceBas);
    */
    pts_Quad.push_back(pt_Quadmin);
    pts_Quad.push_back(pt_Quadmax);

    //construct_aabb(pts_Quad);
    AABB quad_locale = construct_aabb(pts_Quad);
    retrieve_corners(quad_locale);

	return Object::compute_aabb();
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
    //radius -> rayon du cercle
    //half_height -> demi hauteur du cylindre par rapport à l'origine

    //Fonctionne pas pcq pas carré?

    std::vector<double3> pts_Cyl;

    double3 pt_milieu_Cyl = {0, 0, 0};

    //Checker si c'est le fait que ca soit pas carré
    double largeur_cyl;
    double hauteur_cyl;
    if (Cylinder::radius < half_height){
        largeur_cyl = half_height;
        hauteur_cyl = half_height;
    }else{
        largeur_cyl = Cylinder::radius;
        hauteur_cyl = Cylinder::radius;
    }

    double3 pt_Cylmin = {-(pt_milieu_Cyl.x + hauteur_cyl + largeur_cyl), -(pt_milieu_Cyl.y + hauteur_cyl + largeur_cyl), -(pt_milieu_Cyl.z + hauteur_cyl + largeur_cyl)};
    double3 pt_Cylmax = {pt_milieu_Cyl.x + hauteur_cyl + largeur_cyl, pt_milieu_Cyl.y + hauteur_cyl + largeur_cyl, pt_milieu_Cyl.z + hauteur_cyl + largeur_cyl};

    pts_Cyl.push_back(pt_Cylmin);
    pts_Cyl.push_back(pt_Cylmax);

    std::cout <<"cylindre x min :  " << pt_Cylmin.x << "    " <<"cylindre y min :  " << pt_Cylmin.y << "    " <<"cylindre z min :  " << pt_Cylmin.z <<"\n";

    //construct_aabb(pts_Cyl);
    AABB cyl_locale = construct_aabb(pts_Cyl);
    retrieve_corners(cyl_locale);

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
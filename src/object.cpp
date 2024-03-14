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
    //INFORMATION:
    /*
     * -JSP exactement si j'ai appelé retreive_corners et construct_aabb comme il faut
     * -À la fin, on retourne compute_aabb(), mais je semble pas faire grand chose avec les coins
     * -JSP si on doit avoir accès aux coins dans le fichier aabb.cpp, avant je callais retrieve_corners dans le construct_aabb
     * -METTRE la fin de ligne quand on imprime, sinon ca bug
     */


    //DEBUT À IGNORER (relativement) -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*
    //Lignes que je veux pas suprimer, pcq j'ai pt mal fait:

    //PTS externes de la sphère

    //On prend des extrémités.
    //On peut déjà connaitre le maximum et minimum local
    double3 pt_externeX_plus = {pt_milieu_Sphere.x + radius, 0, 0};
    double3 pt_externeX_moins = {-(pt_milieu_Sphere.x + radius), 0, 0};
    double3 pt_externeY_plus = {0, pt_milieu_Sphere.y + radius, 0};
    double3 pt_externeY_moins = {0, -(pt_milieu_Sphere.y + radius), 0};
    double3 pt_externeZ_plus = {0, 0, pt_milieu_Sphere.z + radius};
    double3 pt_externeZ_moins = {0, 0, -(pt_milieu_Sphere.z + radius)};


    //J'avais mis tous les points dans la liste qu'on envoit à construct aab, mais je ne crois pas que ce soit nécessaire puisqu'on a déjà les points min et max

    pts_Sphere.push_back(pt_externeX_plus);
    pts_Sphere.push_back(pt_externeX_moins);
    pts_Sphere.push_back(pt_externeY_plus);
    pts_Sphere.push_back(pt_externeY_moins);
    pts_Sphere.push_back(pt_externeZ_plus);
    pts_Sphere.push_back(pt_externeZ_moins);


    //Calcul initial de points globaux

    //Transformer espace objet vers espace global, on multiplie les pts par la matrice de transformation
    double minX_Shpere = (Sphere::transform[0].x * sphere_locale.min.x) + (Sphere::transform[0].y * sphere_locale.min.y) + (Sphere::transform[0].z * sphere_locale.min.z) + (Sphere::transform[0].w * 1);
    double minY_Shpere = (Sphere::transform[1].x * sphere_locale.min.x) + (Sphere::transform[1].y * sphere_locale.min.y) + (Sphere::transform[1].z * sphere_locale.min.z) + (Sphere::transform[1].w * 1);
    double minZ_Sphere = (Sphere::transform[2].x * sphere_locale.min.x) + (Sphere::transform[2].y * sphere_locale.min.y) + (Sphere::transform[2].z * sphere_locale.min.z) + (Sphere::transform[2].w * 1);
    double3 pt_Shpere_min_globale = { minX_Shpere, minY_Shpere, minZ_Sphere};


    //1ere tentative pt milieu

    //Faudrait transformer le pt_milieu PEUT pas fonctionner pcq c'est tout à 0, à moins que le w transforme comme il faut
    //double3 point_milieu_Sphere_Global;
    //double minX_Shpere = (Sphere::transform[0].x * sphere_locale.min.x) + (Sphere::transform[0].y * sphere_locale.min.y) + (Sphere::transform[0].z * sphere_locale.min.z) + (Sphere::transform[0].w * 1);
    //double minY_Shpere = (Sphere::transform[1].x * sphere_locale.min.x) + (Sphere::transform[1].y * sphere_locale.min.y) + (Sphere::transform[1].z * sphere_locale.min.z) + (Sphere::transform[1].w * 1);
    //double minZ_Sphere = (Sphere::transform[2].x * sphere_locale.min.x) + (Sphere::transform[2].y * sphere_locale.min.y) + (Sphere::transform[2].z * sphere_locale.min.z) + (Sphere::transform[2].w * 1);


    //FIN À IGNORER -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    */


    //DEBUT NOTES --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*
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


    //FACON 1- : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //1-On transforme le point milieu
    double3 point_milieu_Sphere_Global;
    //JSP si c'est mieux de faire "Sphere::transform[0].x" OU "Sphere::transform[0][0]" À VOIR
    double X_Sphere_Global = (Sphere::transform[0].x * pt_milieu_Sphere_locale.x) + (Sphere::transform[0].y * pt_milieu_Sphere_locale.y) + (Sphere::transform[0].z * pt_milieu_Sphere_locale.z) + (Sphere::transform[0].w * 1);
    double Y_Sphere_Global = (Sphere::transform[1].x * pt_milieu_Sphere_locale.x) + (Sphere::transform[1].y * pt_milieu_Sphere_locale.y) + (Sphere::transform[1].z * pt_milieu_Sphere_locale.z) + (Sphere::transform[1].w * 1);
    double Z_Sphere_Global = (Sphere::transform[2].x * pt_milieu_Sphere_locale.x) + (Sphere::transform[2].y * pt_milieu_Sphere_locale.y) + (Sphere::transform[2].z * pt_milieu_Sphere_locale.z) + (Sphere::transform[2].w * 1);
    point_milieu_Sphere_Global = {X_Sphere_Global, Y_Sphere_Global, Z_Sphere_Global};

    //2-On veut la transformation du rayon
    //Au cas où ce n'est pas une transformation uniforme, il faut faire 3 rayons différents
    //On multiplie par la diagonale de la matrice de transformation
    //JSP si c'est mieux de faire "Sphere::transform[0].x" OU "Sphere::transform[0][0]" À VOIR
    double sphere_rayonX = Sphere::radius * Sphere::transform[0].x;
    double sphere_rayonY = Sphere::radius * Sphere::transform[1].y;
    double sphere_rayonZ = Sphere::radius * Sphere::transform[2].z;

    //3-On calcul le min et le max global
    double3 pt_Sphere_min_Global = {-(point_milieu_Sphere_Global.x + sphere_rayonX), -(point_milieu_Sphere_Global.y + sphere_rayonY), -(point_milieu_Sphere_Global.z + sphere_rayonZ)};
    double3 pt_Sphere_max_Global = {point_milieu_Sphere_Global.x + sphere_rayonX, point_milieu_Sphere_Global.y + sphere_rayonY, point_milieu_Sphere_Global.z + sphere_rayonZ};

    //4-On store les pts
    pts_Sphere_Globale.push_back(pt_Sphere_min_Global);
    pts_Sphere_Globale.push_back(pt_Sphere_max_Global);

    //5-On fait le AABB global
    AABB Sphere_Globale = construct_aabb(pts_Sphere_Globale); //avant, je faisais : construct_aabb(pts_Sphere); mais je crois pas que ce soit bon
    Coins_Globaux = retrieve_corners(Sphere_Globale);

    //SI la transformation n'est pas uniforme, il faut

    //FIN NOTES ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    */

    //DEBUT CODE -------------------------------------------------------------------------------------------------------

    //On fait un conteneur pts dans lequel on va mettre les pts
    std::vector<double3> pts_Sphere_Locale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Locaux;

    //On fait la même chose, mais pour les coor globales
    std::vector<double3> pts_Sphere_Globale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Globaux;


    //(1) : Trouver les points locaux


    //On sait que le point milieu c'est (0, 0, 0)
    double3 pt_milieu_Sphere_locale = {0, 0, 0};

    //On trouve le min et le max local, car on a le rayon
    double3 pt_Sphere_Locale_min = {-(pt_milieu_Sphere_locale.x + Sphere::radius), -(pt_milieu_Sphere_locale.y + Sphere::radius), -(pt_milieu_Sphere_locale.z + Sphere::radius)};
    double3 pt_Sphere_Locale_max = {pt_milieu_Sphere_locale.x + Sphere::radius, pt_milieu_Sphere_locale.y + Sphere::radius, pt_milieu_Sphere_locale.z + Sphere::radius};

    //On met le min et le max dans la liste qu'on va envoyer à construct_aabb
    pts_Sphere_Locale.push_back(pt_Sphere_Locale_min);
    pts_Sphere_Locale.push_back(pt_Sphere_Locale_max);

    //On créé un objet spécifique sphere locale qui va calculer les coor min
    // et max locales et on appel les fonctions du fichier aabb
    AABB Sphere_Locale = construct_aabb(pts_Sphere_Locale);
    Coins_Locaux = retrieve_corners(Sphere_Locale);


    //Impression de la matrice de transformation
    //for (int i = 0; i < 4; ++i){
    //    for (int j = 0; j < 4; ++j){
    //        std::cout <<"transform " <<i <<", " <<j <<" : " <<Sphere::transform[i][j] <<"    "; /*<<"    " << "transform y : " <<Sphere::transform[i][j] <<"    " <<"transform z : " << Sphere::transform[i][j] <<"    " <<"transform w : " << Sphere::transform[i][j]<< "\n";*/
    //    }
    //    std::cout << "\n" << std::endl;
    //}


    //(2) : Trouver les points globaux

    //ICI, JAI 2 FACONS:
    //                  1-On transforme le pt milieu, puis on transforme le rayon et on refais les mêmes calculs
    //                      MAIS, chui certain de la boite est allignée avec les axes
    //                      (Je crois pas que c'est bon, car utiliser les fonctions du aabb ne semble pas nécessaire.)
    //                  *** J'ai mis la facon 1- dans les notes ***

    //                  2-On transforme les points de la liste des coins locaux avec la matrice de transformation
    //                      On transforme aussi le rayon, juste pour être sûr.
    //                      (JSP si j'ai fais de la bonne facon)


    //FACON 2- :

    //On fait une boucle sur les coins locaux et on applique les transformations sur ceux-ci
    //On store dans la liste de points globaux
    for (const auto& point : Coins_Locaux){
        //Appllique transformation sur chacun des coins
        double x_Transformation_Globale = (Sphere::transform[0].x * point.x) + (Sphere::transform[0].y * point.y) + (Sphere::transform[0].z * point.z) + (Sphere::transform[0].w * 1);
        double y_Transformation_Globale = (Sphere::transform[1].x * point.x) + (Sphere::transform[1].y * point.y) + (Sphere::transform[1].z * point.z) + (Sphere::transform[1].w * 1);
        double z_Transformation_Globale = (Sphere::transform[2].x * point.x) + (Sphere::transform[2].y * point.y) + (Sphere::transform[2].z * point.z) + (Sphere::transform[2].w * 1);

        //On store le point, car si j'essaie de le mettre directement dans la liste de coins,ca bugg
        double3 pt_Transforme_LtoG = {x_Transformation_Globale, y_Transformation_Globale, z_Transformation_Globale};
        //On store les points dans la liste de coins globaux
        Coins_Globaux.push_back(pt_Transforme_LtoG);
    }

    //IMPRESSION pour vérification
    //for (const auto& point : Coins_Globaux){
    //    std::cout << "coins globaux: (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
    //}

    //On fait le aabb global, mais avec la liste des points globaux, on veut trouver le min et le max
    AABB Sphere_Globale = construct_aabb(Coins_Globaux);

    //Pour s'assurer que ca soit alligné avec les axes, on appel retrieve corners
    Coins_Globaux = retrieve_corners(Sphere_Locale);


	return Object::compute_aabb();
    //return Sphere::compute_aabb();
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
    //Faire carré à l'entour?
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
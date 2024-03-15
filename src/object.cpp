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
    /*
    //                  1-On transforme le pt milieu, puis on transforme le rayon et on refais les mêmes calculs
    //                      MAIS, chui certain de la boite est allignée avec les axes
    //                      (Je crois pas que c'est bon, car utiliser les fonctions du aabb ne semble pas nécessaire.)
    //                  *** J'ai mis la facon 1- dans les notes ***

    //                  2-On transforme les points de la liste des coins locaux avec la matrice de transformation
    //                      On transforme aussi le rayon, juste pour être sûr.
    //                      (JSP si j'ai fais de la bonne facon)
    */

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
        pts_Sphere_Globale.push_back(pt_Transforme_LtoG);
    }

    //IMPRESSION pour vérification
    for (const auto& point : Coins_Globaux){
        std::cout << "coins globaux: (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
    }

    //On fait le aabb global, mais avec la liste des points globaux, on veut trouver le min et le max
    AABB Sphere_Globale = construct_aabb(pts_Sphere_Globale);

    //Pour s'assurer que ca soit alligné avec les axes, on appel retrieve corners
    Coins_Globaux = retrieve_corners(Sphere_Globale);

    for (const auto& point : Coins_Globaux){
        std::cout << "coins globaux 2 : (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
    }

    //return Sphere::compute_aabb();
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
    //NOTES confusion : ------------------------------------------------------------------------------------------------
    /*
    // Espace Local: Quad(Rectangle) centrée à l'origine tel que la normale est Z+ pour une largeur de (2 * half_size) x (2 * half_size)
    //Selon cette phrase (object.h ligne 203), j'assume que le Quad est un carré, mais ici, ca dit que c'est un rectangle.
    //Chui un peu confu pcq ce que je comprend de l'énoncé, on devrait avoir un quadrillatère en général,
    // MAIS on nous a juste fourni une demi largeur, alors forcément, ca doit être un carré.
    //OU c'est un carré dans l'espace local et ensuite, après sa transformation, ca sera autre chose?
    */

    //On fait un conteneur pts dans lequel on va mettre les pts
    std::vector<double3> pts_Quad_Locale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Quad_Locaux;

    //On fait la même chose, mais pour les coor globales
    std::vector<double3> pts_Quad_Globale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Quad_Globaux;

    //On fait des vars au cas où on a un plan après la transformation
    //double3 point_Quad_Global_min_YZ;
    //double3 point_Quad_Global_min_XZ;
    //double3 point_Quad_Global_min_XY;

    //À la place, jvais juste faire un min et un max globale non précis qu'on va changer selon la situation
    double3 pt_Quad_Global_min;
    double3 pt_Quad_Global_max;


    //1-ESPACE LOCAL

    //On cree le point milieu
    double3 pt_milieu_Quad_Local = {0, 0, 0};

    //On trouve le minimum et le maximum
    double3 pt_Quad_Locale_min = {-(pt_milieu_Quad_Local.x + Quad::half_size), -(pt_milieu_Quad_Local.y + Quad::half_size), -(pt_milieu_Quad_Local.z + Quad::half_size)};
    double3 pt_Quad_Locale_max = {pt_milieu_Quad_Local.x + Quad::half_size, pt_milieu_Quad_Local.y + Quad::half_size, pt_milieu_Quad_Local.z + Quad::half_size};

    //On store les points min et max
    pts_Quad_Locale.push_back(pt_Quad_Locale_min);
    pts_Quad_Locale.push_back(pt_Quad_Locale_max);

    //On fait un objet aabb local et on trouve les coins
    AABB Quad_Local = construct_aabb(pts_Quad_Locale);
    Coins_Quad_Locaux = retrieve_corners(Quad_Local);


    //2-ESPACE GLOBAL

    //Je crois pas que ca c'est encore bon
    /*
    //1-On calcule la transformation sur les points min et max en premier
    //  on va voir si le rectangle s'est transformé en plan,
    //      je crois que c'est plus efficace que de faire un aabb et ensuite le modifier ou devoir le suprimer et en refaire un autre(*)
    //(*) : on va peut-être devoir faire ca si on a déjà calculé toutes les transformations et on a fait un aabb
    //jsp trop comment on va faire si on calcule les transformations en premier, va falloir chercher le min et le max quand meme


    pt_Quad_Global_min = {(Quad::transform[0].x * pt_Quad_Locale_min.x) + (Quad::transform[0].y * pt_Quad_Locale_min.y) + (Quad::transform[0].z * pt_Quad_Locale_min.z) + (Quad::transform[0].w * 1),
                          (Quad::transform[1].x * pt_Quad_Locale_min.x) + (Quad::transform[1].y * pt_Quad_Locale_min.y) + (Quad::transform[1].z * pt_Quad_Locale_min.z) + (Quad::transform[1].w * 1),
                          (Quad::transform[2].x * pt_Quad_Locale_min.x) + (Quad::transform[2].y * pt_Quad_Locale_min.y) + (Quad::transform[2].z * pt_Quad_Locale_min.z) + (Quad::transform[2].w * 1)};

    pt_Quad_Global_max = {(Quad::transform[0].x * pt_Quad_Locale_max.x) + (Quad::transform[0].y * pt_Quad_Locale_max.y) + (Quad::transform[0].z * pt_Quad_Locale_max.z) + (Quad::transform[0].w * 1),
                          (Quad::transform[1].x * pt_Quad_Locale_max.x) + (Quad::transform[1].y * pt_Quad_Locale_max.y) + (Quad::transform[1].z * pt_Quad_Locale_max.z) + (Quad::transform[1].w * 1),
                          (Quad::transform[2].x * pt_Quad_Locale_max.x) + (Quad::transform[2].y * pt_Quad_Locale_max.y) + (Quad::transform[2].z * pt_Quad_Locale_max.z) + (Quad::transform[2].w * 1)};

    pts_Quad_Globale.push_back(pt_Quad_Global_min);
    pts_Quad_Globale.push_back(pt_Quad_Global_max);
    */

    //Puisqu'on va faire la transformation du carré, on veut voir si on ne l'a pas transformé en un plan
    //EX : SI x est 0 pour le min et le max, alors le rectangle est un plan YZ
    // On veut alors lui faire un volume.
    //Cependant, je crois qu'il faut le mettre au dessus, à cote ou en dessous => je crois qu'il faut qu'il fasse face au rayon
    if (pt_Quad_Global_min.x == 0 && pt_Quad_Global_max.x == 0){
        std::cout<<"On a un plan" << std::endl;
        // point_Quad_Global_min_YZ

        //pt_Quad_Global_min =
        //pt_Quad_Global_max =
    }
    if (pt_Quad_Global_min.y == 0 && pt_Quad_Global_max.y == 0){
        std::cout<<"On a un plan" << std::endl;
        // point_Quad_Global_min_XZ;

        //pt_Quad_Global_min =
        //pt_Quad_Global_max =
    }
    if (pt_Quad_Global_min.z == 0 && pt_Quad_Global_max.z == 0){
        std::cout<<"On a un plan" << std::endl;
        // point_Quad_Global_min_XY;

        //pt_Quad_Global_min =
        //pt_Quad_Global_max =
    }

    //On va transformer les points
    //Devrait être la même logique que pour la sphère
    for (const auto& point : Coins_Quad_Locaux){
        //Applique transformation sur chacun des coins
        double x_Transformation_Globale = (Quad::transform[0].x * point.x) + (Quad::transform[0].y * point.y) + (Quad::transform[0].z * point.z) + (Quad::transform[0].w * 1);
        double y_Transformation_Globale = (Quad::transform[1].x * point.x) + (Quad::transform[1].y * point.y) + (Quad::transform[1].z * point.z) + (Quad::transform[1].w * 1);
        double z_Transformation_Globale = (Quad::transform[2].x * point.x) + (Quad::transform[2].y * point.y) + (Quad::transform[2].z * point.z) + (Quad::transform[2].w * 1);

        //On store le point, car si j'essaie de le mettre directement dans la liste de coins,ca bugg
        double3 pt_Transforme_LtoG = {x_Transformation_Globale, y_Transformation_Globale, z_Transformation_Globale};
        //On store les points dans la liste de points globaux
        pts_Quad_Globale.push_back(pt_Transforme_LtoG);
    }

    //On fait le aabb
    //Quand on fait le aabb avec la liste des points, le construct va automatiquement trouver le min et le max
    AABB Quad_Global = construct_aabb(pts_Quad_Globale);
    Coins_Quad_Globaux = retrieve_corners(Quad_Global);

    //IMPRESSION pour vérifier
    for (const auto& point : Coins_Quad_Globaux){
        std::cout << "coins globaux: (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
    }
    std::cout <<"min global : \n" <<"x : " <<Quad_Global.min.x <<"    " <<"y : " <<Quad_Global.min.y <<"    " <<"z : " <<Quad_Global.min.z <<"\n" <<std::endl;
    std::cout <<"min global : \n" <<"x : " <<Quad_Global.max.x <<"    " <<"y : " <<Quad_Global.max.y <<"    " <<"z : " <<Quad_Global.max.z <<"\n" <<std::endl;

    //EST_CE QUE ya d'autres cas à faire?

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

    //On fait un conteneur pts dans lequel on va mettre les pts
    std::vector<double3> pts_Cyl_Locale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Cyl_Locaux;

    //On fait la même chose, mais pour les coor globales
    std::vector<double3> pts_Cyl_Globale;
    //On fait une liste qui contient les coins raportés par retrieve corners
    std::vector<double3> Coins_Cyl_Globaux;

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
	return Object::compute_aabb();
}
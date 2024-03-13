#include <iostream>
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
//CODE BING AI
//double Dist = hypot(hypot(aabb.max.x - aabb.min.x, aabb.max.y - aabb.min.y), aabb.max.z - aabb.min.z); //On veut la distance entre les deux points, MAIS chui pas capable de get les x,y,z ...... nvm, c'est comme ca :"aabb.min.x"
std::vector<double3> retrieve_corners(AABB aabb) {
    //On appel cette fonc à partir de construct_aabb
    //On recoit les min et les max

    std::vector<double3> coins;
    //On veut mettre les points uns par dessus les autres

    //PROF ----------------------------------------------------------------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //On doit trouver les matrices de T pour mettre dans l'espace golbal
    //PUIS faire en sorte que la boite reste allignée
    //S'assurer d'un volume = si rectangle est planaire bha ya uncun volume


    //So basically, si on se met exactement au dessu du point minimum, mais à la hauteur du point maximum, ca nous donne le coin droit,devant,haut
    //Gauche/Droite : aabb.min.x/aabb.max.x
    //Devant/Derriere : aabb.min.y/aabb.max.y
    //Bas/Haut : aabb.min.z/aabb.max.z

    //devantBasGauche = point minimum
    //derriereHautDroit = point maximum
    double3 coinGaucheDevantHaut = {aabb.min.x, aabb.min.y, aabb.max.z}; //On a un point au dessus du point minimum, mais à la hauteur du point max
    double3 coinDroitDevantHaut = {aabb.max.x, aabb.min.y, aabb.max.z};
    double3 coinDroitDevantBas = {aabb.max.x, aabb.min.y, aabb.min.z};
    double3 coinGaucheDevantBas = aabb.min; //{aabb.min.x, aabb.min.y, aabb.min.z};

    double3 coinGaucheDerriereHaut = {aabb.min.x, aabb.max.y, aabb.max.z};
    double3 coinDroitDerriereHaut = aabb.max; //{aabb.max.x, aabb.max.y, aabb.max.z};
    double3 coinDroitDerriereBas = {aabb.max.x, aabb.max.y, aabb.min.z};
    double3 coinGaucheDerriereBas = {aabb.min.x, aabb.max.y, aabb.min.z};

    coins.push_back(coinGaucheDevantHaut);
    coins.push_back(coinDroitDevantHaut);
    coins.push_back(coinDroitDevantBas);
    coins.push_back(coinGaucheDevantBas);

    coins.push_back(coinGaucheDerriereHaut);
    coins.push_back(coinDroitDerriereHaut);
    coins.push_back(coinDroitDerriereBas);
    coins.push_back(coinGaucheDerriereBas);

    /*
    //IMPRIMER pour valider
    for (const auto& point : coins) {
        std::cout << "x : " << point.x <<"\n" << "y : " << point.y << "\n" << "z : " << point.z << "\n" << "\n";
    }
    */

    //return std::vector<double3>{};
    return coins;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
    //ICI, on fait boucle sur les pts,
    // on trouves les extrémités,                       (1)
    // puis on appel retrive_corners                    (2)
    //      avec retrive_corners, on get les coins      (3)
    //Une fois qu'on get les coins, on calc le volume.  (4)


    //(1)
    //On définie les extrémités
    //CODE CHAT GPT
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double minZ = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    double maxZ = std::numeric_limits<double>::lowest();

    //On trouve les extrémités parmis la liste des pts
    for (const auto& point : points) {
        minX = std::min(minX, point.x);
        minY = std::min(minY, point.y);
        minZ = std::min(minZ, point.z);
        maxX = std::max(maxX, point.x);
        maxY = std::max(maxY, point.y);
        maxZ = std::max(maxZ, point.z);
    }

    //J'ai rajouté/importé std::count, juste pour imprimer,
    // JSP si on peut changer ou si ca affecte quoi que ce soit d'autre
    std::cout << "minX : " << minX <<"\n";
    std::cout << "minY : " << minY <<"\n";
    std::cout << "minZ : " << minZ <<"\n";
    std::cout << "maxX : " << maxX <<"\n";
    std::cout << "maxY : " << maxY <<"\n";
    std::cout << "maxZ : " << maxZ <<"\n";

    //(2)
    //ENSUITE, on appel la fonction retrieve_corners avec notre nouvel AABB
    //retrieve_corners(AABB{double3{minX,minY,minZ},double3{maxX,maxY,maxZ}});
    //On va appeler cette fonction à partir de object
    //Pour ca, faut retourner le AABB

	//return AABB{double3{-DBL_MAX,-DBL_MAX,-DBL_MAX},double3{DBL_MAX,DBL_MAX,DBL_MAX}};
    return AABB{double3{minX,minY,minZ},double3{maxX,maxY,maxZ}};
    /*
    for (const auto& point : points){
        int x = 8.0;
        if (point.x > x | point.y > x | point.z > x){
            std::cout << "Point: (" << point.x << ", " << point.y << ", " << point.z << ") ETRANGE, > que 8 \n";
        } else{
            std::cout << "Point: (" << point.x << ", " << point.y << ", " << point.z << ")\n";
        }
        std::cout << "Point: (" << point.x << ", " << point.y << ", " << point.z << ")\n";
        // POUR FAIRE CA, FAUT IMPORTER UN AUTRE BHAY DE LA LIBRARY, MAIS JE VX PAS PRENDRE LE RISQUE DE FUCK UP
        //pour le moment, ca fonctionne -> on passe à travers la liste des pts en ordre.

    }*/
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};
#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    if (root == nullptr)
        return false;

    // Create a vector to store nodes to visit (simulating a stack)
    std::vector<BVHNode*> nodes_to_visit;
    nodes_to_visit.push_back(root);

    bool has_intersection = false;
    double closest_intersection = t_max;

    // Traverse the BVH tree using depth-first search
    while (!nodes_to_visit.empty()) {
        BVHNode* node = nodes_to_visit.back();
        nodes_to_visit.pop_back();

        if (node->aabb.intersect(ray, t_min, closest_intersection)) {
            if (node->left == nullptr && node->right == nullptr) {
                Intersection temp_hit;
                if (objects[node->idx]->intersect(ray, t_min, closest_intersection, &temp_hit)) {
                    has_intersection = true;
                    closest_intersection = temp_hit.depth;
                    *hit = temp_hit;
                }
            } else {
                nodes_to_visit.push_back(node->left);
                nodes_to_visit.push_back(node->right);
            }
        }
    }

    return has_intersection;
}


// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.

bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool hasIntersection = false;
    double closestIntersection = t_max;

    for(auto obj : objects) {
        Intersection tempHit;

        AABB boundingBox = obj->compute_aabb();

        if(boundingBox.intersect(ray, t_min, t_max)){
            if(obj->intersect(ray, t_min, t_max, &tempHit)) {
                if(tempHit.depth < closestIntersection) {
                    hasIntersection = true;
                    closestIntersection = tempHit.depth;
                    *hit = tempHit; // Update hit avec les donées de la nouvelle intersection
                }
            }
        }
    }

    return hasIntersection;
}

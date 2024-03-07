#include "raytracer.h"

void Raytracer::render(const Scene& scene, Frame* output)
{
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far;
    }

    // Calcule les vecteurs de base de la caméra.
    double3 cameraPos = scene.camera.position;
    double3 lookAt = normalize(scene.camera.center - cameraPos);
    double3 right = normalize(cross(lookAt, scene.camera.up));
    double3 up = normalize(cross(right, lookAt));
    
    // Calcule les dimensions de la fenêtre de visualisation.
    double fov = scene.camera.fovy;
    double aspect_ratio = scene.camera.aspect;
    double viewport_height = scene.camera.z_near * tan(deg2rad(fov*0.5))*2;
    double viewport_width = viewport_height * aspect_ratio;

    double3 bottomLeftLocal = cameraPos - right * (viewport_width / 2.0) - up * (viewport_height / 2.0) + lookAt * scene.camera.z_near;

    double jittering_radius = scene.jitter_radius;

    // Itère sur tous les pixels de l'image.
    for (int y = 0; y < scene.resolution[1]; y++) {
        if (y % 40){
            std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
        }

        for(int x = 0; x < scene.resolution[0]; x++) {

            // generate random offset for jittering
            int avg_z_depth = 0;
            double3 avg_ray_color{0,0,0};
            
            for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
                // Génère le rayon approprié pour ce pixel.
                Ray ray;
                // Initialise la profondeur de récursivité du rayon.
                int ray_depth = 0;
                // Initialize la couleur du rayon
                double3 ray_color{0,0,0};

                // Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
                // Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius. 
                // Faites la moyenne des différentes couleurs obtenues suite à la récursion.
                double2 randomOffset = random_in_unit_disk() * jittering_radius;
                double ray_depth_out = scene.camera.z_far;
                double deltaX = viewport_width/(double)scene.resolution[0];
                double deltaY = viewport_height/(double)scene.resolution[1];

                double3 viewportPixelCoord = bottomLeftLocal + right * (x + randomOffset.x) * deltaX + up * (y + randomOffset.y) * deltaY;

                ray.origin = cameraPos;
                ray.direction = normalize(viewportPixelCoord - cameraPos);

                trace(scene, ray, ray_depth, &ray_color, &ray_depth_out);

                avg_ray_color += ray_color;
                avg_z_depth += (int)ray_depth_out;
            }

            avg_z_depth = avg_z_depth / scene.samples_per_pixel;
            avg_ray_color = avg_ray_color / scene.samples_per_pixel;

            // Test de profondeur
            if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
                avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
                z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

                // Met à jour la couleur de l'image (et sa profondeur)
                output->set_color_pixel(x, y, avg_ray_color);
                output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / (scene.camera.z_far-scene.camera.z_near));
            }
        }
    }

    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeur.
void Raytracer::trace(const Scene& scene,
					  Ray ray, int ray_depth,
					  double3* out_color, double* out_z_depth)
{
	Intersection hit;
	// Fait appel à l'un des containers spécifiées.
    
    if(scene.container->intersect(ray,EPSILON,*out_z_depth,&hit)) {		
		Material& material = ResourceManager::Instance()->materials[hit.key_material];

		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
		
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		// 
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.

		*out_color = shade(scene, hit);
		*out_z_depth = hit.depth;
	}
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
	// Material& material = ResourceManager::Instance()->materials[hit.key_material]; lorsque vous serez rendu à la partie texture.
	return hit.normal;
}

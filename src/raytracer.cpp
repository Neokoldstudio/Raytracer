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

        double3 reflexion_color = {0,0,0};
        double3 refraction_color = {0,0,0};


        if(ray_depth <= 1){
            // @@@@@@ VOTRE CODE ICI
            // Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
            Ray reflexion;
            reflexion.origin = hit.position;
            reflexion.direction = -normalize(2 * dot(ray.direction, hit.normal) * hit.normal - ray.direction); //formule de reflection mirroir
            double reflexionDepth = scene.camera.z_far;
            trace(scene, reflexion, ray_depth+1, &reflexion_color, &reflexionDepth);
            // @@@@@@ VOTRE CODE ICI
            // Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
            //
            // Assumez que l'extérieur/l'air a un indice de réfraction de 1.
            //
            // Toutes les géométries sont des surfaces et non pas de volumes.
            Ray refraction;
            refraction.origin = hit.position;
            double eta = 1.0/material.refractive_index;
           double cos_theta_i = dot(ray.direction, hit.normal);

            double3 perpendicular = eta * (ray.direction - cos_theta_i * hit.normal);

            double parallel_length = -sqrt(1.0 - pow(length(perpendicular),2));
            double3 parallel = parallel_length * hit.normal;

            refraction.direction = normalize(perpendicular + parallel);
            double refractionDepth = scene.camera.z_far;
            trace(scene, refraction, ray_depth + 1, &refraction_color, &refractionDepth);

        }
		
		*out_color = shade(scene, hit) + reflexion_color * material.k_reflection + refraction_color * material.k_refraction;
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
    Material& material = ResourceManager::Instance()->materials[hit.key_material]; // Assuming you have access to the material data.

    double3 viewDirection = normalize(scene.camera.position - hit.position);

    double3 color = material.color_albedo;

     // si une texture est présente, on récupère la couleur aux coordonées UV
    if (material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) {
        // Calculate UV coordinates
        unsigned int texture_width = material.texture_albedo.width()-1;
        unsigned int texture_height = material.texture_albedo.height()-1;
        unsigned int u = (unsigned int)(hit.uv.x * texture_width);
        unsigned int v = (unsigned int)(hit.uv.y * texture_height);
        // Sample texture color using UV coordinates
        unsigned char red, green, blue;
        material.texture_albedo.get_pixel(u, v, red, green, blue);

        // Convert color components to double values and normalize them
        color = { (double)red / 255.0, (double)green / 255.0, (double)blue / 255.0 };
    }

    double3 lightDirection;
    double3 lightDiffuse = {0.0,0.0,0.0};
    double3 lightSpecular = {0.0,0.0,0.0};
    double3 currentLight = {0.0,0.0,0.0};
    double3 finalLight = {0.0,0.0,0.0};
    double3 ambient = {0.0, 0.0, 0.0};


    for(auto light : scene.lights) {
        lightDirection = light.position - hit.position;
        double lightDistance = length(lightDirection);
        double radius = light.radius;
        double lightContribution = 1.0;
        lightDirection = normalize(lightDirection);

        Ray shadowRay;
        shadowRay.origin = hit.position;
        shadowRay.direction = lightDirection;

        double out_umbra_depth;
        Intersection hitUmbra;

        

// Check if the light source has a radius
        if(light.radius > 0.0) {
            // Sample points within the light source's area for soft shadows
            int numSamples = 5; // Number of samples
            int hitCount = 0;

            for(int i = 0; i < numSamples; i++) {
                // Generate random offset within a disk
                double2 randomOffset = random_in_unit_disk() * light.radius;

                // Compute sample point within the light source's area
                double3 samplePoint = light.position + randomOffset.x * normalize(cross(lightDirection, {0, 1, 0})) +
                                    randomOffset.y * normalize(cross(lightDirection, {1, 0, 0}));

                // Compute distance from hit point to sample point
                double samplePointDistance = length(samplePoint - hit.position);

                // Check if the shadow ray is occluded by geometry
                Ray shadowRay;
                shadowRay.origin = hit.position;
                shadowRay.direction = normalize(samplePoint - hit.position);
                Intersection hitUmbra;

                if(scene.container->intersect(shadowRay, EPSILON, samplePointDistance, &hitUmbra)) {
                    hitCount++;
                }
            }

            // Calculate soft shadow coefficient based on number of unoccluded samples
            lightContribution = 1.0 - double(hitCount) / numSamples;
        }
        else{
            if(scene.container->intersect(shadowRay,EPSILON,lightDistance,&hitUmbra)) lightContribution = 0.0;
        }

        double lambertCoef = std::max(dot(hit.normal, lightDirection),0.0);

        double3 halfwayVec = normalize(viewDirection + lightDirection);
        double specularCoef = pow(std::max(dot(hit.normal, halfwayVec),0.0), material.shininess);

        lightDiffuse = material.k_diffuse * color * lambertCoef;

        lightSpecular = material.k_specular * (material.metallic * color + (1 - material.metallic))* specularCoef;//R

        currentLight = (lightDiffuse + lightSpecular)* light.emission/pow(lightDistance,2);
        finalLight += currentLight * lightContribution;
    }

    ambient = scene.ambient_light * material.k_ambient * color;

    finalLight += ambient;

    return finalLight;
}
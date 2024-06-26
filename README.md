
![splash](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/7d059b3d-11d5-4fab-a3a1-9a5d73ffb55b)

**SpatialRay** is an open sourced CPU based pathtraced rendering engine. Based on the BRDF reflectance model, it produces photolike outputs of complex scenes. 

## Demo Renders

![the big small render](https://github.com/CharlesCowdery/RayTracing/assets/54870004/302b9f06-bb9e-4698-b969-b1657ac4e76d)

![scandianIITheScandaning-saturated](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/5875afe9-f6f2-4608-a954-d9b2a906085e)

## Features
- Support for the open source GLTF file standard for 3d scenes.
- Materials are based on the [OpenPBR BRDF reflectance model](https://academysoftwarefoundation.github.io/OpenPBR/) and [(Heitz)](https://inria.hal.science/hal-00942452v1/document).
  - Specularity modeled by IOR via the GGX Microfacet BRDF model.
  - Supports Full material transmission
  - Utilizes the Oren-Nayar diffuse equation.
- Includes AI denoising via OpenImageDenoise
- Full material texturing support, including normal mapping.
- Vertex smoothing/smooth shading.
- Direct light sampling for all implemented puctual light types.
- Utilizes the OCIO standard to implement color spaces.
- Adaptive sampling allows for noise aware rendering
  - The noise analysis strategy is based roughly on [(Dammertz, Hanika)](https://www.semanticscholar.org/paper/A-Hierarchical-Automatic-Stopping-Condition-for-Dammertz-Hanika/8329759ae51c924557f375707e4989549c6c1b46). 
  - Based on a region priority model that directly relates samples per pass of each region to its noise.

![kitchen](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/c6fe4d82-bf8c-4750-bdb7-9348dbd71569)

![loft_saturated](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/ef6c09f9-3a00-4d08-aa5d-0e64fa53e390)

## Technical Features
- Raycasting is fully multithreaded, and utilizes AVX-256 SIMD hardware acceleration.
- Acceleration is based on an AVX specialized 8 wide QBVH [(Wald, Benthin, Boulos)](https://www.cs.cmu.edu/afs/cs/academic/class/15869-f11/www/readings/wald08_widebvh.pdf).
- Utilizes the visible normal distribution PDF detailed in [(Heitz,d'Eon)](https://inria.hal.science/hal-00996995v1/document#page=11&zoom=100,96,180).
- Russian Roulette based ray termination
- Textures are stored via Morton Tiling for increased cache coherance.
- Multithreaded texture loading.
### Noise Analysis and Adapative Sampling View
![image](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/505a151f-72a0-4809-bb98-f828741650d9)
Red highlighting denotes areas of higher noise, and thus increased sampling priority.

### AI denoising
![image](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/487d462b-4dd9-4acc-b4b0-0993add81ce7)

## Future Improvements
- Subsurface Scattering
- BVH spatial splitting
- Multithreaded QBVH specialized generation
- Improved Monte Carlo Sampling / Multiple Importance Sampling
- Accelerated loading of large textures

## About the project
This is an ongoing educational research project I took on in 2023 to exercise my skills in C++, and to pursue an interest in simulation. As such, this project is composed entirely of my own code, with the exception of the standard library, SFML + Dear ImGui (basic graphics & UI lib) and tinyGLTF (a gltf file parser). 

The primary focus of this is not neccessarily to produce photorealistic outputs, but to instead optimize and iterate on the concept of raytracing. What started as an interest in spatial acceleration structures has long since evolved past a 1-2 month digression. 

I've learned an immeasurable amount about programming from this project, from nitty gritty C++ specifics, to the overall tradeoffs of optimization complexity. It has most importantly taught me how to research. How to understand and breakdown whitepapers, grasping the meaning behind the methodologies. I've come to appreciate the importance of iteration, to push bounds slowly but methodically. I've had enforced upon me strong morals about project management and feature scope. Features cannot be just _tackled_, they must be approached.

# Technical structure
This section will document the loose inner workings of the engine. 
## Core
Spatial Ray at its heart is a branched path-tracer. It does not follow a single paper for its implementation, and is instead a hodgepodge of techniques and intuition.
  ### Data storage
   - All major structures involved directly with ray tracing have a structured and workable form, and a "packaged" form with a smaller memory footprint.
   - Hierachical structures are flattened. Notably the BVH is flattened into an array, and uses array indices instead of pointers.
   - Object inheritance chains have their transforms applied and are added directly into global space.
  ### Tracing a Ray Down the BVH
  Spatialray utilizes a modified QBVH specialized for 8 wide simd. Traversal is sorted by nearest intersection.
  - When a ray is dispatched, its slope, origin, inverse slope, and -position/slope, are loaded into AVX-256 registers (the values are copied into all 8 positions).
  - A stack is initialized with the BVH root index loaded as the first element. An additional stack is created for both intersection distances and tri count.
  - The process of navigation follows:
    - The BVH on the top of the stack is loaded.
    - If a BVH does not have triangles in it, it has BVH children nodes.
      - An AVX version of the branchless Slab method is used to intersect the ray against the current BVH section's children.
      - Intersection distance of the children is sorted nearest first, and then the all intersected children are added onto the stack.
    - If it has tris, an AVX version of the [Moller-Trumbore intersection algorithmn](https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm) is used to determine Tri intersection distance and UV coordinates.
  ### Bounce Processing
  1) After an intersection is found, the material UV of the intersection point is calculated, and the material is sampled at that coordinate.
  2) A normal is calculated via linear interpolation of the vertex normals.
  3) If a material has a tangent space normal map, A TBN matrix is calculated, and used to transfrom the normal pertubation into local space.
  4) Transform matrices are calculated to align vectors from Y up space to normal-up space.
  5) The emission of the material is calculated, and added to the rays return.
  6) The number of branched rays to generate is calculated.
     1) The visible NDF Probability Distribution is sampled to generate a microfacet in the surface.
     2) The specular reflection chance of this facet is calculated via the fresnel equation, and randomly chooses between ray transmission/absorbtion and reflection.
     3) If the ray reflects specularly, the specular BRDF is used to determine the amount of the child rays light is conserved during reflection.
     4) If the ray is absorbed and remitted diffusely, it uses the Oren Nayar diffuse BRDF to determine the amount of light conserved during reflection.
     5) If reflected specularly, the microfacet reflection is used as the child ray direction, otherwise a direction is generated via a cosine weighted hemisphere aligned with the normal.
     6) If the final direction falls below the geometric normal (that of the tris geometry) or the smoothed normal (the theoretical geometry the smoothed/perturbed normal simulates), the ray is rejected and regenerated.
     7) The ray is then probabalistically terminated (rejected without regeneration) via russian roulette, with the probability of termination being that of the percent of light conserved during the bounce.
     8) if not terminated, the ray is emitted and processed, and the light return is added to the parent rays return.


## Important/Useful references
- Zap Andersson, et al. OpenPBR Surface specification, version 0.3, 2024-02-21. - [Link](https://academysoftwarefoundation.github.io/OpenPBR/)
- Eric Heitz, Eugene d’Eon. Importance Sampling Microfacet-Based BSDFs using the Distribution of
Visible Normals. Computer Graphics Forum, 2014. ffhal-00996995v1f - [link](https://inria.hal.science/hal-00996995v1/document)
- Eric Heitz. Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs. \[Research
Report\] RR-8468, 2014. ffhal-00942452v1f - [link](https://inria.hal.science/hal-00942452v1/document)
- Wenzel Jakob, Greg Humphreys. Physically Based Rendering: From Theory To Implementation. 4th edition, 2023. - [link](https://www.pbr-book.org/)
- Naty Hoffman. Physics and Math of Shading. Siggraph 2015. - [Link](https://blog.selfshadow.com/publications/s2015-shading-course/hoffman/s2015_pbs_physics_math_slides.pdf)
- Jakub Boksansky. Crash Course in BRDF Implementation. - [Link](https://boksajak.github.io/files/CrashCourseBRDF.pdf)
- Jacco Bikker. How to Build a BVH - [Link](https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/)
- Dammertz, H., Hanika, J., Keller, A., & Lensch, H.P. (2010). A Hierarchical Automatic Stopping Condition for Monte Carlo Global Illumination. - [Link](https://www.semanticscholar.org/paper/A-Hierarchical-Automatic-Stopping-Condition-for-Dammertz-Hanika/8329759ae51c924557f375707e4989549c6c1b46)




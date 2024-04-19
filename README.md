
![splash](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/7d059b3d-11d5-4fab-a3a1-9a5d73ffb55b)

**SpatialRay** is an open sourced CPU based raytraced rendering engine. Based on the BRDF reflectance model, it produces photolike outputs of complex scenes. 

## Demo Renders

![the big small render](https://github.com/CharlesCowdery/RayTracing/assets/54870004/302b9f06-bb9e-4698-b969-b1657ac4e76d)

![scandianIITheScandaning-saturated](https://github.com/CharlesCowdery/SpatialRay/assets/54870004/5875afe9-f6f2-4608-a954-d9b2a906085e)



## Features
- Support for the open source GLTF file standard for 3d scenes.
- Materials are based on the [OpenPBR BRDF reflectance model](https://academysoftwarefoundation.github.io/OpenPBR/) and [(Heitz)](https://inria.hal.science/hal-00942452v1/document).
  - Specularity modeled by IOR via the GGX Microfacet BRDF model.
  - Utilizes the Oren-Nayar diffuse equation.
- Full material texturing support, including normal mapping.
- Vertex smoothing/smooth shading.
- Direct light sampling for all implemented puctual light types.
- Utilizes the OCIO standard to implement color spaces.
- Adaptive sampling allows for noise aware rendering
  - The noise analysis strategy is based roughly on [(Dammertz, Hanika)](https://www.semanticscholar.org/paper/A-Hierarchical-Automatic-Stopping-Condition-for-Dammertz-Hanika/8329759ae51c924557f375707e4989549c6c1b46). 
  - Based on a region priority model that directly relates samples per pass of each region to its noise.

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

## Future Improvements
- Transmissive materials
- Subsurface Scattering
- BVH spatial splitting
- Multithreaded QBVH specialized generation
- Improved Monte Carlo Sampling / Importance Sampling
- Accelerated loading of large textures

## About the project
This is an ongoing educational research project I took on in 2023 to exercise my skills in C++, and to pursue an interest in simulation. As such, this project is composed entirely of my own code, with the exception of the standard library, SFML + Dear ImGui (basic graphics & UI lib) and tinyGLTF (a gltf file parser). 

The primary focus of this is not neccessarily to produce photorealistic outputs, but to instead optimize and iterate on the concept of raytracing. What started as an interest in spatial acceleration structures has long since evolved past a 1-2 month digression. 

I've learned an immeasurable amount about programming from this project, from nitty gritty C++ specifics, to the overall tradeoffs of optimization complexity. It has most importantly taught me how to research. How to understand and breakdown whitepapers, grasping the meaning behind the methodologies. I've come to appreciate the importance of iteration, to push bounds slowly but methodically. I've had enforced upon me strong morals about project management and feature scope. Features cannot be just _tackled_, they must be approached.

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




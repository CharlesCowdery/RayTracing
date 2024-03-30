
![splash](https://github.com/CharlesCowdery/RayTracing/assets/54870004/03d5a988-0da0-4adf-86bb-f77941aafbc6)

**SpatialRay** is an open sourced CPU based raytraced rendering engine. Based on the BRDF reflectance model, it produces photolike outputs of complex scenes. 

## Demo Renders

![the big small render](https://github.com/CharlesCowdery/RayTracing/assets/54870004/302b9f06-bb9e-4698-b969-b1657ac4e76d)

## Features
- Support for the open source GLTF file standard allows for easy interoperability with other software.
- Materials are based on the industry standard BRDF reflectance model.
- Full material texturing support, including normal mapping.
- Direct light sampling for all puctual light types.
- Utilizes the OCIO standard to implement color spaces

## Technical Features
- Raycasting is fully multithreaded, and utilizes AVX-256 SIMD hardware acceleration.
- Acceleration is based on an AVX specialized 8 wide QBVH [(Wald, Benthin, Boulos)](https://www.cs.cmu.edu/afs/cs/academic/class/15869-f11/www/readings/wald08_widebvh.pdf).
- Textures are stored via Morton Tiling for increased cache coherance.


## About the project
This was an educational project I took on in 2023 to exercise my skills in C++, and to pursue an interest in simulation. As such, this project is composed entirely of my own code, with the exception of the standard library, SFML (a basic graphics lib) and tinyGLTF (a gltf file parser). 

The primary focus of this is not neccessarily to produce photorealistic outputs, but to instead optimize and iterate on the concept of raytracing. I had started this project out of an interest in spatial acceleration structures, and it was meant to be a simple 1-2 month digression. This has since become my most ambitious and challenging project to date, and the mark it has left on my understanding of both programming and research cannot be understated.

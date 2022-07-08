# pore-size
This is a research project for HEMI (JH extreme materials lab). Pore sizes between particles. Begins with 2D, ends with 3D. All the code that has helped me with this.


# daily-log
June 21:
Started everything up, not too much

June 22:
Today I learned basic matrix operations (I have not taken linear algebra) and basic matlab commands. I used these to make a 2D Delaunay triangulation with circles at each of the vertices and find the area of just the triangles in matlab (minus the area of the circles). I started trying to do this in a 3D space - I was able to plot a tetrahedron mesh with spheres on each of the vertices and find the volume of the tetramesh total, but am still figuring out how to find the volume of the tetramesh minus part of the volume of the sphere (I haven't yet learned solid angles and need to learn that to do it).

June 23:
Today I figured out how to use vectors to calculate solid angles, use these angles to calculate surface area on a part of a curved sphere, then use that to find the volume of a sector of a sphere. This allowed me to finish my tetrahedron-area-finding-between-spheres program. I began learning linux and about image processing. Watershed segmentation to map out granular particles. 

June 24:
Began to try watershed image segmentation. Started writing a program that segmented a slideshow of images of particles. This will be necessary to connect to my volume-finding thing (which I will need to binarize he he) 

June 27:
Finished watershed segmentation program. The program now outputs x,y,z of centroids, and also gives mean radii. 

June 28:
Started on making a program to find area of pore sizes in a 2D space. This was the beginning of the approximation program. It uses the watershed image segmentation program. 

June 29:
Continued work on the approximation of pore sizes.

June 30:
Continued work on the approximation of pore sizes.

July 1:
Finished 2D approximation of pore sizes. Began a program to find pore sizes between particles in a 2D space using binarization, which is a far better way to do it. It also uses the watershed code I made earlier. 


July 2:
Finished 2D binarization code. Beginning to debug. We want approximation and binarization areas to be about the same.

July 5:
Continuing to debug code to find pore size. Made some diagrams that may be helpful later (showing pores before impact in a 3D space, including a tetramesh, and pores after impact in a 3D space, with another tetramesh.)

July 6:
Finished debugging pore size code. Approximation and binarization are now within ~400 pixels per area of each other. Made a quick program to show if they are similar–called “checkmatrices.m” Made some diagrams and a quick summary of my 2D approximations, one with spheres and one with binarization. 

July 7:
Finished 3D pore size code (the approximation version). Creates a histogram and everything!! I started on the binarization version (the more accurate one). 

July 8:
Continued binarization code. Optimizing to make it take less time to run (I got it down from ~43 hours to 50 minutes, but… its still gonna need some work.) Also, I’m debugging it… since it doesn’t match the approximation code. Llooll. 


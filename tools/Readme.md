## distance_to_walls_3D.cpp: 

Utility program to calculate distance to walls from binary segmented data. 

A value 0 signify a wall, while every other value is integer distance from the wall, i.e. the minimum distance a ray have to travel to hit the wall from the center of that voxel. Maximum distance of 255 (unsigned char max value).

No dependencies other than the binary file have to be in ../data/


Alternatively, get the precalculated version from http://folk.ntnu.no/mortevas/distance_segmented_castle_512.ubc


## distance_to_walls.cpp: 

Utility program to visualize 2D version of distance_to_walls_3D.cpp.

Press 'D' a few times to grow the white squares (walls), then press 'E' to calculate distances

Left click + drag to move camera, scroll to zoom. 

A black circle will show the minimum distance a ray can move before hitting a wall from the cursor position.

Scroll far enough, and it will print the actual distances. 

Depend on glfw3 and stb_easy_font.h


Use Makefile to compile (tools).
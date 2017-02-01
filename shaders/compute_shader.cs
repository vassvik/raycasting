/*

Compute shader ray tracer of voxels in a uniform grid

Source tracing/stepping along the ray: 
    http://www.cse.chalmers.se/edu/year/2011/course/TDA361/grid.pdf

Source bounding box intersection test: 
    http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection

Others:
    https://github.com/LWJGL/lwjgl3-wiki/wiki/2.6.1.-Ray-tracing-with-OpenGL-Compute-Shaders-(Part-I)
    http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays
    http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-overview/light-transport-ray-tracing-whitted

*/

#version 430

layout (local_size_x = 8, local_size_y = 8) in;                             // size of local work groups
layout(rgba32f) writeonly uniform image2D framebufferImage;    // texture covering window

uniform sampler3D voxelTextureSampler;  // 3D texture for voxe ldata, unsigned char mapped to unit floats

uniform vec3 orig;         // camera position
uniform vec3 f;            // forward camera vector
uniform vec3 r;            // right camera vector
uniform vec3 u;            // up camera vector

uniform ivec3 N;           // number of voxels in each direction
uniform ivec3 chosenVoxel; // chosen voxel from CPU ray tracer from cursor
uniform int chosenFace;    // chosen face of chosen voxel from CPU
uniform int textureMode;

uniform float half_width;  // half the frustum width at 1 depth unit
uniform float half_height; // half the frustum height at 1 depth unit
uniform float fogDistance;

uniform int dist_clip;

// Color of each cube face
uniform vec3 colors[6] = {
    vec3(0.55, 0.9, 0.9), 
    vec3(0.9, 0.55, 0.9), 
    vec3(0.9, 0.9, 0.5), 
    vec3(0.5, 0.9, 0.5), 
    vec3(0.5, 0.5, 0.9), 
    vec3(0.9, 0.6, 0.6)
};

// local right-handed coordinate system per face
uniform ivec3 face_vectors[3][6] = {
    {   // tangents
        ivec3( 0,  0,  1), 
        ivec3( 1,  0,  0), 
        ivec3( 0,  1,  0), 
        ivec3( 0,  1,  0), 
        ivec3( 0,  0,  1), 
        ivec3( 1,  0,  0)
    },
    {   // bitangents
        ivec3( 0,  1,  0), 
        ivec3( 0,  0,  1), 
        ivec3( 1,  0,  0), 
        ivec3( 0,  0,  1), 
        ivec3( 1,  0,  0), 
        ivec3( 0,  1,  0)
    },
    {   // normals
        ivec3(-1,  0,  0), 
        ivec3( 0, -1,  0), 
        ivec3( 0,  0, -1), 
        ivec3( 1,  0,  0), 
        ivec3( 0,  1,  0), 
        ivec3( 0,  0,  1)
    }
};

struct Hit {
    vec3 color;
    vec3 hit_pos;
    bool hit;
    int face;
    float t;
};


#define min3(t) min(min(t.x, t.y), t.z)
#define max3(t) max(max(t.x, t.y), t.z)

//#define testVoxel(pos) (texture( voxelTextureSampler, ((pos + vec3(0.5, 0.5, 0.5))/vec3(N))).x <= 0.001953125)
#define testVoxel(pos) (texture( voxelTextureSampler, ((pos + vec3(0.5, 0.5, 0.5))/vec3(N))).x <= (dist_clip+0.5)/256.0)

/*
bool testVoxel(ivec3 pos) {
    return texture( voxelTextureSampler, ((pos + vec3(0.5, 0.5, 0.5))/vec3(N))).x == 0.0;
}
*/


vec3 colorFromFace2(int face, ivec3 pos, vec3 rest_pos, float t) {
    // Distance to walls calculations, used for outline rendering
    vec3 e = 2*rest_pos - vec3(1.0, 1.0, 1.0);
    float distance_from_edge1 = dot(e, face_vectors[0][face]);
    float distance_from_edge2 = dot(e, face_vectors[1][face]);
    vec2 distance_from_edges = vec2(distance_from_edge1, distance_from_edge2);
    float max_distance_from_edges = 1.0 - max(abs(distance_from_edges[0]), abs(distance_from_edges[1]));
    float min_distance_from_edges = 1.0 - min(abs(distance_from_edges[0]), abs(distance_from_edges[1]));

    // outline width scaling by distance from camera and field of view
    ivec2 size = imageSize(framebufferImage);
    float z = 50*half_height/size.y;
    float w = 60*z*t/1024 + z*(1-t/1024);
    float a = atan(distance_from_edges.y, distance_from_edges.x);
    w = max(w, 0.001)*(1 + 0.5*cos(2*3.141592*min_distance_from_edges + 12*3.14145*a));

    // color voxel under cursor specially
    if (chosenVoxel == pos) {
        if (face == -1 || max_distance_from_edges > w) {
            return vec3(1.0, 0.5, 0.0);
        } else {
            return vec3(0.01, 0.01, 0.01);
        }
    }

    // check if boundary of bounding box
    if (face == -1) 
        return vec3(1.0, 0.0, 0.0);

    // check if inside
    if (max_distance_from_edges > w) {
        return colors[face];
        float A = (pos.x + pos.y + pos.z)*10;

        float cx = distance_from_edge1*cos(A) - 0.0 + 0.5*cos((pos.x - 256)/1.0);
        float cy = distance_from_edge2*cos(A) + 0.5*cos((pos.y - 256)/1.0);


        float tmp1 = cx;
        float tmp2 = cy;

        cx = tmp1*cos(A) - tmp2*sin(A);
        cy = tmp1*sin(A) + tmp2*cos(A);

        float x = 0.0;
        float y = 0.0;

        int iter = 0;
        int maxiter = 50;

        while (x*x + y*y < 4 && iter < maxiter) {
            float tmp = x*x - y*y + cx;
            y = 2*x*y + cy;
            x = tmp;
            iter++;
        }

        if (iter == maxiter) {
            float A = (pos.x + pos.y + pos.z)*10;

            float cx = cos(A);
            float cy = sin(A);


 

            float x = distance_from_edge1*cos(A) - 0.0 + 0.5*cos((pos.x - 256)/1.0);
            float y = distance_from_edge2*cos(A) + 0.5*cos((pos.y - 256)/1.0);

            int iter = 0;
            int maxiter = 50;

            while (x*x + y*y < 4 && iter < maxiter) {
                float tmp = x*x - y*y + cx;
                y = 2*x*y + cy;
                x = tmp;
                iter++;
            }
            if (iter == maxiter)
                return vec3(0.0, 0.0, 0.0);
            else
                return vec3(cos(1.0*iter), cos(2.0*iter), cos(3.0*iter));
        }
        else
            return vec3(cos(1.0*iter), cos(2.0*iter), cos(3.0*iter));


        return colors[face];
    }

    // Check edges
    for (int i = 0; i < 2; i++) {
        for (int j = -1; j <= 1; j += 2) {
            if (j == 1 ? distance_from_edges[i] > (1.0 - w) : distance_from_edges[i] < -(1.0 - w)) {
                if (testVoxel(pos+j*face_vectors[i][face]+face_vectors[2][face])) {
                    return vec3(0.0, 0.0, 0.0);
                } else {
                    if (!testVoxel(pos+j*face_vectors[i][face])) 
                        return vec3(0.0, 0.0, 0.0);
                }
            }  
        }
    }

    // check corners
    for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
            float dx = abs(i - distance_from_edges[0]);
            float dy = abs(j - distance_from_edges[1]);

            if (dx*dx + dy*dy < w*w) {
                if (testVoxel(pos+i*face_vectors[0][face]+j*face_vectors[1][face]+face_vectors[2][face])) {
                    return vec3(0.0, 0.0, 0.0);
                } else {
                    if (!testVoxel(pos+i*face_vectors[0][face]+j*face_vectors[1][face])) {
                        return vec3(0.0, 0.0, 0.0);
                    }
                }
            }  
        }
    }   


    return colors[face];

}

Hit traceWrap(vec3 orig, vec3 dir) {
    // Setup intersection test
    vec3 invdir = 1.0/dir;
    ivec3 sign = ivec3((dir.x < 0 ? 1 : 0), (dir.y < 0 ? 1 : 0), (dir.z < 0 ? 1 : 0));

    // ray hit a voxel (wraps), find out which
    float tmin = 0.0;
    ivec3 pos = ivec3(floor(orig));                       // determine the position of the first hit voxel
    
    int face = -1;              // assume illegal face (boundary)
    bool out_of_bounds = true;  // assume ray goes through voxels without hitting
    float tnew = tmin;          // value for t which hit first/current voxel

    if (testVoxel(pos)) {
        // test if first voxel is filled, early exit
        // could add face detection instead of face == -1
        out_of_bounds = false;
    } else {
        // first voxel is not filled, iterate along the ray until one is found or out of bounds
        vec3 tDelta = abs(invdir);                       // the amount to move along the ray to cross a voxel in each direction
        ivec3 step = -2*sign + 1;
        vec3 t = ((pos + 0.5 + 0.5*step) - orig)*invdir; // the distance to move *to* along the ray for which the *next* voxel is hit, in each direction
        int ctr = 0;
        while (tnew < fogDistance) {
            // while inside bounding box

            if (testVoxel(pos)) { // HIT!
                out_of_bounds = false;
                break;
            }

            // NO HIT! Move to next voxel

            tnew = min3(t);

            if (tnew == t.x) {
                pos.x = pos.x + step.x;
                t.x = t.x + tDelta.x;
                face = 0 + 3*sign.x;
            } else if (tnew == t.y) {
                pos.y = pos.y + step.y;
                t.y = t.y + tDelta.y;
                face = 1 + 3*sign.y;
            } else if (tnew == t.z) {
                pos.z = pos.z + step.z;
                t.z = t.z + tDelta.z;
                face = 2 + 3*sign.z;
            }

        }  
    }

    // COLORING
    vec3 rest_pos = clamp(orig + tnew*dir - pos, 0.0, 1.0); // position relative to cube, cubemap lookup

    float dist = clamp(1.0 - tnew/fogDistance, 0.0, 1.0);
     if (pos == chosenVoxel)
        dist = 1.0;

    Hit ret;
    ret.color = colorFromFace2(face, pos, rest_pos, tnew)*dist + vec3(0.8, 0.8, 0.8)*(1.0 - dist);
    ret.hit_pos = orig + tnew*dir;
    ret.hit = true;
    ret.face = face;
    ret.t = tnew;
    return ret;
}



Hit trace(vec3 orig, vec3 dir) {
    // Setup intersection test
    vec3 invdir = 1.0/dir;
    ivec3 sign = ivec3((dir.x < 0 ? 1 : 0), (dir.y < 0 ? 1 : 0), (dir.z < 0 ? 1 : 0));


    vec3 t_back = (N*sign - orig)*invdir;       // distance to back planes (relative to direction)
    vec3 t_front = (N*(1-sign) - orig)*invdir;   // distance to front planes
    float tmin = max3(t_back);                  // distance to closest plane behind camera
    float tmax = min3(t_front);                  // diistance to closest plane in front camera

    // (tmin > t_front.y || t_back.y > tmax || tmin > t_front.z || t_back.z > tmax) checks for bounding box test 
    // (max(tmin, tmax) < 0) remove intersecting "inverse" when looking away from the voxels
    if (tmin > t_front.y || t_back.y > tmax || tmin > t_front.z || t_back.z > tmax || max(tmin, tmax) < 0) {
        // no intersection of bounding box, ray did not hit a voxel
    } else { 
        // ray hit a voxel, find out which
        tmin = (tmin < 0 && tmax > 0 ? 0.0 : tmin);      // if inside bounding box, set tmin to zero.  

        ivec3 pos = ivec3(orig + tmin*dir);     // determine the position of the first hit voxel
        pos = clamp(pos, ivec3(0,0,0), N-1);    // clamp to valid ranges, in case of loss of floating point precision
        
        int face = -1;              // assume illegal face (boundary)
        bool out_of_bounds = true;  // assume ray goes through voxels without hitting
        float tnew = tmin;          // value for t which hit first/current voxel

        if (testVoxel(pos)) {
            // test if first voxel is filled, early exit
            // could add face detection instead of face == -1
            out_of_bounds = false;
        } else {
            // first voxel is not filled, iterate along the ray until one is found or out of bounds
            vec3 tDelta = abs(invdir);                       // the amount to move along the ray to cross a voxel in each direction
            ivec3 step = -2*sign + 1;
            vec3 t = ((pos + 0.5 + 0.5*step) - orig)*invdir; // the distance to move *to* along the ray for which the *next* voxel is hit, in each direction

            while (pos.x < N.x && pos.y < N.y && pos.z < N.z && pos.x >= 0 && pos.y >= 0 && pos.z >= 0) {
                // while inside bounding box

                if (testVoxel(pos)) { // HIT!
                    out_of_bounds = false;
                    break;
                }

                // NO HIT! Move to next voxel

                tnew = min3(t);

                if (tnew == t.x) {
                    pos.x = pos.x + step.x;
                    t.x = t.x + tDelta.x;
                    face = 0 + 3*sign.x;
                } else if (tnew == t.y) {
                    pos.y = pos.y + step.y;
                    t.y = t.y + tDelta.y;
                    face = 1 + 3*sign.y;
                } else if (tnew == t.z) {
                    pos.z = pos.z + step.z;
                    t.z = t.z + tDelta.z;
                    face = 2 + 3*sign.z;
                }

            }  
        }

        // COLORING
        if (!out_of_bounds) {
            vec3 rest_pos = clamp(orig + tnew*dir - pos, 0.0, 1.0); // position relative to cube, cubemap lookup


            float dist = clamp(1.0 - tnew/fogDistance, 0.0, 1.0); 
            if (pos == chosenVoxel)
                dist = 1.0;

            // blend face color with fog based on distance
            Hit ret;
            ret.color = colorFromFace2(face, pos, rest_pos, tnew)*dist + vec3(0.8, 0.8, 0.8)*(1.0 - dist);
            ret.hit_pos = orig + tnew*dir;
            ret.hit = true;
            ret.face = face;
            ret.t = tnew;
            return ret;
            
        }
    }

    Hit ret;
    ret.color = vec3(0.0, 0.0, 0.0);
    ret.hit_pos = vec3(-1.0, -1.0, -1.0);
    ret.hit = false;
    ret.face = -1;
    ret.t = -1;
    return ret;

    
}



/*
    Compute shader -- ray tracing of voxels
     - Invoked once for every pixel in the framebuffer/window. 
     - Shoots a ray from the camera position until it hits a filled voxel
     - If no voxels are hit, color pixel black
       Otherwise, color it according to which face it hits
     - Color is written to a framebuffer/image/texture
*/
void main() {
    ivec2 pixel = ivec2(gl_GlobalInvocationID.xy);
    ivec2 size = imageSize(framebufferImage);

    if (pixel.x >= size.x || pixel.y >= size.y)
        return; // invocation out of bounds (size not multiple of work size)
 
    // normalized device coordinates each pixel
    float dx = 2.0*(pixel.x+0.5)/float(size.x) - 1.0;
    float dy = 2.0*(pixel.y+0.5)/float(size.y) - 1.0;

    // calculate ray for pixel
    vec3 dir = f + half_height*dy*u + half_width*dx*r;
    vec3 color;
    Hit ret;
    if (textureMode == 0)
        ret = trace(orig,dir);
    else
        ret = traceWrap(orig,dir);

    imageStore(framebufferImage, pixel, vec4(ret.color, 1.0));
    return;

    vec3 lightPos = vec3(20,20, 10);
    vec3 lightDir = normalize(ret.hit_pos - lightPos);
    Hit ret2;

    if (textureMode == 0)
        ret2 = trace(lightPos,lightDir);
    else
        ret2 = traceWrap(lightPos,lightDir);




    float occlusionDist = 1000000;
    if (ret.hit && ret.face != -1) {
        vec3 normalDir = face_vectors[2][ret.face];

        Hit ret3 = trace(ret.hit_pos + 0.01*normalDir, normalDir);
        if (ret3.hit)
            occlusionDist = ret3.t;
    }

    float d = 30;

    if (ret.hit && length(ret2.hit_pos - ret.hit_pos) < 0.01)
        imageStore(framebufferImage, pixel, vec4(ret.color*0.5  *min(occlusionDist/d, 1.0), 1.0));
    else
        imageStore(framebufferImage, pixel, vec4(ret.color*min(occlusionDist/d, 1.0), 1.0));
    // write to frame buffer texture
     //   color = vec3(1,0,0);
    //imageStore(framebufferImage, pixel, vec4(ret.color, 1.0));

}
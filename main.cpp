#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <imgui.h>
#include "imgui_impl_glfw_gl3.h"
#include "imgui_impl_glfw_gl3.cpp"
#define IM_ARRAYSIZE(_ARR)  ((int)(sizeof(_ARR)/sizeof(*_ARR)))

#include <GLFW/glfw3.h>
#include <gl_math.h>


#include <utils.h>


int mod2(double i, int n) {
    return i - n*floor(i/n);
}
int mod2(double i, double n) {
    return i - n*floor(i/n);
}

double time_raytrace = 0;
double time_blocks = 0;
double alpha_time = 0.01;
double time_uploadXY = 0.0;
double time_uploadXZ = 0.0;
double time_uploadYZ = 0.0;


int intersect = 1;
double t = 0;


bool removeBlock = false;
bool addBlock = false;


GLFWwindow* window;
double prevx = -1, prevy = -1;
double resx = 1280,resy = 720;

int clickedButtons = 0;
enum buttonMaps { FIRST_BUTTON=1, SECOND_BUTTON=2, THIRD_BUTTON=4, FOURTH_BUTTON=8, FIFTH_BUTTON=16, NO_BUTTON=0 };
enum modifierMaps { CTRL=2, SHIFT=1, ALT=4, META=8, NO_MODIFIER=0 };

GLuint programID;
GLuint VertexArrayID;

GLuint computeProgramID, computerRenderProgramID;

GLuint texture, texture3D, textureCubemap;

//GLuint textureMode = GL_MIRRORED_REPEAT;
//GLuint textureMode = GL_REPEAT;
GLint textureMode = GL_CLAMP_TO_BORDER;


// vertex position for two triangles (using TRIANGLE_STRIP) covering the whole screen
GLuint vertexbuffer_quad;
static const GLfloat g_vertex_buffer_data[] = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
    -1.0f,  1.0f, 0.0f,
     1.0f,  1.0f, 0.0f,
};

// uv texture coordinates for the vertices of the two triangles
GLuint uvbuffer_quad;
static const GLfloat g_uv_buffer_data[] = {
    0.0f, 0.0f,
    1.0f, 0.0f,
    0.0f, 1.0f,
    1.0f, 1.0f,
};




vec3 cam_pos = vec3(20.0f, 20.0f, 2.0f);
vec3 f, u, r;
float phi = PI/4;
float theta = PI/2;
float fov = 60.0f;

unsigned char dist_clipped = 0;


char *readFile(const char *filename);
GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path);
GLuint LoadComputeShader(const char * compute_file_path);
void initGL();
void Draw();
void cameraStuff();

ivec3 trace(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out);
ivec3 traceWrap2(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out);
ivec3 traceWrap(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out);

void key_callback(GLFWwindow* win, int key, int scancode, int action, int mods);
void mousebutton_callback(GLFWwindow* win, int button, int action, int mods);
void mousepos_callback(GLFWwindow* win, double xpos, double ypos);
void mousewheel_callback(GLFWwindow* win, double xoffset, double yoffset);
void windowsize_callback(GLFWwindow *win, int width, int height);

const int Nx = 512, Ny = Nx, Nz = Nx;
ivec3 N = ivec3(Nx, Ny, Nz);

float fogDistance = 128;

unsigned char*data;

ivec3 pos;
int chosenFace;

// Color of each cube face
vec3 colors[6] = {
    vec3(0.55, 0.9, 0.9), 
    vec3(0.9, 0.55, 0.9), 
    vec3(0.9, 0.9, 0.5), 
    vec3(0.5, 0.9, 0.5), 
    vec3(0.5, 0.5, 0.9), 
    vec3(0.9, 0.6, 0.6)
};


void updateComputeShader() {
    double t1 = glfwGetTime();
    GLuint computeProgramID_tmp = LoadComputeShader("shaders/compute_shader.cs");
    glDeleteProgram(computeProgramID);
    computeProgramID = computeProgramID_tmp;
    double t2 = glfwGetTime();
    printf("update compute, time spent = %fms\n", (t2-t1)*1000); fflush(stdout);
}


void doInput() {
    ImGuiIO& io = ImGui::GetIO();


    if (!io.WantCaptureMouse) {
        if (io.MouseDown[0]) {
            double anglesPerPixel = 0.1 * PI / 180.0;
            phi = fmod(phi - io.MouseDelta.x * anglesPerPixel + 2*PI, 2*PI);   // periodic
            theta = max(1.0*PI/180.0, min(179.0*PI/180.0, theta + io.MouseDelta.y * anglesPerPixel));
        }

        if (io.MouseWheel != 0.0) {
            double zoomFactor = pow(0.95, io.MouseWheel);
            if (io.KeysDown[GLFW_KEY_LEFT_CONTROL]) {
                fogDistance *= zoomFactor;
            } else if (io.KeysDown[GLFW_KEY_LEFT_SHIFT]) {
                dist_clipped = mod2(dist_clipped + io.MouseWheel, 255);
                printf("%d\n", dist_clipped); fflush(stdout);
            } else {
                fov = max(15, min(120.0, fov*zoomFactor));
            }
        }

        double t1 = glfwGetTime();

        double aspect = resx/resy;
        double FOV = fov*PI/180.0;
        double half_height = tan(FOV/2.0);
        double half_width = tan(FOV/2.0)*aspect;

        double screen_posx, screen_posy;
        glfwGetCursorPos(window,&screen_posx, &screen_posy); 
        double posx = 2*CLAMP(screen_posx/resx, 0.0, 1.0)-1;
        double posy = 2*CLAMP((resy - screen_posy)/resy, 0.0, 1.0)-1;

        vec3 dir = f + half_height*posy*u + half_width*posx*r;
        
        bool out_of_bounds;
        if (textureMode == GL_REPEAT)
            pos = traceWrap(cam_pos, dir, out_of_bounds, chosenFace);
        else if (textureMode == GL_MIRRORED_REPEAT)
            pos = traceWrap2(cam_pos, dir, out_of_bounds, chosenFace);
        else
            pos = trace(cam_pos, dir, out_of_bounds, chosenFace);


        double t2 = glfwGetTime();
        time_raytrace = 1e9*(t2-t1);

        if (!out_of_bounds) {
            glBindTexture(GL_TEXTURE_3D, texture3D);
            ivec3 pos2 = pos;
            ivec3 pos3 = pos;

            if (textureMode == GL_MIRRORED_REPEAT || textureMode == GL_REPEAT) {
                pos2.x = mod2(pos.x, Nx);
                pos2.y = mod2(pos.y, Ny);
                pos2.z = mod2(pos.z, Nz);
            }

            if (textureMode == GL_MIRRORED_REPEAT) {
                pos3.x = floor(pos.x/float(Nx));
                pos3.y = floor(pos.y/float(Ny));
                pos3.z = floor(pos.z/float(Nz));

                if (pos3.x % 2 != 0) pos2.x = Nx-1 - pos2.x;
                if (pos3.y % 2 != 0) pos2.y = Ny-1 - pos2.y;
                if (pos3.z % 2 != 0) pos2.z = Nz-1 - pos2.z;
            }
            int X = pos2.x, Y = pos2.y, Z = pos2.z;

            if (removeBlock) {
                unsigned char values1[Ny*Nz], values2[Nx*Nz], values3[Nx*Ny];
                for (int i = 0; i < Ny*Nz; i++) values1[i] = 255;
                for (int i = 0; i < Nx*Nz; i++) values2[i] = 255;
                for (int i = 0; i < Nx*Ny; i++) values3[i] = 255;

                glPixelStorei( GL_UNPACK_ALIGNMENT, 1);

                double tt1 = glfwGetTime();
                if (1) {
                    for (int j = 0; j < Ny; j++) 
                        for (int k = 0; k < Nz; k++) 
                            data[k*Nx*Ny+j*Nx+pos2.x] = 255;
                    glTexSubImage3D(GL_TEXTURE_3D, 0,  pos2.x, 0, 0,  1, Ny, Nz,  GL_RED, GL_UNSIGNED_BYTE, &values1[0]);
                } 
                double tt2 = glfwGetTime();
                if (1) {
                    for (int i = 0; i < Nx; i++) 
                        for (int k = 0; k < Nz; k++) 
                            data[k*Nx*Ny+pos2.y*Nx+i] = 255;
                    glTexSubImage3D(GL_TEXTURE_3D, 0,  0, pos2.y, 0,  Nx, 1, Nz,  GL_RED, GL_UNSIGNED_BYTE, &values2[0]);
                } 
                double tt3 = glfwGetTime();
                if (1) {
                    for (int i = 0; i < Nx; i++) 
                        for (int j = 0; j < Ny; j++) 
                            data[pos2.z*Nx*Ny+j*Nx+i] = 255;
                    glTexSubImage3D(GL_TEXTURE_3D, 0,  0, 0, pos2.z,  Nx, Ny, 1,  GL_RED, GL_UNSIGNED_BYTE, &values3[0]);
                }

                double tt4 = glfwGetTime();

                if (time_uploadXY == 0.0) {
                    time_uploadXY = tt4-tt3;
                    time_uploadXZ = tt3-tt2;
                    time_uploadYZ = tt2-tt1;
                } else {
                    time_uploadXY = alpha_time*(tt4-tt3) + (1.0 - alpha_time)*time_uploadXY;
                    time_uploadXZ = alpha_time*(tt3-tt2) + (1.0 - alpha_time)*time_uploadXZ;
                    time_uploadYZ = alpha_time*(tt2-tt1) + (1.0 - alpha_time)*time_uploadYZ;
                }
                
            } else if (addBlock) {
                if (chosenFace != -1) {
                    if (chosenFace == 0) {
                        X = mod2(pos2.x-1, Nx);
                    } else  if (chosenFace == 1) {
                        Y = mod2(pos2.y-1, Ny);
                    } else  if (chosenFace == 2) {
                        Z = mod2(pos2.z-1, Nz);
                    } else  if (chosenFace == 3) {
                        X = mod2(pos2.x+1, Nx);
                    } else  if (chosenFace == 4) {
                        Y = mod2(pos2.y+1, Ny);
                    } else  if (chosenFace == 5) {
                        Z = mod2(pos2.z+1, Nz);
                    }
                    data[Z*Nx*Ny + Y*Nx + X] = 0;
                    glTexSubImage3D(GL_TEXTURE_3D, 0, X, Y, Z, 1, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &data[Z*Nx*Ny + Y*Nx + X]);  
                }
            }

            double t3 = glfwGetTime();
            if (addBlock || removeBlock)
                time_blocks = 1e3*(t3-t1);
        }
        removeBlock = false;
        addBlock = false;

    }
    if (!io.WantCaptureKeyboard) {
        if (ImGui::IsKeyPressed(GLFW_KEY_TAB, false)) {
            updateComputeShader();
        }
        if (ImGui::IsKeyPressed(GLFW_KEY_F5, false)) {
            printf("screen\n"); fflush(stdout);
        }

        if (ImGui::IsKeyPressed(GLFW_KEY_1, false)) {
            textureMode = GL_CLAMP_TO_BORDER;
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, textureMode);
        }

        if (ImGui::IsKeyPressed(GLFW_KEY_2, false)) {
            textureMode = GL_REPEAT;
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, textureMode);
        }
        if (ImGui::IsKeyPressed(GLFW_KEY_3, false)) {
            textureMode = GL_MIRRORED_REPEAT;
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, textureMode);
        }

        if (ImGui::IsKeyPressed(GLFW_KEY_SPACE, true)) {
            if (io.KeysDown[GLFW_KEY_LEFT_CONTROL]) {
                addBlock = true;
            } else {
                removeBlock = true;
            }
        }

    }

    if (!io.WantCaptureKeyboard) {
        cameraStuff();
    }
}

void doGUI() {
    static bool show_test_window = false;


    // 1. Show a simple window
    // Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appears in a window automatically called "Debug"
    {
        const char* side_labels[] = {"N/A", "-X", "-Y", "-Z", "+X", "+Y", "+Z"};

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        ImGui::Text("res = %dx%d", (int)resx, (int)resy);
        ImGui::Text("%dx%dx%d = %d voxels", Nx, Ny, Nz, Nx*Ny*Nz);
        ImGui::Text("time raytrace = %f ns", time_raytrace);
        ImGui::Text("time adding/removing blocks = %f ms", time_blocks);
        ImGui::Text("distance to block = %f", t);
        ImGui::Text("pos = (%.3f, %.3f, %.3f)", cam_pos.x, cam_pos.y, cam_pos.z);
        ImGui::Text("dir = (%.3f, %.3f, %.3f)", f.x, f.y, f.z);
        ImGui::Text("angles = (%.3f, %.3f)", phi*180/PI, theta*180/PI);
        ImGui::Text("time upload XY slice = %.3f ms", 1000*time_uploadXY);
        ImGui::Text("time upload XZ slice = %.3f ms", 1000*time_uploadXZ);
        ImGui::Text("time upload YZ slice = %.3f ms", 1000*time_uploadYZ);









        double hfov = atan(tan((fov*PI/180.0)/2.0)*resx/resy)*180/PI*2;
        ImGui::Text("hfov = %.3f, vfov = %.3f", hfov, fov);
        // convert chosen voxel position in case of texture wrapping
        ivec3 pos2 = pos;
        if (textureMode == GL_MIRRORED_REPEAT || textureMode == GL_REPEAT) {
            pos2.x = mod2(pos.x, Nx);
            pos2.y = mod2(pos.y, Ny);
            pos2.z = mod2(pos.z, Nz);
        }

        if (textureMode == GL_MIRRORED_REPEAT) {
            ivec3 pos3 = pos;
            pos3.x = floor(pos.x/float(Nx));
            pos3.y = floor(pos.y/float(Ny));
            pos3.z = floor(pos.z/float(Nz));

            if (pos3.x % 2 != 0) pos2.x = Nx-1 - pos2.x;
            if (pos3.y % 2 != 0) pos2.y = Ny-1 - pos2.y;
            if (pos3.z % 2 != 0) pos2.z = Nz-1 - pos2.z;
        }

        if ((pos.x < 0 || pos.y < 0 || pos.z < 0 || pos.x > Nx-1 || pos.y > Ny-1 || pos.z > Nz-1) && (textureMode == GL_CLAMP_TO_BORDER)) {
            ImGui::Text("chosen voxel = OUT_OF_BOUNDS");
            ImGui::Text("chosen face = OUT_OF_BOUNDS");
        } else {
            ImGui::Text("chosen voxel = (%d, %d, %d)", pos2.x, pos2.y, pos2.z);
            ImGui::Text("chosen face = %s", side_labels[chosenFace+1]);
        }

        if (ImGui::CollapsingHeader("Adjustables")) {
            for (int i = 0; i < 6; i++)
                ImGui::ColorEdit3(side_labels[i+1], &colors[i].x);
            
            int dist_clipped_int = dist_clipped;
            ImGui::SliderInt("wall distance", &dist_clipped_int, 0, 20);
            dist_clipped = dist_clipped_int;

            // BUG?: when changing to repeat, clamped?
            const char* texture_mode_labels[] = {"clamped", "repeat", "mirrored repeat"};
            int mode = textureMode == GL_REPEAT ? 1 : textureMode == GL_MIRRORED_REPEAT ? 2 : 0;
            ImGui::Combo("Drawing mode", &mode, texture_mode_labels, IM_ARRAYSIZE(texture_mode_labels));
            textureMode = (mode == 1 ? GL_REPEAT : mode == 2 ? GL_MIRRORED_REPEAT : GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, textureMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, textureMode);

            float fog = fogDistance;
            ImGui::SliderFloat("fog distance", &fog, 0.0, 256.0, "%.3f");
            fogDistance = fog;

            const char* on_off[] = {"OFF", "ON"};
            static int vsync = 1;
            ImGui::Combo("VSYNC", &vsync, on_off, IM_ARRAYSIZE(on_off));
            glfwSwapInterval(vsync);

            if (ImGui::CollapsingHeader("Test")) {
                if (ImGui::Button("Test Window"))    show_test_window ^= 1;
            }
        }   
    }

    if (show_test_window)
    {
        ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
        ImGui::ShowTestWindow(&show_test_window);
    }
}


int main() {

    // If using "segmented_castle_512.ubc" instead of "distance_segmented_castle_512.ubc", no distance estimation
    FILE *fp = fopen("data/distance_segmented_castle_512.ubc", "rb");
    //FILE *fp = fopen("data/segmented_castle_512.ubc", "rb");

    
    
    int size = Nx*Ny*Nz;
    data = new unsigned char[size]; 
    fread(data, sizeof(unsigned char), size, fp);

    
    fclose(fp);


    initGL();

    ImGui_ImplGlfwGL3_Init(window, true);


    while ( !glfwWindowShouldClose(window)) {   
        // Tell OpenGL which program to use (can use multiple)
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        resx = w, resy = h;

        glClearColor(1.0, 1.0, 1.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ImGui_ImplGlfwGL3_NewFrame();



        Draw();
        
        doInput();
        doGUI();

        ImGui::Render();  // DOESNT WORK, WHY? TEXTURE BINDING INTERFERENCE?

        // Swap buffers
        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VertexArrayID);
    glDeleteProgram(programID);

    glfwTerminate();

    return 0;
}


char *readFile(const char *filename) {
    // Read content of "filename" and return it as a c-string.
    printf("Reading %s\n", filename);
    FILE *f = fopen(filename, "rb");

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    printf("Filesize = %d\n", int(fsize));

    char *string = (char*)malloc(fsize + 1);
    fread(string, fsize, 1, f);
    string[fsize] = '\0';
    fclose(f);

    return string;
}


GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path){
    // This function will load the vertex and fragment shader source code from files
    // then compile them and create a program that is run on the GPU every frame. 

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // The vertex shader is called on every vertex, usually transforming it using a projection view matrix
    // MVP = P * V * M
    // Other informaion can be sent to the fragment shader

    // Create the Vertex shader
    GLuint VertexShaderID   = glCreateShader(GL_VERTEX_SHADER);
    char *VertexShaderCode   = readFile(vertex_file_path);

    // Compile Vertex Shader
    printf("Compiling shader : %s\n", vertex_file_path); fflush(stdout);
    glShaderSource(VertexShaderID, 1, (const char**)&VertexShaderCode , NULL);
    glCompileShader(VertexShaderID);

    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ){
        char VertexShaderErrorMessage[9999];
        glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, VertexShaderErrorMessage);
        printf("%s\n", VertexShaderErrorMessage); fflush(stdout);
    }


    // The fragment shader is called on every fragment (pixel/sample), usually for coloring. 

    // Create the Fragment shader
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    char *FragmentShaderCode = readFile(fragment_file_path);

    // Compile Fragment Shader
    printf("Compiling shader : %s\n", fragment_file_path); fflush(stdout);
    glShaderSource(FragmentShaderID, 1, (const char**)&FragmentShaderCode , NULL);
    glCompileShader(FragmentShaderID);

    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
        char FragmentShaderErrorMessage[9999];
        glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, FragmentShaderErrorMessage);
        printf("%s\n", FragmentShaderErrorMessage); fflush(stdout);
    }


    // Create and Link the program
    printf("Linking program\n"); fflush(stdout);
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ){
        GLchar ProgramErrorMessage[9999];
        glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
        printf("%s\n", &ProgramErrorMessage[0]); fflush(stdout);
    }

    // Can delete shaders objects after program has been linked/created
    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);
    free(FragmentShaderCode);
    free(VertexShaderCode);


    return ProgramID;
}

GLuint LoadComputeShader(const char * compute_file_path){

    // Create the shaders
    GLuint ComputeShaderID   = glCreateShader(GL_COMPUTE_SHADER);
    char *ComputeShaderCode   = readFile(compute_file_path);


    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    printf("Compiling shader : %s\n", compute_file_path);
    glShaderSource(ComputeShaderID, 1, (const char**)&ComputeShaderCode , NULL);
    glCompileShader(ComputeShaderID);

    // Check Vertex Shader
    glGetShaderiv(ComputeShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(ComputeShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ){
        char ComputeShaderErrorMessage[9999];
        glGetShaderInfoLog(ComputeShaderID, InfoLogLength, NULL, ComputeShaderErrorMessage);
        printf("%s\n", ComputeShaderErrorMessage);
    }


    // Link the program
    printf("Linking program\n");
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, ComputeShaderID);
    glLinkProgram(ProgramID);

    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ){
        GLchar ProgramErrorMessage[9999];
        glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
        printf("%s\n", &ProgramErrorMessage[0]);
    }

    glDeleteShader(ComputeShaderID);

    glUseProgram(ProgramID);
    
    glUniform1i(glGetUniformLocation(ProgramID, "framebufferImage"), 0);
    glUniform1i(glGetUniformLocation(ProgramID, "voxelTextureSampler"), 1);

    return ProgramID;
}


void cameraStuff() {
    static double t1 = glfwGetTime(); // store previous time
    // "speed" is a misnomer, as it is dependent on the time difference. speed = (t2-t1)*3 means the speed is 3 units per second. 
    double t2 = glfwGetTime();
    double speed = (t2-t1)*2; // units to move since last frame
    
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) {
        speed *= 10;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
        speed *= 10;
    }

    char str[256];
    static double dt = 1.0/60;
    double alpha = 0.05;
    dt = (t2-t1)*alpha + dt*(1-alpha);
    sprintf(str, "fps = %7.1f, time raytrace = %6.1f ns, N = %ld voxels, intersect = %d, t = %f, fov = %f\n", 1.0/dt, time_raytrace, (long int)Nx*Ny*Nz, intersect, t, fov);
    glfwSetWindowTitle(window, str);
    t1 = t2;


    // Movement, detached from callbacks. Only move when WASD or QE is pushed. 
    // W and S is forward/backward
    // A and D is left/right
    // Q and E is down/up (relative to camera direction)

    // forward and backwards
    ///*
    double deltaMoveForward = speed*(glfwGetKey(window, GLFW_KEY_W) - glfwGetKey(window, GLFW_KEY_S));
    double deltaMoveRight   = speed*(glfwGetKey(window, GLFW_KEY_D) - glfwGetKey(window, GLFW_KEY_A));
    double deltaMoveUp      = speed*(glfwGetKey(window, GLFW_KEY_E) - glfwGetKey(window, GLFW_KEY_Q));
    //*/
    

    // Convert camera direction to a right handed coordinate system
    // with one vector pointing in the direction of the camera (y-axis), 
    // one vector poiting to the "right" of this one (x-axis)
    // and one orthogonal to these two, the "up" axis (z-axis). r x f = u
    // 
    // These are automatically normalized. "Easily" derived by hand. Physicist notation is used, i.e. theta is the polar axis. 
    // 
    // These are convenient for moving around and creating the view matrix. 
    double sinp = sin(phi),   cosp = cos(phi);
    double sint = sin(theta), cost = cos(theta);
    f = vec3(cosp*sint, sinp*sint, cost);                   // forward vector, normalized, spherical coordinates
    r = vec3(sinp, -cosp, 0.0f);                            // right vector, relative to forward
    u = vec3(-cosp*cost, -sinp*cost, sint);                 // "up" vector, u = r x f

    // advance position
    cam_pos = cam_pos + deltaMoveForward*f + deltaMoveRight*r + deltaMoveUp*u;

}



ivec3 trace(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out) {
    vec3 invdir = 1.0/dir;
    ivec3 sign = ivec3((dir.x < 0 ? 1 : 0), (dir.y < 0 ? 1 : 0), (dir.z < 0 ? 1 : 0));



    vec3 t_back = (vec3(N*sign) - orig)*invdir;
    vec3 t_front = (vec3(N*(1-sign)) - orig)*invdir;

    float tmin = max3(t_back);
    float tmax = min3(t_front);

    ivec3 pos;
    out = -1;
    
    out_of_bounds = true;  // assume ray goes through voxels without hitting


    if (tmin > t_front.y || t_back.y > tmax || tmin > t_front.z || t_back.z > tmax || max(tmin, tmax) < 0) {
        pos = ivec3(-1, -1, -1);
    } else {
        tmin = (tmin < 0 && tmax > 0 ? 0.0 : tmin); 
        pos = ivec3(orig + tmin*dir);
        pos = clamp(pos, ivec3(0,0,0), -1 + N);

        int face = -1;              // assume illegal face (boundary)
        float tnew = tmin;

        if (data[pos.z*Nx*Ny + pos.y*Nx + pos.x] == 0) {
            out_of_bounds = false;
        } else {
            vec3 tDelta = abs(invdir);                       // the amount to move along the ray to cross a voxel in each direction
            ivec3 step = 1 - 2*sign;
            vec3 t = ((0.5f + 0.5f*vec3(step) + vec3(pos)) - orig)*invdir;
            out = -1;
            while (pos.x < N.x && pos.y < N.y && pos.z < N.z && pos.x >= 0 && pos.y >= 0 && pos.z >= 0) {
                // while inside bounding box

                if (data[pos.z*Nx*Ny + pos.y*Nx + pos.x] == 0) { // HIT!
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

        out = face;
        t = tnew;
    }

    return pos;
}




ivec3 traceWrap2(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out) {
    vec3 invdir = 1.0/dir;
    ivec3 sign = ivec3((dir.x < 0 ? 1 : 0), (dir.y < 0 ? 1 : 0), (dir.z < 0 ? 1 : 0));

    ivec3 pos;
    out = -1;
    
    out_of_bounds = true;  // assume ray goes through voxels without hitting


    float tmin = 0.0;
    pos = ivec3(floor(orig));

    int face = -1;              // assume illegal face (boundary)
    float tnew = tmin;


    ivec3 pos2 = pos;
    pos2.x = mod2(pos.x, Nx);
    pos2.y = mod2(pos.y, Ny);
    pos2.z = mod2(pos.z, Nz);

    ivec3 pos3 = pos;
    pos3.x = pos.x/Nx;
    pos3.y = pos.y/Ny;
    pos3.z = pos.z/Nz;

    if (pos3.x % 2 != 0) pos2.x = Nx-1 - pos2.x;
    if (pos3.y % 2 != 0) pos2.y = Ny-1 - pos2.y;
    if (pos3.z % 2 != 0) pos2.z = Nz-1 - pos2.z;




    if (data[pos2.z*Nx*Ny + pos2.y*Nx + pos2.x] == 0) {
        out_of_bounds = false;
    } else {
        vec3 tDelta = abs(invdir);                       // the amount to move along the ray to cross a voxel in each direction
        ivec3 step = 1 - 2*sign;
        vec3 t = ((0.5f + 0.5f*vec3(step) + vec3(pos)) - orig)*invdir;
        out = -1;
        
        while (tnew < 1024) {
            // while inside bounding box
            pos2 = pos;
            pos2.x = mod2(pos.x, Nx);
            pos2.y = mod2(pos.y, Ny);
            pos2.z = mod2(pos.z, Nz);

            pos3 = pos;
            pos3.x = floor(pos.x/float(Nx));
            pos3.y = floor(pos.y/float(Ny));
            pos3.z = floor(pos.z/float(Nz));

            if (pos3.x % 2 != 0) pos2.x = Nx-1 - pos2.x;
            if (pos3.y % 2 != 0) pos2.y = Ny-1 - pos2.y;
            if (pos3.z % 2 != 0) pos2.z = Nz-1 - pos2.z;

            if (data[pos2.z*Nx*Ny + pos2.y*Nx + pos2.x] == 0) { // HIT!
                out_of_bounds = false;
                break;
            }

            // NO HIT! Move to next voxel

            tnew = min3(t);

            if (tnew == t.x) {
                pos.x = pos.x + step.x;
                t.x = t.x + tDelta.x;
                face = 0 + 3*sign.x;
                if (pos3.x % 2 != 0) face = (face + 3) % 6;
            } else if (tnew == t.y) {
                pos.y = pos.y + step.y;
                t.y = t.y + tDelta.y;
                face = 1 + 3*sign.y;
                if (pos3.y % 2 != 0) face = (face + 3) % 6;
            } else if (tnew == t.z) {
                pos.z = pos.z + step.z;
                t.z = t.z + tDelta.z;
                face = 2 + 3*sign.z;
                if (pos3.z % 2 != 0) face = (face + 3) % 6;
            }

        }  
    }
    t = tnew;
    out = face;

    return pos;
}


ivec3 traceWrap(const vec3 &orig, const vec3 &dir, bool &out_of_bounds, int &out) {
    vec3 invdir = 1.0/dir;
    ivec3 sign = ivec3((dir.x < 0 ? 1 : 0), (dir.y < 0 ? 1 : 0), (dir.z < 0 ? 1 : 0));

    ivec3 pos;
    out = -1;
    
    out_of_bounds = true;  // assume ray goes through voxels without hitting


    float tmin = 0.0;
    pos = ivec3(floor(orig));

    int face = -1;              // assume illegal face (boundary)
    float tnew = tmin;


    ivec3 pos2 = pos;
    pos2.x = mod2(pos.x, Nx);
    pos2.y = mod2(pos.y, Ny);
    pos2.z = mod2(pos.z, Nz);

    ivec3 pos3 = pos;
    pos3.x = pos.x/Nx;
    pos3.y = pos.y/Ny;
    pos3.z = pos.z/Nz;


    if (data[pos2.z*Nx*Ny + pos2.y*Nx + pos2.x] == 0) {
        out_of_bounds = false;
    } else {
        vec3 tDelta = abs(invdir);                       // the amount to move along the ray to cross a voxel in each direction
        ivec3 step = 1 - 2*sign;
        vec3 t = ((0.5f + 0.5f*vec3(step) + vec3(pos)) - orig)*invdir;
        out = -1;
        
        while (tnew < 1024) {
            // while inside bounding box
            pos2 = pos;
            pos2.x = mod2(pos.x, Nx);
            pos2.y = mod2(pos.y, Ny);
            pos2.z = mod2(pos.z, Nz);

            pos3 = pos;
            pos3.x = floor(pos.x/float(Nx));
            pos3.y = floor(pos.y/float(Ny));
            pos3.z = floor(pos.z/float(Nz));

            if (data[pos2.z*Nx*Ny + pos2.y*Nx + pos2.x] == 0) { // HIT!
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
    t = tnew;
    out = face;


    return pos;
}

void Draw() {


    
    double aspect = resx/resy;
    double FOV = fov*PI/180.0;
    double half_height = tan(FOV/2.0);
    double half_width = tan(FOV/2.0)*aspect;




    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, resx, resy, 0, GL_RED, GL_FLOAT, NULL);

    glUseProgram(computeProgramID);
    glUniform3fv(glGetUniformLocation(computeProgramID, "orig"), 1,     &cam_pos.x);
    glUniform3fv(glGetUniformLocation(computeProgramID, "f"),    1,     &f.x);
    glUniform3fv(glGetUniformLocation(computeProgramID, "r"),    1,     &r.x);
    glUniform3fv(glGetUniformLocation(computeProgramID, "u"),    1,     &u.x);
    glUniform1f( glGetUniformLocation(computeProgramID, "half_width"),  half_width);
    glUniform1f( glGetUniformLocation(computeProgramID, "half_height"), half_height);
    glUniform1f( glGetUniformLocation(computeProgramID, "fogDistance"), fogDistance);
    
    glUniform3iv(glGetUniformLocation(computeProgramID, "chosenVoxel"), 1,  &pos.x);
    glUniform3i( glGetUniformLocation(computeProgramID, "N"),           Nx, Ny, Nz);
    glUniform1i( glGetUniformLocation(computeProgramID, "chosenFace"),  chosenFace);
    glUniform1i( glGetUniformLocation(computeProgramID, "dist_clip"),   dist_clipped);
    glUniform1i( glGetUniformLocation(computeProgramID, "textureMode"), textureMode == GL_REPEAT ? 1 : textureMode == GL_MIRRORED_REPEAT ? 2 : 0);

    glUniform3fv(glGetUniformLocation(computeProgramID, "colors"), 6, &colors[0].x);


    glDispatchCompute(ceil(resx/8), ceil(resy/8), 1); // blocks of 8^2


    glUseProgram(computerRenderProgramID);

    glBindVertexArray(VertexArrayID);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer_quad);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuffer_quad);
    glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,0,(void*)0);
    glScissor(0, 0, resx, resy);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);


}
void initGL() {
    printf("Initializing OpenGL/GLFW\n"); 
    if (!glfwInit()) {
        printf("Could not initialize\n");
        exit(-1);
    }
    glfwWindowHint(GLFW_SAMPLES, 4);    // samples, for antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4); // shader version should match these
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // do not use deprecated functionality

    window = glfwCreateWindow(resx, resy, "GLSL template", 0, 0);
    if (!window) {
        printf("Could not open glfw window\n");
        glfwTerminate();
        exit(-2);
    }
    glfwMakeContextCurrent(window); 

    if (!gl_lite_init()) {
        printf("error init gl_lite\n"); fflush(stdout);
    }


    // Create and bind  VBO (Vertex Buffer Object). Vertex and element buffers are attached to this. 
    // Can have multiple VBO's
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);


    computerRenderProgramID = LoadShaders( "shaders/vertex_shader_render.vs", "shaders/fragment_shader_render.fs" );

    glGenBuffers(1, &vertexbuffer_quad);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer_quad);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    glGenBuffers(1, &uvbuffer_quad);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuffer_quad);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_uv_buffer_data), g_uv_buffer_data, GL_STATIC_DRAW);


    glGenTextures(1, &texture);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, resx, resy, 0, GL_RED, GL_FLOAT, NULL);
    glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    computeProgramID = LoadComputeShader("shaders/compute_shader.cs");



    glGenTextures(1, &texture3D);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, texture3D);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, Nx, Ny, Nz, 0, GL_RED, GL_UNSIGNED_BYTE, data);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, textureMode);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, textureMode);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, textureMode);

    // VSync on
    // set background color and enable depth testing
    glfwSwapInterval(0);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

}


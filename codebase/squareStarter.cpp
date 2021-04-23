//HW 0 - Moving Square
//Starter code for the first homework assignment.
//This code assumes SDL2 and OpenGL are both properly installed on your system

#include "glad/glad.h"  //Include order can matter here
#if defined(__APPLE__) || defined(__linux__)
 #include <SDL2/SDL.h>
 #include <SDL2/SDL_opengl.h>
#else
 #include <SDL.h>
 #include <SDL_opengl.h>
#endif
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include "image_lib.h"
#include"raytracer.h"

using namespace std;

float vertices[] = {  //The function updateVertices() changes these values to match p1,p2,p3,p4
//  X     Y     R     G     B     U    V
  1.0f,  1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,  // top right
  1.0f, -1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, // bottom right
  -1.0f,  1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f,  // top left 
  -1.0f, -1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f,  // bottom left
};

//////////////////////////
///  Begin your code here
/////////////////////////

bool fullscreen = false;
bool done = false;

void keyReleased(int scancode) {
    if (scancode == SDLK_f) //If "f" is pressed
        fullscreen = !fullscreen;
    if (scancode == SDLK_ESCAPE)
        done = true; //Exit event loop
}

// Shader sources
const GLchar* vertexSource =
   "#version 430 core\n"
   "in vec2 position;"
   "in vec3 inColor;"
   "in vec2 inTexcoord;"
   "out vec3 Color;"
   "out vec2 texcoord;"
   "void main() {"
   "   Color = inColor;"
   "   gl_Position = vec4(position, 0.0, 1.0);"
   "   texcoord = inTexcoord;"
   "}";
    
const GLchar* fragmentSource =
   "#version 430 core\n"
   "uniform sampler2D tex0;"
   "in vec2 texcoord;"
   "out vec3 outColor;"
   "void main() {"
   "   outColor = texture(tex0, texcoord).xyz;"
   "}";
    
struct MaterialGL {
    float ka[3];
    float padding1;
    float kd[3];
    float padding2;
    float ks[3];
    float padding3;
    float kt[3];
    float padding4;
    // 3 vec4
    float ns;
    float ior;
    float padding5;
    float padding6;
    // 2 vec
    MaterialGL() {
        ka[0] = 1.0f;
        ka[1] = 1.0f;
        ka[2] = 1.0f;
        kd[0] = 1.0f;
        kd[1] = 1.0f;
        kd[2] = 1.0f;
        ks[0] = 0.0f;
        ks[1] = 0.0f;
        ks[2] = 0.0f;
        kt[0] = 0.0f;
        kt[1] = 0.0f;
        kt[2] = 0.0f;
        ns = 5.0f;
        ior = 1.0f;
    };
};

struct TriangleGL {
    float p1[3]; // 1 vec3
    float padding1;
    float p2[3]; // 1 vec3
    float padding2;
    float p3[3]; // 1 vec3
    float padding3;
    float n1[3]; // 1 vec3
    float padding4;
    float n2[3]; // 1 vec3
    float padding5;
    float n3[3]; // 1 vec3
    float padding6;
    MaterialGL mat;
};

struct LightGL {
    float pos[4];
    float dir[4];
    float clr[4];
    // assume point light,
    // try adding an int later
};

int max_tri = 0;
MaterialGL* mats = nullptr;
TriangleGL* tris = nullptr;
LightGL* lights = nullptr;
int num_lights = 0;

size_t width = 1080;
size_t height = 720;
float half_fov = 45.0;
float eye[3] = { 0.0 ,0.0, 0.0 };
float fwd[3] = { 0.0, 0.0, -1.0 };
float cam_r[3] = { 1.0, 0.0, 0.0 };
float up[3] = { 0.0, 1.0, 0.0 };
float b_clr[3] = { 0.0, 0.0, 0.0 };
void loadFromFile(string input_file_name) {
    FILE* in_file;

    in_file = fopen(input_file_name.c_str(), "r");
    if (in_file == NULL) {
        std::cerr << "Couldn't open file: " << input_file_name << std::endl;
        exit(1);
    }
    // The file was opened, so we are good to parse data now!
    MaterialGL cur_mat;
    // we allow at most 100 materials in the shader
    int num_mats = 1.0;
    mats = new MaterialGL[100];
    lights = new LightGL[100];
    mats[0] = cur_mat;

    float* verts = nullptr;
    float* norms = nullptr;
    char arg[1024];
    Point3D min_pt(INFINITY, INFINITY, INFINITY), max_pt(-INFINITY, -INFINITY, -INFINITY);
    int max_vert(-1), max_norm(-1), nindex(0), vindex(0), tindex(0);

    while (fgets(arg, 1024, in_file)) {
        // Check if it is a comment, if so skip!
        if (arg[0] == '#') {
            continue;
        }
        // This contains valid data, so extract the command
        char command[100];
        int fieldsread = sscanf(arg, "%s ", command);
        std::string commandstr = command;
        if (fieldsread < 1) {
            continue;
        }
        if (commandstr == "camera_fwd:") {
            sscanf(arg, "camera_fwd: %f %f %f", &fwd[0], &fwd[1], &fwd[2]);
            printf("found forward\n");
        }
        else if (commandstr == "camera_up:") {
            sscanf(arg, "camera_up: %f %f %f", &up[0], &up[1], &up[2]);
        }
        else if (commandstr == "camera_pos:") {
            sscanf(arg, "camera_pos: %f %f %f", &eye[0], &eye[1], &eye[2]);
        }
        else if (commandstr == "material:") {
            float a_r, a_g, a_b, d_r, d_g, d_b, s_r, s_g, s_b, t_r, t_g, t_b, ns, ior;
            sscanf(arg, "material: %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &a_r, &a_g, &a_b,
                &d_r, &d_g, &d_b, &s_r, &s_g, &s_b, &ns, &t_r, &t_g, &t_b, &ior);
            cur_mat = MaterialGL();
            cur_mat.ka[0] = a_r;
            cur_mat.ka[1] = a_g;
            cur_mat.ka[2] = a_b;

            cur_mat.kd[0] = d_r;
            cur_mat.kd[1] = d_g;
            cur_mat.kd[2] = d_b;

            cur_mat.ks[0] = s_r;
            cur_mat.ks[1] = s_g;
            cur_mat.ks[2] = s_b;

            cur_mat.kt[0] = t_r;
            cur_mat.kt[1] = t_g;
            cur_mat.kt[2] = t_b;

            cur_mat.ns = ns;
            cur_mat.ior = ior;
            mats[num_mats] = cur_mat;
            num_mats++;
        }
        else if (commandstr == "max_vertices:") {
            sscanf(arg, "max_vertices: %d", &max_vert);
            verts = new float[3 * max_vert];
        }
        else if (commandstr == "max_normals:") {
            sscanf(arg, "max_normals: %d", &max_norm);
            norms = new float[3 * max_norm];
        }
        else if (commandstr == "max_triangles:") {
            sscanf(arg, "max_triangles: %d", &max_tri);
            tris = new TriangleGL[max_tri];
        }
        else if (commandstr == "vertex:") {
            sscanf(arg, "vertex: %f %f %f", &verts[vindex], &verts[vindex + 1], &verts[vindex + 2]);
            vindex += 3;
        }
        else if (commandstr == "normal:") {
            sscanf(arg, "normal: %f %f %f", &norms[nindex], &norms[nindex + 1], &norms[nindex + 2]);
            nindex += 3;
        }
        else if (commandstr == "triangle:") {
            if (max_vert < 0) {
                printf("ERROR: NUMBER OF VERTICES NOT SPECIFIED. SKIPPING\n");
                continue;
            }
            int p1, p2, p3;
            TriangleGL t;
            sscanf(arg, "triangle: %d %d %d", &p1, &p2, &p3);
            memcpy(t.p1, verts + p1*3, 3 * sizeof(float));
            memcpy(t.p2, verts + p2*3, 3 * sizeof(float));
            memcpy(t.p3, verts + p3*3, 3 * sizeof(float));
            // now get the face normal
            float norm[3];
            Dir3D face = cross(Point3D(t.p2[0], t.p2[1], t.p2[2]) - Point3D(t.p1[0], t.p1[1], t.p1[2]), Point3D(t.p3[0], t.p3[1], t.p3[2]) - Point3D(t.p1[0], t.p1[1], t.p1[2]));
            face = face.normalized();
            norm[0] = face.x;
            norm[1] = face.y;
            norm[2] = face.z;
            memcpy(t.n1, norm, 3 * sizeof(float));
            memcpy(t.n2, norm, 3 * sizeof(float));
            memcpy(t.n3, norm, 3 * sizeof(float));
            t.mat = cur_mat;
            tris[tindex++] = t;
        }
        else if (commandstr == "normal_triangle:") {
            if (max_vert < 0 || max_norm < 0) {
                printf("ERROR: NUMBER OF VERTICES/NORMALS NOT SPECIFIED. SKIPPING\n");
                continue;
            }
            int p1, p2, p3, n1, n2, n3;
            TriangleGL t;
            sscanf(arg, "normal_triangle: %d %d %d %d %d %d", &p1, &p2, &p3, &n1, &n2, &n3);
            memcpy(t.p1, verts + p1*3, 3 * sizeof(float));
            memcpy(t.p2, verts + p2*3, 3 * sizeof(float));
            memcpy(t.p3, verts + p3*3, 3 * sizeof(float));

            memcpy(t.n1, norms + n1*3, 3 * sizeof(float));
            memcpy(t.n2, norms + n2*3, 3 * sizeof(float));
            memcpy(t.n3, norms + n3*3, 3 * sizeof(float));
            t.mat = cur_mat;
            tris[tindex++] = t;
        }
        else if (commandstr == "background:") {
            int val = sscanf(arg, "background: %f %f %f", &b_clr[0], &b_clr[1], &b_clr[2]);
        }
        else if (commandstr == "point_light:") {
            LightGL light;
            sscanf(arg, "point_light: %f %f %f %f %f %f", &light.clr[0], &light.clr[1], &light.clr[2], &light.pos[0], &light.pos[1], &light.pos[2]);
            lights[num_lights++] = light;
        }
        else if (commandstr == "camera_fov_ha:") {
            sscanf(arg, "camera_fov_ha: %f", &half_fov);
        }
    }
    fclose(in_file);
    printf("Loaded %d triangles\n", max_tri);
    // orthogonalize the camera basis
    Dir3D f(fwd[0], fwd[1], fwd[2]);
    Dir3D u(up[0], up[1], up[2]);
    Dir3D r = cross(u, f).normalized();
    u = cross(f, r).normalized();
    f = f.normalized();

    fwd[0] = f.x;
    fwd[1] = f.y;
    fwd[2] = f.z;

    cam_r[0] = r.x;
    cam_r[1] = r.y;
    cam_r[2] = r.z;

    up[0] = u.x;
    up[1] = u.y;
    up[2] = u.z;
    printf("Background Color: %f %f %f\n", b_clr[0], b_clr[1], b_clr[2]);
    printf("Camera eye: %f %f %f\n", eye[0], eye[1], eye[2]);
    printf("Camera forward: %f %f %f\n", fwd[0], fwd[1], fwd[2]);
    printf("Camera right: %f %f %f\n", cam_r[0], cam_r[1], cam_r[2]);
    printf("Camera up: %f %f %f\n", up[0], up[1], up[2]);
}

int main(int argc, char *argv[]){
   //Image image(width, height);
   string filename;
   cin >> filename;
   float *img_data = new float[(size_t)4*width*height];
   memset(img_data, 1, (size_t)4 * width * height);
   loadFromFile(filename);
   SDL_Init(SDL_INIT_VIDEO);  //Initialize Graphics (for OpenGL)

   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);

   
   SDL_Window* window = SDL_CreateWindow("My OpenGL Program", 100, 100, width, height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
   //Create a context to draw in
   SDL_GLContext context = SDL_GL_CreateContext(window);

   //OpenGL functions using glad library
   if (gladLoadGLLoader(SDL_GL_GetProcAddress)) {
       printf("OpenGL loaded\n");
       printf("Vendor:   %s\n", glGetString(GL_VENDOR));
       printf("Renderer: %s\n", glGetString(GL_RENDERER));
       printf("Version:  %s\n", glGetString(GL_VERSION));
   }
   else {
       printf("ERROR: Failed to initialize OpenGL context.\n");
       return -1;
   }
   //glViewport(0, 0, width, height);
   //// Allocate Texture 0 (Created in Load Image) ///////
   GLuint texture;
   glGenTextures(1, &texture);
   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, texture);

   //What to do outside 0-1 range
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); //GL_LINEAR
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); //GL_LINEAR
   //TODO: TEST your understanding: Try GL_LINEAR instead of GL_NEAREST on the 4x4 test image. What is happening?

   //Load the texture into memory
   glBindTexture(GL_TEXTURE_2D, texture);
   glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, img_data);
   glGenerateMipmap(GL_TEXTURE_2D);
   //// End Allocate Texture ///////


   //Build a Vertex Array Object. This stores the VBO and attribute mappings in one object
   GLuint vao;
   glGenVertexArrays(1, &vao); //Create a VAO
   glBindVertexArray(vao); //Bind the above created VAO to the current context


   //Allocate memory on the graphics card to store geometry (vertex buffer object)
   GLuint vbo;
   glGenBuffers(1, &vbo);  //Create 1 buffer called vbo
   glBindBuffer(GL_ARRAY_BUFFER, vbo); //Set the vbo as the active array buffer (Only one buffer can be active at a time)
   glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW); //upload vertices to vbo
   //GL_STATIC_DRAW means we won't change the geometry, GL_DYNAMIC_DRAW = geometry changes infrequently
   //GL_STREAM_DRAW = geom. changes frequently.  This effects which types of GPU memory is used
   GLuint rayTracer, computeShader;
   // Load the compute shader
   ifstream computeFile("..\\..\\rayTrace_compute.glsl");
   string computeSource((istreambuf_iterator<char>(computeFile)), istreambuf_iterator<char>());;
   computeShader = glCreateShader(GL_COMPUTE_SHADER);
   const char* src = computeSource.c_str();
   glShaderSource(computeShader, 1, &src, NULL);
   glCompileShader(computeShader);
   GLint status;
   glGetShaderiv(computeShader, GL_COMPILE_STATUS, &status);
   if (!status) {
       char buffer[512];
       glGetShaderInfoLog(computeShader, 512, NULL, buffer);
       printf("Compute Shader Compile Failed. Info:\n\n%s\n", buffer);
   }

   rayTracer = glCreateProgram();
   glAttachShader(rayTracer, computeShader);
   glLinkProgram(rayTracer);
   glUseProgram(rayTracer);

   // set all the compute shader uniforms
   // 0 since we are using texture 0
   glUniform1i(glGetUniformLocation(rayTracer, "result"), 0);
   // we won't change the camera or anything, so let's set those too
   GLenum err = glGetError();
   if (err != GL_NO_ERROR) {
       cout << err << endl;
   }

   float d = (height * .5) / tanf(half_fov * (M_PI / 180.0));
   glUniform1i(glGetUniformLocation(rayTracer, "num_lights"), num_lights);
   glUniform1i(glGetUniformLocation(rayTracer, "num_tris"), max_tri);
   glUniform1f(glGetUniformLocation(rayTracer, "width"), width);
   glUniform1f(glGetUniformLocation(rayTracer, "half_width"), width * .5);
   glUniform1f(glGetUniformLocation(rayTracer, "height"), height);
   glUniform1f(glGetUniformLocation(rayTracer, "half_height"), height * .5);
   glUniform1f(glGetUniformLocation(rayTracer, "d"), d);
   glUniform3fv(glGetUniformLocation(rayTracer, "background_clr"), 1, b_clr);
   glUniform3fv(glGetUniformLocation(rayTracer, "eye"), 1, eye);
   glUniform3fv(glGetUniformLocation(rayTracer, "forward"), 1, fwd);
   glUniform3fv(glGetUniformLocation(rayTracer, "right"), 1, cam_r);
   glUniform3fv(glGetUniformLocation(rayTracer, "up"), 1, up);
   // create the triangle buffer
   GLuint tri_ssbo, light_ssbo;
   glGenBuffers(1, &tri_ssbo);
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, tri_ssbo);
   glBufferData(GL_SHADER_STORAGE_BUFFER, max_tri * sizeof(TriangleGL), tris, GL_STREAM_READ);
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, tri_ssbo);
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
   err = glGetError();
   if (err != GL_NO_ERROR) {
       int val;
       cout << err << endl;
       cin >> val;
   }
   glGenBuffers(1, &light_ssbo);
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, light_ssbo);
   glBufferData(GL_SHADER_STORAGE_BUFFER, num_lights * sizeof(LightGL), lights, GL_STREAM_READ);
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, light_ssbo);
   glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // reset the bound buffer
   err = glGetError();
   if (err != GL_NO_ERROR) {
       cout << err << endl;
   }
   //Load the vertex Shader
   GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
   glShaderSource(vertexShader, 1, &vertexSource, NULL);
   glCompileShader(vertexShader);

   //Let's double check the shader compiled 
   glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
   if (!status) {
       char buffer[512];
       glGetShaderInfoLog(vertexShader, 512, NULL, buffer);
       printf("Vertex Shader Compile Failed. Info:\n\n%s\n", buffer);
   }

   GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
   glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
   glCompileShader(fragmentShader);

   //Double check the shader compiled
   glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
   if (!status) {
       char buffer[512];
       glGetShaderInfoLog(fragmentShader, 512, NULL, buffer);
       printf("Fragment Shader Compile Failed. Info:\n\n%s\n", buffer);
   }

   //Join the vertex and fragment shaders together into one program
   GLuint shaderProgram = glCreateProgram();
   glAttachShader(shaderProgram, vertexShader);
   glAttachShader(shaderProgram, fragmentShader);
   glBindFragDataLocation(shaderProgram, 0, "outColor"); // set output
   glLinkProgram(shaderProgram); //run the linker

   glUseProgram(shaderProgram); //Set the active shader (only one can be used at a time)


   //Tell OpenGL how to set fragment shader input 
   glBindBuffer(GL_ARRAY_BUFFER, vbo);
   GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
   //               Attribute, vals/attrib., type, normalized?, stride, offset
   glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(float), 0);
   glEnableVertexAttribArray(posAttrib); //Binds the VBO current GL_ARRAY_BUFFER 

   GLint colAttrib = glGetAttribLocation(shaderProgram, "inColor");
   glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(2 * sizeof(float)));
   glEnableVertexAttribArray(colAttrib);

   GLint texAttrib = glGetAttribLocation(shaderProgram, "inTexcoord");
   glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(5 * sizeof(float)));
   glEnableVertexAttribArray(texAttrib);
   //Event Loop (Loop forever processing each event as fast as possible)
   glUseProgram(rayTracer);
   auto s = chrono::high_resolution_clock::now();
   glDispatchCompute(width/10, height/10, 1);
   glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
   auto e = chrono::high_resolution_clock::now();
   auto dur = e - s;
   auto ms = std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count();
   cout << "time elapsed: " << ms << endl;
   glUseProgram(shaderProgram);
   SDL_Event windowEvent;
   int t = 0;
   while (!done) {
       while (SDL_PollEvent(&windowEvent)) {  //Process input events (e.g., mouse & keyboard)
           if (windowEvent.type == SDL_QUIT) 
               done = true;
           if (windowEvent.type == SDL_KEYUP)
               keyReleased(windowEvent.key.keysym.sym);
       }
       t = (t++) % 200;
       float amt = 4 + 8 * .01 * (-.01 * (t - 100) * (t - 100) + 100);
       lights[0].clr[0] = 
       lights[0].clr[1] = 
       lights[0].clr[2] = amt;

       glBindBuffer(GL_SHADER_STORAGE_BUFFER, light_ssbo);
       glBufferData(GL_SHADER_STORAGE_BUFFER, num_lights * sizeof(LightGL), lights, GL_STREAM_READ);
       glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, light_ssbo);
       glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // reset the bound buffer

       glUseProgram(rayTracer);
       glDispatchCompute(width / 10, height / 10, 1);
       glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
       glUseProgram(shaderProgram);
       
       glDrawArrays(GL_TRIANGLE_STRIP, 0, 4); //Draw the two triangles (4 vertices) making up the square
       SDL_GL_SwapWindow(window); //Double buffering
   }
   delete[] img_data;
   if (mats != nullptr)
       delete[] mats;
   if (tris != nullptr)
       delete[] tris;
   if (mats != nullptr)
       delete[] mats;
   if (lights != nullptr)
       delete[] lights;

   glDeleteProgram(rayTracer);
   glDeleteShader(computeShader);
   glDeleteProgram(shaderProgram);
   glDeleteShader(fragmentShader);
   glDeleteShader(vertexShader);

   glDeleteBuffers(1, &vbo);
   glDeleteVertexArrays(1, &vao);


   //Clean Up
   SDL_GL_DeleteContext(context);
   SDL_Quit();
   return 0; 
}
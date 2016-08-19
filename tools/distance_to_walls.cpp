#include <stdio.h> // printf, sprintf, fflush, stdout
#include <math.h>  // atan, pow, ceil (for std_easy_font.h)

double PI = 4.0*atan(1.0);


// OpenGL/GLFW stuff
#include <GLFW/glfw3.h>

GLFWwindow* window;
double resx = 1440/2, resy = 1150/2;
double prevx = -1, prevy = -1;     // to track how much the mouse moved between frames
double cam_x = -0.1, cam_y = -0.1; // world coordinates of lower-left corner of the window
double cam_height = 1.2;
double cam_width = cam_height*resx/double(resy); 
int clickedButtons = 0;

enum buttonMaps { FIRST_BUTTON=1, SECOND_BUTTON=2, THIRD_BUTTON=4, FOURTH_BUTTON=8, FIFTH_BUTTON=16, NO_BUTTON=0 };
enum modifierMaps { CTRL=2, SHIFT=1, ALT=4, META=8, NO_MODIFIER=0 };

void initGL(void);
void windowsize_callback(GLFWwindow*, int, int);
void key_callback(GLFWwindow*, int, int, int, int);
void mousebutton_callback(GLFWwindow*, int, int, int);
void mousepos_callback(GLFWwindow*, double, double);
void mousewheel_callback(GLFWwindow*, double, double);
void Draw(void);


#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define MAX(a,b)   ((a) > (b) ? (a) : (b))


int mod2(double i, int n) {
    return i - n*floor(i/n);
}
int mod2(double i, double n) {
    return i - n*floor(i/n);
}

// Text stuff
#include <stb_easy_font.h>

enum TextAlignment {ALIGN_MID, ALIGN_LEFT, ALIGN_RIGHT};

void print_string(float, float, char*, float, float, float);
void print_string_world(float, float, char, float, float, float, float, int);
void print_string_screen(float, float, char*, float, float, float, float);


const int N = 1000;
unsigned char *grid, *grid1, *grid2;


int testValue = 0;

struct Dist {
	short dx, dy;
	float dist;

	bool operator<(const Dist &rhs) {
		return this->dist < rhs.dist;
	}
	Dist& operator=(const Dist &rhs) {
		this->dx = rhs.dx;
		this->dy = rhs.dy;
		this->dist = rhs.dist;

		return *this;
	}
};

Dist *distances;


#define SWAP(TYPE,a,b)  \
    do { TYPE stb__t; stb__t = (a); (a) = (b); (b) = stb__t; } while (0)


void mergeArrays(Dist *arr1, Dist *arr2, int N1, int N2, Dist *out) {
	int i = 0, j = 0, k = 0;
	while (k < N1 + N2) {
		//if ( ((i < N1 && j < N2)  && arr2[j] < arr1[i]) || j == N2 )
		if ( ((i < N1 && j < N2)  &&  arr1[i] < arr2[j]) || j == N2 )
			out[k++] = arr1[i++];
		else
			out[k++] = arr2[j++];
	}
}

void mergeSort(Dist **arr, int N0) {
	Dist *arr1 = *arr;
	Dist *arr2 = (Dist*)malloc(sizeof(Dist)*N0);

	int M = log2(N0);
	int m = 1;
	int a = N0/2;
	for (int i = 0; i <= M; i++) {
		int b = a << (i+1);

		for (int j = 0; j < a; j++) {
			int k = 2*j*m;
			mergeArrays(&arr1[k], &arr1[k+m], m, m, &arr2[k]);
		}
		if (N0 - b > m) {
			int k = 2*a*m;
			mergeArrays(&arr1[k], &arr1[k+m], m, (N0-b)&(m-1), &arr2[k]); 
		} else {
			for (int k = 0; k < N0-b; k++) 
				arr2[N0-k-1] = arr1[N0-k-1];
		}
		SWAP(Dist*, arr1,arr2);
		m <<= 1;
		a >>= 1;
	}

	*arr = arr1;

	free(arr2);
}

int isSorted(Dist *arr, int n) {
	for (int i = 0; i < n-1; i++) {
		if (arr[i+1] < arr[i])
			return 0;
	}

	return 1;
}

void reset() {
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			grid[j*N+i] = 255;
		}
	}
	for (int i = 0; i < N; i++) {
		int j = int((0.5 + 0.5*sin(2*PI*i/100.0))*N);
		j = MAX(0, MIN(j, N-1));
		grid[j*N+i] = 0;
	}
	for (int i = 0; i < N; i++) {
		int j = int((0.5 + 0.5*sin(2*PI*i/100.0))*N);
		j = MAX(0, MIN(j, N-1));
		grid[i*N+j] = 0;
	}
	testValue = 0;
}


int main() {
	cam_height = MIN(500,N);
	cam_width = cam_height*resx/double(resy);

	distances = new Dist[N*N];

	for (int k = 0; k < N*N; k++) {
		int dx = k % N;
		int dy = k / N;
		if (dx > N/2.0) dx -= N;
		if (dx <= -N/2.0) dx += N;
		if (dy > N/2.0) dy -= N;
		if (dy <= -N/2.0) dy += N;

		distances[k].dx = dx;
		distances[k].dy = dy;
		distances[k].dist = sqrt(dx*dx+dy*dy);
	}



	mergeSort(&distances, N*N);

	grid1 = new unsigned char[N*N];
	grid2 = new unsigned char[N*N];
	grid = grid1;


	
	reset();
	



	printf("sorted = %s\n", isSorted(distances, N*N) ? "true" : "false"); fflush(stdout);


	initGL();



	// Main loop
	while ( !glfwWindowShouldClose(window)) {	
		Draw();
		glfwPollEvents();
	}

	return 0;
}

void distanceEstimate() {
	double t1 = glfwGetTime(); 
	for (int k = 0; k < N*N; k++) {
		int i = k % N;
		int j = k / N;

		if (grid[k] == 0) 
			continue;

		double min_dist = 100000000;
		for (int K = 0; K < N*N; K++) {
			int x = i+distances[K].dx;
			if (x < 0) x += N;
			if (x >= N) x -= N;
			int y = j+distances[K].dy;
			if (y < 0) y += N;
			if (y >= N) y -= N;
			int z = y*N+x;

			if (z == k)
				continue;

			if (grid[z] == 0) {
				min_dist = distances[K].dist;
				break;
			}
		}
		grid[k] = MIN(255, int(min_dist));

		if (k % (100*N) == 0) {
			double t2 = glfwGetTime();
			printf("%d/%d: %f\n",k, N*N-1, t2-t1); fflush(stdout);
		}
			
	}
	double t2 = glfwGetTime();
	printf("%f\n",t2-t1); fflush(stdout);
}


void dilate2() {
	printf("dilate\n"); fflush(stdout);

	unsigned char *othergrid;
	if (grid == grid1)
		othergrid = grid2;
	else
		othergrid = grid1;

	for (int k = 0; k < N*N; k++) {
		int x = k % N;
		int y = k / N;

		int count = 0;

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				//if (!(i==0 && j == 0) && (abs(i)+abs(j) == 1)) {
				if (!(i==0 && j == 0)) {
					int X = (x + i + N) % N;
					int Y = (y + j + N) % N;
					count += (grid[Y*N+X] == testValue);
				}
			}
		}

		if (grid[y*N+x] == 255 && count > 0) {
			othergrid[y*N+x] = 0;
		} else {
			othergrid[y*N+x] = grid[y*N+x];
		}
	}


	grid = othergrid;

}

void dilate() {
	printf("dilate\n"); fflush(stdout);
	if (testValue == 255)
		return;

	unsigned char *othergrid;
	if (grid == grid1)
		othergrid = grid2;
	else
		othergrid = grid1;

	for (int k = 0; k < N*N; k++) {
		int x = k % N;
		int y = k / N;

		int count = 0;

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				//if (!(i==0 && j == 0) && (abs(i)+abs(j) == 1)) {
				if (!(i==0 && j == 0)) {
					int X = (x + i + N) % N;
					int Y = (y + j + N) % N;
					count += (grid[Y*N+X] == testValue);
				}
			}
		}

		if (grid[y*N+x] == 255 && count > 0) {
			othergrid[y*N+x] = testValue+1;
		} else {
			othergrid[y*N+x] = grid[y*N+x];
		}
	}

	testValue++;

	grid = othergrid;

}

// Printing text using an 8 pixel bitmap font (amiga clone?)
// Assume ortogonal projection and screen coordinates
// Draws at depth 0
// from https://github.com/nothings/stb/
void print_string(float x, float y, char *text, float r, float g, float b) {
	static char buffer[99999]; // ~500 chars, enough to fill screen
	int num_quads;
	num_quads = stb_easy_font_print(x, y, text, NULL, buffer, sizeof(buffer));
	glColor3f(r,g,b);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_FLOAT, 16, buffer);
	glDrawArrays(GL_QUADS, 0, num_quads*4);
	glDisableClientState(GL_VERTEX_ARRAY);
}

// Draw text in screen coordinates
// Saves matrix states
void print_string_screen(float x, float y, char *text, float r, float g, float b, float scale = 1.0) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, resx, resy, 0.0, -1.0, 1.0); // top to bottom for text

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(0.0, 0.0, 0.1);    // Move text layer up 0.1 units to be above 0-level
	glScalef(scale, scale, scale);

	print_string(x, y, text, r, g, b);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

// Draw text in world coordinates
// Saves matrix states
// Uses an arbitrary scale
// ALIGN_MID (default), ALIGN_LEFT or ALIGN_RIGHT to justify/align the text relative to position
void print_string_world(float x, float y, char *text, float r, float g, float b, float scale, int alignment=ALIGN_MID) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, resx, resy, 0.0, -1000, 1000);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	int text_width = stb_easy_font_width(text);
	scale *= 0.01/(cam_width/resy); // convert scale from world to screen coordinates using camera information
	
	if (alignment == ALIGN_MID) {
		glTranslatef((x - cam_x)/cam_width * resx - 0.5*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	} else if (alignment == ALIGN_RIGHT) {
		glTranslatef((x - cam_x)/cam_width * resx - 1.0*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	} else if (alignment == ALIGN_LEFT) {
		glTranslatef((x - cam_x)/cam_width * resx - 0.0*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	}
	glScalef(scale, scale, scale);
	print_string(0.0, 0.0, text, r, g, b);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}



// Main drawing function
void Draw() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw non-text stuff
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glLineWidth(1.0);
	glColor3f(0,0,0);
	glBegin(GL_LINE_LOOP);
	glVertex3f(0, 0, 0);
	glVertex3f(N, 0, 0);
	glVertex3f(N, N, 0);
	glVertex3f(0, N, 0);
	glEnd();

	double xpos, ypos;

	glfwGetCursorPos(window, &xpos, &ypos);
	double xend = cam_x + xpos*(cam_width)/resx;
	double yend = cam_y + (resy - ypos)*cam_height/resy;	
	glLineWidth(4.0);
	glColor3f(0,0,0);
	glBegin(GL_LINE_LOOP);
	float r = grid[mod2(int(yend),N)*N+mod2(int(xend),N)];
	r -= 1.5;
	if (r <= 0)
		r = 0;
	for (int i = 0; i < 32; i++) {
		glVertex3f(xend + r*cos(2*PI*i/32), yend + r*sin(2*PI*i/32), 0.5);	
	}
	glEnd();

	for (int j = cam_y; j < cam_y+cam_height; j++) {
		for (int i = cam_x; i < cam_x+cam_width; i++) {
			int i2 = mod2(i, N);
			int j2 = mod2(j, N);
			glBegin(GL_QUADS);
			int c = grid[j2*N+i2];
			float R = 0.5 + 0.5*cos(7*c/10.0);
			float G = 0.5 + 0.5*cos(13*c/10.0);
			float B = 0.5 + 0.5*cos(23*c/10.0);
			glColor3f(R, G, B);
			glVertex3f(i,   j, -0.1);
			glVertex3f(i+1, j, -0.1);
			glVertex3f(i+1, j+1, -0.1);
			glVertex3f(i,   j+1, -0.1);
			glEnd();
			if (resx/cam_width > 64) {	
				if (i > cam_x-1 && i < cam_x + cam_width && j > cam_y-1 && j < cam_y + cam_height) {
					char str[32];
					sprintf(str, "%d", grid[j2*N+i2]);
					print_string_world(i+0.5, j+0.5, str, 1, 0, 0, 2.0);
				}
			}
		}	
	}


	// Draw Text
	static char str_text[256];
	static int frame = 0;
	int pos = 0;
	sprintf(str_text, "This is frame #%d", frame++);
	print_string_screen(3, 3 + 10*pos++, str_text, 0.0, 0.0, 0.0, 1.5);
	sprintf(str_text, "This is on the next line. Camera position = (%.4f, %.4f), (%.4f, %.4f)", cam_x, cam_y, cam_x + cam_width, cam_y + cam_height);
	print_string_screen(3, 3 + 10*pos++, str_text, 0.0, 0.0, 0.0, 1.5);

	// done, possibly wait for vsync
	glfwSwapBuffers(window);
}

// Initialize GLFW and OpenGL:
// Craetes a window, make it the current content
// Setup orthogonal coordinates
// Set background to white
// Tries to enable VSync
// Set up callback functions
void initGL() {
	printf("Initializing OpenGL/GLFW\n"); 
	if (!glfwInit()) {
		printf("Could not initialize\n");
		exit(-1);
	}

	glfwWindowHint(GLFW_SAMPLES , 4);
	window = glfwCreateWindow(resx, resy, "window title", 0, 0);
	if (!window) {
		printf("Could not open glfw window\n");
		glfwTerminate();
		exit(-2);
	}
	glfwMakeContextCurrent(window);

	glLoadIdentity();
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1.0, 1.0);	
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glfwSwapInterval(1);
	glEnable(GL_DEPTH_TEST);

	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mousebutton_callback);
	glfwSetScrollCallback(window, mousewheel_callback);
	glfwSetCursorPosCallback(window, mousepos_callback);
	glfwSetWindowSizeCallback(window, windowsize_callback);
}

// Callback function called every time the window size change
// Adjusts the camera width and heigh so that the scale stays the same
// Resets projection matrix
void windowsize_callback(GLFWwindow *win, int width, int height) {
	double distance_per_pixel = cam_height/resy; // assuming cam_height/resy == cam_width/resx

	resx = width;
	resy = height;
	cam_width = distance_per_pixel*resx;
	cam_height = distance_per_pixel*resy;

	glLoadIdentity();
	glViewport(0, 0, resx, resy);
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1, 1);
}

// Callback function called every time a keyboard key is pressed, released or held down
void key_callback(GLFWwindow* win, int key, int scancode, int action, int mods) {
	 printf("key = %d, scancode = %d, action = %d, mods = %d\n", key, scancode, action, mods); fflush(stdout);

	// Close window if escape is released
	if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) {
		glfwSetWindowShouldClose(win, GL_TRUE);
	}

	if (key == GLFW_KEY_SPACE && action) {
		dilate();
	}

	if (key == GLFW_KEY_D && action) {
		dilate2();
	}

	if (key == GLFW_KEY_E && action) {
		distanceEstimate();
	}

	if (key == GLFW_KEY_R && action) {
		reset();
	}
}

// Callback function called every time a mouse button pressed or released
void mousebutton_callback(GLFWwindow* win, int button, int action, int mods) {
	// get current cursor position, convert to world coordinates
	glfwGetCursorPos(win, &prevx, &prevy);
	double xend = cam_x + prevx*(cam_width)/resx;
	double yend = cam_y + (resy - prevy)*cam_height/resy;	

	printf("button = %d, action = %d, mods = %d at (%f %f)\n", button, action, mods, xend, yend); fflush(stdout);

	// To track the state of buttons outside this function
	if (action == 1)
		clickedButtons |= (1 << button);
	else
		clickedButtons &= ~(1 << button);


	// Test each button
	if (clickedButtons&FIRST_BUTTON) {
		
	} else if (clickedButtons&SECOND_BUTTON) {

	} else if (clickedButtons&THIRD_BUTTON) {

	} else if (clickedButtons&FOURTH_BUTTON) {

	} else if (clickedButtons&FIFTH_BUTTON) {

	}
}

// Callback function called every time a the mouse is moved
void mousepos_callback(GLFWwindow* win, double xpos, double ypos) {
	// move the camera if the first mouse button is held down
	// the cursor will stay on at the same location relative to world coordinates after movement
	if (clickedButtons&FIRST_BUTTON) {
		cam_x -= (xpos-prevx)/resx*cam_width;
		cam_y += (ypos-prevy)/resy*cam_height;

		prevx = xpos;
		prevy = ypos;
	} else if (clickedButtons&SECOND_BUTTON) {

	} else if (clickedButtons&THIRD_BUTTON) {

	} else if (clickedButtons&FOURTH_BUTTON) {

	} else if (clickedButtons&FIFTH_BUTTON) {

	}
}

// Callback function called every time a the mouse scroll wheel is moved
// yoffset = up-down
// xoffset = left-right
void mousewheel_callback(GLFWwindow* win, double xoffset, double yoffset) {
	double zoomFactor = pow(0.95,yoffset);

	glfwGetCursorPos(win, &prevx, &prevy);                  // get mouse position in window
	double xend = cam_x + prevx*(cam_width)/resx;			// convert screen position to world cooridnates
	double yend = cam_y + (resy - prevy)*cam_height/resy;	
	cam_x = (1.0-zoomFactor)*xend + zoomFactor*cam_x;		// update lower left corner
	cam_y = (1.0-zoomFactor)*yend + zoomFactor*cam_y;

	cam_width *= zoomFactor;
	cam_height *= zoomFactor;

	// zoom the camera by keeping the center constant and move the edges accordingly
	//cam_x = cam_x + cam_width*(1.0 - zoomFactor)/2.0;
	//cam_y = cam_y + cam_height*(1.0 - zoomFactor)/2.0;

	// cam_width *= zoomFactor;
	//cam_height *= zoomFactor;
}

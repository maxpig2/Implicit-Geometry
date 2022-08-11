
// std
#include <iostream>
#include <string>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"


using namespace std;
using namespace cgra;
using namespace glm;


void basic_model::draw(const glm::mat4 &view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;
	
	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

	//drawCylinder();
	//drawSphere();
	mesh.draw(); // draw
}

Application::Application(GLFWwindow *window) : m_window(window) {
	
	shader_builder sb;
    sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//oren_nayar_frag.glsl"));
	

	GLuint shader = sb.build();

	m_texture = rgba_image(CGRA_SRCDIR + std::string("/res//textures//Texture.png")).uploadTexture(GL_RGBA8,0);
	m_normal = rgba_image(CGRA_SRCDIR + std::string("/res//textures//NormalMap.png")).uploadTexture(GL_RGBA8,0);
	

	m_model.shader = shader;
	m_model.mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//teapot.obj")).build();
	m_model.color = vec3(1, 0, 0);
	m_cam_pos = vec2( 0, 0 );
}


void Application::render() {
	
	// retrieve the window height
	int width, height;
	glfwGetFramebufferSize(m_window, &width, &height); 

	m_windowsize = vec2(width, height); // update window size
	glViewport(0, 0, width, height); // set the viewport to draw to the entire window

	// clear the back-buffer
	glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	// enable flags for normal/forward rendering
	glEnable(GL_DEPTH_TEST); 
	glDepthFunc(GL_LESS);

	// projection matrix
	mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	mat4 view = translate(mat4(1), vec3( m_cam_pos.x, m_cam_pos.y, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw,   vec3(0, 1, 0));


	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

	
	//Textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_texture);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, m_normal);

	
	//Shadow Mapping
	GLuint depthMapFramebuff = 0;
	glGenFramebuffers(1, &depthMapFramebuff);
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFramebuff);
	//DepthTexture
	GLuint depthTexture;
	glGenTextures(1, &depthTexture);
	glBindTexture(GL_TEXTURE_2D, depthTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT16, 1024,1024,0 ,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glFramebufferTexture(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,depthTexture,0);
	glDrawBuffer(GL_NONE);

	vec3 lightInvDir = vec3(0.5f,2,2);
	//Compute MVP Matrix
	float near_plane = -10.0f;
	float far_plane = 20.0f;
	mat4 depthProjectionMatrix = ortho(-10.0f,10.0f,-10.0f,10.0f,near_plane,far_plane);

	mat4 depthViewMatrix = lookAt(lightInvDir, vec3(0,0,0), vec3(0,1,0));
	mat4 depthModelMatrix = mat4(1.0);
	mat4 depthMVP = depthProjectionMatrix * depthViewMatrix * depthModelMatrix;
	//glUniformMatrix4fv(depthMatrixID, 1, GL_FALSE, &depthMVP[0][0]);



	glBindFramebuffer(GL_FRAMEBUFFER,0);
	

	// draw the model
	m_model.draw(view, proj);
}


//Shape Attributes
int spherelatD = 15;
int spherelongD = 15;
int cubeSphereDivisions = 20;
int cubeSphereRadius = 5;
int torusLatD = 15;
int torusLongD = 15;
float torusRadiusIn = 2.0f;
float torusRadiusOut = 5.0f;

void Application::renderGUI() {

	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiSetCond_Once);
	ImGui::Begin("Options", 0);

	// display current camera parameters
	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);

	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	ImGui::Separator();
	ImGui::Text("Mesh");
	// example of how to use input boxes
	if (ImGui::Button("Teapot")){
		m_model.mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//teapot.obj")).build();
	}
	ImGui::SameLine();
	if (ImGui::Button("Sphere")){
		m_model.mesh = sphere_latlong(spherelatD,spherelongD,2);
	}
	ImGui::SameLine();
	if (ImGui::Button("Cube")){
		m_model.mesh = sphere_from_cube(cubeSphereDivisions, cubeSphereRadius);
	}
	ImGui::SameLine();
	if (ImGui::Button("Torus")){
		m_model.mesh = torus(torusLatD,torusLongD, torusRadiusIn, torusRadiusOut);
	}

	ImGui::Separator();
	ImGui::Text("Shaders");
	if (ImGui::Button("Simple Shader")){
		shader_builder sb;
   		 sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
		sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
		GLuint shader = sb.build();
		m_model.shader = shader;
	}
	ImGui::SameLine();
	if (ImGui::Button("Cook Torrance")){
		shader_builder sb;
   		 sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
		sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//cook_torrance_frag.glsl"));
		GLuint shader = sb.build();
		m_model.shader = shader;
	}
	ImGui::SameLine();
	if (ImGui::Button("Oren Nayar")){
		shader_builder sb;
   		 sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
		sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//oren_nayar_frag.glsl"));
		GLuint shader = sb.build();
		m_model.shader = shader;
	}
	ImGui::SameLine();
	if (ImGui::Button("Textures")){
		shader_builder sb;
   		 sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//texture_vert.glsl"));
		sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//texture_frag.glsl"));
		GLuint shader = sb.build();
		m_model.shader = shader;
	}


	ImGui::Separator();
	ImGui::Text("Sphere Attributes");
	if (ImGui::InputInt("Sphere Latitude Divisions", &spherelatD)) {m_model.mesh = sphere_latlong(spherelatD,spherelongD,2);}
	if (ImGui::InputInt("Sphere Longitude Divisions", &spherelongD)) {m_model.mesh = sphere_latlong(spherelatD,spherelongD,2);}

	ImGui::Separator();
	ImGui::Text("Cube Attributes");
	if (ImGui::InputInt("Cube Divisions", &cubeSphereDivisions)) {m_model.mesh = sphere_from_cube(cubeSphereDivisions, cubeSphereRadius);}
	if (ImGui::InputInt("Cube Radius", &cubeSphereRadius)) {m_model.mesh = sphere_from_cube(cubeSphereDivisions, cubeSphereRadius);}

	ImGui::Separator();	
	ImGui::Text("Torus Attributes");
	if (ImGui::InputInt("Torus Latitude Divisions", &torusLatD)) {m_model.mesh = torus(torusLatD,torusLongD, torusRadiusIn, torusRadiusOut);}
	if (ImGui::InputInt("Torus Longitude Divisions", &torusLongD)) {m_model.mesh = torus(torusLatD,torusLongD, torusRadiusIn, torusRadiusOut);}
	if (ImGui::InputFloat("Torus Inner Radius", &torusRadiusIn)) {m_model.mesh = torus(torusLatD,torusLongD, torusRadiusIn, torusRadiusOut);}
	if (ImGui::InputFloat("Torus Outer Radius", &torusRadiusOut)) {m_model.mesh = torus(torusLatD,torusLongD, torusRadiusIn, torusRadiusOut);}

	// finish creating window
	ImGui::End();
}


void Application::cursorPosCallback(double xpos, double ypos) {
	vec2 whsize = m_windowsize / 2.0f;

	double y0 = glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f);
	double y = glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f);
	double dy = -( y - y0 );

	double x0 = glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f);
	double x = glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f);
	double dx = x - x0;

	if (m_leftMouseDown) {
		// clamp the pitch to [-pi/2, pi/2]
		m_pitch += float( acos(y0) - acos(y) );
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

		// wrap the yaw to [-pi, pi]
		m_yaw += float( acos(x0) - acos(x) );
		if (m_yaw > pi<float>()) 
			m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) 
			m_yaw += float(2 * pi<float>());
	} else if ( m_rightMouseDown ) {
		m_distance += dy * 10;
	} else if ( m_middleMouseDown ) {
		m_cam_pos += vec2( dx, dy ) * 10.f;
	}

	// updated mouse position
	m_mousePosition = vec2(xpos, ypos);
}


void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods; // currently un-used

	// capture is left-mouse down
	if (button == GLFW_MOUSE_BUTTON_LEFT)
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
	else if (button == GLFW_MOUSE_BUTTON_RIGHT)
		m_rightMouseDown = (action == GLFW_PRESS);
	else if (button == GLFW_MOUSE_BUTTON_MIDDLE)
		m_middleMouseDown = (action == GLFW_PRESS);
}


void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset; // currently un-used
	m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
	(void)key, (void)scancode, (void)action, (void)mods; // currently un-used
}


void Application::charCallback(unsigned int c) {
	(void)c; // currently un-used
}



gl_mesh Application::sphere_latlong(int latitudeDivisions, int longitudeDivisions, float radius){
	mesh_builder sphere;
	//Object Information
	vector<vec3> positions;
	vector<vec3> normals;
	vector<vec2> uvs;
	vector<vec3> tangents;
	vector<int> indices;
	//Steps
	float latitudeStep = (2*PI)/latitudeDivisions;	//Azimuth Direction
	float longitudeStep = (PI)/longitudeDivisions; //Elevation Direction
	//Angles
	float azimuthalAngle;
	float elevationAngle;
	//Iterate in a circle. 
	for (int i = 0; i <= longitudeDivisions; i++){
		elevationAngle = i*longitudeStep;
		for (int j = 0; j <= latitudeDivisions; j++) {
			azimuthalAngle = j*latitudeStep;
			//Calculate Vertex Positions
			float x = radius * (cosf(azimuthalAngle) * sinf(elevationAngle));
			float y = radius * (sinf(azimuthalAngle) * sinf(elevationAngle));
			float z = radius * (cosf(elevationAngle));
			positions.push_back(vec3(x,y,z));
			//Calculate Normals
			float oneOverRadius = 1.0f/radius;
			vec3 nnnn = vec3(x*oneOverRadius,y*oneOverRadius,z*oneOverRadius);
			normals.push_back(vec3(x*oneOverRadius,y*oneOverRadius,z*oneOverRadius));
			//Calculate UVs
			uvs.push_back(vec2((float)i/latitudeDivisions, (float)j/longitudeDivisions));
			tangents.push_back(normalize(vec3(cross(vec3(0,1,0),vec3(x,y,z)))));
		}
	}
	//Triangles!
	for (int i = 0; i < longitudeDivisions; i++){
		int m = i * (latitudeDivisions + 1);
		int n = m + latitudeDivisions + 1;
		for (int j = 0; j < latitudeDivisions; j++, m++, n++) {
            	//Define the bottom triangle
                 indices.push_back(m);
                 indices.push_back(n);
                 indices.push_back(m + 1);
                  
				//Define the top triangle
                indices.push_back(m + 1);
                indices.push_back(n);
                indices.push_back(n + 1);
		}
	}
	for (int k = 0; k < indices.size(); ++k) {
             sphere.push_index(k);
             sphere.push_vertex(mesh_vertex{
                 positions[indices[k]],
                 normals[indices[k]],
                 uvs[indices[k]],
				 tangents[indices[k]]
             });
         }

	//We made a sphere, lets return it!
	return sphere.build();
}

gl_mesh Application::sphere_from_cube(float subdivisions, float radius){
	mesh_builder cube;
	//Object Information
	vector<vec3> positions;
	vector<vec3> normals;
	vector<vec2> uvs;
	vector<int> indices;

	//Front
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = ((i-(subdivisions/2))/(subdivisions/2));
			float y = ((j-(subdivisions/2))/(subdivisions/2));
			float z = 1;
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
	//Left
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = -1;
			float y = ((j-(subdivisions/2))/(subdivisions/2));
			float z = ((i-(subdivisions/2))/(subdivisions/2));
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
	//Back
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = -((i-(subdivisions/2))/(subdivisions/2));
			float y = ((j-(subdivisions/2))/(subdivisions/2));
			float z = -1;
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
	//Right
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = 1;
			float y = ((j-(subdivisions/2))/(subdivisions/2));
			float z = ((i-(subdivisions/2))/(subdivisions/2));
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
	//Top
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = -((i-(subdivisions/2))/(subdivisions/2));
			float y = 1;
			float z = ((j-(subdivisions/2))/(subdivisions/2));
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
	//Bottom
	for (int i = 0; i <= subdivisions; i++) {
		for (int j = 0; j <= subdivisions; j++) {
			float x  = -((i-(subdivisions/2))/(subdivisions/2));
			float y = -1;
			float z = ((j-(subdivisions/2))/(subdivisions/2));
			float xSquared = x*x;
			float ySquared = y*y;
			float zSquared = z*z;
			x = x * sqrt(1-(ySquared/2)-(zSquared/2)+((ySquared*zSquared)/3));
			y = y * sqrt(1-(zSquared/2) - (xSquared/2) + ((zSquared*xSquared)/3));
			z = z * sqrt(1-(xSquared/2) - (ySquared/2) + ((xSquared*ySquared)/3));	
			x *= radius;
			y *= radius;
			z *= radius;
			positions.push_back(vec3(x,y,z));
			normals.push_back(vec3(0,1,0));
		}
	}
for (int i = 0; i <= subdivisions * 6 + 4; i++){
	if (i != subdivisions && i != subdivisions * 2 + 1 && i != subdivisions * 3 + 2 && i != subdivisions * 4 + 3 && i != subdivisions * 5 + 4) {
		int m = i * (subdivisions + 1);
		int n = m + subdivisions + 1;
		for (int j = 0; j < subdivisions; j++, m++, n++) {
            	//Define the bottom triangle
                 indices.push_back(m);
                 indices.push_back(n);
                 indices.push_back(m + 1); 
				//Define the top triangle
                indices.push_back(m + 1);
                indices.push_back(n);
                indices.push_back(n + 1);
		}
	}
}
	for (int k = 0; k < indices.size(); ++k) {
             cube.push_index(k);
             cube.push_vertex(mesh_vertex{
                 positions[indices[k]],
                 normals[indices[k]]
             });
         }
	return cube.build();
}

gl_mesh Application::torus(int latitudeDivisions, int longitudeDivisions, float radiusInner, float radiusOuter){

	mesh_builder torus;
	//Object Information
	vector<vec3> positions;
	vector<vec3> normals;
	vector<vec2> uvs;
	vector<int> indices;
	//Variables
	float R = radiusOuter;
	float r = radiusInner;
	float a = 0.5*(R-r);
	float c = 0.5*(R+r);
	//Steps
	float latitudeStep = (2*PI)/latitudeDivisions;	//Azimuth Direction
	float longitudeStep = (2*PI)/longitudeDivisions; //Elevation Direction
	//Angles
	float azimuthalAngle;
	float elevationAngle;

	//Iterate in a circle. 
	for (int i = 0; i <= longitudeDivisions; i++){
		elevationAngle = i*longitudeStep;
		for (int j = 0; j <= latitudeDivisions; j++) {
			azimuthalAngle = j*latitudeStep;
			//Calculate Vertex Positions
			float x =  ((c + a*(cosf(elevationAngle))) * cosf(azimuthalAngle));
			float y = ((c + a*(cosf(elevationAngle))) * sinf(azimuthalAngle));
			float z = a * (sinf(elevationAngle));
			positions.push_back(vec3(x,y,z));
			//Calculate Normals
			float normalizer = 1.0f/radiusInner;
			normals.push_back(vec3(x*normalizer,y*normalizer,z*normalizer));
			//Calculate UVs
			uvs.push_back(vec2((float)j/latitudeDivisions, (float)i/longitudeDivisions));
		}
	}
	//Triangles!
	for (int i = 0; i < longitudeDivisions; i++){
		int m = i * (latitudeDivisions + 1);
		int n = m + latitudeDivisions + 1;
		for (int j = 0; j < latitudeDivisions; j++, m++, n++) {
            	//Define the bottom triangle
                 indices.push_back(m);
                 indices.push_back(n);
                 indices.push_back(m + 1);
                  
				//Define the top triangle
                indices.push_back(m + 1);
                indices.push_back(n);
                indices.push_back(n + 1);
		}
	}
	for (int k = 0; k < indices.size(); ++k) {
             torus.push_index(k);
             torus.push_vertex(mesh_vertex{
                 positions[indices[k]],
                 normals[indices[k]]
             });
         }

	//We made a Torus, lets return it!
	return torus.build();
}

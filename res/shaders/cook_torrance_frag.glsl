#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform sampler2D uColourTexture;
uniform sampler2D uNormalTexture;

// viewspace data (this must match the output of the fragment shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
	mat3 tbn;
} f_in;

// framebuffer output
out vec4 fb_color;


void main() {
	vec3 eye = normalize(-f_in.position);
	vec3 light = vec3(0.5,0.5,-1);
	
	float PI = 3.1415;
	//Colours
	vec3 objectColour = vec3(1,0,0);
	vec3 specularColour = vec3(1,1,1);
	vec3 lightColour = vec3(1,1,1);

	vec3 N = normalize(f_in.normal); //Unit surface normal

	vec3 V = eye; //Unit vector in direction of the viewer
	vec3 L = normalize(-light); //Unit vector in the direction of a light
	vec3 H = normalize(V+L); //Unit Angular Bisector of V and L

	float G = min(1,min(((2*dot(N,H)*dot(N,V))/dot(V,H)),((2*dot(N,H)*dot(N,L))/dot(V,H)))); //Geometrical Attenuation Factor

	float a = acos(dot(N,H)); //Angle between N and H
	float m = 0.1; //Root mean square slope of facets
	float D = (exp((-pow(tan(a),2)/ pow(m,2)))) / (PI*pow(m,2)*pow(cos(a),4)) ; //Facet slope distribution function

	float c = dot(V,H);
	float n = 0.1;//indexOfRefraction 
	float F0 = pow(n-1,2)/pow(n+1,2);
	float F = F0 + (1 -F0) * pow(1.0 - c, 5);//Reflectance of a perfectly smooth surface

	float Rs = (F/PI) * ((D*G)/(dot(N,L)*dot(N,V))); //Specular Bidirectional Reflectance
	Rs = (D*G*F)/(4*dot(N,L)*dot(N,V));

	float diffuse = max(dot(N,L),0.0);
	float k = 0.0; // Extinction Coefficient

	if (diffuse == 0.0) {lightColour = vec3(0,0,0);}
	vec3 fragColour = (diffuse * lightColour * objectColour) + (lightColour * specularColour * (k + Rs * (1.0 - k))); // Final Result

	// output to the frambuffer
	fb_color = vec4(fragColour, 1);
}











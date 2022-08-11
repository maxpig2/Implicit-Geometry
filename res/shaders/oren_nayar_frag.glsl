#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

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
	float PI = 3.1415;
	// calculate lighting (hack)
	vec3 eye = normalize(-f_in.position);
	vec3 light = vec3(0.5,0.5,-1);
	//Variables
	vec3 L = normalize(-light); //Light Source Direction
	vec3 V = eye; //Viewing Direction
	vec3 N = normalize(f_in.normal); //Unit surface normal

	//Calculate S and T
	float s = dot(L,V) - (dot(N,L)*dot(N,V));
	float t = 0.0;
	if (s <= 0) {
		t = 1.0;
	} else {
		t = max(dot(N,L),dot(N,V));
	}
	
	float sigma = 0.1f;
	float sigmaSquared = sigma*sigma;
	float p = 3.0;
	float A = (1/PI)* ( 1 - 0.5 * (sigmaSquared/(sigmaSquared+0.33) + (0.17*p) * (sigmaSquared/(sigmaSquared+0.13))  )       );
	float B = (1/PI) * (0.45 * (sigmaSquared / (sigmaSquared + 0.09)));

	//"Improved" Oren Nayar Implementation
	//float A = 1/(PI + ((PI/2.0) - (2.0/3.0)) * sigma);
	//float B = sigma/(PI + ((PI/2.0) - (2.0/3.0)) * sigma);

	float diffuse = (p * dot(N,L))*(A + (B * (s/t)));

	float specularity  = 0.0;
	float shininess = 15;
	if (diffuse > 0){
	vec3 H = normalize(L+V);
	float specularAngle = max(dot(H,N),0.0);
	specularity = pow(specularAngle,shininess);
	}

	vec3 objectColour = vec3(1,0,0);
	vec3 fragColour = (diffuse * objectColour) + specularity;

	// output to the frambuffer
	fb_color = vec4(fragColour, 1);
}
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

	//Normal Data
	vec3 N = texture(uNormalTexture, f_in.textureCoord).rgb;
	N = normalize((2.0 * N) - 1.0);
	vec3 colour = texture(uColourTexture, f_in.textureCoord).rgb;

	light = normalize(f_in.tbn * -light);
	float diffuse = max(dot(N, light), 0.0);

	float specularity = 0.0;
	
	if (diffuse > 0.0){
		vec3 V = normalize(-f_in.position);
		vec3 H = normalize(light + V);
		float specularAngle = max(dot(H,N), 0.0);
		specularity = pow(specularAngle, 30);
	}
	
	float ambient = 0.1;
	vec3 fragColour =  (diffuse * (colour)) + (ambient * colour) + specularity;

	// output to the frambuffer
	fb_color = vec4(fragColour, 1);
}


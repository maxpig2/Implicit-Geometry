#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

// mesh data
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;
layout(location = 3) in vec3 aTangent;

// model data (this must match the input of the vertex shader)
out VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
	mat3 tbn;
} v_out;

void main() {	



	//Calculate TBN
	//vec3 T = normalize((uModelViewMatrix * vec4(aTangent,1.0))).xyz;
	//vec3 N = normalize((uModelViewMatrix * vec4(aNormal,1.0))).xyz;
	
	//https://learnopengl.com/Advanced-Lighting/Normal-Mapping
	vec3 T = normalize(vec3(uModelViewMatrix * vec4(aTangent,1.0)));
	vec3 N = normalize(vec3(uModelViewMatrix * vec4(aNormal,1.0)));
	vec3 B = normalize(vec3(uModelViewMatrix * vec4(cross(N,T),1.0)));

	//vec3 B = cross(N,T).xyz;
	mat3 TBN = mat3(T,B,N);
	v_out.tbn = transpose(TBN);

	v_out.position = (uModelViewMatrix * vec4(aPosition, 1)).xyz;
	v_out.normal = normalize((uModelViewMatrix * vec4(aNormal, 0)).xyz);
	v_out.textureCoord = aTexCoord;

	// set the screenspace position (needed for converting to fragment data)
	gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aPosition, 1);
}
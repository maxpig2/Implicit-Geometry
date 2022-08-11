#version 330 core

// mesh data
layout(location = 0) in vec3 aPosition;

//Uniforms
uniform mat4 depthMVP;

// model data (this must match the input of the vertex shader)
out VertexData {
} v_out;

void main() {
		gl_Position = depthMVP * vec4(aPosition, 1);
}
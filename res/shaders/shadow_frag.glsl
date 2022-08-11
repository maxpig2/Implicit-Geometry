#version 330 core


// viewspace data (this must match the output of the fragment shader)
in VertexData {
} f_in;

// framebuffer output
layout (location = 0) out float fragmentdepth;

void main() {
	fragmentdepth = gl_FragCoord.z
}
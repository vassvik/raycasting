#version 430 core

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec2 uvCoord;

out vec2 uv;

void main(){
	uv = uvCoord;
    gl_Position = vec4(vertexPosition, 1.0);
}


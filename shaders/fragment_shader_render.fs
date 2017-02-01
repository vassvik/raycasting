#version 430 core

in vec2 uv;

uniform sampler2D framebufferTexture;

out vec3 color;

void main()
{
    color = texture(framebufferTexture, uv).xyz;;
}
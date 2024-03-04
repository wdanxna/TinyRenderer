#pragma once
void rasterize(Matrix &projection, Matrix &view, Matrix &modelworld, Matrix &vp, float zbuffer[250000], TGAImage &framebuffer);

void rasterize(Model &model, Matrix &projection, Matrix &view, Matrix &modelworld, Matrix &vp, float zbuffer[250000], TGAImage &framebuffer);

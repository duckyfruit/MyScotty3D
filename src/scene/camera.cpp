
#include "camera.h"
#include "../gui/manager.h"
#include "../pathtracer/samplers.h"
#include "../test.h"

std::pair<Ray, float> Camera::sample_ray(RNG &rng, uint32_t px, uint32_t py) {
	//A3T1 - step 1 - camera rays

	//Sample a ray that starts at the origin and passes through pixel (px,py) + random offset on the sensor plane.
	//
	//Because cameras look down the -z axis, the "sensor plane" is
	// the rectangle from (-w/2,-h/2,-1) to (w/2,h/2,-1)
	// where:
	//  h is such that the angle made by (0,-h/2,-1) - (0,0,0) - (0,h/2,-1) is `vertical_fov`
	//  and w / h is given by `aspect_ratio`.
	//
	//The relationship between sensor pixels and the sensor plane is such that
	//  sensor pixel location (0,0) maps to (-w/2,-h/2,-1),
	//  and sensor pixel location (film.width,film.height) maps to (w/2,h/2,-1).

	//Compute the position on the sensor (in pixels):
	Samplers::Rect s;
	Vec2 offset = s.sample(rng);
	float offset_pdf = s.pdf(offset);
	Vec2 sensor_pixel = Vec2(float(px), float(py)) + offset;


	//TODO: Transform from sensor pixels into world position on the sensor plane
	//map from one place to another --> put px to somewhere between (-w/2, -h/2 , -1) and (w/2, h/2 ,-1)
	//sensor plane is w long and h high, which means we can divide sensor pixel location by w and add w/2
	//(void)sensor_pixel;
	float newh = tan(Radians(vertical_fov/2.0f)) * 2.0f; //use fov to get to new height
	float neww = aspect_ratio * newh; //get new width using aspect ratio
	float w = float(film.width);
	float h = float(film.height); 
	float mapX = ((sensor_pixel.x / w) - 0.5f) * neww; //get the new direction from mapping
	float mapY = ((sensor_pixel.y / h) - 0.5f) * newh;

	
	//Build ray:
	Ray ray;
	ray.point = Vec3(); //ray should start at the origin
	ray.dir = Vec3(mapX,mapY,-1.0f); //TODO: compute from sensor plane position
	ray.depth = film.max_ray_depth; //rays should, by default, go as deep as the max depth parameter allows
	//printf("%f %f \n %f %f \n %f %f\n %f %f\n", float(px), float(py), float(film.width), float(film.height), mapX, mapY, ray.dir.x, ray.dir.y);
   	return {ray, offset_pdf};
}


Mat4 Camera::projection() const {
	return Mat4::perspective(vertical_fov, aspect_ratio, near_plane);
}

GL::Lines Camera::to_gl() const {

	GL::Lines cage;

	// float ap = near_plane * std::tan(Radians(vertical_fov) / 2.0f);
	float h = 2.0f * std::tan(Radians(vertical_fov) / 2.0f);
	float w = aspect_ratio * h;

	Vec3 tr = Vec3(0.5f * w, 0.5f * h, -1.0f) * focal_dist;
	Vec3 tl = Vec3(-0.5f * w, 0.5f * h, -1.0f) * focal_dist;
	Vec3 br = Vec3(0.5f * w, -0.5f * h, -1.0f) * focal_dist;
	Vec3 bl = Vec3(-0.5f * w, -0.5f * h, -1.0f) * focal_dist;

	Vec3 ftr = Vec3(0.5f * aperture_size, 0.5f * aperture_size, -near_plane);
	Vec3 ftl = Vec3(-0.5f * aperture_size, 0.5f * aperture_size, -near_plane);
	Vec3 fbr = Vec3(0.5f * aperture_size, -0.5f * aperture_size, -near_plane);
	Vec3 fbl = Vec3(-0.5f * aperture_size, -0.5f * aperture_size, -near_plane);

	cage.add(ftl, ftr, Gui::Color::black);
	cage.add(ftr, fbr, Gui::Color::black);
	cage.add(fbr, fbl, Gui::Color::black);
	cage.add(fbl, ftl, Gui::Color::black);

	cage.add(ftr, tr, Gui::Color::black);
	cage.add(ftl, tl, Gui::Color::black);
	cage.add(fbr, br, Gui::Color::black);
	cage.add(fbl, bl, Gui::Color::black);

	cage.add(bl, tl, Gui::Color::black);
	cage.add(tl, tr, Gui::Color::black);
	cage.add(tr, br, Gui::Color::black);
	cage.add(br, bl, Gui::Color::black);

	return cage;
}

bool operator!=(const Camera& a, const Camera& b) {
	return a.vertical_fov != b.vertical_fov || a.aspect_ratio != b.aspect_ratio || a.near_plane != b.near_plane
		   || a.aperture_shape != b.aperture_shape || a.aperture_size != b.aperture_size || a.focal_dist != b.focal_dist
	       || a.film.width != b.film.width || a.film.height != b.film.height
	       || a.film.samples != b.film.samples || a.film.max_ray_depth != b.film.max_ray_depth
	       || a.film.sample_pattern != b.film.sample_pattern
	;
}

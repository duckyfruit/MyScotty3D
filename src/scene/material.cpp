
#include "material.h"
#include "../util/rand.h"

namespace Materials {

Vec3 reflect(Vec3 dir) {
	
	//A3T5 Materials - reflect helper

    // Return direction to incoming light that would be
	// reflected out in direction dir from surface
	// with normal (0,1,0)
	//Vec3 normal = Vec3{0.f, 1.f, 0.f};
	//Vec3 ref = -1.f * dir + (2.f* dir.y) * normal;
    return Vec3{-dir.x, dir.y, -dir.z};
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {
	//A3T5 Materials - refract helper

	// Use Snell's Law to refract out_dir through the surface.
	// Return the refracted direction. Set was_internal to true if
	// refraction does not occur due to total internal reflection,
	// and false otherwise.
	Vec3 in_dir;
	Vec3 normal = Vec3{0.f, 1.f, 0.f};
	float cosThetaT = dot(out_dir,normal);

	float ni;
	float nt;
	if(out_dir.y < 0.f)
	{
		ni = 1.f;
		nt = index_of_refraction;
		normal = -normal;
		cosThetaT = -cosThetaT;
	}
		
	else
	{
		nt = 1.f;
		ni = index_of_refraction;
	}

	//Vec3 unit = (out_dir / sqrt(pow(out_dir.x,2.0f) + pow(out_dir.y,2.0f) + pow(out_dir.z,2.0f))); //they're all unit
	
	
	float sin2ThetaT =1.0f - pow(cosThetaT,2.0f);
	float sin2ThetaI = pow((nt/ni), 2.0f) * sin2ThetaT;
	float cosThetaI = sqrt(1.0f - sin2ThetaI);

	if(sin2ThetaI >= 1.0f)
	was_internal = true;
	else
	was_internal = false;

	in_dir = (nt/ni) * -out_dir + ((nt/ni) * cosThetaT - cosThetaI) * normal;
	return in_dir;

	/*float fetat = acos(1.f/sqrt(pow(out_dir.unit().x,2.f) + pow(out_dir.unit().y,2.f) + pow(out_dir.unit().z,2.f)));

	float fetai = asin((nt * sin(fetat)) / ni);
	
	if(sin((fetai)) > 1)
	was_internal = true;
	else
	was_internal = false;

	//if(out_dir.y < 0.f)
	//in_dir = ((nt/ni) * (-1.f *out_dir) + ((nt/ni )* cos(fetat) - cos(fetai)) * (-1.f * normal)); //two cases
	//else
	in_dir = ((nt/ni) * (-1.f *out_dir) + ((nt/ni )* cos(fetat) - cos(fetai)) * normal); //two cases
	//enter and exit
	// The surface normal is (0,1,0) */

	return in_dir;
}

float schlick(Vec3 in_dir, float index_of_refraction) {
	//A3T5 Materials - Schlick's approximation helper

	// Implement Schlick's approximation of the Fresnel reflection factor.
	float ni;
	float nt;
	Vec3 normal = Vec3{0.f, 1.f, 0.f};
	float cosThetaT = dot(in_dir,normal);
	if(in_dir.y < 0.f)
	{
		ni = 1.f;
		nt = index_of_refraction;
		normal = normal * -1;
		cosThetaT = cosThetaT * -1;
	}
		
	else
	{
		nt = 1.f;
		ni = index_of_refraction;
	}

	//Vec3 unit = (in_dir / sqrt(pow(in_dir.x,2.0f) + pow(in_dir.y,2.0f) + pow(in_dir.z,2.0f)));
	
	//float fetat = acos(1.f/sqrt(pow(unit.x,2.f) + pow(unit.y,2.f) + pow(unit.z,2.f)));

	//float fetai = asin((nt * sin(fetat)) / ni);
	float sin2ThetaT =1.0f - pow(cosThetaT,2.0f);
	float sin2ThetaI = pow((nt/ni), 2.0f) * sin2ThetaT;
	float cosThetaI = sqrt(1.0f - sin2ThetaI);

	float rparallel = (nt * cosThetaI - ni * cosThetaT)/ (nt * cosThetaI + ni * cosThetaT);
	float rperpen =  (ni * cosThetaI - nt * cosThetaT)/ (ni * cosThetaI + nt * cosThetaT);

	//float rparallel =  (nt * cos(fetai) - ni * cos(fetat))/ (nt * cos(fetai) + ni * cos(fetat));
	
	//float rperpen =  (ni * cos(fetai) - nt * cos(fetat))/ (ni * cos(fetai) + nt * cos(fetat));



	return (0.5f * (pow(rparallel, 2.f) + pow(rperpen, 2.f)));
}

Spectrum Lambertian::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF evaluation

    // Compute the ratio of outgoing/incoming radiance when light from in_dir
    // is reflected through out_dir: (albedo / PI_F) * cos(theta).
    // Note that for Scotty3D, y is the 'up' direction.

	//y is up, we got a radiance ray coming in, and a radiance ray going out, and uv 

	//we need BRDF function
	//we need incident radiance
	//computes the ratio of outgoing to incoming radiance given a pair of directions
	//how do I get incident radiance?

	//theta = uv in arctan?

	float cos = dot(Vec3{0.f,1.f,0.f}, in);
	//just the y component :>>
	//float mag = sqrt(pow(in.x, 2.0f) + pow(in.y, 2.0f) + pow(in.z, 2.0f)); 
	Spectrum scatterF = ((((albedo.lock() -> evaluate(uv))/PI_F)) * std::max(0.0f,cos)); //scatter function at current point
	//Spectrum ratio = Spectrum((scatterF[0]), (scatterF[1]), (scatterF[2]));


    return scatterF;
}

Scatter Lambertian::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF scattering
	//Select a scattered light direction at random from the Lambertian BSDF

	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this will be useful

	Scatter ret;
	//TODO: sample the direction the light was scatter from from a cosine-weighted hemisphere distribution:
	ret.direction = sampler.sample(rng);
	
	//TODO: compute the attenuation of the light using Lambertian::evaluate():
	ret.attenuation = evaluate(out,ret.direction,uv);

	return ret;
}

float Lambertian::pdf(Vec3 out, Vec3 in) const {
	//A3T4: Materials - Lambertian BSDF probability density function
    // Compute the PDF for sampling in_dir from the cosine-weighted hemisphere distribution.
	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this might be handy!
	
    return sampler.pdf(in);;
}

Spectrum Lambertian::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Lambertian::display() const {
	return albedo;
}

void Lambertian::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(albedo);
}

Spectrum Mirror::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Mirror::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5: mirror

	// Use reflect to compute the new direction
	// Don't forget that this is a discrete material!
	// Similar to albedo, reflectance represents the ratio of incoming light to reflected light

    Scatter ret;
    ret.direction = reflect(out);
	//float pi = 2*acos(0.0f);
    ret.attenuation = reflectance.lock() -> evaluate(uv);
    return ret;
}

float Mirror::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Mirror::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Mirror::display() const {
	return reflectance;
}

void Mirror::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
}

Spectrum Refract::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Refract::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - refract

	// Use refract to determine the new direction - what happens in the total internal reflection case?
    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
	// Don't forget that this is a discrete material!
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction
	bool wasInt;
    Scatter ret;
	
    ret.direction = refract(out,ior,wasInt);
	if(wasInt)
	{
		ret.direction = reflect(out);

		ret.attenuation = Spectrum(1.f);

	}
	else
	{
		if(out.y < 0.f)
		ret.attenuation = transmittance.lock()->evaluate(uv) * 1.f/ior;
		else 
		ret.attenuation = transmittance.lock()->evaluate(uv) * ior;
	}
	


	return ret;


}

float Refract::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Refract::emission(Vec2 uv) const {
	return {};
}

bool Refract::is_emissive() const {
	return false;
}

bool Refract::is_specular() const {
	return true;
}

bool Refract::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Refract::display() const {
	return transmittance;
}

void Refract::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(transmittance);
}

Spectrum Glass::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Glass::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - glass

    // (1) Compute Fresnel coefficient. Tip: Schlick's approximation.
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
    // What happens upon total internal reflection?
    // When debugging Glass, it may be useful to compare to a pure-refraction BSDF
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

    Scatter ret;
	bool wasInt;

	
	ret.direction = refract(out,ior,wasInt);
	if(wasInt)
	{
		ret.direction = reflect(out);
		ret.attenuation = reflectance.lock() -> evaluate(uv);
	}
	else
	{
		float Fresnal = schlick(ret.direction,ior);
		if(rng.coin_flip(Fresnal))
		{
			ret.direction = reflect(out);
			ret.attenuation = reflectance.lock() -> evaluate(uv);
			
		}
		else
		{
			if(out.y < 0.f)
			ret.attenuation = transmittance.lock()->evaluate(uv) * 1.f/ior;
			else 
			ret.attenuation = transmittance.lock()->evaluate(uv) * ior;
		}
	

		
	}
	
	return ret;

 
}

float Glass::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Glass::emission(Vec2 uv) const {
	return {};
}

bool Glass::is_emissive() const {
	return false;
}

bool Glass::is_specular() const {
	return true;
}

bool Glass::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Glass::display() const {
	return transmittance;
}

void Glass::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
	f(transmittance);
}

Spectrum Emissive::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Emissive::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	Scatter ret;
	ret.direction = {};
	ret.attenuation = {};
	return ret;
}

float Emissive::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Emissive::emission(Vec2 uv) const {
	return emissive.lock()->evaluate(uv);
}

bool Emissive::is_emissive() const {
	return true;
}

bool Emissive::is_specular() const {
	return true;
}

bool Emissive::is_sided() const {
	return false;
}

std::weak_ptr<Texture> Emissive::display() const {
	return emissive;
}

void Emissive::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(emissive);
}

} // namespace Materials

bool operator!=(const Materials::Lambertian& a, const Materials::Lambertian& b) {
	return a.albedo.lock() != b.albedo.lock();
}

bool operator!=(const Materials::Mirror& a, const Materials::Mirror& b) {
	return a.reflectance.lock() != b.reflectance.lock();
}

bool operator!=(const Materials::Refract& a, const Materials::Refract& b) {
	return a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Glass& a, const Materials::Glass& b) {
	return a.reflectance.lock() != b.reflectance.lock() ||
	       a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Emissive& a, const Materials::Emissive& b) {
	return a.emissive.lock() != b.emissive.lock();
}

bool operator!=(const Material& a, const Material& b) {
	if (a.material.index() != b.material.index()) return false;
	return std::visit(
		[&](const auto& material) {
			return material != std::get<std::decay_t<decltype(material)>>(b.material);
		},
		a.material);
}

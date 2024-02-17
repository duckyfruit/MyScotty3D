
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	// A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'

	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f); //at 1.5 3.0

	//printf("%f %f\n" , x, y);

	int floorx = (int)(floor(x-0.5f)); //x'
	int floory = (int)(floor(y-0.5f)); //y'
	
	int floorxincr = (int)(floorx + 1);

	int flooryincr = (int)(floory + 1);

	float deltax = abs(x - (floorx + 0.5f)); // change in x
	float deltay = abs(y - (floory + 0.5f)); // change in y


	if((floorxincr == (int)(image.w) )|| (flooryincr == (int)(image.h)) )
	{
		if((floorxincr == (int)(image.w)))
		{ 
			floorxincr = floorx;
		}
		if((flooryincr == (int)(image.h)))
		{ 
			flooryincr = floory; 
		}
			
	}
	if(floorx < 0 || floory < 0)
	{
		if((floorx < 0))
		{ 
			floorx = 0;
		}
		if((floory < 0))
		{ 
			floory = 0;
		}

	}

	auto t1 = image.at(floorx,floory); //texture lookups of x, y x+1, y  x, y+1 x+1, y+1
	auto t2 = image.at(floorxincr,floory); //texture lookups of x, y x+1, y  x, y+1 x+1, y+1
	auto t3 = image.at(floorx,flooryincr); //texture lookups of x, y x+1, y  x, y+1 x+1, y+1
	auto t4 = image.at(floorxincr,flooryincr); //texture lookups of x, y x+1, y  x, y+1 x+1, y+1
	//printf("%d %d\n" , floorx, floory);

	//printf("%d %d\n" , floorxincr, flooryincr);

	

	//printf("%f %f\n" , deltax, deltay);

	
	//auto t2 = image.at(floorxincr, floory);
	//auto t3 = image.at(floorx , flooryincr);
	//auto t4 = image.at(floorxincr, flooryincr);

	std::string test = to_string(t1);
	std::string test2 = to_string(t2);
	std::string test3 = to_string(t3);
	std::string test4 = to_string(t4);
	/*printf(" t1: %s\n",test.c_str());
	printf(" t2: %s\n",test2.c_str());
	printf(" t3: %s\n",test3.c_str());
	printf(" t4: %s\n",test4.c_str()); */

	Spectrum tx;
	Spectrum ty;


	tx = (1.0f - deltax) * t1 + deltax * t2;
	ty = (1.0f - deltax) * t3 + deltax * t4;
 

	std::string testx = to_string(tx);
	std::string testy = to_string(ty);

	//printf(" tx: %s\n",testx.c_str());
	//printf(" ty: %s\n",testy.c_str());
	return (1.0f - deltay) * tx + deltay * ty;
}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	// A1T6: sample_trilinear
	//TODO: implement trilinear sampling strategy on using mip-map 'levels'
	int floord = (int)(floor(lod));
	int floordincr = (int)(floord + 1);
	float deltad = abs(lod - floord);
	Spectrum td;
	Spectrum td1;


	//printf("%d\n" , floord);
	//printf("%d\n" , floordincr);

	//printf("%d\n" , (int)(levels.size()));


	if((floord) >= (int)(levels.size()))
	{
		floord = (int)(levels.size());
		floordincr = floord;
	}

	if(floord <= 0) //use base case
	{
		floord = 0;
		td = sample_bilinear(base,uv);
		td1 = sample_bilinear(levels[floord],uv);
		return ((1.0f - deltad) * td + deltad * td1);

	}
	else
	{
		td = sample_bilinear(levels[floord - 1],uv);
		td1 = sample_bilinear(levels[floordincr - 1],uv);
		return ((1.0f - deltad) * td + deltad * td1);

	}

	//printf("%f\n" , uv.x);
	//printf("%f\n" , uv.y);

	//printf("%d\n" , floord);

	//printf("%f\n" , deltad);

	
	//return sample_nearest(base, uv); //placeholder so image doesn't look blank
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);
		Spectrum x1;
		Spectrum y1;
		Spectrum x2;
		Spectrum y2;

		for(int j = 0; j < (int)(src.h); j +=2 )
		{
			for(int i = 0; i < (int)(src.w); i+=2)
			{
				if((i/2 < (int)(dst.w)) && (j/2 < (int)(dst.h)))
				{
					x1 = src.at(i,j);
					x2 = src.at(i+1,j);
					y1 = src.at(i,j+1);
					y2 = src.at(i+1,j+1);
					dst.at(i/2, j/2) = (x1 + x2 + y1 + y2)/4.0f;
				}

				/*std::string test = to_string(x1);
				std::string test2 = to_string(x2);
				std::string test3 = to_string(y1);
				std::string test4 = to_string(y2);

				printf(" t2: %s\n",test.c_str());
				printf(" t2: %s\n",test2.c_str());
				printf(" t2: %s\n",test3.c_str());
				printf(" t2: %s\n",test4.c_str()); */

				
			}
		}


		// A1T6: generate
		//TODO: Write code to fill the levels of the mipmap hierarchy by downsampling

		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.

	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (image.w == 0 && image.h == 0) return Spectrum();
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}

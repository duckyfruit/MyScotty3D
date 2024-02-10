// clang-format off
#include "pipeline.h"

#include <iostream>

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include "framebuffer.h"
#include "sample_pattern.h"
template<PrimitiveType primitive_type, class Program, uint32_t flags>
void Pipeline<primitive_type, Program, flags>::run(std::vector<Vertex> const& vertices,
                                                   typename Program::Parameters const& parameters,
                                                   Framebuffer* framebuffer_) {
	// Framebuffer must be non-null:
	assert(framebuffer_);
	auto& framebuffer = *framebuffer_;

	// A1T7: sample loop
	// TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	//  	 This will probably involve inserting a loop of the form:
	// 		 	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//      	for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//   	 around some subset of the code.
	// 		 You will also need to transform the input and output of the rasterize_* functions to
	// 	     account for the fact they deal with pixels centered at (0.5,0.5).

	std::vector<ShadedVertex> shaded_vertices;
	shaded_vertices.reserve(vertices.size());

	//--------------------------
	// shade vertices:
	for (auto const& v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	// assemble + clip + homogeneous divide vertices:
	std::vector<ClippedVertex> clipped_vertices;

	// reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		// clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		// clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}
	// clang-format off

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	// helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const& sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	// actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// rasterize primitives:

	std::vector<Fragment> fragments;

	// helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const& f) { fragments.emplace_back(f); };

	// actually do rasterization:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
			rasterize_line(clipped_vertices[i], clipped_vertices[i + 1], emit_fragment);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
			rasterize_triangle(clipped_vertices[i], clipped_vertices[i + 1], clipped_vertices[i + 2], emit_fragment);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// depth test + shade + blend fragments:
	uint32_t out_of_range = 0; // check if rasterization produced fragments outside framebuffer 
							   // (indicates something is wrong with clipping)
	for (auto const& f : fragments) {

		// fragment location (in pixels):
		int32_t x = (int32_t)std::floor(f.fb_position.x);
		int32_t y = (int32_t)std::floor(f.fb_position.y);

		// if clipping is working properly, this condition shouldn't be needed;
		// however, it prevents crashes while you are working on your clipping functions,
		// so we suggest leaving it in place:
		if (x < 0 || (uint32_t)x >= framebuffer.width || 
		    y < 0 || (uint32_t)y >= framebuffer.height) {
			++out_of_range;
			continue;
		}

		// local names that refer to destination sample in framebuffer:
		float& fb_depth = framebuffer.depth_at(x, y, 0);
		Spectrum& fb_color = framebuffer.color_at(x, y, 0);


		// depth test:
		if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
			// "Always" means the depth test always passes.
			//fb_depth = f.fb_position.z;
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
			// "Never" means the depth test never passes.
			continue; //discard this fragment
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
			fb_depth = f.fb_position.z;

			// "Less" means the depth test passes when the new fragment has depth less than the stored depth.
			// A1T4: Depth_Less
			// TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less".
		} else {
			static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
		}

		// if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
		if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
			fb_depth = f.fb_position.z;
		}

		// shade fragment:
		ShadedFragment sf;
		sf.fb_position = f.fb_position;
		Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

		// write color to framebuffer if color writes aren't disabled:
		if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {
			// blend fragment:
			if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
				fb_color = sf.color;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
				// A1T4: Blend_Add
				// TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
				fb_color = sf.color + fb_color; //<-- replace this line
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
				// A1T4: Blend_Over
				// TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
				// 		 You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
				fb_color = sf.color + fb_color * (1 - sf.opacity);//<-- replace this line
			} else {
				static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
			}
		}
	}
	if (out_of_range > 0) {
		if constexpr (primitive_type == PrimitiveType::Lines) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely "
			     "wrong with the clip_line function.",
			     out_of_range);
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely "
			     "wrong with the clip_triangle function.",
			     out_of_range);
		}
	}
}

// -------------------------------------------------------------------------
// clipping functions

// helper to interpolate between vertices:
template<PrimitiveType p, class P, uint32_t F>
auto Pipeline<p, P, F>::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  	va, vb: endpoints of line
 *  	emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_line(ShadedVertex const& va, ShadedVertex const& vb,
                                      std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// Determine portion of line over which:
	// 		pt = (b-a) * t + a
	//  	-pt.w <= pt.x <= pt.w
	//  	-pt.w <= pt.y <= pt.w
	//  	-pt.w <= pt.z <= pt.w
	// ... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		// restrict range such that:
		// l + t * dl <= r + t * dr
		// re-arranging:
		//  l - r <= t * (dr - dl)
		if (dr == dl) {
			// want: l - r <= 0
			if (l - r > 0.0f) {
				// works for none of range, so make range empty:
				min_t = 1.0f;
				max_t = 0.0f;
			}
		} else if (dr > dl) {
			// since dr - dl is positive:
			// want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { // dr < dl
			// since dr - dl is negative:
			// want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	// local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
	}
}

/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  	va, vb, vc: vertices of triangle
 *  	emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// A1EC: clip_triangle
	// TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

// -------------------------------------------------------------------------
// rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *   (x,y+0.5) -  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result. 
 *    We will not explicitly check for 45 degree lines along the diamond edges (this will be extra credit),
 *    but you should be able to handle 45 degree lines in every other case (such as starting from pixel centers)
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function<void(Fragment const&)> const& emit_fragment) {
		//inputs include emit_fragment function, clippedvertex va & vb where va and vb are endpoints
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	// A1T2: rasterize_line

	// TODO: Check out the block comment above this function for more information on how to fill in
	// this function!
	// The OpenGL specification section 3.5 may also come in handy.


	//goal: draw the entire line from va to vb
	//first: pick a major axis

	/* --------  VARIABLES -------*/
	bool majorAxisX = true; //major axis variable
 
	float vaX = va.fb_position.x; //va.x
	float vaY = va.fb_position.y; //va.y
	float vaZ = va.fb_position.z; //va.z

	float vbX = vb.fb_position.x; //vb.x
	float vbY = vb.fb_position.y; //vb.y
	float vbZ = vb.fb_position.z; //vb.z

	float vaXO = va.fb_position.x; //va.x original
	float vaYO = va.fb_position.y; //va.y original
	float vbXO = vb.fb_position.x; //vb.x original
	float vbYO = vb.fb_position.y; //vb.y original

	float s1 = (vbY - vaY) / (vbX - vaX);  //s1 of the line between x and y
	float s2 = (vbY - vaY) / (vbZ - vaZ); //s2 of the line between y and z
	float s3 =  (vbX - vaX) / (vbZ - vaZ); //s3 of the line between x and z

	float b1 = vaY - s1 * vaX; //intercept of the line for x and y
	float b2 = vaY - s2 * vaZ; //intercept of the line for y and z
	float b3 = vaX - s3 * vaZ; //intercept of the line for x and z

	double w; // for calculations later
	double v; //for calculations later
	Fragment shade; //the shade fragment

	/* --------  VARIABLES -------*/



	/* ----------- HELPER FUNCTIONS --------------*/

   auto shadePixel = [&] (float x, float y, float z)
    {
		shade.fb_position.x = floor(x) + 0.5f; //the shaded pixel x position
        shade.fb_position.y = floor(y) + 0.5f; //the shaded pixel y position
		shade.fb_position.z = floor(z) + 0.5f; //the shaded pixel z position
		shade.attributes = va.attributes; //take the attributes of va 
		shade.derivatives.fill(Vec2(0.0f, 0.0f)); //and do derivatives fill?
		emit_fragment(shade); //and call the function
    };

	auto endDiamondPoint = [&] (float x, float y, float z, bool majorAxis, float i, float j)
	{
		float xPcoor = x - floor(x); //calculate the pixel x corr from 0 - 1
		float yPcoor = y - floor(y); //calculate the pixel y corr from 0 - 1


		if(abs(x - i) <= 1.0f && abs(y - j) <= 1.0f) //if the other two coordinates are within the pixel, check whether they are in the other quadrants
		{	return; }

		
			if(xPcoor == 0.5f && yPcoor == 0.0f) //Case 1: bottom diamond point
			{ return; } //shade the pixel 
			else if (xPcoor == 0.0f && yPcoor == 0.5f) //Case 2:left diamond point
			{ return; } //shade the pixel 
			else if(xPcoor <= 0.5f && yPcoor <= 0.5f)//Case 3: if it is in the 3rd quadrant
			{ 
				if(yPcoor >= (-xPcoor + 0.5f)) // if it is inside the diamond
				{ return; } //shade pixel
				if(majorAxis)
				{
					if(vbX < vaX) //if the endpoint is less than the startpoint
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY < vaY) //if the edpoint is less than the startpoints
					{shadePixel(x, y, z); } //shade pixel
				}
				
			} //shade the pixel  
			else if(xPcoor <= 0.5f && yPcoor >= 0.5f) //Case 4: if it is in the 2nd quadrant
			{
				if(yPcoor <= (xPcoor + 0.5f)) // if it is inside the diamond
				{ return; } //shade the pixel  
				if(majorAxis)
				{
					if(vbX < vaX)
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY > vaY)
					{shadePixel(x, y, z); } //shade pixel
				}

			}
			else if(xPcoor >= 0.5f && yPcoor >= 0.5f) //Case 5: if it is in the 1st quadrant
			{
				if(yPcoor <= (-xPcoor + 1.5f)) // if it is inside the diamond
				{ return; } //do not shade the pixel  
				if(majorAxis)
				{
					if(vbX > vaX)
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY > vaY)
					{shadePixel(x, y, z); } //shade pixel
				}
			}
			else if(xPcoor >= 0.5f && yPcoor <= 0.5f) //Case 6: if it is in the 4th quadrant
			{
				if(yPcoor >= (xPcoor - 0.5f)) // if it is inside the diamond
				{ return; } //do not shade the pixel  
				if(majorAxis)
				{ 
					//check whether the line enters the diamond based on the slope & intercept of the line
					if(vbX > vaX)
					{shadePixel(x, y, z); } //shade pixel
				}  
				else
				{
					if(vbY < vaY)
					{shadePixel(x, y, z); } //shade pixel
				}
			}
	};


   auto beginDiamondPoint = [&] (float x, float y, float z, bool majorAxis, float i, float j)
    {
		float xPcoor = x - floor(x); //calculate the pixel x corr from 0 - 1
		float yPcoor = y - floor(y); //calculate the pixel y corr from 0 - 1

		float iPcoor = i - floor(i); //calculate the pixel x corr from 0 - 1
		float jPcoor = j - floor(j); //calculate the pixel y corr from 0 - 1

									//check if the other two coordinates are within the pixel 
		if(((abs(x - i) <= 1.0f) && (floor(x) == floor(i))) && ((abs(y - j) <= 1.0f) && (floor(y) == floor(j)))) //if the other two coordinates are within the pixel, check whether they are in the other quadrants
		{		
			if(xPcoor <= 0.5f && yPcoor <= 0.5f)//Case 3: if it is in the 3rd quadrant
			{ 
				if(iPcoor <= 0.5f && jPcoor <= 0.5f) //if the other coordinates are in the same quadrant
				{
					if(jPcoor >= (-iPcoor + 0.5f)) //if endpoint is inside diamond
					{return;}

					if(yPcoor >= (-xPcoor + 0.5f)) // if it is inside the diamond
					{  shadePixel(x, y, z); }
				}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor <= (-iPcoor + 1.5f)) && (iPcoor >= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor >= (iPcoor - 0.5f)) && (iPcoor >= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else
				{shadePixel(x, y, z);} //shade pixel
			} //shade the pixel  
			else if(xPcoor <= 0.5f && yPcoor >= 0.5f) //Case 4: if it is in the 2nd quadrant
			{

				if(iPcoor <= 0.5f && jPcoor >= 0.5f) //if the other coordinates are in the same quadrant
				{
					if(jPcoor <= (iPcoor + 0.5f)) //if endpoint is inside diamond
					{return;}

					if(yPcoor <= (xPcoor + 0.5f)) // if it is inside the diamond
					{  shadePixel(x, y, z); }
					
				}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor <= (-iPcoor + 1.5f)) && (iPcoor >= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor >= (iPcoor - 0.5f)) && (iPcoor >= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else
				{shadePixel(x, y, z);} //shade pixel
			}
			else if(xPcoor >= 0.5f && yPcoor >= 0.5f) //Case 5: if it is in the 1st quadrant
			{
				if(iPcoor >= 0.5f && jPcoor >= 0.5f) //if the other coordinates are in the same quadrant
				{
					if(jPcoor <= (-iPcoor + 1.5f)) //if endpoint is inside diamond
					{return;}

					if(yPcoor <= (-xPcoor + 1.5f)) // if it is inside the diamond
					{  shadePixel(x, y, z); }
				}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor >= (iPcoor - 0.5f)) && (iPcoor >= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else
				{shadePixel(x, y, z);} //shade pixel
			}
			else if(xPcoor >= 0.5f && yPcoor <= 0.5f) //Case 6: if it is in the 4th quadrant
			{
				if(iPcoor >= 0.5f && jPcoor <= 0.5f) //if the other coordinates are in the same quadrant
				{
					if(jPcoor >= (iPcoor - 0.5f)) //if endpoint is inside diamond
					{return;}

					if(yPcoor >= (xPcoor - 0.5f)) // if it is inside the diamond
					{  shadePixel(x, y, z); }
				}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor <= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor <= (iPcoor + 0.5f)) && (iPcoor <= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else if((jPcoor <= (-iPcoor + 1.5f)) && (iPcoor >= 0.5f && jPcoor >= 0.5f)) // if it is inside the diamond
				{return;}
				else
				{shadePixel(x, y, z);} //shade pixel
			}
		}
		else
		{
			if(xPcoor == 0.5f && yPcoor == 0.0f) //Case 1: bottom diamond point
			{ shadePixel(x, y, z); } //shade the pixel 
			else if (xPcoor == 0.0f && yPcoor == 0.5f) //Case 2:left diamond point
			{ shadePixel(x, y, z); } //shade the pixel 
			else if(xPcoor <= 0.5f && yPcoor <= 0.5f)//Case 3: if it is in the 3rd quadrant
			{ 
				if(yPcoor >= (-xPcoor + 0.5f)) // if it is inside the diamond
				{shadePixel(x, y, z); } //shade pixel
				//check whether the line enters the diamond based on the slope & intercept of the line
				if(majorAxis)
				{
					if(vbX > vaX) //if the endpoint is less than the startpoint
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY > vaY) //if the edpoint is less than the startpoints
					{shadePixel(x, y, z); } //shade pixel
				}
				
			} //shade the pixel  
			else if(xPcoor <= 0.5f && yPcoor >= 0.5f) //Case 4: if it is in the 2nd quadrant
			{
				if(yPcoor <= (xPcoor + 0.5f)) // if it is inside the diamond
				{ shadePixel(x, y, z); } //shade the pixel  
				if(majorAxis)
				{
					if(vbX > vaX)
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY < vaY)
					{shadePixel(x, y, z); } //shade pixel
				}

			}
			else if(xPcoor >= 0.5f && yPcoor >= 0.5f) //Case 5: if it is in the 1st quadrant
			{
				if(yPcoor <= (-xPcoor + 1.5f)) // if it is inside the diamond
				{ shadePixel(x, y, z); } //shade the pixel  
				if(majorAxis)
				{
					if(vbX < vaX)
					{shadePixel(x, y, z); } //shade pixel		
				} //shade pixel
				else
				{
					if(vbY < vaY)
					{shadePixel(x, y, z); } //shade pixel
				}
			}
			else if(xPcoor >= 0.5f && yPcoor <= 0.5f) //Case 6: if it is in the 4th quadrant
			{
				if(yPcoor >= (xPcoor - 0.5f)) // if it is inside the diamond
				{ shadePixel(x, y, z); } //shade the pixel  
				if(majorAxis)
				{ 
					//check whether the line enters the diamond based on the slope & intercept of the line
					if(vbX < vaX)
					{shadePixel(x, y, z); } //shade pixel
				}  
				else
				{
					if(vbY > vaY)
					{shadePixel(x, y, z); } //shade pixel
				}
			}
		}

    };


	/* ----------- HELPER FUNCTIONS --------------*/


	float deltaX = abs(vbX - vaX); // change in X
	float deltaY = abs(vbY - vaY); //change in Y

	if(deltaX <= deltaY){ majorAxisX = false; } // Y is the Major Axis
											  // else, X remains the Major Axis

	printf("%d ", majorAxisX);
	beginDiamondPoint(vaX, vaY, vaZ, majorAxisX, vbX, vbY); //figure out whether to shade the 1st pixel or not
	endDiamondPoint(vbX, vbY, vbZ, majorAxisX, vaX, vaY); //figure out whether to shade the last pixel or not

	if( majorAxisX == false )
	{ //if the major Axis was determined to be Y

		if(vaY > vbY) // check if va y is > vb y
		{ //if so, switch their positions!
			vaX = vb.fb_position.x; 
			vaY = vb.fb_position.y;
			vaZ = vb.fb_position.z;

			vbX = va.fb_position.x;
			vbY = va.fb_position.y;
			vbZ = va.fb_position.z;
			// this one's hard 
			vaXO = vb.fb_position.x; 
			vaYO = vb.fb_position.y;
			vbXO = va.fb_position.x;
			vbYO = va.fb_position.y;

			//s1 also needs to be updated!
			
			s1 = (vbY-vaY) / (vbX - vaX);  //s1 of the line between x and y
			s2 = (vbY - vaY) / (vbZ - vaZ); //s2 of the line between y and z

			b1 = vaY - s1 * vaX; //intercept of the line for x and y 
			b2 = vaY - s2 * vaZ; //intercept of the line for y and z
			
		}
		vaY = floor(vaY);
		vbY = floor(vbY);
		vaY += 1.0f; //increment to the next pixel Y
		vaZ = ((vaY - b2 ) / s2); // get the next pixel z
		//next, draw the line!

		while(vaY < vbY) //while we are still going from the start
																  //to the end of the line
		{
			w = (((double(vaY) + 0.5) - double(vaYO) ) / (double(vbYO) - double(vaYO))); 
			v = (w * (double(vbXO) - double(vaXO)) + double(vaXO));
			//calculate w and v with the equations given in the slides

			shadePixel(float(floor(v)), vaY, vaZ); // shade the pixel

			vaY += 1.0f; //increment to the next pixel Y
			vaZ = ((vaY - b2 ) / s2); // get the next pixel z
		}  

		//figure out if the endpoint should be shaded or not

	}
	else 
	{ //if not, the major Axis is X
		
		if(vaX > vbX) // check if va x is > vb x
		{ //if so, switch their positions!
			vaX = vb.fb_position.x; 
			vaY = vb.fb_position.y;
			vaZ = vb.fb_position.z;

			vbX = va.fb_position.x;
			vbY = va.fb_position.y;
			vbZ = va.fb_position.z;

			vaXO = vb.fb_position.x; 
			vaYO = vb.fb_position.y;
			vbXO = va.fb_position.x;
			vbYO = va.fb_position.y;
			//s1 also needs to be updated!
			s1 = (vbY-vaY) / (vbX - vaX);  //s1 of the line
			s3 =  (vbX - vaX) / (vbZ - vaZ); //s2 for z 
			
			b1 = vaY - s1 * vaX; //intercept of the line
			b3 = vaX - s3 * vaZ; //intercept of the line for z
		}
		vaX = floor(vaX);
		vbX = floor(vbX);
		vaX += 1.0f; //increment to the next pixel X
		vaZ = (s3 * vaX + b3); // get the next pixel z
		//next, draw the line!

		while(vaX < vbX) //while we are still going from the start
																  //to the end of the line
		{
			w = (((double(vaX) + 0.5) - double(vaXO)) / (double(vbXO) - double(vaXO))); 
			v = (w * (double(vbYO) - double(vaYO)) + double(vaYO));
			//calculate w and v with the equations given in the slides
			shadePixel(vaX, float(floor(v)), vaZ);

			vaX += 1.0f; //increment to the next pixel X
			vaZ = (s3 * vaX + b3); // get the next pixel z
		} 
		
	} 

}

/*
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  	(x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Smooth: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 * 	The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function<void(Fragment const&)> const& emit_fragment) {
	// NOTE: it is okay to restructure this function to allow these tasks to use the
	//  same code paths. Be aware, however, that all of them need to remain working!
	//  (e.g., if you break Flat while implementing Correct, you won't get points
	//   for Flat.)
	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {

		//first, make a square that is va-vb x va-vc in area

		//check the points in that area

		//if it is, shade the point

		//if not 
		/* --------  VARIABLES -------*/
 
	float vaX = va.fb_position.x; //va.x
	float vaY = va.fb_position.y; //va.y
	float vaZ = va.fb_position.z; //va.z

	float vbX = vb.fb_position.x; //vb.x
	float vbY = vb.fb_position.y; //vb.y
	float vbZ = vb.fb_position.z; //vb.z


	float vcX = vc.fb_position.x; //vc.x
	float vcY = vc.fb_position.y; //vc.y
	float vcZ = vc.fb_position.z; //vc.z

	float x1 = 0.0f; //horz values for translating bounding box
	float x2 = 0.0f;

	float y1 = 0.0f; //vert values for translating bounding box
	float y2 = 0.0f;

	Fragment shade; //the shade fragment

	/* --------  VARIABLES -------*/

		/* --------  HELPER FUNCTIONS -------*/

		auto dotProduct = [&](float x, float y) //calculate the dot product
		{ return float(x*y); };

		auto crossProduct = [&](float x, float y, float a, float b) //calculate the cross product
		{ return float(x*b - y*a); };
		auto shadePixel = [&](float x, float y, float z)
		{
			shade.fb_position.x = floor(x) + 0.5f; //the shaded pixel x position
			shade.fb_position.y = floor(y) + 0.5f; //the shaded pixel y position
			shade.fb_position.z = floor(z) + 0.5f; //the shaded pixel z position
			shade.attributes = va.attributes; //take the attributes of va 
			shade.derivatives.fill(Vec2(0.0f, 0.0f)); //and do derivatives fill?
			emit_fragment(shade); //and call the function
		};

		auto topRule = [&] (float Ay, float By, float Cy, float x, float y)
		{
			if((Ay == By) && (Ay > Cy))
			{shadePixel(x, y, (vaZ + vbZ + vcZ)/ 3.0f);}
		};

		auto leftRule = [&] (float Ay, float Bx, float By, float Cx, float x, float y)
		{
			if((Ay < By) && (Bx <= Cx))
			{shadePixel(x, y, (vaZ + vbZ + vcZ)/ 3.0f);}

		};

		auto pointInside = [&](float Ax, float Ay, float Bx, float By, float Cx, float Cy, float x, float y) //check if a point is inside the triangle
		{
			float check1 = dotProduct(crossProduct(Cx- Ax, Cy - Ay, Bx - Ax, By - Ay ) , crossProduct(Cx - Ax, Cy - Ay, x - Ax, y - Ay)); //ac ab and aq
			float check2 = dotProduct(crossProduct(Bx- Cx, By - Cy, Ax - Cx, Ay - Cy ) , crossProduct(Bx - Cx, By - Cy, x - Cx, y - Cy));//cb ca and cq
			float check3 = dotProduct(crossProduct(Ax- Bx, Ay - By, Cx - Bx, Cy - By ) , crossProduct(Ax - Bx, Ay - By, x - Bx, y - By));//cb ca and cq
			if((check1 >= 0) && (check2 >= 0) && (check3 >= 0))
			{	
				if((check1 == 0))
				{
					leftRule(Ay, Cx, Cy, Bx, x, y);
					topRule(Ay, Cy, By, x, y);
				}
				else if(check2 == 0)
				{
					leftRule(Cy, Bx, By, Ax, x, y);
					topRule(Cy, By, Ay, x, y);
				}
				else if(check3 == 0)
				{
					leftRule(By, Ax, Ay, Cx, x, y);
					topRule(By, Ay, Cy, x, y);
				}
				else
				{shadePixel(x, y, (vaZ + vbZ + vcZ)/ 3.0f);}
			}
		};		

		auto longestSide = [&] (float Ax, float Bx, float Cx) //find the side with the largest delta x
		{
			float a = abs(Ax - Bx);
			float b = abs(Bx - Cx);
			float c = abs(Ax - Cx);

			if(a >= c && a >= b)
			{
				if(Ax < Bx)
				{
					x1 = Ax;
					x2 = Bx;
				}
				else
				{
					x1 = Bx;
					x2 = Ax;
				}
			}
			else if(b >= c && b >= a)
			{
				if(Bx < Cx)
				{
					x1 = Bx;
					x2 = Cx;
				}
				else
				{
					x1 = Cx;
					x2 = Bx;
				}

			}
			else 
			{
				if(Ax < Cx)
				{
					x1 = Ax;
					x2 = Cx;
				}
				else
				{
					x1 = Cx;
					x2 = Ax;
				}
			}
			
		};

		auto tallestSide = [&] (float Ay, float By, float Cy) //find the side with the largest delta y
		{
			float a = abs(Ay - By);
			float b = abs(By - Cy);
			float c = abs(Ay - Cy);

			if(a >= c && a >= b)
			{
				if(Ay < By)
				{
					y1 = Ay;
					y2 = By;
				}
				else
				{
					y1 = By;
					y2 = Ay;
				}
			}
			else if(b >= c && b >= a)
			{
				if(By < Cy)
				{
					y1 = By;
					y2 = Cy;
				}
				else
				{
					y1 = Cy;
					y2 = By;
				}
			}
			else 
			{
				if(Cy < Ay)
				{
					y1 = Cy;
					y2 = Ay;
				}
				else
				{
					y1 = Ay;
					y2 = Cy;
				}
			}
		};

		auto orientation = [&] (float Ax, float Ay, float Bx, float By, float Cx, float Cy) //returns true for clockwise false for counter clockwise
		{
			if(Ax < Bx)
			{
				if(By < Cy)
				{return false;} //counter clockwise
				else
				{ return true;} //clockwise
			}
			else
			{
				if(By >= Cy)
				{return false;} //counter clockwise
				else
				{ return true;} //clockwise
			}

		};
		/* --------  HELPER FUNCTIONS -------*/

	longestSide(vaX, vbX, vcX); //get the longest side
	tallestSide(vaY, vbY, vcY); //get the tallest side

	//find longest and tallest side of triangle, make a box with the x and y variables from those sides
	for(float i = floor(x1); i <= floor(x2); i+= 1.0f) //go through from x1 -> x2
	{
		for(float j = floor(y1); j <= floor(y2); j+= 1.0f) //go through from y1 -> y2
		{
			if(orientation(vaX, vaY, vbX, vbY, vcX, vcY)) //if the triangle is going clockwise
			{pointInside(vaX, vaY, vcX, vcY, vbX, vbY, floor(i) + 0.5f, floor(j) + 0.5f);}
			else //if the triangle is going counter clockwise
			{pointInside(vaX, vaY, vbX, vbY, vcX, vcY, floor(i) + 0.5f, floor(j) + 0.5f);}
		}
	}



		// A1T3: flat triangles
		// TODO: rasterize triangle (see block comment above this function).

		// As a placeholder, here's code that draws some lines:
		//(remove this and replace it with a real solution)
		/*Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(va, vb, emit_fragment);
		Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(vb, vc, emit_fragment);
		Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(vc, va, emit_fragment); */
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Smooth) {
		// A1T5: screen-space smooth triangles
		// TODO: rasterize triangle (see block comment above this function).

		// As a placeholder, here's code that calls the Flat interpolation version of the function:
		//(remove this and replace it with a real solution)
		Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Flat>::rasterize_triangle(va, vb, vc, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		// A1T5: perspective correct triangles
		// TODO: rasterize triangle (block comment above this function).

		// As a placeholder, here's code that calls the Screen-space interpolation function:
		//(remove this and replace it with a real solution)
		Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Smooth>::rasterize_triangle(va, vb, vc, emit_fragment);
	}
}

//-------------------------------------------------------------------------
// compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
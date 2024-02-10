#include "test.h"
#include "rasterizer/pipeline.h"
#include "rasterizer/programs.h"

#include <limits>
#include <iomanip>
#include <algorithm>
#include <unordered_set>

//Failing the 0.999999 test case? Use doubles over floats -Nico

using TestPipeline = Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;

namespace std {
	template< >
	struct hash< Vec2 > {
		size_t operator()(const Vec2 &v) const {
			static hash< float > hf;
			size_t x = hf(v.x);
			size_t y = hf(v.y);
			return x ^ (y << 16) ^ (y >> (sizeof(y)*8-16));
		}
	};
}

//check that line produces exactly the listed fragments:
void check_line_covers(std::string const &desc, std::vector< Vec2 > const &line_strip, std::unordered_set< Vec2 > const &expected) {

	std::unordered_set< Vec2 > got;
	for (uint32_t i = 0; i + 1 < line_strip.size(); ++i) {
		TestPipeline::ClippedVertex a,b;
		a.fb_position = Vec3(line_strip[i].x, line_strip[i].y, 0.25f);
		a.inv_w = 1.0f;
		a.attributes.fill(1.0f);
		b.fb_position = Vec3(line_strip[i+1].x, line_strip[i+1].y, 0.75f);
		b.inv_w = 1.0f;
		b.attributes.fill(2.0f);
		TestPipeline::rasterize_line(a, b, [&](TestPipeline::Fragment const &frag){
			got.emplace(frag.fb_position.x, frag.fb_position.y);
		});
	}

	std::vector< std::string > raster;

	raster.emplace_back(".");

	uint32_t out_of_raster = 0;

	auto draw = [&raster,&out_of_raster](Vec2 const &px, char c) {
		int32_t x = int32_t(std::floor(px.x));
		int32_t y = int32_t(std::floor(px.y));

		if (x < 0 || y < 0 || x > 10 || y > 10) {
			++out_of_raster;
			return;
		}

		if (uint32_t(y) >= raster.size()) {
			raster.resize(y+1, "");
		}
		if (uint32_t(x) >= raster[y].size()) {
			raster[y].resize(x+1, '.');
		}
		raster[y][x] = c;
	};

	uint32_t matched = 0;
	uint32_t missed = 0;
	uint32_t extra = 0;

	for (auto const &f : got) {
		if ((f.x - std::floor(f.x) != 0.5f) || (f.y - std::floor(f.y) != 0.5f)) {
			throw Test::error("Rasterizing '" + desc + "', got fragment at (" + std::to_string(f.x) + ", " + std::to_string(f.y) + "), which isn't at a pixel center.");
		}
		if (expected.count(f)) {
			draw(f, '#');
			++matched;
		} else {
			draw(f, '!');
			++extra;
		}
	}
	for (auto const &f : expected) {
		if (!got.count(f)) {
			draw(f, '?');
			++missed;
		}
	}

	if (extra > 0 || missed > 0) {
		//failed!
		std::string info = "Example '" + desc + "' missed " + std::to_string(missed) + " ('?'); had " + std::to_string(extra) + " extra ('!'); and matched " + std::to_string(matched) + " ('#') fragments:";

		//square up the raster:
		size_t width = 0;
		for (auto const &row : raster) {
			width = std::max(width, row.size());
		}
		for (auto &row : raster) {
			row.resize(width, '.');
		}

		for (uint32_t y = static_cast<uint32_t>(raster.size()) - 1; y < static_cast<uint32_t>(raster.size()); --y) {
			info += "\n    " + raster[y];
		}

		if (out_of_raster) info += "\n    (" + std::to_string(out_of_raster) + " out-of-range fragments not plotted.)";

		puts(""); //because "test..."
		info("%s", info.c_str());

		throw Test::error("Example '" + desc + "' didn't match expected.");
	}

	//if nothing extra and nothing missed, success!
	assert(matched == expected.size());
}

//check that line produces exactly the fragments drawn in a fancy picture:
void check_line_covers(std::string const &desc, std::initializer_list< Vec2 > const &line_strip, std::initializer_list< std::string > const &raster_) {
	//convert raster to set of points ( with lower-left being (0,0) ):
	std::vector< std::string > raster(raster_);
	std::unordered_set< Vec2 > expected;
	for (uint32_t y = 0; y < raster.size(); ++y) {
		std::string const &row = raster[raster.size()-1-y];
		for (uint32_t x = 0; x < row.size(); ++x) {
			if (row[x] != '.') {
				expected.emplace(x + 0.5f, y + 0.5f);
			}
		}
	}
	//use list-of-points version:
	check_line_covers(desc, line_strip, expected);
}

//--------------------------------------------------
//entering/exiting diamond at (1,1):
// only lines that *exit* the diamond should produce a fragment.


Test test_a1_task2_diamond_inside("a1.task2.diamond.inside", []() {
	check_line_covers(
		"line inside diamond (1,1)",
		{ Vec2(1.5f, 1.25f), Vec2(1.25f, 1.5f) },
		{"...",
		 "...",
		 "..."}
	);
});


Test test_a1_task2_diamond_outside("a1.task2.diamond.outside", []() {
	check_line_covers(
		"line outside diamond (1,1)",
		{ Vec2(1.125f, 1.25f), Vec2(1.25f, 1.125f) },
		{"...",
		 "...",
		 "..."}
	);
});


//----------------------------
//simple horizontal and vertical lines (set up so that no enter/exit logic needed):

Test test_a1_task2_simple_horizontal("a1.task2.simple.horizontal", []() {
	check_line_covers(
		"horizontal line from (1.125, 1.125) to (4.875, 1.125)",
		{ Vec2(1.125f, 1.125f), Vec2(4.875f, 1.125f) },
		{"......",
		 ".####.",
		 "......"}
	);
});


Test test_a1_task2_simple_vertical("a1.task2.simple.vertical", []() {
	check_line_covers(
		"vertical line from (1.125, 1.125) to (1.125, 4.875)",
		{ Vec2(1.125f, 1.125f), Vec2(1.125f, 4.875f) },
		{"...",
		 ".#.",
		 ".#.",
		 ".#.",
		 ".#.",
		 "..."}
	);
});

Test test_a1_task2_single_valid_pos("a1.task2.single.valid.positive", []() {
 check_line_covers(
  "single pixel emitted with positive slope",
  { Vec2(1.1f, 1.15f), Vec2(1.9f, 1.2f) },
  {"...",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_single_valid_neg("a1.task2.single.valid.negative", []() {
 check_line_covers(
  "single pixel emitted with negative slope",
  { Vec2(1.125f, 1.675f), Vec2(1.875, 1.125f) },
  {"...",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_single_valid_middle("a1.task2.single.valid.middle", []() {
 check_line_covers(
  "single pixel emitted but start point is from middle",
  { Vec2(1.5f, 1.5f), Vec2(1.5f, 2.0f) },
  {"...",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_single_valid_endright("a1.task2.single.valid.endright", []() {
 check_line_covers(
  "single pixel emitted but end point is right corner of diamond",
  { Vec2(1.125f, 1.875f), Vec2(2.0f, 1.5f) },
  {"...",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_single_invalid_left("a1.task2.single.invalid.left", []() {
 check_line_covers(
  "end point on left corner = line doesn't 'exit' diamond; also need to swap points for bresenham alg",
  { Vec2(1.875f, 1.875f), Vec2(1.0f, 1.5f) },
  {"...",
   "..."}
 );
});
// edge cases
Test test_a1_task2_edge_bottom("a1.task2.edge.bottom", []() {
 check_line_covers(
  "just crossing bottom corner = enter and exit diamond",
  { Vec2(1.0f, 1.0f), Vec2(2.0f, 1.0f) },
  {"...",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_edge_left("a1.task2.edge.left", []() {
 check_line_covers(
  "just crossing left corner = enter and exit diamond",
  { Vec2(1.0f, 1.0f), Vec2(1.0f, 3.0f) },
  {"...",
   ".#.",
   ".#.",
   "..."}
 );
});
Test test_a1_task2_edge_invalid("a1.task2.edge.invalid", []() {
 check_line_covers(
  "line crosses 2 pixels but emits none",
  { Vec2(1.1f, 1.2f), Vec2(1.25f, 0.8f) },
  {"...",
   "..."}
 );
});

// breadth cases
Test test_a1_task2_breadth_pos_1("a1.task2.breadth.pos.1", []() {
 check_line_covers(
  "line with slope m<=1",
  { Vec2(1.25f, 1.875f), Vec2(6.0f, 6.625f) },
  {".....#",
   "....#.",
   "...#..",
   "..#...",
   ".#....",
   "......",
   "......"}
 );
});
Test test_a1_task2_breadth_pos_small("a1.task2.breadth.pos.small", []() {
 check_line_covers(
  "line with slope 0<=m<=1",
  { Vec2(0.5f, 0.5f), Vec2(8.0f, 2.0f) },
  {"...#####",
   "###....."}
 );
});
Test test_a1_task2_breadth_pos_large("a1.task2.breadth.pos.large", []() {
 check_line_covers(
  "line with slope m>1 and swap",
  { Vec2(2.0f, 4.5f), Vec2(0.1f, 0.25f) },
  {"...",
   "..#",
   ".#.",
   ".#.",
   "#..",
   "#.."}
 );
});
Test test_a1_task2_breadth_neg_small("a1.task2.breadth.neg.small", []() {
 check_line_covers(
  "line with slope -1 <= m <= 0",
  { Vec2(0.5f, 3.5f), Vec2(6.5f, 0.5f) },
  {"##....",
   "..##..",
   "....##",
   "......"}
 );
});
Test test_a1_task2_breadth_neg_large("a1.task2.breadth.neg.large", []() {
 check_line_covers(
  "line with slope -1<m",
  { Vec2(1.1f, 5.9f), Vec2(2.75f, 2.1f) },
  {".#.",
   ".#.",
   "..#",
   "..#",
   "...",
   "..."}
 );
});

// More edge cases based on order because it matters apparently

Test test_a1_task2_order_valid1("a1.task2.order.valid1", []() {
 check_line_covers(
  "direct opposite of a1.task2.single.invalid.left, but DOES emit pixel",
  { Vec2(1.0f, 1.5f), Vec2(1.875f, 1.875f) },
  {".#.",
   "..."}
 );
});

Test test_a1_task2_order_valid2("a1.task2.order.valid2", []() {
 check_line_covers(
  "similar to a1.task2.order.valid1 but starting point is on the right half to make sure you don't have weird heuristics",
  { Vec2(1.6f, 1.45f), Vec2(1.875f, 1.875f) },
  {".#.",
   "..."}
 );
});

Test test_a1_task2_order_valid_incl_incl("a1.task2.order.valid.inclusive-inclusive", []() {
 check_line_covers(
  "inclusive start and end point - line exits the start diamond (even though it starts inside) and end diamond",
  { Vec2(0.6f, 0.5f), Vec2(4.875f, 2.875f) },
  {"...##",
   "..#..",
   "##..."} 
 );
});

Test test_a1_task2_order_valid_incl_excl("a1.task2.order.valid.inclusive-exclusive", []() {
 check_line_covers(
  "swapped endpts of above - inclusive start and exclusive end point - line exits the start diamond but doesn't exit the end diamond",
  { Vec2(4.875f, 2.875f), Vec2(0.6f, 0.5f) },
  {"...##",
   "..#..",
   ".#..."}
 );
});

// start and end at corners (no 45 deg tests though)
// a1.task2.breadth.pos.large already implicitly checks for starting at left corner
Test test_a1_task2_corner_left_end("a1.task2.corner.left.end", []() {
 check_line_covers(
  "ending at left corner = no emit",
  { Vec2(0.0f, 1.0f), Vec2(0.0f, 0.5f) },
  {"..."}
 );
});

Test test_a1_task2_corner_top_start("a1.task2.corner.top.start", []() {
 check_line_covers(
  "starting 'outside' = no emit where you think it should but it does emit",
  { Vec2(0.5f, 1.0f), Vec2(1.0f, 1.0f) },
  {"#..",
   "..."}
 );
});
// ending at top corner depends on start case sooo I'll exclude this case for now
Test test_a1_task2_corner_right_start("a1.task2.corner.right.start", []() {
 check_line_covers(
  "starting 'outside' = no emit where you think it should but it does emit",
  { Vec2(1.0f, 0.5f), Vec2(1.0f, 1.0f) },
  {".#."}
 );
});
// ending at top corner depends on start case sooo I'll exclude this case for now
Test test_a1_task2_corner_bottom_start("a1.task2.corner.bottom.start", []() {
 check_line_covers(
  "starting at bottom corner => emit because it immediately leaves (depending on endpt tho)",
  { Vec2(0.5f, 0.0f), Vec2(0.0f, 0.0f) },
  {"#.."}
 );
});
Test test_a1_task2_corner_bottom_start_invalid("a1.task2.corner.bottom.start.invalid", []() {
 check_line_covers(
  "starting at bottom corner => no emit because endpt is within diamond",
  { Vec2(0.5f, 0.0f), Vec2(0.3f, 0.6f) },
  {"..."}
 );
});
Test test_a1_task2_corner_bottom_end("a1.task2.corner.bottom.end", []() {
 check_line_covers(
  "ending at bottom corner = no emit",
  { Vec2(0.0f, 1.0f), Vec2(0.5f, 0.0f) },
  {"..."}
 );
}); 

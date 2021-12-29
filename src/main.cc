#include <chrono>          // timing utilities
#include <iostream>       // text i/o
#include <iomanip>
#include <stdio.h>      // printf if needed
#include <vector>      // std::vector
#include <random>     // prng
#include <string>    // std::string
#include <sstream>    // std::stringstream
#include <algorithm> // clamp
#include <atomic>   // atomic_llong
#include <thread>  // threads
#include <memory> // shared_ptr

using std::cerr, std::cin, std::cout, std::endl, std::flush;
using std::chrono::high_resolution_clock, std::chrono::duration_cast, std::chrono::milliseconds;

// my vector library
#include "lib/vector.h"

// image input
#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"

// image output
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

// default types
#define baseType double
using vec2 = vector2< baseType >;
using vec3 = vector3< baseType >;

// render parameters
constexpr long long X_IMAGE_DIM = 1920;
constexpr long long Y_IMAGE_DIM = 1080;
constexpr long long TILESIZE_XY = 8;
constexpr long long MAX_BOUNCES = 10;
constexpr long long NUM_SAMPLES = 32;
constexpr long long NUM_THREADS = 5;
constexpr baseType  IMAGE_GAMMA = 2.2;
constexpr baseType  HIT_EPSILON = baseType( std::numeric_limits< baseType >::epsilon() );
constexpr baseType  DMAX_TRAVEL = baseType( std::numeric_limits< baseType >::max() ) / 10.0;
constexpr baseType  FIELD_OF_VIEW = 0.35;
constexpr baseType  PALETTE_SCALAR = 16.18;
constexpr baseType  BRIGHTNESS_SCALAR = 16.18;
constexpr long long REPORT_DELAY = 16; 				// reporter thread sleep duration, in ms
constexpr long long NUM_PRIMITIVES = 69;
constexpr long long PROGRESS_INDICATOR_STOPS = 69; // cli spaces to take up

// ray representation (origin+direction)
struct ray {
	vec3 origin;
	vec3 direction;
};

// represents a ray hit and the associated information
struct hitrecord {									// hit record
	vec3 position;										// position
	vec3 normal;											// normal
	baseType dtransit = DMAX_TRAVEL;	// how far the ray traveled, initially very large
	int materialID = -1;							// material (indexed into scene list)
	int primitiveID = 0;							// so you can refer to the values for this primitive later
	vec2 uv;													// used for triangles, barycentric coords
	bool front;												// hit on frontfacing side
};

// inline uint32_t wang_hash( uint32_t x ){
// 	x = ( x ^ 12345391 ) * 2654435769;
// 	x ^= ( x << 6 ) ^ ( x >> 26 );
// 	x *= 2654435769;
// 	x += ( x << 5 ) ^ ( x >> 12 );
// 	return x;
// }

// Random Utilities
baseType rng( std::shared_ptr< std::mt19937_64 > gen ) { // gives a value in the range 0.-1.
	std::uniform_real_distribution< baseType > distribution( 0., 1. );
	return distribution( *gen );
}
vec3 randomVector( std::shared_ptr< std::mt19937_64 > gen ) { // random vector centered around 0.
	return vec3( rng( gen ), rng( gen ), rng( gen ) ) - vec3( 0.5 );
}
vec3 randomUnitVector( std::shared_ptr< std::mt19937_64 > gen ) { // random direction vector (unit length)
	baseType z = rng( gen ) * 2.0 - 1.0;
	baseType a = rng( gen ) * 2.0 * pi;
	baseType r = sqrt( 1.0 - z * z );
	baseType x = r * cos( a );
	baseType y = r * sin( a );
	return vec3( x, y, z );
}
vec3 random_in_unit_disk( std::shared_ptr< std::mt19937_64 > gen ) { // random in unit disk (xy plane)
	vec3 val = randomUnitVector( gen );
	return vec3( val.values[ 0 ], val.values[ 1 ], 0.0 );
}

// iq style palette
vec3 palette( baseType t,
	vec3 a = vec3( 0.50, 0.50, 0.50 ),
	vec3 b = vec3( 0.50, 0.50, 0.50 ),
	vec3 c = vec3( 1.00, 1.00, 1.00 ),
	vec3 d = vec3( 0.00, 0.33, 0.67 ) ) {
	vec3 temp = ( c * t + d ) * 2.0 * pi;
	return a + b * vec3( cos( temp.values[ 0 ] ), cos( temp.values[ 1 ] ), cos( temp.values[ 2 ] ) );
}

class primitive { // base class for primitives
public:
	virtual hitrecord intersect( ray r ) const = 0; // pure virtual, base definition dne
	int materialID; // indexes into scene material list
};

// sphere
class sphere : public primitive {
public:
	sphere( vec3 c, baseType r, int m ) : center( c ), radius( r ) { materialID = m; }
	hitrecord intersect( ray r ) const override {
		hitrecord h; h.dtransit = DMAX_TRAVEL;
		vec3 disp = r.origin - center;
		baseType b = dot( r.direction, disp );
		baseType c = dot( disp, disp ) - radius * radius;
		baseType des = b * b - c; // b squared minus c - discriminant of the quadratic
		if( des >= 0.0 ){ // hit at either one or two points
			baseType d = std::min( std::max( -b + std::sqrt( des ), 0.0 ), std::max( -b - std::sqrt( des ), 0.0 ) );
			if( d >= 0.0 ){ // make sure at least one intersection point is in front of the camera before continuing
				h.dtransit = d;
				h.materialID = materialID;
				h.position = r.origin + h.dtransit * r.direction;
				h.normal = normalize( h.position - center );
				h.front = dot( h.normal, r.direction ) < 0.0 ? true : false;
			}
		}
		return h;
	}
	private:  // geometry parameters
		vec3 center;
		baseType radius;
};

// triangle
class triangle : public primitive {
public:
	triangle( vec3 p0, vec3 p1, vec3 p2, int m ) : points{ p0, p1, p2 } { materialID = m; }
	hitrecord intersect ( ray r ) const override { // Möller–Trumbore intersection algorithm
		hitrecord hit;  hit.dtransit = DMAX_TRAVEL;
		const vec3 edge1 = points[ 1 ] - points[ 0 ];
		const vec3 edge2 = points[ 2 ] - points[ 0 ];
		const vec3 pvec = cross( r.direction, edge2 );
		const baseType determinant = dot( edge1, pvec );
		if ( determinant > -HIT_EPSILON && determinant < HIT_EPSILON )
			return hit; // no hit, return

		const baseType inverseDeterminant = 1.0f / determinant;
		const vec3 tvec = r.origin - points[ 0 ];
		hit.uv.values[ 0 ] = dot( tvec, pvec ) * inverseDeterminant; // u value
		if ( hit.uv.values[ 0 ] < 0.0f || hit.uv.values[ 0 ] > 1.0f )
			return hit; // no hit, return

		const vec3 qvec = cross( tvec, edge1 );
		hit.uv.values[ 1 ] = dot( r.direction, qvec ) * inverseDeterminant; // v value
		if ( hit.uv.values[ 1 ] < 0.0f || hit.uv.values[ 0 ] + hit.uv.values[ 1 ] > 1.0f )
			return hit; // no hit, return

		hit.dtransit = dot( edge2, qvec ) * inverseDeterminant; // distance term to hit
		hit.position = points[ 0 ] + hit.uv.values[ 0 ] * edge2 + hit.uv.values[ 1 ] * edge1;
		hit.normal = cross( edge1, edge2 );
		hit.materialID = materialID;
		hit.front = dot( hit.normal, r.direction ) < 0.0 ? true : false; // determine front or back

		return hit; // return true result with all relevant info
	}
private:  // geometry parameters
	vec3 points[ 3 ];
};


// axis aligned bounding box
class aabb : public primitive {
public:
	aabb( vec3 mins, vec3 maxs, int m ) : bounds{ mins, maxs } { materialID = m; }

	//	Alexander Majercik, Cyril Crassin, Peter Shirley, Morgan McGuire
	//	A Ray-Box Intersection Algorithm and Efficient Dynamic Voxel Rendering
	//	Journal of Computer Graphics Techniques Vol. 7, No. 3, 2018
	// bool slabs( vec3 corner0, vec3 corner1, vec3 rayOrigin, vec3 rayDir ) {
	// 	vec3 invRaydir = vec3( 1.0 ) / rayDir;
	// 	vec3 t0 = ( corner0 - rayOrigin ) * invRaydir;
	// 	vec3 t1 = ( corner1 - rayOrigin ) * invRaydir;
	// 	vec3 tmin = min( t0, t1 ); // min on each x,y,z
	// 	vec3 tmax = max( t0, t1 ); // max on each x,y,z
	// 	return max_component(tmin) <= min_component(tmax); // boolean hit
	// }


	//	Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
	//	"An Efficient and Robust Ray-Box Intersection Algorithm"
	//	Journal of graphics tools, 10(1):49-54, 2005
	hitrecord intersect ( ray r ) const override {
		hitrecord hit; // initially at max distance

		baseType t0 = 0.0, t1 = DMAX_TRAVEL;
		baseType tmin, tmax, tymin, tymax, tzmin, tzmax;
		if (r.direction.x() >= 0) {
			tmin = (bounds[0].x() - r.origin.x()) / r.direction.x();
			tmax = (bounds[1].x() - r.origin.x()) / r.direction.x();
		}	else {
			tmin = (bounds[1].x() - r.origin.x()) / r.direction.x();
			tmax = (bounds[0].x() - r.origin.x()) / r.direction.x();
		}
		if (r.direction.y() >= 0) {
			tymin = (bounds[0].y() - r.origin.y()) / r.direction.y();
			tymax = (bounds[1].y() - r.origin.y()) / r.direction.y();
		}
		else {
			tymin = (bounds[1].y() - r.origin.y()) / r.direction.y();
			tymax = (bounds[0].y() - r.origin.y()) / r.direction.y();
		}
		if ( (tmin > tymax) || (tymin > tmax) ){
			tmin = DMAX_TRAVEL;
			goto end;
		}
		if (tymin > tmin)
			tmin = tymin;
		if (tymax < tmax)
			tmax = tymax;
		if (r.direction.z() >= 0) {
			tzmin = (bounds[0].z() - r.origin.z()) / r.direction.z();
			tzmax = (bounds[1].z() - r.origin.z()) / r.direction.z();
		} else {
			tzmin = (bounds[1].z() - r.origin.z()) / r.direction.z();
			tzmax = (bounds[0].z() - r.origin.z()) / r.direction.z();
		}
		if ( (tmin > tzmax) || (tzmin > tmax) ){
			tmin = DMAX_TRAVEL;
			goto end;
		}
		if (tzmin > tmin)
			tmin = tzmin;
		if (tzmax < tmax)
			tmax = tzmax;
		if ( (tmin > t1) || (tmax < t0) ) {
			tmin = DMAX_TRAVEL;
			goto end;
		}


	end:
		hit.dtransit = std::min( DMAX_TRAVEL, tmin );
		hit.position = r.origin + hit.dtransit * r.direction;

		// WIP - selection based on hit face
		baseType dMX = abs( bounds[ 0 ].values[ 0 ] - hit.position.values[ 0 ] );
		baseType dPX = abs( bounds[ 1 ].values[ 0 ] - hit.position.values[ 0 ] );
		baseType dMY = abs( bounds[ 0 ].values[ 1 ] - hit.position.values[ 1 ] );
		baseType dPY = abs( bounds[ 1 ].values[ 1 ] - hit.position.values[ 1 ] );
		baseType dMZ = abs( bounds[ 0 ].values[ 2 ] - hit.position.values[ 2 ] );
		baseType dPZ = abs( bounds[ 1 ].values[ 2 ] - hit.position.values[ 2 ] );
		baseType dmin = std::min( std::min( dMX, dPX ), std::min( std::min( dMY, dPY ), std::min( dMZ, dPZ ) ) );
		if ( dmin == dMX ) {
			hit.normal = vec3( -1.0,  0.0,  0.0 );
		} else if ( dmin == dPX ) {
			hit.normal = vec3(  1.0,  0.0,  0.0 );
		} else if ( dmin == dMY ) {
			hit.normal = vec3(  0.0, -1.0,  0.0 );
		} else if ( dmin == dPY ) {
			hit.normal = vec3(  0.0,  1.0,  0.0 );
		} else if ( dmin == dMZ ) {
			hit.normal = vec3(  0.0,  0.0, -1.0 );
		} else if ( dmin == dPZ ) {
			hit.normal = vec3(  0.0,  0.0,  1.0 );
		}

		hit.materialID = materialID;
		hit.front = dot( hit.normal, r.direction ) < 0.0 ? true : false;

		// send back the hit record
		return hit;
	}

private: // geometry parameters : min ( index 0 ), max ( index 1 ) on x,y,z axis
	vec3 bounds[ 2 ];
};


//   todo : material handling


class camera{ // camera class generates view vectors from a set of basis vectors
public:
	camera(){}
	void lookat( const vec3 from, const vec3 at, const vec3 up ){
		position = from;
		bz = normalize( at - from );
		bx = normalize( cross( up, bz ) );
		by = normalize( cross( bx, bz ) );
	}
	ray sample( const vec2 p ) const { // argument is pixel location - assumes any desired jitter is applied at call site
		ray r;
		r.origin = position;
		// remap [0, dimension] indexing to [-dimension/2., dimension/2.]
		baseType lx = ( p.values[ 0 ] - baseType( x / 2.0 ) ) / baseType( x / 2.0 );
		baseType ly = ( p.values[ 1 ] - baseType( y / 2.0 ) ) / baseType( y / 2.0 );
		baseType aspect_ratio = baseType( x ) / baseType( y );            // calculate pixel offset
		r.direction = normalize( aspect_ratio * lx * bx + ly * by + ( 1.0 / FoV ) * bz ); // construct from basis
		return r;
	}
private:
	vec3 position;    // location of viewer
	vec3 bx, by, bz;  // basis vectors for sample calcs
	int x = X_IMAGE_DIM, y = Y_IMAGE_DIM; // overal dimensions of the screen
	baseType FoV = FIELD_OF_VIEW; // field of view
};


class scene{ // scene as primitive list + material list container
public:
	scene() { }
	void clear() { contents.clear(); }
	void populate(){
		std::random_device r;
		std::seed_seq s{ r(), r(), r(), r(), r(), r(), r(), r(), r() };
		auto gen = std::make_shared< std::mt19937_64 >( s );
		// for (int i = 0; i < 7; i++)
			// contents.push_back(std::make_shared<sphere>(0.8*randomVector(gen), 0.03*rng(gen), 0));
		// for ( int i = 0; i < NUM_PRIMITIVES; i++ ){
		//   baseType yval = ( ( baseType( i ) / NUM_PRIMITIVES ) - 0.5 ) * 2.0;
		//   vec3 p1 = vec3( std::cos( yval * 6.5 ), std::sin( yval * 14.0 ), yval * pi );
		//   vec3 p2 = vec3( std::cos( yval * 9.7 ) + 0.1, std::sin( yval * 16.4 ) + 0.7, yval);
		//   vec3 p3 = vec3( std::cos( yval * 15.8 ) - 0.3, std::sin( yval * 19.2 ) + 0.1, yval + 0.3 * rng( gen ) ) + randomVector( gen ) * 0.04;
		//   contents.push_back( std::make_shared< triangle >( p1, p2, rng( gen ) < 0.1 ? randomVector( gen ) : p3, rng( gen ) < 0.1 ? 1 : 3 ) );
		//   contents.push_back( std::make_shared< sphere >( randomVector( gen ), 0.4 * rng( gen ), rng( gen ) < 0.4 ? 0 : 2 ) );
		// }
		//
		// contents.push_back( std::make_shared< sphere >( vec3( 1.5, 0.0, 0.0 ), 0.10, 2 ) );
		// contents.push_back( std::make_shared< sphere >( vec3( 0.0, 1.5, 0.0 ), 0.15, 2 ) );
		// contents.push_back( std::make_shared< sphere >( vec3( 0.0, 0.0, 1.5 ), 0.20, 2 ) );


		for( int i = 0; i < 200; i++ ) {
			float x = rng( gen ) - 0.5, y = rng( gen ) - 0.5, z = rng( gen ) - 0.5;
			if ( abs( x ) < 0.25 && abs( y ) < 0.25 && abs( z ) < 0.25 )
				contents.push_back( std::make_shared< sphere >( vec3( x, y, z ), 0.125 * rng( gen ), 2 ) );



			int select = int( floor( rng(gen) * 3.0 ) );
			vec3 v0 = randomVector( gen ) * 0.12577;
			vec3 v1 = randomVector( gen ) * 0.12577;

			vec3 base = vec3( 1.0 );
			for( int j = 0; j < 3; j++ )
				if( j != select )
					base.values[ j ] *= 0.0;

			vec3 base2 = base;
			base.values[select] *= rng( gen ) ? 1.0 : -1.0;

			contents.push_back( std::make_shared< aabb >( base + vec3( std::min( v0.values[ 0 ], v1.values[ 0 ] ), std::min( v0.values[ 1 ], v1.values[ 1 ] ), std::min( v0.values[ 2 ], v1.values[ 2 ] ) ),
																											base2 + vec3( std::max( v0.values[ 0 ], v1.values[ 0 ] ), std::max( v0.values[ 1 ], v1.values[ 1 ] ), std::max( v0.values[ 2 ], v1.values[ 2 ] ) ), 1 ) );
		}

	}
	hitrecord rayQuery( ray r ) const {
		hitrecord h; // iterate through primitives and check for nearest intersection
		baseType currentMin = DMAX_TRAVEL; // initially 'a big number'
		for( int i = 0; i < contents.size(); i++ ) {
			hitrecord temp = contents[ i ]->intersect( r ); temp.primitiveID = i;
			if( temp.dtransit < DMAX_TRAVEL && temp.dtransit > 0. && temp.dtransit < currentMin ) {
				currentMin = temp.dtransit;
				h = temp;
			}
		}
		return h;
	}
	std::vector< std::shared_ptr< primitive > > contents; // list of primitives making up the scene
	// std::vector<std::shared_ptr<material>> materials; // list of materials present in the scene
};


class renderer{
public:
	std::atomic< unsigned long long > tileIndexCounter { 0 }; // used to get new tiles
	std::atomic< unsigned long long > tileFinishCounter{ 0 }; // used for status reporting
	const unsigned long long totalTileCount = std::ceil( X_IMAGE_DIM / TILESIZE_XY ) * ( std::ceil( Y_IMAGE_DIM / TILESIZE_XY ) + 1 );

	renderer() { bytes.resize( xdim * ydim * 4, 0 ); s.populate(); rng_seed(); }
	void renderAndSaveTo( std::string filename ) {
		// c.lookat( vec3( 0.0, 0.0, 2.0 ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );
		// c.lookat( vec3( 5.0, 5.0, 5.0 ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );
		c.lookat( randomUnitVector( gen[ 0 ] ) * ( 2.2 + rng( gen[ 0 ] ) ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );
		std::thread threads[ NUM_THREADS + 1 ];                 // create thread pool
		for ( int id = 0; id <= NUM_THREADS; id++ ){         // do work
			threads[ id ] = ( id == NUM_THREADS ) ? std::thread(// reporter thread
				[ this, id ] () {
					const auto tstart = high_resolution_clock::now();
					while ( true ) { // report timing
						// show status - break on 100% completion
						cout << "\r\033[K";
						const baseType frac = baseType( tileFinishCounter ) / baseType( totalTileCount );

						cout << "["; //  [=====....................] where equals shows progress
						for( int i = 0; i <= PROGRESS_INDICATOR_STOPS * frac; i++ ) cout << "=";
						for( int i = 0; i < PROGRESS_INDICATOR_STOPS * ( 1 - frac); i++ ) cout << ".";
						cout << "]" << std::flush;

						// const int tile_width_char = std::ceil(std::log10(totalTileCount));
						// cout << " (" << std::setw(tile_width_char) << tileFinishCounter << " / " << std::setw(tile_width_char) << totalTileCount << ") " << std::flush;
						cout << "[" << std::setw( 3 ) << 100.0 * frac << "% " << std::flush;

						cout << std::setw( 7 ) << std::showpoint << duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0 << " sec]" << std::flush;

						if( tileFinishCounter >= totalTileCount ){
							const float seconds = duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0;
							const long long total_rays = X_IMAGE_DIM * Y_IMAGE_DIM * NUM_SAMPLES * MAX_BOUNCES;

							cout << "\r\033[K[" << std::string( PROGRESS_INDICATOR_STOPS + 1, '=' ) << "] " << seconds << " sec - total rays: " << total_rays << " (" << total_rays / seconds << "rays/sec)" << endl;
							break;
						}

						// sleep for some amount of time before showing again
						std::this_thread::sleep_for( milliseconds( REPORT_DELAY ) );
					}
				}
			) : std::thread(
				[ this, id ] () {
					// now tile based
					while ( true ) {
						// solve for x and y from the index
						unsigned long long index = tileIndexCounter.fetch_add( 1 );
						if ( index >= totalTileCount ) break;

						constexpr int numTilesX = int( std::ceil( float( X_IMAGE_DIM ) / float( TILESIZE_XY ) ) );
						constexpr int numTilesY = int( std::ceil( float( Y_IMAGE_DIM ) / float( TILESIZE_XY ) ) );

						const int tile_x_index = index % numTilesX;
						const int tile_y_index = ( index / numTilesX );

						const int tile_base_x = tile_x_index * TILESIZE_XY;
						const int tile_base_y = tile_y_index * TILESIZE_XY;

						for ( int y = tile_base_y; y < tile_base_y + TILESIZE_XY; y++ )
						for ( int x = tile_base_x; x < tile_base_x + TILESIZE_XY; x++ ) {
							vec3 running_color = vec3( 0. );								// initially zero, averages sample data
							for ( int s = 0; s < nsamples; s++ )						// get sample data (n samples)
								running_color += pathtraceSample( x, y, id );	// bounce your ray around
							running_color /= baseType( nsamples );					// sample averaging
							tonemapAndGamma( running_color );								// tonemapping + gamma
							write( running_color, vec2( x, y ) );						// write final output values
						}
						tileFinishCounter.fetch_add( 1 );
					}
				}
			);
		}
		for ( int id = 0; id <= NUM_THREADS; id++ )
			threads[ id ].join();
		cout << "Writing \'" << filename << "\'";
		const auto tistart = high_resolution_clock::now();
		stbi_write_png( filename.c_str(), xdim, ydim, 4, &bytes[0], xdim * 4 );
		cout << " - " << duration_cast< milliseconds >( high_resolution_clock::now() - tistart ).count() / 1000.0 << " seconds to write to disk" << endl;
	}
private:
	camera c; 		// generates view rays
	scene s; 			// holds all scene geometry + their associated materials
	int xdim = X_IMAGE_DIM, ydim = Y_IMAGE_DIM, nsamples = NUM_SAMPLES, bmax = MAX_BOUNCES;
	std::vector< unsigned char > bytes;											// image buffer for stb_image_write
	std::vector< std::shared_ptr< std::mt19937_64 > > gen;	// PRNG states per thread
	void rng_seed() {
		std::random_device r;
		for ( int i = 0; i < NUM_THREADS; i++ ) {
			std::seed_seq s{ r(), r(), r(), r(), r(), r(), r(), r(), r() };
			gen.push_back( std::make_shared< std::mt19937_64 >( s ) );
		}
	}

	// todo - add switch handling, for some different tonemap curves
	void tonemapAndGamma( vec3& in ){
		in *= 0.6f; // function to tonemap color value in place
		baseType a = 2.51;
		baseType b = 0.03;
		baseType c = 2.43;
		baseType d = 0.59;
		baseType e = 0.14;
		in = ( in * ( a * in + vec3( b ) ) ) / ( in * ( c * in + vec3( d ) ) + e ); // tonemap
		in.values[ 0 ] = std::pow( std::clamp( in.values[ 0 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA ); // gamma correct
		in.values[ 1 ] = std::pow( std::clamp( in.values[ 1 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA );
		in.values[ 2 ] = std::pow( std::clamp( in.values[ 2 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA );
  }

	vec3 pathtraceSample( const int x, const int y, const int id ){
		// throughput's initial value of 1. in each channel indicates that it is initially
		// capable of carrying all of the light intensity possible (100%), and it is reduced
		vec3 throughput	= vec3( 1.0 );	// by the albedo of the material on each bounce
		vec3 current		= vec3( 0.0 );	// init to zero - initially no light present
		vec3 old_ro;										// old_ro holds previous hit location, unitialized

		// get initial ray origin + ray direction from camera
		ray r = c.sample( vec2( x + rng( gen[ id ] ), y + rng( gen[ id ] ) ) );
		for ( int bounce = 0; bounce < MAX_BOUNCES; bounce++ ) {
			old_ro = r.origin; 							// cache old hit location
			hitrecord h = s.rayQuery( r );	// get a new hit location (scene query)

			if ( h.dtransit == DMAX_TRAVEL ) break;
			r.origin = r.origin + h.dtransit * r.direction + h.normal * HIT_EPSILON;
			r.direction = normalize( ( 1.0 + HIT_EPSILON ) * h.normal + randomUnitVector( gen[ id ] ) ); // diffuse reflection

			// the form is:
				// current    += throughput*current_emission // emission term
				// throughput *= albedo                     // diffuse absorption term

			// if ( h.materialID == 0 ){ // the spheres
			// 	// current += throughput * vec3(1.3, 1.2, 1.1);
			// 	// throughput *= vec3(0.999);
			// 	throughput *= palette( h.primitiveID * PALETTE_SCALAR );
			// } else if ( h.materialID == 1 && h.front ) {
			// 	current += throughput * BRIGHTNESS_SCALAR * vec3( h.uv.values[ 0 ], h.uv.values[ 1 ], 1 - h.uv.values[ 0 ] - h.uv.values[ 1 ] );
			// } else if ( h.materialID == 1 && !h.front ) {
			// 	current += throughput * ( palette( h.primitiveID * PALETTE_SCALAR ) * BRIGHTNESS_SCALAR );
			// } else if ( h.materialID == 2 ) {
			if ( h.materialID == 2 ) {
				r.direction = reflect( r.origin - old_ro, h.normal );
				throughput *= vec3( 0.89 );
			} else if ( h.materialID == 1 ) {
				current += throughput * abs( h.normal );
			}


			// } else if ( h.materialID == 3 ) {
			// 	throughput *= vec3( 0.999 );
			// } else if ( h.dtransit == DMAX_TRAVEL ) {
			// // current += throughput * 0.1 * vec3(0.918, 0.75, 0.6); // sky color and escape
			// 	break; // escape
			// }


			baseType p = std::max( throughput.values[ 0 ], std::max( throughput.values[ 1 ], throughput.values[ 2 ] ) );
			if ( rng( gen[ id ] ) > p ) break;	// russian roulette termination check
			throughput *= 1.0 / p; 							// russian roulette compensation term
		}
		return current;
	}
	void write( vec3 col, vec2 loc ){ // writes to image buffer
		if ( loc.values[ 0 ] < 0 || loc.values[ 0 ] >= X_IMAGE_DIM ) return;
		if ( loc.values[ 1 ] < 0 || loc.values[ 1 ] >= Y_IMAGE_DIM ) return;
		const int index = 4 * ( loc.values[ 1 ] * xdim + loc.values[ 0 ] );
		for( int c = 0; c < 4; c++ )
			bytes[ index + c ] = ( c == 3 ) ? 255 : col.values[ c ] * 255.0;
	}
};

int main ( int argc, char const *argv[] ) {
	// std::string filename = std::string(argv[1]); // from CLI
	const auto tstart = high_resolution_clock::now();

	for ( size_t i = 0; i <= 3; i++ ) {
		std::stringstream s; s << "outputs/out" << std::to_string( i ) << ".png";
		renderer r;
		r.renderAndSaveTo( s.str() );
	}

	cout << "Total Render Time: " <<
		duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0
			<< " seconds" << endl; return 0;
}

#include <chrono>				// timing utilities
#include <iostream>			// text i/o
#include <iomanip>			// cout formatting
#include <stdio.h>			// printf if needed
#include <vector>				// std::vector
#include <random>				// prng
#include <string>				// std::string
#include <sstream>			// std::stringstream
#include <algorithm>		// clamp
#include <atomic>				// atomic_llong
#include <thread>				// threads
#include <memory>				// shared_ptr

using std::cerr, std::cin, std::cout, std::endl, std::flush;
using std::chrono::high_resolution_clock, std::chrono::duration_cast, std::chrono::milliseconds;
using std::make_shared;

// my vector library
#include "lib/vector.h"

// primitives
#include "lib/primitives.h"

// config
#include "lib/defines.h"

// image input
#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"

// image output
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

class wang {
public:
	wang( uint32_t s ) : seed( s ) {}
	uint32_t seed;
	void hash() {
		seed = ( seed ^ 12345391 ) * 2654435769;
		seed ^= ( seed << 6 ) ^ ( seed >> 26 );
		seed *= 2654435769;
		seed += ( seed << 5 ) ^ ( seed >> 12 );
//--
		// seed += ~(seed << 15);
		// seed ^= (seed >> 10);
		// seed += (seed << 3);
		// seed ^= (seed >> 6);
		// seed += ~(seed << 11);
		// seed ^= (seed >> 16);
	}

	baseType getNum() {
		hash();
		return baseType( seed ) / 4294967295.0;
	}
};


// std::random Utilities
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



// iq style palettes
int paletteToUse = 0;
vec3 palette( baseType t ) {
	vec3 a, b, c, d;
	switch ( paletteToUse ) {
		case 0:
			a = vec3( 0.50, 0.50, 0.50 );
			b = vec3( 0.50, 0.50, 0.50 );
			c = vec3( 1.00, 1.00, 1.00 );
			d = vec3( 0.00, 0.10, 0.20 );
			break;
		case 1:
			a = vec3( 0.50, 0.50, 0.50 );
			b = vec3( 0.50, 0.50, 0.50 );
			c = vec3( 1.00, 1.00, 0.50 );
			d = vec3( 0.80, 0.90, 0.30 );
			break;
		case 2:
			a = vec3( 0.50, 0.50, 0.50 );
			b = vec3( 0.50, 0.50, 0.50 );
			c = vec3( 1.00, 0.70, 0.40 );
			d = vec3( 0.00, 0.15, 0.20 );
			break;

		default: break;
	}
	vec3 temp = ( c * t + d ) * 2.0 * pi;
	return a + b * vec3( cos( temp.values[ 0 ] ), cos( temp.values[ 1 ] ), cos( temp.values[ 2 ] ) );
}



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

	void recursiveSplit( vec3 min, vec3 max /*, int previousAxisPick */ ) {
		if( length( min - max ) < 0.1 ||
				abs( max.values[ 0 ] - min.values[ 0 ] ) < 0.01 ||
				abs( max.values[ 1 ] - min.values[ 1 ] ) < 0.01 ||
				abs( max.values[ 2 ] - min.values[ 2 ] ) < 0.01 )
			return;

		if( rng( gen ) < 0.1 ) {
			contents.push_back( make_shared< aabb > ( min + vec3( 0.005 ), max - vec3( 0.005 ), rng( gen ) < 0.3 ? rng( gen ) < 0.2 ? 2 : 3 : 1 ) ); // shrink slightly, to create gaps
		} else {
			// pick an axis

			int axisPick;
			// while( previousAxisPick == ( axisPick = int( floor( rng( gen ) * 3.0 ) ) ) );
			axisPick = int( floor( rng( gen ) * 3.0 ) );
			// cout << "picked axis " << axisPick << endl << std::flush;

			// get a value 0-1 to split by
			baseType axisSplit = ( rng( gen ) * 0.75 ) + 0.15;
			baseType midpointLoc = ( max.values[ axisPick ] - min.values[ axisPick ] ) * axisSplit + min.values[ axisPick ];

			// cout << "input mins: " << min.values[ 0 ] << " " << min.values[ 1 ] << " " << min.values[ 2 ] << std::endl << std::flush;
			// cout << "input maxs: " << max.values[ 0 ] << " " << max.values[ 1 ] << " " << max.values[ 2 ] << std::endl << std::flush;


			// generate the aabb dimensions for the two split boxes
			vec3 minMiddle = min;
			minMiddle.values[ axisPick ] = midpointLoc;
			// cout << "computed minMiddle: " << minMiddle.values[ 0 ] << " " << minMiddle.values[ 1 ] << " " << minMiddle.values[ 2 ] << std::endl << std::flush;

			vec3 maxMiddle = max;
			maxMiddle.values[ axisPick ] = midpointLoc;
			// cout << "computed maxMiddle: " << maxMiddle.values[ 0 ] << " " << maxMiddle.values[ 1 ] << " " << maxMiddle.values[ 2 ] << std::endl << std::flush;
			// cout << endl << endl;

			// recurse, branch in 2 split on that axis
			recursiveSplit( min, maxMiddle /* , axisPick */ );
			recursiveSplit( minMiddle, max /* , axisPick */ );
		}
	}

	void recursiveWangSplit( vec3 min, vec3 max, int previousAxisPick, wang myWang, int depth ) {
		if( length( min - max ) < BOX_EPSILON ||
				abs( max.values[ 0 ] - min.values[ 0 ] ) < BOX_EPSILON / 10.0 ||
				abs( max.values[ 1 ] - min.values[ 1 ] ) < BOX_EPSILON / 10.0 ||
				abs( max.values[ 2 ] - min.values[ 2 ] ) < BOX_EPSILON / 10.0 )
			return;

		// for( int i = 0; i < depth; i++ ) myWang.getNum(); // cycle it a few times

		// if( myWang.getNum() < 0.01 && rng( gen ) < 0.5 ) {
		if( myWang.getNum() * depth > 0.45 * depth  ) {
			contents.push_back( make_shared< aabb > ( min + vec3( 0.005 ), max - vec3( 0.005 ), myWang.getNum() < 0.3 ? 3 : 1 ) ); // shrink slightly, to create gaps
		} else {
			int axisPick; // pick an axis, make sure it's different than the previousy picked one
			while( previousAxisPick == ( axisPick = int( floor( myWang.getNum() * 3.0 ) ) ) );

			// baseType axisSplit = ( myWang.getNum() ) + 0.15; // get a value 0-1 to split by
			baseType axisSplit = 0.5;
			baseType midpointLoc = ( max.values[ axisPick ] - min.values[ axisPick ] ) * axisSplit + min.values[ axisPick ];

			// generate the aabb dimensions for the two split boxes
			vec3 minMiddle = min; minMiddle.values[ axisPick ] = midpointLoc;
			vec3 maxMiddle = max; maxMiddle.values[ axisPick ] = midpointLoc;

			// recurse, branch in 2 split on that axis
			recursiveWangSplit( min, maxMiddle, axisPick, myWang, depth + 1 );
			recursiveWangSplit( minMiddle, max, axisPick, myWang, depth + 1 );
		}
	}


	void recursiveMultiSplit( vec3 min, vec3 max, int previousAxisPick ) {
	// five options:
		// epsilon condition and break
		// draw box and break
		// draw no box and break
		// split box( regular splits ) and recurse
		// split box( irregular splits ) and recurse
		if( length( min - max ) < BOX_EPSILON || // break on tiny box
			abs( max.values[ 0 ] - min.values[ 0 ] ) < BOX_EPSILON / 10.0 ||
			abs( max.values[ 1 ] - min.values[ 1 ] ) < BOX_EPSILON / 10.0 ||
			abs( max.values[ 2 ] - min.values[ 2 ] ) < BOX_EPSILON / 10.0 ) {
			return;
		} else if ( rng( gen ) < 0.1 ) { // draw box and break out
			contents.push_back( make_shared< aabb > ( min + vec3( 0.005 ), max - vec3( 0.005 ), rng( gen ) < 0.3 ? 3 : 1 ) ); // shrink slightly, to create gaps
			return;
		} else { // continue down the tree
			int axisPick;
			while( previousAxisPick == ( axisPick = int( floor( rng( gen ) * 3.0 ) ) ) );

			baseType axisLength = max.values[ axisPick ] - min.values[ axisPick ];
			int numSplits = int( floor( rng( gen ) * 6.0 ) );
			baseType sectionWidth = axisLength / baseType( numSplits + 1 );	// 1 cuts in 2, 2 cuts in 3...

			vec3 cutMin = min, cutMax = max;
			cutMax.values[ axisPick ] = cutMin.values[ axisPick ] + sectionWidth;
			for ( int i = 0; i <= numSplits; i++ ) {
				// increment to walk along the selected axis
				cutMin.values[ axisPick ] += sectionWidth;
				cutMax.values[ axisPick ] += sectionWidth;
				// recursive call
				recursiveMultiSplit( cutMin, cutMax, axisPick );
			}
		}
	}

	void recursiveWangMultiSplit( vec3 min, vec3 max, int previousAxisPick, wang myWang ) {
	// five options:
		// epsilon condition and break
		// draw box and break
		// draw no box and break
		// split box( regular splits ) and recurse
		// split box( irregular splits ) and recurse

		if( length( min - max ) < BOX_EPSILON || // break on tiny box
			abs( max.values[ 0 ] - min.values[ 0 ] ) < BOX_EPSILON / 10.0 ||
			abs( max.values[ 1 ] - min.values[ 1 ] ) < BOX_EPSILON / 10.0 ||
			abs( max.values[ 2 ] - min.values[ 2 ] ) < BOX_EPSILON / 10.0 ) {
				// cout << " breaking on epsilon " << endl;
			return;
		} else if ( myWang.getNum() < 0.1618 ) { // draw box and break out
			contents.push_back( make_shared< aabb > ( min + vec3( 0.005 ), max - vec3( 0.005 ), myWang.getNum() < 0.3 ? 3 : 1 ) ); // shrink slightly, to create gaps
			// cout << "drawing box " << contents.size();
			return;
		} else { // continue down the tree
			// cout << " branching " << endl;
			int axisPick;
			// while( previousAxisPick == ( axisPick = int( floor( myWang.getNum() * 3.0 ) ) ) );
			axisPick = int( floor( myWang.getNum() * 3.0 ) );

			baseType axisLength = max.values[ axisPick ] - min.values[ axisPick ];
			int numSplits = int( floor( myWang.getNum() * 6.0 ) );	// todo: replace with wang hash rng
			baseType sectionWidth = axisLength / baseType( numSplits + 1 );	// 1 cuts in 2, 2 cuts in 3...

			vec3 cutMin = min, cutMax = max;
			cutMax.values[ axisPick ] = cutMin.values[ axisPick ] + sectionWidth;
			for ( int i = 0; i <= numSplits; i++ ) {
				// increment to walk along the selected axis
				cutMin.values[ axisPick ] += sectionWidth;
				cutMax.values[ axisPick ] += sectionWidth;
				// recursive call
				recursiveWangMultiSplit( cutMin, cutMax, axisPick, myWang );
			}
		}
	}

	// validating deterministic wang rng logic
	void recursive( wang myWang, int depth ){
		if( depth == 5 ) return;
		cout << "seeing " << myWang.getNum() << " at depth " << depth << endl;
		recursive( myWang, depth + 1 );
		recursive( myWang, depth + 1 );
	}

	void populate(){
		std::random_device r;
		std::seed_seq s{ r(), r(), r(), r(), r(), r(), r(), r(), r() };
		gen = make_shared< std::mt19937_64 >( s );
		// for (int i = 0; i < 7; i++)
			// contents.push_back(make_shared<sphere>(0.8*randomVector(gen), 0.03*rng(gen), 0));
		// for ( int i = 0; i < NUM_PRIMITIVES; i++ ){
		//   baseType yval = ( ( baseType( i ) / NUM_PRIMITIVES ) - 0.5 ) * 2.0;
		//   vec3 p1 = vec3( std::cos( yval * 6.5 ), std::sin( yval * 14.0 ), yval * pi );
		//   vec3 p2 = vec3( std::cos( yval * 9.7 ) + 0.1, std::sin( yval * 16.4 ) + 0.7, yval);
		//   vec3 p3 = vec3( std::cos( yval * 15.8 ) - 0.3, std::sin( yval * 19.2 ) + 0.1, yval + 0.3 * rng( gen ) ) + randomVector( gen ) * 0.04;
		//   contents.push_back( make_shared< triangle >( p1, p2, rng( gen ) < 0.1 ? randomVector( gen ) : p3, rng( gen ) < 0.1 ? 1 : 3 ) );
		//   contents.push_back( make_shared< sphere >( randomVector( gen ), 0.4 * rng( gen ), rng( gen ) < 0.4 ? 0 : 2 ) );
		// }
		//
		// contents.push_back( make_shared< sphere >( vec3( 1.5, 0.0, 0.0 ), 0.10, 2 ) );
		// contents.push_back( make_shared< sphere >( vec3( 0.0, 1.5, 0.0 ), 0.15, 2 ) );
		// contents.push_back( make_shared< sphere >( vec3( 0.0, 0.0, 1.5 ), 0.20, 2 ) );


		vec3 v1, v2, v3;
		v1 = randomUnitVector( gen );
		v2 = randomUnitVector( gen );
		v3 = randomUnitVector( gen );

		for( int i = 0; i < NUM_PRIMITIVES; i++ ) {

			// shell
			contents.push_back( make_shared< sphere >( randomUnitVector( gen ) * 1.4, 0.1 * rng( gen ), rng( gen ) < 0.1618 ? 3 : 2 ) );

			// segments
			contents.push_back( make_shared< sphere >( v1 + baseType( i ) * ( ( v2 - v1 ) / NUM_PRIMITIVES ), 0.1 * rng( gen ), rng( gen ) < 0.1618 ? 3 : 2 ) );
			contents.push_back( make_shared< sphere >( v2 + baseType( i ) * ( ( v3 - v2 ) / NUM_PRIMITIVES ), 0.1 * rng( gen ), rng( gen ) < 0.1618 ? 3 : 2 ) );
			contents.push_back( make_shared< sphere >( v3 + baseType( i ) * ( ( v1 - v3 ) / NUM_PRIMITIVES ), 0.1 * rng( gen ), rng( gen ) < 0.1618 ? 3 : 2 ) );



			// int select = int( floor( rng( gen ) * 3.0 ) );
			// vec3 v0 = randomVector( gen );
			// vec3 v1 = randomVector( gen );

			// contents.push_back( make_shared< aabb >( vec3( std::min( v0.values[ 0 ], v1.values[ 0 ] ), std::min( v0.values[ 1 ], v1.values[ 1 ] ), std::min( v0.values[ 2 ], v1.values[ 2 ] ) ),
				//																	vec3( std::max( v0.values[ 0 ], v1.values[ 0 ] ), std::max( v0.values[ 1 ], v1.values[ 1 ] ), std::max( v0.values[ 2 ], v1.values[ 2 ] ) ), rng( gen ) < 0.5 ? 3 : 1 ) );
		}



		wang myWang( 69420 * rng( gen ) * 10000 );


		// recursiveSplit( vec3( -1.0 ), vec3( 1.0 ), 0 );
		// recursiveSplit( vec3( -1.0, -0.25, -0.5 ), vec3( 1.0, 0.25, 0.5 ) );

		float x, y, z;
		x = rng( gen ) * 0.25 + 0.25;
		y = rng( gen ) * 0.15 + 0.50;
		z = rng( gen ) * 0.33 + 0.60;

		// recursiveWangSplit( vec3( -1.0 * x, -1.0 * y, -1.0 * z ), vec3( x, y, z ), int( floor( rng( gen ) * 3.0 ) ), myWang, 0 );
		recursiveMultiSplit( vec3( -1.0 * x, -1.0 * y, -1.0 * z ), vec3( x, y, z ), int( floor( rng( gen ) * 3.0 ) ) );


		// recursiveWangMultiSplit( vec3( -1.0, -0.25, -0.5 ), vec3( 1.0, 0.25, 0.5 ), 1, 69420 * rng( gen ) );
		// recursiveWangMultiSplit( vec3( -1.0 ), vec3( 1.0 ), 1, 69420 * rng( gen ) );
		// cout << "drawing with " << contents.size() << " primitives " << endl;

		//recursive( myWang, 0 );
		//cout << endl;


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

	std::shared_ptr< std::mt19937_64 > gen;
};


class renderer{
public:
	std::atomic< unsigned long long > tileIndexCounter { 0 }; // used to get new tiles
	std::atomic< unsigned long long > tileFinishCounter{ 0 }; // used for status reporting
	const unsigned long long totalTileCount = std::ceil( X_IMAGE_DIM / TILESIZE_XY ) * ( std::ceil( Y_IMAGE_DIM / TILESIZE_XY ) + 1 );

	renderer() { bytes.resize( xdim * ydim * 4, 0 ); rng_seed(); }
	void renderAndSaveTo( std::string filename ) {
		// c.lookat( vec3( 0.0, 0.0, 2.0 ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );
		// c.lookat( vec3( 5.0, 5.0, 5.0 ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );

		while ( s.contents.size() < 10 + NUM_PRIMITIVES * 4 ) {
			s.clear();
			s.populate();
		}

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
			gen.push_back( make_shared< std::mt19937_64 >( s ) );
		}
	}

	// todo - add switch handling, for some different tonemap curves
	void tonemapAndGamma( vec3& in ){
		// in *= 0.6f; // function to tonemap color value in place
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
				// current += throughput * abs( h.normal );
				current += throughput * palette( ( 2.2 / NUM_PRIMITIVES ) * h.primitiveID * PALETTE_SCALAR );
			} else if ( h.materialID == 3 ) {
				throughput *= vec3( 0.65 );
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

	int numImages = IMAGE_SEQUENCE_END_INDEX - IMAGE_SEQUENCE_START_INDEX;
	cout << "Drawing " << numImages << "x " << X_IMAGE_DIM << "x" << Y_IMAGE_DIM <<  " images at " << NUM_SAMPLES << "spp with " << MAX_BOUNCES << " max bounces" << endl;
	cout << "Rendering with " << NUM_THREADS-1 << " worker threads starting at " << IMAGE_SEQUENCE_START_INDEX << endl;


	for ( size_t i = IMAGE_SEQUENCE_START_INDEX; i <= IMAGE_SEQUENCE_END_INDEX; i++ ) {
		paletteToUse = i%3;
		std::stringstream s; s << "outputs/out" << std::to_string( i ) << ".png";
		renderer r;
		r.renderAndSaveTo( s.str() );
	}

	cout << "Total Render Time: " <<
		duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0
			<< " seconds" << endl; return 0;
}

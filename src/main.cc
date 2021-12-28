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

// todo:
  // cleanup the handling of thread state
  // more extensive performance report on completion (samples/second, resolution, bounces, etc)

// my vector library
#include "lib/vector.h"

// image input
#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"

// image output
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

// default types
#define base_type double
using vec2 = vector2< base_type >;
using vec3 = vector3< base_type >;

// render parameters
constexpr long long X_IMAGE_DIM = 1920 / 2;
constexpr long long Y_IMAGE_DIM = 1080 / 2;
constexpr long long TILESIZE_XY = 8;
constexpr long long MAX_BOUNCES = 69;
constexpr long long NUM_SAMPLES = 420;
constexpr long long NUM_THREADS = 4;
constexpr base_type IMAGE_GAMMA = 2.2;
constexpr base_type HIT_EPSILON = base_type( std::numeric_limits< base_type >::epsilon() );
constexpr base_type DMAX_TRAVEL = base_type( std::numeric_limits< base_type >::max() ) / 10.0;
constexpr base_type FIELD_OF_VIEW = 0.69420;
constexpr base_type PALETTE_SCALAR = 16.18;
constexpr base_type BRIGHTNESS_SCALAR = 16.18;
constexpr long long REPORT_DELAY = 618; // reporter thread sleep duration, in ms
constexpr long long NUM_PRIMITIVES = 69;
constexpr long long PROGRESS_INDICATOR_STOPS = 69; // cli spaces to take up

// ray representation (origin+direction)
struct ray {
  vec3 origin;
  vec3 direction;
};

// represents a ray hit and the associated information
struct hitrecord {                // hit record
    vec3 position;               // position
    vec3 normal;                // normal
    base_type dtransit = DMAX_TRAVEL; // how far the ray traveled, initially very large
    int material_index = -1;         // material (indexed into scene list)
    int primitive_index = 0;        // so you can refer to the values for this primitive later
    vec2 uv;                       // used for triangles, barycentric coords
    bool front;                   // hit on frontfacing side
};

inline uint32_t wang_hash( uint32_t x ){
    x = ( x ^ 12345391 ) * 2654435769;
    x ^= ( x << 6 ) ^ ( x >> 26 );
    x *= 2654435769;
    x += ( x << 5 ) ^ ( x >> 12 );
    return x;
}

// Random Utilities
base_type rng( std::shared_ptr< std::mt19937_64 > gen ) { // gives a value in the range 0.-1.
  std::uniform_real_distribution< base_type > distribution( 0., 1. );
  return distribution( *gen );
}
vec3 random_vector( std::shared_ptr< std::mt19937_64 > gen ) { // random vector centered around 0.
  return vec3( rng( gen ), rng( gen ), rng( gen ) ) - vec3( 0.5 );
}
vec3 random_unit_vector( std::shared_ptr< std::mt19937_64 > gen ) { // random direction vector (unit length)
  base_type z = rng( gen ) * 2.0 - 1.0;
  base_type a = rng( gen ) * 2.0 * pi;
  base_type r = sqrt( 1.0 - z * z );
  base_type x = r * cos( a );
  base_type y = r * sin( a );
  return vec3( x, y, z );
}
vec3 random_in_unit_disk( std::shared_ptr< std::mt19937_64 > gen ) { // random in unit disk (xy plane)
  vec3 val = random_unit_vector( gen );
  return vec3( val.values[ 0 ], val.values[ 1 ], 0.0 );
}

// iq style palette
vec3 palette( base_type t, vec3 a = vec3( 0.5, 0.5, 0.5 ), vec3 b = vec3( 0.5, 0.5, 0.5 ), vec3 c = vec3( 1.0, 1.0, 1.0 ), vec3 d = vec3( 0.00, 0.33, 0.67 ) ) {
  vec3 temp = ( c * t + d ) * 2.0 * pi;
  return a + b * vec3( cos( temp.values[ 0 ] ), cos( temp.values[ 1 ] ), cos( temp.values[ 2 ] ) );
}

class primitive { // base class for primitives
public:
  virtual hitrecord intersect( ray r ) const = 0; // pure virtual, base definition dne
  int material_index; // indexes into scene material list
};

// sphere
class sphere : public primitive {
public:
  sphere( vec3 c, base_type r, int m ) : center( c ), radius( r ) { material_index = m; }
  hitrecord intersect( ray r ) const override {
    hitrecord h; h.dtransit = DMAX_TRAVEL;
    vec3 disp = r.origin - center;
    base_type b = dot( r.direction, disp );
    base_type c = dot(disp, disp) - radius*radius;
    base_type des = b * b - c; // b squared minus c - discriminant of the quadratic
    if( des >= 0.0 ){ // hit at either one or two points
      base_type d = std::min( std::max( -b + std::sqrt( des ), 0.0 ), std::max( -b - std::sqrt( des ), 0.0 ) );
      if( d >= 0.0 ){ // make sure at least one intersection point is in front of the camera before continuing
        h.dtransit = d;
        h.material_index = material_index;
        h.position = r.origin + h.dtransit * r.direction;
        h.normal = normalize( h.position - center );
        h.front = dot( h.normal, r.direction ) < 0.0 ? true : false;
      }
    }
    return h;
  }
private:  // geometry parameters
  vec3 center;
  base_type radius;
};

// triangle
class triangle : public primitive {
public:
  triangle( vec3 p0, vec3 p1, vec3 p2, int m ) : points{ p0, p1, p2 } { material_index = m; }
  hitrecord intersect( ray r ) const override { // Möller–Trumbore intersection algorithm
    hitrecord hit;  hit.dtransit = DMAX_TRAVEL;
    const vec3 edge1 = points[ 1 ] - points[ 0 ];
    const vec3 edge2 = points[ 2 ] - points[ 0 ];
    const vec3 pvec = cross( r.direction, edge2 );
    const base_type det = dot( edge1, pvec );
    if ( det > -HIT_EPSILON && det < HIT_EPSILON )
      return hit; // no hit, return

    const base_type invDet = 1.0f / det;
    const vec3 tvec = r.origin - points[ 0 ];
    hit.uv.values[ 0 ] = dot( tvec, pvec ) * invDet; // u value
    if ( hit.uv.values[ 0 ] < 0.0f || hit.uv.values[ 0 ] > 1.0f )
      return hit; // no hit, return

    const vec3 qvec = cross( tvec, edge1 );
    hit.uv.values[ 1 ] = dot( r.direction, qvec ) * invDet; // v value
    if ( hit.uv.values[ 1 ] < 0.0f || hit.uv.values[ 0 ] + hit.uv.values[ 1 ] > 1.0f )
      return hit; // no hit, return

    hit.dtransit = dot( edge2, qvec ) * invDet; // distance term to hit
    hit.position = points[ 0 ] + hit.uv.values[ 0 ] * edge2 + hit.uv.values[ 1 ] * edge1;
    hit.normal = cross( edge1, edge2 );
    hit.material_index = material_index;
    hit.front = dot( hit.normal, r.direction ) < 0.0 ? true : false; // determine front or back

    return hit; // return true result with all relevant info
  }
private:  // geometry parameters
  vec3 points[ 3 ];
};


// axis aligned bounding box
class aabb : public primitive {

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
    base_type lx = ( p.values[ 0 ] - base_type( x / 2.0 ) ) / base_type( x / 2.0 );
    base_type ly = ( p.values[ 1 ] - base_type( y / 2.0 ) ) / base_type( y / 2.0 );
    base_type aspect_ratio = base_type( x ) / base_type( y );            // calculate pixel offset
    r.direction = normalize( aspect_ratio * lx * bx + ly * by + ( 1.0 / FoV ) * bz ); // construct from basis
    return r;
  }
private:
  vec3 position;    // location of viewer
  vec3 bx, by, bz;  // basis vectors for sample calcs
  int x = X_IMAGE_DIM, y = Y_IMAGE_DIM; // overal dimensions of the screen
  base_type FoV = FIELD_OF_VIEW; // field of view
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
      // contents.push_back(std::make_shared<sphere>(0.8*random_vector(gen), 0.03*rng(gen), 0));
    for ( int i = 0; i < NUM_PRIMITIVES; i++ ){
      base_type yval = ( ( base_type( i ) / NUM_PRIMITIVES ) - 0.5 ) * 2.0;
      vec3 p1 = vec3( std::cos( yval * 6.5 ), std::sin( yval * 14.0 ), yval * pi );
      vec3 p2 = vec3( std::cos( yval * 9.7 ) + 0.1, std::sin( yval * 16.4 ) + 0.7, yval);
      vec3 p3 = vec3( std::cos( yval * 15.8 ) - 0.3, std::sin( yval * 19.2 ) + 0.1, yval + 0.3 * rng( gen ) ) + random_vector( gen ) * 0.04;
      contents.push_back( std::make_shared< triangle >( p1, p2, rng( gen ) < 0.1 ? random_vector( gen ) : p3, rng( gen ) < 0.1 ? 1 : 3 ) );
      contents.push_back( std::make_shared< sphere >( random_vector( gen ), 0.4 * rng( gen ), rng( gen ) < 0.4 ? 0 : 2 ) );
    }
  }
  hitrecord ray_query( ray r ) const {
    hitrecord h; // iterate through primitives and check for nearest intersection
    base_type current_min = DMAX_TRAVEL; // initially 'a big number'
    for( int i = 0; i < contents.size(); i++ ) {
      hitrecord temp = contents[ i ]->intersect( r ); temp.primitive_index = i;
      if( temp.dtransit < DMAX_TRAVEL && temp.dtransit > 0. && temp.dtransit < current_min ) {
        current_min = temp.dtransit;
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
  std::atomic< unsigned long long > tile_index_counter{ 0 }; // used to get new tiles
  std::atomic< unsigned long long > tile_finish_counter{ 0 }; // used for status reporting
  const unsigned long long total_tile_count = std::ceil( X_IMAGE_DIM / TILESIZE_XY ) * ( std::ceil( Y_IMAGE_DIM / TILESIZE_XY ) + 1 );

  renderer() { bytes.resize( xdim * ydim * 4, 0 ); s.populate(); rng_seed(); }
  void render_and_save_to( std::string filename ) {
    // c.lookat(vec3(0., 0., 2.), vec3(0.), vec3(0.,1.,0.));
    c.lookat( random_unit_vector( gen[ 0 ] ) * ( 2.2 + rng( gen[ 0 ] ) ), vec3( 0.0 ), vec3( 0.0, 1.0, 0.0 ) );
    std::thread threads[ NUM_THREADS + 1 ];                 // create thread pool
    for ( int id = 0; id <= NUM_THREADS; id++ ){         // do work
      threads[ id ] = ( id == NUM_THREADS ) ? std::thread(// reporter thread
        [ this, id ] () {
          const auto tstart = high_resolution_clock::now();
          while ( true ) { // report timing
            // show status - break on 100% completion
            cout << "\r\033[K";
            const base_type frac = base_type( tile_finish_counter ) / base_type( total_tile_count );

            cout << "["; //  [=====....................] where equals shows progress
            for( int i = 0; i <= PROGRESS_INDICATOR_STOPS * frac; i++ ) cout << "=";
            for( int i = 0; i < PROGRESS_INDICATOR_STOPS * ( 1 - frac); i++ ) cout << ".";
            cout << "]" << std::flush;

            // const int tile_width_char = std::ceil(std::log10(total_tile_count));
            // cout << " (" << std::setw(tile_width_char) << tile_finish_counter << " / " << std::setw(tile_width_char) << total_tile_count << ") " << std::flush;
            cout << "[" << std::setw( 3 ) << 100.0 * frac << "% " << std::flush;

            cout << std::setw( 7 ) << std::showpoint << duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0
                  << " sec]" << std::flush;

            if( tile_finish_counter >= total_tile_count ){
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
            unsigned long long index = tile_index_counter.fetch_add( 1 );
            if ( index >= total_tile_count ) break;

            constexpr int num_tiles_x = int( std::ceil( float( X_IMAGE_DIM ) / float( TILESIZE_XY ) ) );
            constexpr int num_tiles_y = int( std::ceil( float( Y_IMAGE_DIM ) / float( TILESIZE_XY ) ) );

            const int tile_x_index = index % num_tiles_x;
            const int tile_y_index = ( index / num_tiles_x );

            const int tile_base_x = tile_x_index * TILESIZE_XY;
            const int tile_base_y = tile_y_index * TILESIZE_XY;

            for ( int y = tile_base_y; y < tile_base_y + TILESIZE_XY; y++ )
            for ( int x = tile_base_x; x < tile_base_x + TILESIZE_XY; x++ ) {
              vec3 running_color = vec3( 0. );      // initially zero, averages sample data
              for ( int s = 0; s < nsamples; s++ ) // get sample data (n samples)
                running_color += get_pathtrace_color_sample( x, y, id );
              running_color /= base_type( nsamples );  // sample averaging
              tonemap_and_gamma( running_color );     // tonemapping + gamma
              write( running_color, vec2( x, y ) );  // write final output values
            }
            tile_finish_counter.fetch_add( 1 );
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
  camera c; // generates view rays
  scene s; // holds all scene geometry + their associated materials
  int xdim = X_IMAGE_DIM, ydim = Y_IMAGE_DIM, nsamples = NUM_SAMPLES, bmax = MAX_BOUNCES;
  std::vector< unsigned char > bytes;        // image buffer for stb_image_write
  std::vector< std::shared_ptr< std::mt19937_64 > > gen; // PRNG states per thread
  void rng_seed() {
    std::random_device r;
    for ( int i = 0; i < NUM_THREADS; i++ ) {
      std::seed_seq s{ r(), r(), r(), r(), r(), r(), r(), r(), r() };
      gen.push_back( std::make_shared< std::mt19937_64 >( s ) );
    }
  }
  void tonemap_and_gamma( vec3& in ){
    in *= 0.6f; // function to tonemap color value in place

    // todo - add switch handling, for some different tonemap curves

    base_type a = 2.51;
    base_type b = 0.03;
    base_type c = 2.43;
    base_type d = 0.59;
    base_type e = 0.14;
    in = ( in * ( a * in + vec3( b ) ) ) / ( in * ( c * in + vec3( d ) ) + e ); // tonemap
    in.values[ 0 ] = std::pow( std::clamp( in.values[ 0 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA ); // gamma correct
    in.values[ 1 ] = std::pow( std::clamp( in.values[ 1 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA );
    in.values[ 2 ] = std::pow( std::clamp( in.values[ 2 ], 0.0, 1.0 ), 1.0 / IMAGE_GAMMA );
  }
  vec3 get_pathtrace_color_sample( const int x, const int y, const int id ){
    // throughput's initial value of 1. in each channel indicates that it is initially
    // capable of carrying all of the light intensity possible (100%), and it is reduced
    vec3 throughput = vec3( 1.0 ); // by the albedo of the material on each bounce
    vec3 current    = vec3( 0.0 ); // init to zero - initially no light present
    vec3 old_ro; // old_ro holds previous hit location, unitialized

    // get initial ray origin + ray direction from camera
    ray r = c.sample( vec2( x + rng( gen[ id ] ), y + rng( gen[ id ] ) ) );
    for ( int bounce = 0; bounce < MAX_BOUNCES; bounce++ ) {
      old_ro = r.origin; // cache old hit location
      hitrecord h = s.ray_query( r ); // get a new hit location (scene query)

      r.origin = r.origin + h.dtransit * r.direction + h.normal * HIT_EPSILON;
      r.direction = normalize( ( 1.0 + HIT_EPSILON ) * h.normal + random_unit_vector( gen[ id ] ) ); // diffuse reflection

      // the form is:
        // current    += throughput*current_emission // emission term
        // throughput *= albedo                     // diffuse absorption term

      if ( h.material_index == 0 ){ // the spheres
        // current += throughput * vec3(1.3, 1.2, 1.1);
        // throughput *= vec3(0.999);
        throughput *= palette( h.primitive_index * PALETTE_SCALAR );
      } else if ( h.material_index == 1 && h.front ) {
        current += throughput * BRIGHTNESS_SCALAR * vec3( h.uv.values[ 0 ], h.uv.values[ 1 ], 1 - h.uv.values[ 0 ] - h.uv.values[ 1 ] );
      } else if ( h.material_index == 1 && !h.front ) {
        current += throughput * ( palette( h.primitive_index * PALETTE_SCALAR ) * BRIGHTNESS_SCALAR );
      } else if ( h.material_index == 2 ) {
        r.direction = reflect( r.origin - old_ro, h.normal );
        throughput *= vec3( 0.89 );
      } else if ( h.material_index == 3 ) {
        throughput *= vec3( 0.999 );
      } else if ( h.dtransit == DMAX_TRAVEL ) {
        // current += throughput * 0.1 * vec3(0.918, 0.75, 0.6); // sky color and escape
        break; // escape
      }

      base_type p = std::max( throughput.values[ 0 ], std::max( throughput.values[ 1 ], throughput.values[ 2 ] ) );
      if ( rng( gen[ id ] ) > p ) // russian roulette termination check
        break;

      throughput *= 1.0 / p; // russian roulette compensation term

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

  // renderer r; r.render_and_save_to(filename);

  for ( size_t i = 72; i <= 100; i++ ) {
    std::stringstream s; s << "outputs/out" << std::to_string( i ) << ".png";
    renderer r; r.render_and_save_to( s.str() );
  }

  cout << "Total Render Time: " <<
    duration_cast< milliseconds >( high_resolution_clock::now() - tstart ).count() / 1000.0
        << " seconds" << endl; return 0;
}

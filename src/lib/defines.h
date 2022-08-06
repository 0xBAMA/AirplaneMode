#ifndef DEFINES
#define DEFINES

#include "vector.h"

// default types
#define baseType double
using vec2 = vector2< baseType >;
using vec3 = vector3< baseType >;

// render parameters
constexpr long long X_IMAGE_DIM = 1000;
constexpr long long Y_IMAGE_DIM = 1000;
constexpr long long TILESIZE_XY = 100;
constexpr long long MAX_BOUNCES = 15;
constexpr long long NUM_SAMPLES = 100;
constexpr long long NUM_THREADS = 17;

constexpr long long IMAGE_SEQUENCE_START_INDEX = 0;
constexpr long long IMAGE_SEQUENCE_END_INDEX = 10;

constexpr baseType  IMAGE_GAMMA = 2.2;
constexpr baseType  HIT_EPSILON = baseType( std::numeric_limits< baseType >::epsilon() * 10.0 );
constexpr baseType  DMAX_TRAVEL = baseType( std::numeric_limits< baseType >::max() ) / 10.0;
constexpr baseType  FIELD_OF_VIEW = 0.3;
constexpr baseType  PALETTE_SCALAR = 0.1618;
constexpr baseType  BRIGHTNESS_SCALAR = 16.18;
constexpr baseType  BOX_EPSILON = 0.1618;
constexpr long long REPORT_DELAY = 32; 				// reporter thread sleep duration, in ms
constexpr long long NUM_PRIMITIVES = 69;
constexpr long long PROGRESS_INDICATOR_STOPS = 69; // cli spaces to take up

#endif

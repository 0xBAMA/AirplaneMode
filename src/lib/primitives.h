#ifndef PRIMITIVES
#define PRIMITIVES

#include "vector.h"
#include "defines.h"

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
	aabb( vec3 mins, vec3 maxs, int m ) : bounds{ mins, maxs } {
		materialID = m;
		vec3 tempMins = vec3( std::min( bounds[ 0 ].values[ 0 ], bounds[ 1 ].values[ 0 ] ), std::min( bounds[ 0 ].values[ 1 ], bounds[ 1 ].values[ 1 ] ), std::min( bounds[ 0 ].values[ 2 ], bounds[ 1 ].values[ 2 ] ) );
		vec3 tempMaxs = vec3( std::max( bounds[ 0 ].values[ 0 ], bounds[ 1 ].values[ 0 ] ), std::max( bounds[ 0 ].values[ 1 ], bounds[ 1 ].values[ 1 ] ), std::max( bounds[ 0 ].values[ 2 ], bounds[ 1 ].values[ 2 ] ) );
		bounds[ 0 ] = tempMins;
		bounds[ 1 ] = tempMaxs;
  }

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

#endif

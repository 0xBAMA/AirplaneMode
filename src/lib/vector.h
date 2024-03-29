#ifndef VECTOR
#define VECTOR

// by jb

#include <cmath>
#include <random>
#include <chrono>

constexpr double pi = 3.1415926535897932384626433832795;

template < class T >
class vector2 {
public:
	T values [ 2 ];

	vector2 ()						{ values[ 0 ] = T( 0 );	values[ 1 ] = T( 0 );	}
	vector2 ( T val )			{ values[ 0 ] = val;		values[ 1 ] = val;		}
	vector2 ( T x, T y )	{ values[ 0 ] = x;			values[ 1 ] = y;			}

	// wip - I don't like the semantics of this
	const T x() const { return values[ 0 ]; }
	const T y() const { return values[ 1 ]; }

	//  +,- operators
	const vector2< T > operator+( const vector2< T >& other ) const { return vector2( this->values[ 0 ] + other.values[ 0 ], this->values[ 1 ] + other.values[ 1 ]); }
	const vector2< T >& operator+=( const vector2< T >& other ) { this->values[ 0 ] += other.values[ 0 ], this->values[ 1 ] += other.values[ 1 ]; }

	const vector2< T > operator-( const vector2< T >& other ) const { return vector2( this->values[ 0 ] - other.values[ 0 ], this->values[ 1 ] - other.values[ 1 ]); }
	const vector2< T >& operator-=( const vector2< T >& other ) { this->values[ 0 ] -= other.values[ 0 ], this->values[ 1 ] -= other.values[ 1 ]; }

	//  multiplication/division by a scalar of the same type - scale each element
	const vector2< T > operator*( const T&  scalar ) const { return vector2( this->values[ 0 ] * scalar, this->values[ 1 ] * scalar ); }
	const vector2< T >& operator*=( const T&  scalar ) { this->values[ 0 ] *= scalar, this->values[ 1 ] *= scalar; return *this; }

	const vector2< T > operator/( const T&  scalar ) const { return vector2( this->values[ 0 ] / scalar, this->values[ 1 ] / scalar ); }
	const vector2< T >& operator/=( const T&  scalar ) { this->values[ 0 ] /= scalar, this->values[ 1 ] /= scalar; return *this; }

	//  multiplication by a vector of the same type - elementwise multiply (hadamard product)
	const vector2< T > operator*( const vector2< T >& other ) const { return vector2( this->values[ 0 ] * other.values[ 0 ], this->values[ 1 ] * other.values[ 1 ] ); }
	const vector2< T >& operator*=( const vector2< T >& other ) { this->values[ 0 ] *= other.values[ 0 ], this->values[ 1 ] *= other.values[ 1 ]; return *this; }

	//  division by a vector of the same type - elementwise divide
	const vector2< T > operator/( const vector2< T >& other ) const { return vector2( this->values[ 0 ] / other.values[ 0 ], this->values[ 1 ] / other.values[ 1 ] ); }
	const vector2< T >& operator/=( const vector2< T >& other ) { this->values[ 0 ] /= other.values[ 0 ], this->values[ 1 ] /= other.values[ 1 ]; return *this; }
};

// addition/subtraction in the other order
template < class T >
const vector2< T > operator+( const T& scalar, const vector2< T >& vec ) { return vector2< T >(vec.values[ 0 ] + scalar, vec.values[ 1 ] + scalar ); }
template < class T >
const vector2< T > operator-( const T& scalar, const vector2< T >& vec ) { return vector2< T >(vec.values[ 0 ] - scalar, vec.values[ 1 ] - scalar ); }

// multiplication/division in the other order
template < class T >
const vector2< T > operator*( const T& scalar, const vector2< T >& vec ) { return vector2< T >(vec.values[ 0 ] * scalar, vec.values[ 1 ] * scalar ); }
template < class T >
const vector2< T > operator/( const T& scalar, const vector2< T >& vec ) { return vector2< T >(vec.values[ 0 ] / scalar, vec.values[ 1 ] / scalar ); }

template < class T >
T dot( vector2< T > v1, vector2< T > v2 ) { return v1.values[ 0 ] * v2.values[ 0 ] + v1.values[ 1 ] * v2.values[ 1 ]; }

template < class T > // squared length
T lengthSquared( vector2< T > v ) { return dot( v, v ); }

template < class T > // vector length
T length( vector2< T > v ) { return sqrt( lengthSquared( v ) ); }

template < class T > // return unit length colinear vector
vector2< T > normalize( vector2< T > in ) { T len = length( in ); return in / len; }

template < class T > // absolute value of all vector elements
vector2< T > abs( vector2< T > v ) { return vector2< T >( abs( v.values[ 0 ] ), abs( v.values[ 1 ] ) ); }

// rotate the 2d vector, by specified number degrees
template < class T >
vector2< T > rotate2D( vector2< T > in, double amt ) {
	double radians = ( pi / 180.0 ) * amt;
	return vector2< T >( in.values[ 0 ] * std::cos( radians ) - in.values[ 1 ] * std::sin( radians ), in.values[ 0 ] * std::sin( radians ) + in.values[ 1 ] * std::cos( radians ) );
}

// --------
// --------
// --------

template < class T >
class vector3 {
public:
	T values [ 3 ];

	vector3 ()								{ values[ 0 ]= T( 0 );	values[ 1 ] = T( 0 );	values[ 2 ]= T( 0 );	}
	vector3 ( T val )					{ values[ 0 ]= val;			values[ 1 ] = val;		values[ 2 ]= val;			}
	vector3 ( T x, T y, T z )	{ values[ 0 ]= x;				values[ 1 ] = y;			values[ 2 ]= z;				}

	// wip - I don't like the semantics of this
	const T x() const { return values[ 0 ]; }
	const T y() const { return values[ 1 ]; }
	const T z() const { return values[ 2 ]; }

	//  +,- operators
	const vector3< T > operator+( const vector3< T >& other ) const { return vector3( this->values[ 0 ] + other.values[ 0 ], this->values[ 1 ] + other.values[ 1 ], this->values[ 2 ] + other.values[ 2 ] ); }
	const vector3< T >& operator+=( const vector3< T >& other ) { this->values[ 0 ] += other.values[ 0 ], this->values[ 1 ] += other.values[ 1 ], this->values[ 2 ] += other.values[ 2 ]; return *this; }

	const vector3< T > operator-( const vector3< T >& other ) const { return vector3( this->values[ 0 ] - other.values[ 0 ], this->values[ 1 ] - other.values[ 1 ], this->values[ 2 ] - other.values[ 2 ]); }
	const vector3< T >& operator-=( const vector3< T >& other ) { this->values[ 0 ] -= other.values[ 0 ], this->values[ 1 ] -= other.values[ 1 ], this->values[ 2 ] -= other.values[ 2 ]; return *this; }

	//  multiplication/division by a scalar of the same type - scale each element
	const vector3< T > operator*( const T& scalar ) const { return vector3( this->values[ 0 ] * scalar, this->values[ 1 ] * scalar, this->values[ 2 ] * scalar ); }
	const vector3< T >& operator*=( const T& scalar ) { this->values[ 0 ] *= scalar, this->values[ 1 ] *= scalar, this->values[ 2 ] *= scalar; return *this; }

	const vector3< T > operator/( const T& scalar ) const { return vector3( this->values[ 0 ] / scalar, this->values[ 1 ] / scalar, this->values[ 2 ] / scalar ); }
	const vector3< T >& operator/=( const T& scalar ) { this->values[ 0 ] /= scalar, this->values[ 1 ] /= scalar, this->values[ 2 ] /= scalar; return *this; }

	//  multiplication by a vector of the same type - elementwise multiply (hadamard product)
	const vector3< T > operator*( const vector3< T >& other ) const { return vector3( this->values[ 0 ] * other.values[ 0 ], this->values[ 1 ] * other.values[ 1 ], this->values[ 2 ] * other.values[ 2 ] ); }
	const vector3< T >& operator*=( const vector3< T >& other ) { this->values[ 0 ] *= other.values[ 0 ], this->values[ 1 ] *= other.values[ 1 ], this->values[ 2 ] *= other.values[ 2 ]; return *this; }

	//  division by a vector of the same type - elementwise divide
	const vector3< T > operator/( const vector3< T >& other ) const { return vector3( this->values[ 0 ] / other.values[ 0 ], this->values[ 1 ] / other.values[ 1 ], this->values[ 2 ] / other.values[ 2 ]); }
	const vector3< T >& operator/=( const vector3< T >& other ) { this->values[ 0 ] /= other.values[ 0 ], this->values[ 1 ] /= other.values[ 1 ], this->values[ 2 ] /= other.values[ 2 ]; return *this; }
};

// addition/subtraction in the other order
template < class T >
const vector3< T > operator+( const T& scalar, const vector3< T >& vec ) { return vector3< T >( vec.values[ 0 ] + scalar, vec.values[ 1 ] + scalar, vec.values[ 2 ] + scalar ); }
template < class T >
const vector3< T > operator-( const T& scalar, const vector3< T >& vec ) { return vector3< T >( vec.values[ 0 ] - scalar, vec.values[ 1 ] - scalar, vec.values[ 2 ] - scalar ); }

// multiplication/division in the other order
template < class T >
const vector3< T > operator*( const T& scalar, const vector3< T >& vec ) { return vector3< T >( vec.values[ 0 ] * scalar, vec.values[ 1 ] * scalar, vec.values[ 2 ] * scalar ); }
template < class T >
const vector3< T > operator/( const T& scalar, const vector3< T >& vec ) { return vector3< T >( vec.values[ 0 ] / scalar, vec.values[ 1 ] / scalar, vec.values[ 2 ] / scalar ); }

template < class T >
const T dot( const vector3< T > v1, const vector3< T > v2 ) { return v1.values[ 0 ] * v2.values[ 0 ] + v1.values[ 1 ] * v2.values[ 1 ] + v1.values[ 2 ] * v2.values[ 2 ]; }

template < class T > // squared length
const T lengthSquared( const vector3< T > v ) { return dot( v, v ); }

template < class T > // vector length
const T length( const vector3< T > v ) { return sqrt( lengthSquared( v ) ); }

template < class T > // return unit length colinear vector
const vector3< T > normalize( const vector3< T > in ) { T len = length( in ); return in / len; }

template < class T > // elementwise absolute value
vector3< T > abs( const vector3< T > v ) { return vector3< T >( abs( v.values[ 0 ] ), abs( v.values[ 1 ] ), abs( v.values[ 2 ] ) ); }

template < class T > //  cross product
const vector3< T > cross( const vector3< T > a, const vector3< T > b ) {
	vector3< T > product;
	product.values[ 0 ] =    a.values[ 1 ] * b.values[ 2 ] - a.values[ 2 ] * b.values[ 1 ];
	product.values[ 1 ] = -( a.values[ 0 ] * b.values[ 2 ] - a.values[ 2 ] * b.values[ 0 ] );
	product.values[ 2 ] =    a.values[ 0 ] * b.values[ 1 ] - a.values[ 1 ] * b.values[ 0 ];
	return product;
}

template < class T > // reflect function
const vector3< T > reflect( const vector3< T > i, const vector3< T > n ){
	return i - 2.0 * dot( n, i ) * n;
}

template < class T >
const vector3< T > mix( const vector3< T > x, const vector3< T > y, const T a ) {
	return x * ( 1.0 - a ) + y * a;
}

template < class T > // Rodrigues rotation formula via https://suricrasia.online/demoscene/functions/
const vector3< T > erot( const vector3< T > point, vector3< T > axis, const T amount ) {
	axis = normalize( axis );
	return mix( dot( axis, point ) * axis, point, cos( amount ) ) + cross( axis, point ) * sin( amount );
}


// refract function would make sense to add
// some of the other demoscene type rotation implementations might make sense as well


// --------
// --------
// --------

template < class T >
class vector4 {
public:
	T values [ 4 ];

	vector4 ()						{ values[ 0 ]= T( 0 );	values[ 1 ] = T( 0 );	values[ 2 ]= T( 0 );	values[ 3 ]= T( 0 ); }
	vector4 ( T val )				{ values[ 0 ]= val;		values[ 1 ] = val;		values[ 2 ]= val;		values[ 3 ]= val;	}
	vector4 ( T x, T y, T z, T w )	{ values[ 0 ]= x;		values[ 1 ] = y;		values[ 2 ]= z;			values[ 3 ]= w;	}

	// wip - I don't like the semantics of this
	const T x() const { return values[ 0 ]; }
	const T y() const { return values[ 1 ]; }
	const T z() const { return values[ 2 ]; }
	const T w() const { return values[ 3 ]; }

	//  +,- operators
	const vector4< T > operator+( const vector4< T >& other ) const { return vector4( this->values[ 0 ] + other.values[ 0 ], this->values[ 1 ] + other.values[ 1 ], this->values[ 2 ] + other.values[ 2 ], this->values[ 3 ] + other.values[ 3 ] ); }
	const vector4< T >& operator+=( const vector4< T >& other ) { this->values[ 0 ] += other.values[ 0 ], this->values[ 1 ] += other.values[ 1 ], this->values[ 2 ] += other.values[ 2 ], this->values[ 3 ] += other.values[ 3 ]; return *this; }

	const vector4< T > operator-( const vector4< T >& other ) const { return vector4( this->values[ 0 ] - other.values[ 0 ], this->values[ 1 ] - other.values[ 1 ], this->values[ 2 ] - other.values[ 2 ], this->values[ 3 ] - other.values[ 3 ] ); }
	const vector4< T >& operator-=( const vector4< T >& other ) { this->values[ 0 ] -= other.values[ 0 ], this->values[ 1 ] -= other.values[ 1 ], this->values[ 2 ] -= other.values[ 2 ], this->values[ 3 ] -= other.values[ 3 ]; return *this; }

	//  multiplication/division by a scalar of the same type - scale each element
	const vector4< T > operator*( const T& scalar ) const { return vector4( this->values[ 0 ] * scalar, this->values[ 1 ] * scalar, this->values[ 2 ] * scalar, this->values[ 3 ] * scalar ); }
	const vector4< T >& operator*=( const T& scalar ) { this->values[ 0 ] *= scalar, this->values[ 1 ] *= scalar, this->values[ 2 ] *= scalar, this->values[ 3 ] *= scalar; return *this; }

	const vector4< T > operator/( const T& scalar ) const { return vector4( this->values[ 0 ] / scalar, this->values[ 1 ] / scalar, this->values[ 2 ] / scalar ); }
	const vector4< T >& operator/=( const T& scalar ) { this->values[ 0 ] /= scalar, this->values[ 1 ] /= scalar, this->values[ 2 ] /= scalar, this->values[ 3 ] /= scalar; return *this; }

	//  multiplication by a vector of the same type - elementwise multiply (hadamard product)
	const vector4< T > operator*( const vector4< T >& other ) const { return vector4( this->values[ 0 ] * other.values[ 0 ], this->values[ 1 ] * other.values[ 1 ], this->values[ 2 ] * other.values[ 2 ], this->values[ 3 ] * other.values[ 3 ] ); }
	const vector4< T >& operator*=( const vector4< T >& other ) { this->values[ 0 ] *= other.values[ 0 ], this->values[ 1 ] *= other.values[ 1 ], this->values[ 2 ] *= other.values[ 2 ], this->values[ 3 ] *= other.values[ 3 ]; return *this; }

	//  division by a vector of the same type - elementwise divide
	const vector4< T > operator/( const vector4< T >& other ) const { return vector4( this->values[ 0 ] / other.values[ 0 ], this->values[ 1 ] / other.values[ 1 ], this->values[ 2 ] / other.values[ 2 ], this->values[ 3 ] / other.values[ 3 ]); }
	const vector4< T >& operator/=( const vector4< T >& other ) { this->values[ 0 ] /= other.values[ 0 ], this->values[ 1 ] /= other.values[ 1 ], this->values[ 2 ] /= other.values[ 2 ], this->values[ 3 ] /= other.values[ 3 ]; return *this; }
};

// addition/subtraction in the other order
template < class T >
const vector4< T > operator+( const T& scalar, const vector4< T >& vec ) { return vector4< T >( vec.values[ 0 ] + scalar, vec.values[ 1 ] + scalar, vec.values[ 2 ] + scalar, vec.values[ 3 ] + scalar ); }
template < class T >
const vector4< T > operator-( const T& scalar, const vector4< T >& vec ) { return vector4< T >( vec.values[ 0 ] - scalar, vec.values[ 1 ] - scalar, vec.values[ 2 ] - scalar, vec.values[ 3 ] - scalar ); }

// multiplication/division in the other order
template < class T >
const vector4< T > operator*( const T& scalar, const vector4< T >& vec ) { return vector4< T >( vec.values[ 0 ] * scalar, vec.values[ 1 ] * scalar, vec.values[ 2 ] * scalar ); }
template < class T >
const vector4< T > operator/( const T& scalar, const vector4< T >& vec ) { return vector4< T >( vec.values[ 0 ] / scalar, vec.values[ 1 ] / scalar, vec.values[ 2 ] / scalar, vec.values[ 3 ] / scalar ); }

template < class T >
const T dot( const vector4< T > v1, const vector4< T > v2 ) { return v1.values[ 0 ] * v2.values[ 0 ] + v1.values[ 1 ] * v2.values[ 1 ] + v1.values[ 2 ] * v2.values[ 2 ] + v1.values[ 3 ] * v2.values[ 3 ]; }

template < class T > // elementwise absolute value
vector4< T > abs( const vector4< T > v ) { return vector4< T >( abs( v.values[ 0 ] ), abs( v.values[ 1 ] ), abs( v.values[ 2 ] ), abs( v.values[ 3 ] ) ); }

template < class T >
const vector4< T > mix( const vector4< T > x, const vector4< T > y, const T a ) {
	return x * ( 1.0 - a ) + y * a;
}

template < class T > // Rodrigues rotation formula via https://suricrasia.online/demoscene/functions/
const vector4< T > erot( const vector4< T > point, vector4< T > axis, const T amount ) {
	axis = normalize( axis );
	return mix( dot( axis, point ) * axis, point, cos( amount ) ) + cross( axis, point ) * sin( amount );
}




#endif

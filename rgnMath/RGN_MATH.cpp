#include "RGN_MATH.h"

float Dot(float4 &first, float4 &second)
{
	return (first.vector.x * second.vector.x +
		first.vector.y * second.vector.y +
		first.vector.z * second.vector.z);
}

float4 Cross(float4 &first, float4 &second)
{
	float4 res;

	res.vector.x = first.vector.y * second.vector.z - first.vector.z * second.vector.y;
	res.vector.y = first.vector.z * second.vector.x - first.vector.x * second.vector.z;
	res.vector.z = first.vector.x * second.vector.y - first.vector.y * second.vector.x;

	return res;
}

void TranslateToMat(float3 &pos, float4x4 &mat)
{
	mat.Identity();

	/*

	00		01		02		03
	10		11		12		13
	20		21		22		23
	pos.x	pos.y	pos.z	33

	*/

	mat.m30 = pos.x;
	mat.m31 = pos.y;
	mat.m32 = pos.z;
}

void ScaleToMat(float3 &scale, float4x4 &mat)
{
	mat.Identity();

	/*

	scale.x	01		02		03
	10		scale.y	12		13
	20		21		scale.z	23
	30		31		32		33

	*/

	mat.m00 = scale.x;
	mat.m11 = scale.y;
	mat.m22 = scale.z;
}

void RotToMat(float3 axis, float angle, float4x4 &mat)
{
	// need to normalize axis before rotation

	float cosA		= cos(angle);
	float invCosA	= 1.0f - cosA;

	float sinA		= sin(angle);

	float x = axis.x;
	float y = axis.y;
	float z = axis.z;

	//-> setting
	mat.Set(
		cosA + invCosA * (x * x),	 invCosA * x * y + sinA * z, invCosA * x * z - sinA * y,	0,
		invCosA * x * y - sinA*z,	 cosA + invCosA * (y * y),	 invCosA * y * z + sinA * x,	0,
		invCosA * x * z + sinA * y,	 invCosA * y * z - sinA * x, cosA + invCosA * (z * z),		0,
		0, 0, 0, 1);
}

#pragma region "###	FLOAT4		###"

/*
	FLOAT4X4 is main class for working with vectors and points

	all operations with vectors use vector3 component
	for transformations ( multiplication on marix ) it use W component
	which is equal 1.0f for points and 0.0f for vectors

*/

float4::float4()
{
	vector.x = vector.y = vector.z = w = 0.0f;
}

float4::float4(float iX, float iY, float iZ,
	float iW)
{
	vector.x = iX;
	vector.y = iY;
	vector.z = iZ;
	w		 = iW;
}

float4::float4(float3 &vec,
	float iW)
{
	vector.x = vec.x;
	vector.y = vec.y;
	vector.z = vec.z;
	w		 = iW;
}

void float4::Set(float iX, float iY, float iZ, 
	float iW)
{
	vector.x = iX;
	vector.y = iY;
	vector.z = iZ;
	w		 = iW;
}

float  float4::Length()
{
	return sqrt(
		vector.x * vector.x +
		vector.y * vector.y +
		vector.z * vector.z
		);
}

void float4::Normalize()
{
	float L = 1.0f / Length();

	vector.x *= L;
	vector.y *= L;
	vector.z *= L;
}

float float4::Dot(const float4 &vec)
{
	return (
		vector.x * vec.vector.x +
		vector.y * vec.vector.y +
		vector.z * vec.vector.z
		);
}

float4 float4::Cross(const float4 &vec)
{
	float4 res;

	res.vector.x = vector.y * vec.vector.z - vector.z * vec.vector.y;
	res.vector.y = vector.z * vec.vector.x - vector.x * vec.vector.z;
	res.vector.z = vector.x * vec.vector.y - vector.y * vec.vector.x;

	return res;
}

void float4::operator+=(const float4 &vec)
{
	vector.x += vec.vector.x;
	vector.y += vec.vector.y;
	vector.z += vec.vector.z;
}

float4 float4::operator+(const float4 &vec)
{
	float4 res;

	res.vector.x = vector.x + vec.vector.x;
	res.vector.y = vector.y + vec.vector.y;
	res.vector.z = vector.z + vec.vector.z;

	return res;
}

void float4::operator-=(const float4 &vec)
{
	vector.x -= vec.vector.x;
	vector.y -= vec.vector.y;
	vector.z -= vec.vector.z;
}

float4 float4::operator-(const float4 &vec)
{
	float4 res;

	res.vector.x = vector.x - vec.vector.x;
	res.vector.y = vector.y - vec.vector.y;
	res.vector.z = vector.z - vec.vector.z;

	return res;
}

void float4::operator*=(float v)
{
	vector.x *= v;
	vector.y *= v;
	vector.z *= v;
}

float4 float4::operator*(float v)
{
	float4 res;

	res.vector.x = vector.x * v;
	res.vector.y = vector.y * v;
	res.vector.z = vector.z * v;

	return res;
}

void float4::operator*=(const float4x4 &mat)
{
	float4 res = *this;

	res.vector.x = vector.x * mat.m00 + vector.y * mat.m10 + vector.z * mat.m20 + w * mat.m30;
	res.vector.y = vector.x * mat.m01 + vector.y * mat.m11 + vector.z * mat.m21 + w * mat.m31;
	res.vector.z = vector.x * mat.m02 + vector.y * mat.m12 + vector.z * mat.m22 + w * mat.m32;
	res.w		 = vector.x * mat.m03 + vector.y * mat.m13 + vector.z * mat.m23 + w * mat.m33;

	*this = res;
}

float4 float4::operator*(const float4x4 &mat)
{
	float4 res;

	res.vector.x = vector.x * mat.m00 + vector.y * mat.m10 + vector.z * mat.m20 + w * mat.m30;
	res.vector.y = vector.x * mat.m01 + vector.y * mat.m11 + vector.z * mat.m21 + w * mat.m31;
	res.vector.z = vector.x * mat.m02 + vector.y * mat.m12 + vector.z * mat.m22 + w * mat.m32;
	res.w		 = vector.x * mat.m03 + vector.y * mat.m13 + vector.z * mat.m23 + w * mat.m33;

	return res;
}

bool float4::operator==(const float4 &vec)
{
	if ( ( (vector.x - vec.vector.x) < F_EPSILON ) &&
		 ( (vector.y - vec.vector.y) < F_EPSILON ) &&
		 ( (vector.z - vec.vector.z) < F_EPSILON ))
		return true;

	return false;
}

bool float4::operator!=(const float4 &vec)
{
	if ( ( (vector.x - vec.vector.x) < F_EPSILON ) &&
		 ( (vector.y - vec.vector.y) < F_EPSILON ) &&
		 ( (vector.z - vec.vector.z) < F_EPSILON ))
		return false;

	return true;
}

#pragma endregion

#pragma region "###	FLOAT4X4	###"

// FLOAT2x2

float2x2::float2x2()
{
	m00 = m01 = m10 = m11 = 0.0f;
}

float2x2::~float2x2()
{
}

void float2x2::Set(
	float iM00, float iM01,
	float iM10, float iM11
	)
{
	m00 = iM00;
	m01 = iM01;
	m10 = iM10;
	m11 = iM11;
}

float float2x2::Det()
{
	return m00 * m11 - m01 * m10;
}

// FLOAT3x3

float3x3::float3x3()
{
	m00 = m01 = m02 = 0.0f;
	m10 = m11 = m12 = 0.0f;
	m20 = m21 = m22 = 0.0f;
}

float3x3::~float3x3()
{
}

void float3x3::Set(
	float iM00, float iM01, float iM02,
	float iM10, float iM11, float iM12,
	float iM20, float iM21, float iM22
	)
{
	m00 = iM00;
	m01 = iM01;
	m02 = iM02;

	m10 = iM10;
	m11 = iM11;
	m12 = iM12;

	m20 = iM20;
	m21 = iM21;
	m22 = iM22;
}

float float3x3::Det()
{
	float2x2 temp;			// for computing of a minor

	float xD1, xD2, xD3;	// values of first row
	float d1, d2, d3;		// value of computing minors

	xD1 = m00;
	xD2 = -m01;				// change a sign ( as at the school =) )
	xD3 = m02;

	// set a minor ( 00 )
	temp.m00 = m11;
	temp.m01 = m12;
	temp.m10 = m21;
	temp.m11 = m22;

	// compute a minor
	d1 = temp.Det();

	// minor 01
	temp.m00 = m10;
	temp.m01 = m12;
	temp.m10 = m20;
	temp.m11 = m22;

	d2 = temp.Det();

	// minor 02
	temp.m00 = m10;
	temp.m01 = m11;
	temp.m10 = m20;
	temp.m11 = m21;

	d3 = temp.Det();

	// compute a determinant 3x3
	return xD1 * d1 + xD2 * d2 + xD3 * d3;
}

// FLOAT4x4

float4x4::float4x4()
{
	m01 = m02 = m03 = 0.0f; // row 0
	m10 = m12 = m13 = 0.0f; // 1 
	m20 = m21 = m23 = 0.0f; // 2 
	m30 = m31 = m32 = 0.0f; // 3

	m00 = m11 = m22 = m33 = 0.0f; // main diagonal
}

float4x4::float4x4(
	float iM00, float iM01, float iM02, float iM03,
	float iM10, float iM11, float iM12, float iM13,
	float iM20, float iM21, float iM22, float iM23,
	float iM30, float iM31, float iM32, float iM33
	)
{
	m00 = iM00;
	m01 = iM01;
	m02 = iM02;
	m03 = iM03;

	m10 = iM10;
	m11 = iM11;
	m12 = iM12;
	m13 = iM13;

	m20 = iM20;
	m21 = iM21;
	m22 = iM22;
	m23 = iM23;

	m30 = iM30;
	m31 = iM31;
	m32 = iM32;
	m33 = iM33;
}

void float4x4::Identity()
{
	m01 = m02 = m03 = 0.0f;
	m10 = m12 = m13 = 0.0f;
	m20 = m21 = m23 = 0.0f;
	m30 = m31 = m32 = 0.0f;

	m00 = m11 = m22 = m33 = 1.0f; // main diagonal
}

void float4x4::Zero()
{
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			mat[row][col] = 0.0f;
}

void float4x4::Set(
	float iM00, float iM01, float iM02, float iM03,
	float iM10, float iM11, float iM12, float iM13,
	float iM20, float iM21, float iM22, float iM23,
	float iM30, float iM31, float iM32, float iM33
	)
{
	m00 = iM00;
	m01 = iM01;
	m02 = iM02;
	m03 = iM03;

	m10 = iM10;
	m11 = iM11;
	m12 = iM12;
	m13 = iM13;

	m20 = iM20;
	m21 = iM21;
	m22 = iM22;
	m23 = iM23;

	m30 = iM30;
	m31 = iM31;
	m32 = iM32;
	m33 = iM33;
}

float4x4 float4x4::Transpose()
{
	float4x4 res;

	res.m01 = m10;
	res.m02 = m20;
	res.m03 = m30;

	res.m10 = m01;
	res.m12 = m21;
	res.m13 = m31;

	res.m20 = m02;
	res.m21 = m12;
	res.m23 = m32;

	res.m30 = m03;
	res.m31 = m13;
	res.m32 = m23;

	// just copy the main diagonal

	res.m00 = m00;
	res.m11 = m11;
	res.m22 = m22;
	res.m33 = m33;

	return res;
}

bool float4x4::operator==(const float4x4 &iMat)
{
	// test each elements of matrix
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			if ((mat[row][col] - iMat.mat[row][col]) > F_EPSILON) // with considering an error
				return false;

	return true;
}

bool float4x4::operator!=(const float4x4 &iMat)
{
	// test each elements of matrix
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			if ((mat[row][col] - iMat.mat[row][col]) > F_EPSILON) // with considering an error
				return true;

	return false;
}

float4x4 float4x4::operator+(const float4x4 &iMat)
{
	float4x4 res;

	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			res.mat[col][row] = mat[col][row] + iMat.mat[col][row];

	return res;
}

void float4x4::operator+=(const float4x4 &iMat)
{
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			mat[col][row] += iMat.mat[col][row];
}

float4x4 float4x4::operator-(const float4x4 &iMat)
{
	float4x4 res;

	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			res.mat[col][row] = mat[col][row] - iMat.mat[col][row];

	return res;
}

void float4x4::operator-=(const float4x4 &iMat)
{
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			mat[col][row] -= iMat.mat[col][row];
}

float4x4 float4x4::operator*(float k)
{
	float4x4 res;

	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			res.mat[col][row] = mat[col][row] * k;

	return res;
}

void float4x4::operator*=(float k)
{
	for (int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			mat[col][row] *= k;
}

float4x4 float4x4::operator*(const float4x4 &iMat)
{
	float4x4 res;

	// cross product each row of this->matrix with each column inputted matrix

	// a little confusedly ...  
	res.m00 = this->m00 * iMat.m00 + this->m01 * iMat.m10 + this->m02 * iMat.m20 + this->m03 * iMat.m30;
	res.m01 = this->m00 * iMat.m01 + this->m01 * iMat.m11 + this->m02 * iMat.m21 + this->m03 * iMat.m31;
	res.m02 = this->m00 * iMat.m02 + this->m01 * iMat.m12 + this->m02 * iMat.m22 + this->m03 * iMat.m32;
	res.m03 = this->m00 * iMat.m03 + this->m01 * iMat.m13 + this->m02 * iMat.m23 + this->m03 * iMat.m33;

	res.m10 = this->m10 * iMat.m00 + this->m11 * iMat.m10 + this->m12 * iMat.m20 + this->m13 * iMat.m30;
	res.m11 = this->m10 * iMat.m01 + this->m11 * iMat.m11 + this->m12 * iMat.m21 + this->m13 * iMat.m31;
	res.m12 = this->m10 * iMat.m02 + this->m11 * iMat.m12 + this->m12 * iMat.m22 + this->m13 * iMat.m32;
	res.m13 = this->m10 * iMat.m03 + this->m11 * iMat.m13 + this->m12 * iMat.m23 + this->m13 * iMat.m33;

	res.m20 = this->m20 * iMat.m00 + this->m21 * iMat.m10 + this->m22 * iMat.m20 + this->m23 * iMat.m30;
	res.m21 = this->m20 * iMat.m01 + this->m21 * iMat.m11 + this->m22 * iMat.m21 + this->m23 * iMat.m31;
	res.m22 = this->m20 * iMat.m02 + this->m21 * iMat.m12 + this->m22 * iMat.m22 + this->m23 * iMat.m32;
	res.m23 = this->m20 * iMat.m03 + this->m21 * iMat.m13 + this->m22 * iMat.m23 + this->m23 * iMat.m33;

	res.m30 = this->m30 * iMat.m00 + this->m31 * iMat.m10 + this->m32 * iMat.m20 + this->m33 * iMat.m30;
	res.m31 = this->m30 * iMat.m01 + this->m31 * iMat.m11 + this->m32 * iMat.m21 + this->m33 * iMat.m31;
	res.m32 = this->m30 * iMat.m02 + this->m31 * iMat.m12 + this->m32 * iMat.m22 + this->m33 * iMat.m32;
	res.m33 = this->m30 * iMat.m03 + this->m31 * iMat.m13 + this->m32 * iMat.m23 + this->m33 * iMat.m33;

	return res;
}

void float4x4::operator*=(const float4x4 &iMat)
{
	float4x4 temp;

	// cross product each row of this->matrix with each column inputted matrix

	temp.m00 = this->m00 * iMat.m00 + this->m01 * iMat.m10 + this->m02 * iMat.m20 + this->m03 * iMat.m30;
	temp.m01 = this->m00 * iMat.m01 + this->m01 * iMat.m11 + this->m02 * iMat.m21 + this->m03 * iMat.m31;
	temp.m02 = this->m00 * iMat.m02 + this->m01 * iMat.m12 + this->m02 * iMat.m22 + this->m03 * iMat.m32;
	temp.m03 = this->m00 * iMat.m03 + this->m01 * iMat.m13 + this->m02 * iMat.m23 + this->m03 * iMat.m33;

	temp.m10 = this->m10 * iMat.m00 + this->m11 * iMat.m10 + this->m12 * iMat.m20 + this->m13 * iMat.m30;
	temp.m11 = this->m10 * iMat.m01 + this->m11 * iMat.m11 + this->m12 * iMat.m21 + this->m13 * iMat.m31;
	temp.m12 = this->m10 * iMat.m02 + this->m11 * iMat.m12 + this->m12 * iMat.m22 + this->m13 * iMat.m32;
	temp.m13 = this->m10 * iMat.m03 + this->m11 * iMat.m13 + this->m12 * iMat.m23 + this->m13 * iMat.m33;

	temp.m20 = this->m20 * iMat.m00 + this->m21 * iMat.m10 + this->m22 * iMat.m20 + this->m23 * iMat.m30;
	temp.m21 = this->m20 * iMat.m01 + this->m21 * iMat.m11 + this->m22 * iMat.m21 + this->m23 * iMat.m31;
	temp.m22 = this->m20 * iMat.m02 + this->m21 * iMat.m12 + this->m22 * iMat.m22 + this->m23 * iMat.m32;
	temp.m23 = this->m20 * iMat.m03 + this->m21 * iMat.m13 + this->m22 * iMat.m23 + this->m23 * iMat.m33;

	temp.m30 = this->m30 * iMat.m00 + this->m31 * iMat.m10 + this->m32 * iMat.m20 + this->m33 * iMat.m30;
	temp.m31 = this->m30 * iMat.m01 + this->m31 * iMat.m11 + this->m32 * iMat.m21 + this->m33 * iMat.m31;
	temp.m32 = this->m30 * iMat.m02 + this->m31 * iMat.m12 + this->m32 * iMat.m22 + this->m33 * iMat.m32;
	temp.m33 = this->m30 * iMat.m03 + this->m31 * iMat.m13 + this->m32 * iMat.m23 + this->m33 * iMat.m33;

	*this = temp;
}

float float4x4::Det()
{
	float3x3 temp; // for computing of a minor

	float xD1, xD2, xD3, xD4;	// value of elements first row
	float d1, d2, d3, d4;		// value of minors

	xD1 = m00;
	xD2 = -m01; // change a sign
	xD3 = m02;
	xD4 = -m03; // here too

	// set a minor3x3 ( 00 )
	temp.m00 = m11;
	temp.m01 = m12;
	temp.m02 = m13;
	
	temp.m10 = m21;
	temp.m11 = m22;
	temp.m12 = m23;

	temp.m20 = m31;
	temp.m21 = m32;
	temp.m22 = m33;


	// compute minor
	d1 = temp.Det();

	// minor 01
	temp.m00 = m10;
	temp.m01 = m12;
	temp.m02 = m13;

	temp.m10 = m20;
	temp.m11 = m22;
	temp.m12 = m23;

	temp.m20 = m30;
	temp.m21 = m32;
	temp.m22 = m33;

	d2 = temp.Det();

	// minor 02
	temp.m00 = m10;
	temp.m01 = m11;
	temp.m02 = m13;

	temp.m10 = m20;
	temp.m11 = m21;
	temp.m12 = m23;

	temp.m20 = m30;
	temp.m21 = m31;
	temp.m22 = m33;

	d3 = temp.Det();

	// minor 03
	temp.m00 = m10;
	temp.m01 = m11;
	temp.m02 = m12;

	temp.m10 = m20;
	temp.m11 = m21;
	temp.m12 = m22;

	temp.m20 = m30;
	temp.m21 = m31;
	temp.m22 = m32;

	d4 = temp.Det();

	// compute a determinant 4x4
	return xD1 * d1 + xD2 * d2 + xD3 * d3 + xD4 * d4;
}

float4x4 float4x4::Inverse()
{
	float4x4 mins;	// matrix of value minors for each element of matrix
	float3x3 temp;	// for computing of minors
	float4x4 res;	// result matrix

	float det = Det(); // found a determinant

	if (!det) // if he is equal zero then inverse matrix is not defined
		return res;

	// form minors for each element

	// row 0
	// 0
	temp.m00 = m11; temp.m01 = m12; temp.m02 = m13;
	temp.m10 = m21; temp.m11 = m22; temp.m12 = m23;
	temp.m20 = m31; temp.m21 = m32; temp.m22 = m33;

	mins.m00 = temp.Det();

	// 1
	temp.m00 = m10; temp.m01 = m12; temp.m02 = m13;
	temp.m10 = m20; temp.m11 = m22; temp.m12 = m23;
	temp.m20 = m30; temp.m21 = m32; temp.m22 = m33;

	mins.m01 = temp.Det();

	// 2
	temp.m00 = m10; temp.m01 = m11; temp.m02 = m13;
	temp.m10 = m20; temp.m11 = m21; temp.m12 = m23;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m33;

	mins.m02 = temp.Det();

	// 3
	temp.m00 = m10; temp.m01 = m11; temp.m02 = m12;
	temp.m10 = m20; temp.m11 = m21; temp.m12 = m22;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m32;

	mins.m03 = temp.Det();

	// row 1
	// 0
	temp.m00 = m01; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m21; temp.m11 = m22; temp.m12 = m23;
	temp.m20 = m31; temp.m21 = m32; temp.m22 = m33;

	mins.m10 = temp.Det();

	// 1
	temp.m00 = m00; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m20; temp.m11 = m22; temp.m12 = m23;
	temp.m20 = m30; temp.m21 = m32; temp.m22 = m33;

	mins.m11 = temp.Det();

	// 2
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m03;
	temp.m10 = m20; temp.m11 = m21; temp.m12 = m23;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m33;

	mins.m12 = temp.Det();

	// 3
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m02;
	temp.m10 = m20; temp.m11 = m21; temp.m12 = m22;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m32;

	mins.m13 = temp.Det();

	// row 2
	// 0
	temp.m00 = m01; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m11; temp.m11 = m12; temp.m12 = m13;
	temp.m20 = m31; temp.m21 = m32; temp.m22 = m33;

	mins.m20 = temp.Det();

	// 1
	temp.m00 = m00; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m10; temp.m11 = m12; temp.m12 = m13;
	temp.m20 = m30; temp.m21 = m32; temp.m22 = m33;

	mins.m21 = temp.Det();

	// 2
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m03;
	temp.m10 = m10; temp.m11 = m11; temp.m12 = m13;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m33;

	mins.m22 = temp.Det();

	// 3
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m02;
	temp.m10 = m10; temp.m11 = m11; temp.m12 = m12;
	temp.m20 = m30; temp.m21 = m31; temp.m22 = m32;

	mins.m23 = temp.Det();

	// row 3
	// 0
	temp.m00 = m01; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m11; temp.m11 = m12; temp.m12 = m13;
	temp.m20 = m21; temp.m21 = m22; temp.m22 = m23;

	mins.m30 = temp.Det();

	// 1
	temp.m00 = m00; temp.m01 = m02; temp.m02 = m03;
	temp.m10 = m10; temp.m11 = m12; temp.m12 = m13;
	temp.m20 = m20; temp.m21 = m22; temp.m22 = m23;

	mins.m31 = temp.Det();

	// 2
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m03;
	temp.m10 = m10; temp.m11 = m11; temp.m12 = m13;
	temp.m20 = m20; temp.m21 = m21; temp.m22 = m23;

	mins.m32 = temp.Det();

	// 3
	temp.m00 = m00; temp.m01 = m01; temp.m02 = m02;
	temp.m10 = m10; temp.m11 = m11; temp.m12 = m12;
	temp.m20 = m20; temp.m21 = m21; temp.m22 = m22;

	mins.m33 = temp.Det();

	// seek cofactors
	mins.m01 = -mins.m01;
	mins.m03 = -mins.m03;

	mins.m10 = -mins.m10;
	mins.m12 = -mins.m12;

	mins.m21 = -mins.m21;
	mins.m23 = -mins.m23;

	mins.m30 = -mins.m30;
	mins.m32 = -mins.m32;

	//
	mins = mins.Transpose();

	// result
	res = mins * (1 / det);

	return res;
}

#pragma endregion

#pragma region "###	QUATERNION	###"

Quaternion::Quaternion() 
{
	vector.x = vector.y = vector.z = 0.0f;
	w = 1.0f;
}

void Quaternion::Identity()
{
	vector.x = vector.y = vector.z = 0.0f;
	w = 1.0f;
}

void Quaternion::Set(float iW, float iX, float iY, float iZ)
{
	vector.x = iX;
	vector.y = iY;
	vector.z = iZ;
	w = iW;
}

void Quaternion::Set(float3 &axis, float angle)
{
	// angle in radian's
	float sinA = sin(angle / 2);

	// quat = ( sin(a/2) * axis.x, sin(a/2) * axis.y, sin(a/2) * axis.z, cos(a/2) )
	w = cos(angle / 2);
	vector.x = sinA * axis.x,
	vector.y = sinA * axis.y,
	vector.z = sinA * axis.z;
}

void Quaternion::Set(float4 &from, float4 &to)
{
	// not tested
	from.Normalize();
	to.Normalize();

	float4 cross = from.Cross(to);
	float dot = from.Dot(to);

	Set(cross.vector.x, cross.vector.y, cross.vector.z, dot);

	w += 1.0f;

	if (w <= F_EPSILON)
	{
		if ((from.vector.z * from.vector.z) > (from.vector.x * from.vector.x))
			Set(0.0f, from.vector.z, -from.vector.y, w);
		else
			Set(from.vector.y, -from.vector.x, 0.0f, w);
	}
}

void Quaternion::Normalize()
{
	float m = Magnitude();

	if (m > F_EPSILON)
	{
		float s = 1.0f / Magnitude();
		*this *= s;
	}
	else
		Identity();
}

void Quaternion::operator+=(const Quaternion &iQuat)
{
	vector.x += iQuat.vector.x;
	vector.y += iQuat.vector.y;
	vector.z += iQuat.vector.z;
	w		 += iQuat.w;
}

Quaternion Quaternion::operator+(const Quaternion &iQuat)
{
	Quaternion temp;

	temp.vector.x = vector.x + iQuat.vector.x;
	temp.vector.y = vector.y + iQuat.vector.y;
	temp.vector.z = vector.z + iQuat.vector.z;
	temp.w		  = w + iQuat.w;

	return temp;
}

void Quaternion::operator-=(const Quaternion &iQuat)
{
	vector.x -= iQuat.vector.x;
	vector.y -= iQuat.vector.y;
	vector.z -= iQuat.vector.z;
	w		 -= iQuat.w;
}

Quaternion Quaternion::operator-(const Quaternion &iQuat)
{
	Quaternion temp;

	temp.vector.x = vector.x - iQuat.vector.x;
	temp.vector.y = vector.y - iQuat.vector.y;
	temp.vector.z = vector.z - iQuat.vector.z;
	temp.w		  = w - iQuat.w;

	return temp;
}

void Quaternion::operator*=(float k)
{
	vector.x *= k;
	vector.y *= k;
	vector.z *= k;
	w		 *= k;
}

Quaternion Quaternion::operator*(float k)
{
	Quaternion temp;

	temp.vector.x = vector.x * k;
	temp.vector.y = vector.y * k;
	temp.vector.z = vector.z * k;
	temp.w		  = w * k;

	return temp;
}

void Quaternion::operator*=(const Quaternion &iQuat)
{
	float4 qVec, qIVec;

	// 'convert' Quaternion to float4 for compting Dot and Cross product with float4 methods
	qVec.vector		= vector;
	qIVec.vector	= iQuat.vector;

	float4 cProd	= qVec.Cross(qIVec);
	float dProd		= qVec.Dot(qIVec);

	float4 resVec	= cProd + (qIVec * w) + (qVec * iQuat.w);	// VV' + wV' + w'V
	float resW		= w * iQuat.w - dProd;						// ww' - V*V

	// restore to Quaternion
	vector	= resVec.vector;
	w		= resW;
}

Quaternion Quaternion::operator*(const Quaternion &iQuat)
{
	Quaternion temp;

	float4 qVec, qIVec;

	// 'convert' to float4 for compting Dot and Cross product with float4 methods
	qVec.vector		= vector;
	qIVec.vector	= iQuat.vector;

	float4 cProd	= qVec.Cross(qIVec);
	float dProd		= qVec.Dot(qIVec);

	float4 resVec	= cProd + (qIVec * w) + (qVec * iQuat.w);	// VV' + wV' + w'V
	float resW		= w * iQuat.w - dProd;						// ww' - V*V

	// restore to Quaternion
	temp.vector = resVec.vector;
	temp.w		= resW;

	return temp;
}

float Quaternion::Norm()
{
	return (
		vector.x * vector.x +
		vector.y * vector.y +
		vector.z * vector.z +
		w * w 
		);
}

float Quaternion::Magnitude()
{
	return sqrt(Norm());
}

Quaternion Quaternion::Conjugate()
{
	Quaternion temp;

	temp.vector.x = -vector.x;
	temp.vector.y = -vector.y;
	temp.vector.z = -vector.z;
	temp.w		  = w;

	return temp;
}

Quaternion Quaternion::Inverse()
{
	Quaternion temp;

	float n = 1.0f / Norm();

	temp.vector.x = -vector.x * n;
	temp.vector.y = -vector.y * n;
	temp.vector.z = -vector.z * n;
	temp.w		  = w * n;

	return temp;
}

float Quaternion::Inner(const Quaternion &quat)
{
	return ( 
		vector.x * quat.vector.x + 
		vector.y * quat.vector.y +
		vector.z * quat.vector.z + 
		w * quat.w );
}

float4 Quaternion::RotFloat(float4 &vector)
{
	float4 res;

	Quaternion temp;	// interim result in Quaternion
	Quaternion tempV;	// converted float4
	Quaternion inverse; // inverted Quaternion

	// 'convert' float4 to Quaternion as (x, y, z, 0) 
	tempV.vector = vector.vector;
	tempV.w = 0.0f;

	inverse = this->Inverse();

	temp = *this*tempV*inverse; // V' = qV~q

	// 'convert' to float4
	res.vector = temp.vector;
	res.w = 0.0f;

	return res;
}

void Quaternion::ToMatrix(float4x4 &outMat)
{
	float wx,	 wy,	 wz;
	float xx,	 yy,	 zz;
	float xy,	 xz,	 yz;
	float x2,	 y2,	 z2;

	float s = 2.0f / Norm();

	x2 = vector.x *s;	 y2 = vector.y * s;	 z2 = vector.z * s;
	xx = vector.x * x2;	 yy = vector.y * y2; zz = vector.z * z2;
	xy = vector.x * y2;	 xz = vector.x * z2; yz = vector.y * z2;
	wx = w * x2;		 wy = w * y2;		 wz = w * z2;

	outMat.Set(
		1.0f - (yy + zz),	 xy + wz,			 xz - wy,			 0,
		xy - wz,			 1.0f - (xx + zz),	 yz + wx,			 0,
		xz + wy,			 yz - wx,			 1.0f - (xx + yy),	 0,
		0, 0, 0, 1.0f
		);
}

void Quaternion::ToAxisAngle(float4 &axis, float &angle)
{
	float vecL = sqrt(
		vector.x * vector.x +
		vector.y * vector.y +
		vector.z * vector.z
		);

	if (vecL > F_EPSILON)
	{
		float iVecL = 1.0f / vecL;

		axis.Set(vector.x * iVecL, vector.y * iVecL, vector.z * iVecL);

		if (w < 0)
			angle = 2.0f * atan2(-vecL, -w);
		else
			angle = 2.0f * atan2(vecL, w);
	}
	else
	{
		axis.Set(0, 0, 0);
		angle = 0.0f;
	}
}

void qLerp::Setup(Quaternion &iFrom, Quaternion &iTo)
{
	from	= iFrom;
	to		= iTo;

	from.Normalize();
	to.Normalize();

	float inner = from.Inner(to);

	if (inner < 0)
		to *= -1.0f;

	to -= from;
}

void qSlerp::Setup(Quaternion &iFrom, Quaternion &iTo)
{
	from = iFrom;
	to	 = iTo;

	from.Normalize();
	to.Normalize();

	float cos_omega = from.Inner(to);

	if (cos_omega < 0)
	{
		cos_omega = -cos_omega;
		to *= -1.0f;
	}
	if (cos_omega > 0.9999f)
		cos_omega = 0.9999f;

	omega = acosf(cos_omega);

	float invSinO = 1.0f / sin(omega);

	from *= invSinO;
	to *= invSinO;
}

#pragma endregion

#pragma once
#include <math.h>

#define F_EPSILON 0.0000005

#define F_PI	3.14159265359
#define F_2PI	6.28318530718

#define F_POINT	 1.0f
#define F_VECTOR 0.0f

// FLOAT4x4
// order: RowMajor

struct float2x2
{
	float m00, m01;
	float m10, m11;

	float2x2();
	~float2x2();

	void Set(
		float iM00, float iM01, 
		float iM10, float iM11
		);

	float Det();
};

struct float3x3
{
	float m00, m01, m02;
	float m10, m11, m12;
	float m20, m21, m22;

	float3x3();
	~float3x3();

	void Set(
		float iM00, float iM01, float iM02, 
		float iM10, float iM11, float iM12, 
		float iM20, float iM21, float iM22
		);

	float Det();
};

class float4x4
{
public:
	union {
		struct {
			float m00, m01, m02, m03; // row 0
			float m10, m11, m12, m13; // 
			float m20, m21, m22, m23; // 
			float m30, m31, m32, m33; // 
		};
		float mat[4][4];
	};

public:
	float4x4();
	float4x4(
		float iM00, float iM01, float iM02, float iM03,
		float iM10, float iM11, float iM12, float iM13,
		float iM20, float iM21, float iM22, float iM23,
		float iM30, float iM31, float iM32, float iM33
		);

	~float4x4() {};

	void Identity();
	void Zero();
	void Set(
		float iM00, float iM01, float iM02, float iM03,
		float iM10, float iM11, float iM12, float iM13,
		float iM20, float iM21, float iM22, float iM23,
		float iM30, float iM31, float iM32, float iM33
		);

	float Det();			//
	float4x4 Inverse();		//
	float4x4 Transpose();	//

	bool operator==(const float4x4 &iMat);		//
	bool operator!=(const float4x4 &iMat);		//

	float4x4 operator+(const float4x4 &iMat);	//
	void operator+=(const float4x4 &iMat);		//
	float4x4 operator-(const float4x4 &iMat);	//
	void operator-=(const float4x4 &iMat);		//

	float4x4 operator*(float k);				//
	void operator*=(float k);					//

	float4x4 operator*(const float4x4 &iMat);	//
	void operator*=(const float4x4 &iMat);		//
};

// FLOAT4

struct float3
{
	float x;
	float y;
	float z;
};

class float4
{
public:

	float3 vector;
	float w;

public:

	float4();
	float4(float iX, float iY, float iZ,
		float iW = F_VECTOR);
	float4(float3 &vec,
		float iW = F_VECTOR);

	~float4() {};

	void Set(float iX, float iY, float iZ,
		float iW = F_VECTOR); 

	// all operations executes only for float3 element
	// * save equality operations and with float4x4

	float Length();					//
	void Normalize();				//

	float Dot(const float4 &vec);			//
	float4 Cross(const float4 &vec);		//

	float4 operator+ (const float4& vec);	//
	void operator+= (const float4 &vec);	//

	float4 operator- (const float4& vec);	//
	void operator-= (const float4 &vec);	//

	float4 operator* (float k);				//
	void operator*= (float k);				//

	float4 operator* (const float4x4 &mat); //
	void operator*= (const float4x4 &mat);	//

	bool operator== (const float4& vec);	//
	bool operator!= (const float4& vec);	//
};

// QUATERNION

class Quaternion
{
public:

	float w;
	float3 vector;

public:

	Quaternion();
	~Quaternion() {};

	void Set(float iW, float iX, float iY, float iZ);
	void Set(float3 &axis, float angle);	// set quaternion via axisAngle 
	void Set(float4 &from, float4 &to);	// set from 2 directions : NOT tested
	void Identity();
	void Normalize();

	float Norm();
	float Magnitude();
	float Inner(const Quaternion &quat);

	Quaternion Conjugate();
	Quaternion Inverse();

	float4 RotFloat(float4 &vector);

	void ToAxisAngle(float4 &axis, float &angle);
	void ToMatrix(float4x4 &outMat);

	Quaternion operator+(const Quaternion &iQuat);
	void operator+=(const Quaternion &iQuat);

	Quaternion operator-(const Quaternion &iQuat);
	void operator-=(const Quaternion &iQuat);

	Quaternion operator*(const Quaternion &iQuat);
	void operator*=(const Quaternion &iQuat);

	Quaternion operator*(float k);
	void operator*=(float k);
};

// LERP & SLERP of QUATERNION

class qLerp
{
public:
	Quaternion from;
	Quaternion to;

public:
	void Setup(Quaternion &iFrom, Quaternion &iTo);
	
	void Interpolate(float t, Quaternion &out) {
		out = from + to * t;
	}
};

class qSlerp
{
public:
	Quaternion from;
	Quaternion to;
	float omega;

public:
	void Setup(Quaternion &iFrom, Quaternion &iTo);

	void Interpolate(float t, Quaternion &out) {
		out = from * sin((1.0f - t) * omega) + to*sin(t * omega);
	}
};

// GLOBAL FUNCTIONS

float Dot(float4 &first, float4 &second);
float4 Cross(float4 &first, float4 &second);

// TRANSFORMATION FUNCTIONS

void TranslateToMat(float3 &pos, float4x4 &mat);
void ScaleToMat(float3 &scale, float4x4 &mat);
void RotToMat(float3 axis, float angle, float4x4 &mat);
/* Copyright 2013 by trilliondimention@163.com */


#include "3Dmath.h"

void math_2DVector::normalize()
{
	float denom = length();
	u/=denom;
	v/=denom;
}

void math_3DVector::normalize()
{
	float denom = length();
	x /= denom;
	y /= denom;
	z /= denom;
}

void math_4DVector::normalize()
{
	float t = length();
	x /= t;
	y /= t;
	z /= t;
	w /= t;
}

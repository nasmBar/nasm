/* Copyright 2013 by trilliondimention@163.com */

#ifndef I3DMATH_H
#define I3DMATH_H

#include <math.h>
#define PI 3.14159265358979323846f

/* 2D vector*/
struct math_2DVector
{
	float u,v;

	math_2DVector()
	{
		u = 0.0;
		v = 0.0;
	}
	math_2DVector(const float u, const float v)
	{
		this->u = u;
		this->v = v;
	}	
	void set(const float u, const float v)
	{
		this->u = u;
		this->v = v;
	}
	
	math_2DVector operator+(const math_2DVector& other)
	{
		math_2DVector ve(u+other.u, v + other.v);
		return ve;
	}
	math_2DVector operator-(const math_2DVector& other)
	{
		math_2DVector ve(u - other.u, v - other.v);
		return ve;
	}
	math_2DVector operator-()
	{
		math_2DVector ve(-u, -v);
		return ve;
	}
	float operator*(const math_2DVector& other)
	{
		return (u * other.u + v * other.v);
	}
	math_2DVector operator*(const float s)
	{
		math_2DVector ve(u * s, v * s);
		return ve;	
	}
	math_2DVector operator/(const float s)
	{
		math_2DVector ve(u/s, v/s);
		return ve;
	}
	float length()
	{
		return (sqrt(u*u + v*v));
	}
	void normalize();
};

/* 3D vector*/
struct math_3DVector
{
	float x,y,z;
	
	math_3DVector()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}	
	math_3DVector(const float x, const float y, const float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	void set(const float x, const float y, const float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	math_3DVector operator+(const math_3DVector& other)
	{
		math_3DVector ve(x + other.x, y + other.y, z + other.z);
		return ve;
	}
	math_3DVector operator-(const math_3DVector& other)
	{
		math_3DVector ve(x - other.x, y - other.y, z - other.z);
		return ve;
	}
	math_3DVector operator-()
	{
		math_3DVector ve(-x, -y, -z);
		return ve;
	}

	math_3DVector operator*(float s)
	{
		math_3DVector ve(x * s, y * s, z * s);
		return ve;
	}

	float operator*(const math_3DVector& other)
	{
		return (x * other.x + y * other.y + z + other.z);
	}

	math_3DVector cross(const math_3DVector& other)
	{
		math_3DVector ve(y * other.z - z * other.y, z * other.x - x*other.z, x * other.y - y * other.x);
	}

	float length()
	{
		return (sqrt(x * x + y * y + z * z));
	}
	float sqrtLength()
	{
		return (x * x + y * y + z * z);
	}	

	void normalize();

	float max()
	{
		float max = x;
		if(max < y )
			max = y;
		if(max < z)
			max = z;
		return max;
	}
	
	float min()
	{
		float min = x;
		if(min >y)
			min = y;
		if(min > z)
			min = z;
		return min;
	}
};

/* 4D vector*/
struct math_4DVector
{
	float x,y,z,w;

	math_4DVector()
	{
		x = 0; y = 0; z = 0; w = 0;
	}	

	math_4DVector(const float x, const float y, const float z, const float w)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	void set(const float x, const float y, const float z, const float w)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	math_4DVector operator+(const math_4DVector& other)
	{
		math_4DVector ve(x + other.x, y + other.y, z + other.z, w + other.w);
		return ve;
	}
	math_4DVector operator-(const math_4DVector& other)
	{
		math_4DVector ve(x - other.x, y - other.y, z - other.z, w - other.w);
		return ve;
	}
	math_4DVector operator-()
	{
		math_4DVector ve(-x, -y, -z, -w);
		return ve;
	}
	
	float operator*(const math_4DVector& other)
	{
		return (x * other.x + y * other.y + z * other.z + w*other.w);
	}

	math_4DVector operator*(float s)
	{
		math_4DVector ve(x * s, y * s, z * s, w * s);
		return ve;
	}

	math_4DVector operator/(float s)
	{
		math_4DVector ve(x/s, y/s, z/s, w/s);
		return ve;
	}	
	
	float length()
	{
		return (sqrt(x * x + y * y + z * z + w * w));
	}

	float sqrtLength()
	{
		return (x * x + y * y + z * z + w * w);
	}

	void normalize();
};

/* 3D matrix */
struct math_3DMatrix
{
	float m[3][3];
	
	math_3DMatrix()
	{
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++)
				m[j][i] = 0.0;
	}
	
	math_3DMatrix(float diagonal)
	{
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++);
			{
				if(i == j)
					m[i][j] = diagonal;
				else
					m[i][j] = 0.0;
			}
	}

	void setDiagonal(const math_3DVector v)
	{
		m[0][0] = v.x;
		m[1][1] = v.y;
		m[2][2] = v.z;
	}

	void makeIdentity()
	{
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++);
			{
				if(i == j)
					m[i][j] = 1.0;
				else
					m[i][j] = 0.0;
			}
	}

	math_3DMatrix operator+(const math_3DMatrix& other)
	{
		math_3DMatrix mat;
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++)
			{
				mat.m[i][j] = m[i][j] + other.m[i][j];
			}
		return mat;
	}

	math_3DMatrix operator-(const math_3DMatrix& other)
	{
		math_3DMatrix mat;
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++)
			{
				mat.m[i][j] = m[i][j] - other.m[i][j];
			}
		return mat;
	}

	math_3DMatrix operator*(const float s)
	{
		math_3DMatrix mat;
		int i,j;
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++)
			{
				mat.m[i][j] = m[i][j] * s;
			}
		return mat;
	}

	void setRow(const int i, const math_3DVector& v)
	{
		m[i][0] = v.x;
		m[i][1] = v.y;
		m[i][2] = v.z;
	}	
	math_3DVector getRow(const int i)
	{
		math_3DVector v(m[i][0], m[i][1], m[i][2]);
		return v;
	}

	void setColumn(const int i, const math_3DVector& v)
	{
		m[0][i] = v.x;
		m[1][i] = v.y;
		m[2][i] = v.z;
	}
	math_3DVector getColumn(const int i)
	{
		math_3DVector v(m[0][i], m[1][i], m[2][i]);
		return v;
	}	

	math_3DVector operator*(const math_3DVector& v)
	{
		math_3DVector res;
		res.x = getRow(0) * v;
		res.y = getRow(1) * v;
		res.z = getRow(2) * v;
	
		return res;
	}

	math_2DVector operator*(const math_2DVector& v)
	{
		math_3DVector t(v.u, v.v, 1.0);
		math_3DVector res;
		res.x = getRow(0) * t;
		res.y = getRow(1) * t;
		res.z = getRow(2) * t;
		
		math_2DVector last(res.x, res.y);
		return last;
	}

	math_3DMatrix operator*(const math_3DMatrix& mat)
	{
		math_3DMatrix res;
		
		int i,j,k;
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				for(k = 0; k < 3; k++)
				{
					res.m[i][j] += m[i][k] * mat.m[k][j]; 
				}
			}
		}	
		return res;
	}

	void setMatrix(const math_3DMatrix other)
	{
		int i,j;
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				m[i][j] = other.m[i][j];
			}
		}
	}

	float det()
	{
		float res;
		res = m[0][0]*m[1][1]*m[2][2]
  		+m[1][0]*m[2][1]*m[0][2]
  		+m[2][0]*m[0][1]*m[1][2]

  		-m[2][0]*m[1][1]*m[0][2]
  		-m[1][0]*m[0][1]*m[2][2]
  		-m[0][0]*m[2][1]*m[1][2];
		
		return res;
	}

	math_3DMatrix transpose()
	{
		math_3DMatrix res;
		int i,j;
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3;j++)
			{
				res.m[j][i] = m[i][j];
			}
		}
		return res;
	}

	math_3DMatrix inverse()
	{
		float d = det();
		math_3DMatrix res;
		if(d * d < 0.00000000001)
		{
			return res;
		}
		float invdet = (float)1.0/d;
		res.m[0][0] =  (m[1][1]*m[2][2]-m[2][1]*m[1][2])*invdet;
		res.m[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1])*invdet;
		res.m[0][2] =  (m[0][1]*m[1][2]-m[0][2]*m[1][1])*invdet;
		res.m[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0])*invdet;
		res.m[1][1] =  (m[0][0]*m[2][2]-m[0][2]*m[2][0])*invdet;
		res.m[1][2] = -(m[0][0]*m[1][2]-m[1][0]*m[0][2])*invdet;
		res.m[2][0] =  (m[1][0]*m[2][1]-m[2][0]*m[1][1])*invdet;
		res.m[2][1] = -(m[0][0]*m[2][1]-m[2][0]*m[0][1])*invdet;
		res.m[2][2] =  (m[0][0]*m[1][1]-m[1][0]*m[0][1])*invdet;	
		return res;
	}
		
	void makeRotateX(const float angle)
	{
		float s = sin(angle);
		float c = cos(angle);
	
		m[1][1] = c;
		m[2][1] = -s;
		m[1][2] = s;
		m[2][2] = c;
	}

	void makeRotateY(const float angle)
	{
		float s = sin(angle);
		float c = cos(angle);

		m[0][0] = c;
		m[0][2] = s;
		m[2][0] = -s;
		m[2][2] = c;
	}

	void makeRotateZ(const float angle)
	{
		float s = sin(angle);
		float c = cos(angle);
	
		m[0][0] = c;
		m[0][1] = -s;
		m[1][0] = s;
		m[1][1] = c;
	}
	/* 3D Math Primer for Graphics and Game Development */
	void makeRotateVector(const float angle, const math_3DVector axis)
	{
		float c = cos(angle);
		float s = sin(angle);
		float oneMinusCos = 1 - c;
		m[0][0] = c + axis.x * axis.x * oneMinusCos;
		m[1][0] = axis.z * s + axis.y * axis.x * oneMinusCos;
		m[2][0] = axis.y * s + axis.z * axis.x * oneMinusCos;
		m[0][1] = axis.z * s + axis.x * axis.y * oneMinusCos;
	
		m[1][1] = c + axis.y * axis.y * oneMinusCos;
		m[2][1] = axis.x * s + axis.z * axis.y * oneMinusCos; 
		m[0][2] = axis.y * s + axis.x * axis.z * oneMinusCos;
		m[1][2] = axis.x * s + axis.y * axis.z * oneMinusCos;
		m[2][2] = c + axis.z * axis.z * oneMinusCos;
	}

	void setPitchYawRoll(math_3DVector pyr)
	{
		float cx,sx,cy,sy,cz,sz;
  		cx=cos(pyr.x); sx=sin(pyr.x);
		cy=cos(pyr.y); sy=sin(pyr.y);
   		cz=cos(pyr.z); sz=sin(pyr.z);
   		m[0][0]=cz*cy+sz*sy*sx;
   		m[1][0]=cx*sz;
   		m[2][0]=-sy*cz+cy*sx*sz;

   		m[0][1]=-sz*cy+cz*sy*sx;
   		m[1][1]=cx*cz;
   		m[2][1]=sy*sz+cz*cy*sx;

   		m[0][2]=sy*cx;
   		m[1][2]=-sx;
   		m[2][2]=cx*cy;
	}

	math_3DVector getPitchYawRoll()
	{	
	 math_3DVector newv;
	 newv.x=atan2(-m[1][2],sqrt(m[0][2]*m[0][2]+m[2][2]*m[2][2]));
	 if ((newv.x-PI/2)*(newv.x-PI/2)<0.00000000001)
  	{ newv.z=0;
    		newv.y=atan2(m[0][1],m[2][1]);
   		return newv;
  	}
  	if ((newv.x+PI/2)*(newv.x+PI/2)<0.00000000001)
  	{
   		newv.z=0;
   		newv.y=-atan2(m[0][1],m[2][1]);
   		if (newv.y>3.14159 || newv.y<-3.14159) newv.y=0; //Just a small hack
     //cout << "SPECIAL A  " << m[0+1*3]<< "  " << m[2+1*3] << endl;
   		return newv;
	}
	float ic=((float)1.0)/cos(newv.x);
	newv.y=atan2(m[0][2]*ic,m[2][2]*ic);
	newv.z=atan2(m[1][0]*ic,m[1][1]*ic);
	return newv;
	}
};

struct math_4DMatrix
{
	float m[4][4];

	math_4DMatrix()
	{
		int i,j;
		for(i = 0; i < 4; i++)
		{
			for(j = 0; j < 4; j++)
			{
				m[j][i] = 0.0;
			}
		}
	}

	math_3DMatrix(float diagonal)	
	{
		int i,j;
		for(i = 0; i < 4;i++)
		{
			for(j = 0; j < 4; j++)
			{
				if(i == j)
					m[i][j] = diagonal;
				else
					m[i][j] = 0;
			}
		}
	}

	
};


#endif

#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>
#include <random>
#include <iostream>
using namespace std;
const double PI = acos(-1.0);
struct WxzRandom
{
	WxzRandom()
	{
		generator=default_random_engine(rd());
		dis=uniform_real_distribution<double>(0.0,1.0);
	}
	double C11Rand()
	{
		return dis(generator);
	}
	default_random_engine generator;
	random_device rd;
	uniform_real_distribution<double> dis; 
};
WxzRandom myRandom;
//double C11Rand()
//{
//	default_random_engine generator(time(nullptr));
//	uniform_real_distribution<double> dis(0.0,1.0); 
//	return dis(generator);
//	srand(time(NULL));
//	return (double)rand()/RAND_MAX;
//}
struct Vec3
{
	double x, y, z;
	Vec3()
	{
		x = 1;
		y = 1;
		z = 1;
	}
	Vec3(double X, double Y, double Z)
	{
		x = X;
		y = Y;
		z = Z;
	}
	double length()
	{
		return sqrt(x*x + y*y + z*z);
	}
	double sqrlength()
	{
		return x*x + y*y + z*z;
	}
	Vec3 normalize()
	{
		double inv = (1 / (this->length()));
		return Vec3(inv*x, inv*y, inv*z);
	}
	Vec3 negate()
	{
		return Vec3(-x, -y, -z);
	}
	Vec3 operator + (const Vec3 b)
	{
		return Vec3(x + b.x, y + b.y, z + b.z);
	}
	Vec3 operator - (const Vec3 b)
	{
		return Vec3(x - b.x, y - b.y, z - b.z);
	}
	friend Vec3 operator * (double a, Vec3 my)
	{
		return Vec3(a*my.x, a*my.y, a*my.z);
	}
	friend Vec3 operator * (Vec3 my, double a)
	{
		return Vec3(a*my.x, a*my.y, a*my.z);
	}
	Vec3 operator /(const double a)
	{
		return Vec3(x / a, y / a, z / a);
	}
	Vec3 operator ^ (Vec3 v)
	{
		return Vec3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
	}
	double dot(Vec3 v)
	{
		return x*v.x + y*v.y + z*v.z;
	}
};
struct Ray
{
	Vec3 origin, direction;
	Ray()
	{
		origin = Vec3(0, 0, 0);
		direction = Vec3(0, 0, 1);
	}
	Ray(Vec3 o, Vec3 dir)
	{
		origin = o;
		direction = dir;
	}
	Vec3 getpoint(double t)
	{
		return origin + direction*t;
	}
};

struct Color
{
	double r, g, b;
	Color() { r = 0; g = 0; b = 0; };
	Color(double R, double G, double B)
	{
		r = R;
		g = G;
		b = B;
	}
	Color operator +(Color v)
	{
		return Color(r + v.r, g + v.g, b + v.b);
	}
	friend Color operator *(Color old, double t)
	{
		return Color(old.r*t, old.g*t, old.b*t);
	}
	friend Color operator *(double t, Color old)
	{
		return Color(old.r*t, old.g*t, old.b*t);
	}
	Color operator *(Color v)	//调制，即三个颜色通道相乘。
	{
		return Color(r*v.r, g*v.g, b*v.b);
	}
	Color operator *(Vec3 v)
	{
		return Color(r*v.x, g*v.y, b*v.z);
	}

};

struct IntersectResult;

class Material
{
public:
	virtual bool Scatter(Ray r_in,IntersectResult result,Vec3 &attenuation,Ray &Scattered)=0;
	virtual Color emitted()  
	{
        return Color(0,0,0); 
	}
};


class Geometry
{
public:
	Material *material;
	virtual IntersectResult intersect(Ray ray) = 0;
};

struct IntersectResult
{
	Geometry* geometry;
	double dist;
	Vec3 pos, normal;
	IntersectResult()
	{
		geometry = NULL;
		dist = 0;
		pos = normal = Vec3(0, 0, 0);
	}
};



struct Sphere :public Geometry
{
	Vec3 center;
	double sqrRadius, radius;
	Sphere() { center = Vec3(0, 0, 0); radius = sqrRadius = 1; }
	Sphere(Vec3 C, double R)
	{
		center = C;
		radius = R;
		sqrRadius = R*R;
	}
	IntersectResult intersect(Ray ray)
	{
		IntersectResult result;
		Vec3 v = ray.origin - center;//光源到球心的方向向量
		double a0 = v.sqrlength() - sqrRadius;//这TM是啥……
		double DdotV = ray.direction.dot(v);//光线方向向量点乘光源到球心的方向向量
		if (DdotV <= 0)
		{
			double discr = DdotV*DdotV - a0;
			if (discr >= 0)
			{
				result.geometry = this;
				result.dist = -DdotV - sqrt(discr);
				result.pos = ray.getpoint(result.dist);
				result.normal = (result.pos - center).normalize();
				return result;
			}
		}
		return result;
	}

};
struct Plane :Geometry
{
	Vec3 normal, pos;
	double d;
	Plane() { normal = Vec3(0, 0, 1); pos = Vec3(0, 0, 0); d = 1; }
	Plane(Vec3 mynormal, double myd)
	{
		normal = mynormal;
		d = myd;
		pos = normal*d;
	}
	IntersectResult intersect(Ray ray)
	{
		IntersectResult result;
		double a = ray.direction.dot(normal);
		if (a >= 0) return result;
		double b = normal.dot(ray.origin - pos);
		result.geometry = this;
		result.dist = -b / a;
		result.pos = ray.getpoint(result.dist);
		result.normal = normal;
		return result;
	}
};
struct Heart :Geometry
{
	Vec3 normal,pos;
	double width;
	Heart() { normal = Vec3(0, 0, 1); pos = Vec3(0, 0, 0); width=0.1; }
	Heart(Vec3 mynormal, Vec3 mypos,double w)
	{
		normal = mynormal;
		pos = mypos;
		width=w;
	}
	IntersectResult intersect(Ray ray)
	{
		IntersectResult result,tmp;
		double a = ray.direction.dot(normal);
		if (a >= 0) return result;
		double b = normal.dot(ray.origin - pos);
		tmp.dist = -b / a;
		tmp.pos = ray.getpoint(tmp.dist);
		double xx=tmp.pos.x-pos.x;
		double yy=tmp.pos.y-pos.y;
//		cout<<"Debug "<<tmp.dist
		double delta=pow((pow(xx,2)+pow(yy,2)-1),3)-pow(xx,2)*pow(yy,3);
		if (delta>width || delta<0) return result; 
		result.geometry = this;
		result.dist=tmp.dist;
		result.pos=tmp.pos; 
		result.normal = normal;
		return result;
	}
};

Vec3 Random_in_Unit_Disk()
{
	Vec3 p;
	do
	{
		p=2.0*Vec3(myRandom.C11Rand(),myRandom.C11Rand(),0)-Vec3(1,1,0);
	}while(p.dot(p)>1.0);
	return p;
}

struct P_Camera //透视摄像机
{
	Vec3 eye, front, refUp, right, up;
	double fov, fovScale,lensRadius;
	P_Camera() {}
	P_Camera(Vec3 Eye, Vec3 Front, Vec3 RefUp, double Fov,double lens)
	{
		eye = Eye;
		front = Front;
		refUp = RefUp;
		fov = Fov;
		right = front^refUp;
		up = right^front;
		fovScale = tan(fov*0.5*PI / 180) * 2;
		lensRadius=lens;
	}
	Ray generateRay(double x, double y)
	{
		Vec3 r = right*(x - 0.5)*fovScale;
		Vec3 u = up*(y - 0.5)*fovScale;
		Vec3 rd=lensRadius*Random_in_Unit_Disk();
		Vec3 offset=up*rd.x+right*rd.y;
		return Ray(eye+offset, (front + r + u).normalize());
	}
};

//struct Camera
//{
//	Vec3 origin;
//	Vec3 lower_left_corner;
//	Vec3 horizontal;
//	Vec3 vertical;
//	Camera()
//	{
//		lower_left_corner=Vec3(-2.0,-1.0,-1.0);
//		horizontal=Vec3(4.0,0.0,0.0);
//		vertical=Vec3(0.0,2.0,0.0);
//		origin=Vec3(0.0,0.0,0.0);
//	}
//	Ray get_ray(double u,double v)
//	{
//		return Ray(origin,lower_left_corner+u*horizontal+v*vertical-origin);
//	}
//};
struct Scene
{
	vector<Geometry*> geometries;
	Scene()
	{
		geometries.clear();
	}
	void Add(Geometry* mygeometry)
	{
		geometries.push_back(mygeometry);
	}
	IntersectResult intersect(Ray ray)
	{
		double minDistance = 1e9;
		IntersectResult result;
		for (int i = 0; i < geometries.size(); i++)
		{
			IntersectResult tmpresult = geometries[i]->intersect(ray);
			if (tmpresult.geometry != NULL && tmpresult.dist < minDistance)
			{
				minDistance = tmpresult.dist;
				result = tmpresult;
			}
		}
		return result;
	}
};


//struct CheckerMaterial :public Material
//{
//	double scale;
//	Color color;
//	CheckerMaterial(double S, double R,Color mycolor)
//	{
//		scale = S;
//		reflectiveness = R;
//		color = mycolor;
//	}
//	Color sample(Ray ray, Vec3 pos, Vec3 normal,Vec3 lightDir, Color lightColor)
//	{
//		double t = abs(floor(pos.y*0.1)+floor(pos.z*scale));
//		int int_t = (int)(t + 0.01);
//		return  (int_t % 2) < 1 ? color : Color(0.9, 0.9, 0.9);
//	}
//};
//
////Vec3 lightDir = Vec3(5, -4, 1).normalize();
////Color lightColor = Color(1, 1, 1);
//
//struct PhongMaterial :public Material
//{
//	Color diffuse, specular;  //漫反射和镜射颜色
//	double shininess;
//	PhongMaterial(Color diff, Color spec, double shin, double reflect)
//	{
//		diffuse = diff;
//		specular = spec;
//		shininess = shin;
//		reflectiveness = reflect;
//	}
//	Color sample(Ray ray, Vec3 pos, Vec3 normal,Vec3 lightDir, Color lightColor)
//	{
//		double NdotL = normal.dot(lightDir);
//		Vec3 H = (lightDir - ray.direction).normalize();
//		double NdotH = normal.dot(H);
//		Color diffuseTerm = diffuse*max(NdotL, (double)0);
//		Color specularTerm = specular*pow(max(NdotH, (double)0), shininess);
//		return lightColor*(diffuseTerm + specularTerm);
//	}
//};
Vec3 Random_in_Unit_Sphere()
{
	Vec3 p;
	do
	{
		p=2.0*Vec3(myRandom.C11Rand(),myRandom.C11Rand(),myRandom.C11Rand())-Vec3(1,1,1);
//		cout<<"Debug "<<p.dot(p)<<endl; 
	}
	while (p.dot(p)>=1.0);
	return p;
}

struct Lambertian:public Material
{
	Lambertian(const Vec3 a):albedo(a){}
	virtual bool Scatter(Ray r_in,IntersectResult result,Vec3& attenuation,Ray& Scattered)
	{
		Vec3 Target=result.pos+result.normal;
		Target=Target+Random_in_Unit_Sphere();
		Scattered=Ray(result.pos,Target-result.pos);
		attenuation=albedo;
		return true;
	}
	Vec3 albedo;
};

Vec3 reflect(Vec3 v,Vec3 n)
{
	return v-2*v.dot(n)*n;
}

struct Metal:public Material
{
	public:
		Metal(const Vec3 a,float f):albedo(a){if (f<1) fuzz=f;else fuzz=1;}
		virtual bool Scatter(Ray r_in,IntersectResult result,Vec3 &attenuation,Ray& Scattered)
		{
			Vec3 reflected=reflect(r_in.direction.normalize(),result.normal);
			Scattered=Ray(result.pos,reflected+fuzz*Random_in_Unit_Sphere());
			attenuation=albedo;
			return (Scattered.direction.dot(result.normal)>0);
		}
		Vec3 albedo;
		double fuzz;
};

bool refract(Vec3 v,Vec3 n,float ni_over_nt,Vec3 &Refracted)
{
	Vec3 uv=v.normalize();
	float dt=uv.dot(n);
	float discriminant=1.0-ni_over_nt*ni_over_nt*(1-dt*dt);
	if (discriminant>0)
	{
		Refracted=ni_over_nt*(v-n*dt)-n*sqrt(discriminant);
		return true;
	}
	else return false;
}

float schlick(float cosine,float ref_idx)
{
	float r0=(1-ref_idx)/(1+ref_idx);
	r0*=r0;
	return r0+(1-r0)*pow((1-cosine),5);
}
struct dielectric:public Material
{
	public:
		dielectric(float ri):ref_idx(ri){}
		virtual bool Scatter(Ray r_in,IntersectResult result,Vec3& attenuation,Ray &Scattered)
		{
			Vec3 outward_normal;
			Vec3 reflected=reflect(r_in.direction,result.normal);
			float ni_over_nt;
			attenuation=Vec3(1,1,1);
			Vec3 refracted;
			float reflect_prob;
			float cosine;
			if (r_in.direction.dot(result.normal)>0)
			{
				outward_normal=-1*result.normal;
				ni_over_nt=ref_idx;
				cosine=ref_idx*r_in.direction.dot(result.normal)/r_in.direction.length();
			}
			else
			{
				outward_normal=result.normal;
				ni_over_nt=1.0/ref_idx;
				cosine=-1*(r_in.direction.dot(result.normal))/r_in.direction.length();
			}
			if (refract(r_in.direction,result.normal,ni_over_nt,refracted))	reflect_prob=schlick(cosine,ref_idx);
			else reflect_prob=1.0;
			if (myRandom.C11Rand()<reflect_prob)
			{
				Scattered=Ray(result.pos,reflected);
			}
			else
			{
				Scattered=Ray(result.pos,refracted);
			}
			return true;
		}
		float ref_idx;
};

struct light : Material  {
    public:
        light(Color a) : emit(a) {}
        virtual bool Scatter(Ray r_in,IntersectResult result,Vec3& attenuation,Ray &Scattered)
        {
        	return false;
		}
        virtual Color emitted() 
		{ 
			return emit;
		}
        Color emit;
};

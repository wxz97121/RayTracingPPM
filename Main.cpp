#include <iostream>
#include <fstream>
#include <vector>
#include "Wxzmath.h"
using namespace std;
int LoopTimes;
Color CalcLoop(Ray r,Scene* world,int depth)
{
	IntersectResult result=world->intersect(r);
	if (result.geometry!=nullptr)
	{
		Ray Scattered;
		Vec3 Attenuation;
		Color Emitted=result.geometry->material->emitted();
		if (depth<LoopTimes && result.geometry->material->Scatter(r,result,Attenuation,Scattered))
		{
			return Emitted+CalcLoop(Scattered,world,depth+1)*Attenuation;
		}
		else
		{
			return Emitted;
		}
	}
	else
	{
//		return Color(0,0,0);
		Vec3 unit=r.direction.normalize();
		double t=0.5*(unit.y+1.0);
		return (1-t)*Color(1,1,1)+t*Color(0.5,0.7,1.0);
	}
}
void InitScene(Scene &myWorld)
{
	auto Sphere1=new Sphere(Vec3(0,1,0),1);
	Sphere1->material=new Lambertian(Vec3(0.1,0.2,0.5));
	myWorld.Add(Sphere1);
	
	auto Sphere2=new Sphere(Vec3(0,-100,0),100);
	Sphere2->material=new Lambertian(Vec3(0.8,0.8,0));
	myWorld.Add(Sphere2);
	
	auto Sphere3=new Sphere(Vec3(2,1,-1),1);
	Sphere3->material=new Metal(Vec3(0.8,0.6,0.2),0.2);
	myWorld.Add(Sphere3);
		
	auto Sphere4=new Sphere(Vec3(-2,1,-1),1);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere4);
	
	auto Sphere5=new Sphere(Vec3(-2,1,-1),-0.95);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere5);
	
	auto myHeart=new Heart(Vec3(0,0,-1),Vec3(0,1.5,-2),0.02);
	myHeart->material=new Lambertian(Vec3(0.9,0.05,0.05));
	myWorld.Add(myHeart);
	

}
void InitSceneWithLight(Scene &myWorld)
{
	auto Sphere1=new Sphere(Vec3(0,1,-1),1);
	Sphere1->material=new Lambertian(Vec3(0.1,0.2,0.5));
	myWorld.Add(Sphere1);
	
//	auto Bottom=new Plane(Vec3(0,1,0),-1);
//	Bottom->material=new Lambertian(Vec3(0.8,0.8,0));
//	myWorld.Add(Bottom);

	auto Sphere2=new Sphere(Vec3(0,-1000,0),1000);
	Sphere2->material=new Lambertian(Vec3(0.8,0.8,0));
	myWorld.Add(Sphere2);
	
	auto Left=new Plane(Vec3(1,0,0),-4);
//	Left->material=new Metal(Vec3(0.9,0.9,0.9),0.1);
	Left->material=new Lambertian(Vec3(0,0.8,0.8));
	myWorld.Add(Left);
	
	auto Right=new Plane(Vec3(-1,0,0),-4);
//	Right->material=new Metal(Vec3(0.9,0.9,0.9),0.1);
	Right->material=new Lambertian(Vec3(0,0.8,0.8));
	myWorld.Add(Right);
	
	auto Sphere3=new Sphere(Vec3(2,1,0),1);
	Sphere3->material=new Metal(Vec3(0.8,0.6,0.2),0.2);
	myWorld.Add(Sphere3);
		
	auto Sphere4=new Sphere(Vec3(-2,1,0),1);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere4);
	
	auto Sphere5=new Sphere(Vec3(-2,1,2),-0.95);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere5);
	
//	auto myHeart=new Heart(Vec3(0,0,-1),Vec3(0,1.5,-2),0.02);
//	myHeart->material=new Lambertian(Vec3(0.9,0.05,0.05));
//	myWorld.Add(myHeart);
	
	auto myTop=new Plane(Vec3(0,-1,0),-5);
	myTop->material=new Lambertian(Vec3(0.7,0.3,0.45));	
	myWorld.Add(myTop);
	
	auto myLight=new Sphere(Vec3(0,6,0),1.5);
	myLight->material=new light(Color(1,1,1));	
	myWorld.Add(myLight);
	
	auto myPlane=new Plane(Vec3(0,0,-1),-6); 
//	myPlane->material=new Metal(Vec3(0.9,0.9,0.9),0.1);
	myPlane->material=new Lambertian(Vec3(0.8,0.1,0.1));
	myWorld.Add(myPlane);
}
void RandomScene(Scene &myWorld)
{
    for (int a = -5; a < 5; a++) {
        for (int b = -2.5; b < 5; b++) {
            float choose_mat =myRandom.C11Rand();
            Vec3 center((double)a+0.9*myRandom.C11Rand(),0.2,(double)b+0.9*myRandom.C11Rand()); 
            if ((center-Vec3(2,1,-1)).length() > 1 && (center-Vec3(-2,1,-1)).length()>1)  
			{ 
			    auto tmpSphere = new Sphere(center, 0.2);
				tmpSphere->material=new Lambertian(Vec3(myRandom.C11Rand()*myRandom.C11Rand(), myRandom.C11Rand()*myRandom.C11Rand(),myRandom.C11Rand()*myRandom.C11Rand()));
                if (choose_mat < 0.8)
					tmpSphere->material=new Lambertian(Vec3(myRandom.C11Rand()*myRandom.C11Rand(), myRandom.C11Rand()*myRandom.C11Rand(),myRandom.C11Rand()*myRandom.C11Rand()));
                else if (choose_mat < 0.95)
                    tmpSphere->material=new Metal(Vec3(0.5*(1 + myRandom.C11Rand()), 0.5*(1 + myRandom.C11Rand()), 0.5*(1 + myRandom.C11Rand())),  0.5*myRandom.C11Rand());
                else tmpSphere->material=new dielectric(1.5);
                myWorld.Add(tmpSphere);
            }
        }
    }
//    auto Sphere4=new Sphere(Vec3(0,1,0),1);
//	Sphere4->material=new dielectric(1.5);
//	myWorld.Add(Sphere4);
//	
//	auto Sphere3=new Sphere(Vec3(2,1,0),1);
//	Sphere3->material=new Metal(Vec3(0.5,0.3,0.2),0.2);
//	myWorld.Add(Sphere3);
//	
//	auto Sphere1=new Sphere(Vec3(0,1,-2),1);
//	Sphere1->material=new Lambertian(Vec3(0.5,0.7,0.9));
//	myWorld.Add(Sphere1);
	
	auto Sphere2=new Sphere(Vec3(0,-100,0),100);
	Sphere2->material=new Lambertian(Vec3(0.8,0.8,0));
	myWorld.Add(Sphere2);
	
	auto Sphere3=new Sphere(Vec3(2,1,-1),1);
	Sphere3->material=new Metal(Vec3(0.8,0.6,0.2),0.2);
	myWorld.Add(Sphere3);
		
	auto Sphere4=new Sphere(Vec3(-2,1,-1),1);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere4);
	
	auto Sphere5=new Sphere(Vec3(-2,1,-1),-0.95);
	Sphere4->material=new dielectric(1.5);
	myWorld.Add(Sphere5);

}
int main()
{
	int nx=1000;
	int ny=1000;
	ios::sync_with_stdio(false); 
	cout<<"Please Input the loop times of each Ray."<<endl;
	cin>>LoopTimes;
	cout<<endl; 
	if (LoopTimes<1) LoopTimes=1;
	if (LoopTimes>1000) LoopTimes=1000;
	
	cout<<"Please Input the width of the image."<<endl;
	cin>>nx;
	cout<<endl;
	if (nx>=4000 || nx<=0) nx=1000;
	
	cout<<"Please Input the height of the image."<<endl;
	cin>>ny;
	cout<<endl;
	if (ny>=4000 || ny<=0) ny=1000;
	
	int SelectScene;
	cout<<"Please Input A INTEGET to choose a Scene."<<endl;
	cin>>SelectScene;
	Scene myWorld;
	if (SelectScene==0) InitScene(myWorld);
	else if(SelectScene==1) InitSceneWithLight(myWorld);
	else RandomScene(myWorld);
	cout<<endl;
	
	int Input=0;
	double Blur=0;
	cout<<"Please Input 1 if you want a image with Defocus Blur"<<endl;
	cin>>Input;
	if (Input==1) Blur=0.02;else Blur=0; 

	cout<<"Calculating..."<<endl;
	ofstream out("Output.PPM");
//	freopen("Output.PPM","w",stdout);

	int ns=200;
	out<<"P3"<<endl<<nx<<" "<<ny<<endl<<255<<endl;

	P_Camera myCamera(Vec3(0,1,-4),Vec3(0,0,1),Vec3(0,1,0),100,Blur);
//	Camera* myCamera=new Camera();
	int percent=0;
	for(int j=ny-1;j>=0;j--)
	{
		for(int i=0;i<nx;i++)
		{
			int jj=ny-1-j;
			if ((int)((jj*nx+i)/((double)nx*ny)*100)>percent)
			{
				percent=(int)((jj*nx+i)/((double)nx*ny)*100);
				cout<<percent<<"%..."<<endl;
			}
			Color myColor[1005];
			#pragma omp parallel for
			for(int s=0;s<ns;s++)
			{
				double u=double(i+myRandom.C11Rand())/double(nx);
				double v=double(j+myRandom.C11Rand())/double(ny);
//				Ray r=myCamera->get_ray(u,v);
				Ray r=myCamera.generateRay(u,v); 
				//IntersectResult TmpResult=myWorld.intersect(r);
				myColor[s]=CalcLoop(r,&myWorld,0);
			}
			Color FinColor(0,0,0);
			for(int s=0;s<ns;s++)
				FinColor=FinColor+myColor[s];
			FinColor=FinColor*(1/double(ns));
			int ir=int(255.99*sqrt(FinColor.r));
			int ig=int(255.99*sqrt(FinColor.g));
			int ib=int(255.99*sqrt(FinColor.b));
			out<<ir<<" "<<ig<<" "<<ib<<endl;
		}
	}
	return 0;
}

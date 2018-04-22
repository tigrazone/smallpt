#define _USE_MATH_DEFINES
#define erand48(dummy) (double(rand()) / RAND_MAX)

#define SPHERES 150

#include <vector>
using std::vector;

//grid element for object list
typedef struct {
	unsigned int idx;
	unsigned int objID;
} gridARRel;

#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
  double x, y, z;                  // position, also color (r,g,b)
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
  Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
  double rad;       // radius
  Vec p, e, c;      // position, emission, color
  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
  Sphere() {}
  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
    if (det<0) return 0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};
vector<Sphere> spheres(SPHERES);

//one grid element and array of all objects conneced to cells
gridARRel gridARRelem;
vector<gridARRel> gridARR;
vector<Vec> bbxmin(SPHERES);
vector<Vec> bbxmax(SPHERES);

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
  double n=SPHERES, d, inf=t=1e20;
  for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
  return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi){
  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return Vec(); // if miss, return black
  const Sphere obj = spheres[id];        // the hit object
  Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
  double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
  if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return obj.e; //R.R.
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
    Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
    return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
  } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
  bool into = n.dot(nl)>0;                // Ray from outside going in?
  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + f.mult(radiance(reflRay,depth,Xi));
  Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
  return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
int main(int argc, char *argv[]){
  int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
  Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
  Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];
  
  //build scene
  
        spheres[0] = Sphere(1e5, Vec(1e5+1,40.8,81.6),   Vec(),Vec(.75,.25,.25),DIFF); // Left
        spheres[1] = Sphere(1e5, Vec(-1e5+99,40.8,81.6), Vec(),Vec(.25,.25,.75),DIFF); // Right
        spheres[2] = Sphere(1e5, Vec(50,40.8, 1e5),      Vec(),Vec(.75,.75,.75),DIFF); // Back
        spheres[3] = Sphere(1e5, Vec(50,40.8,-1e5+170),  Vec(),Vec(),           DIFF); // Front
        spheres[4] = Sphere(1e5, Vec(50, 1e5, 81.6),     Vec(),Vec(.75,.75,.75),DIFF); // Bottom
        spheres[5] = Sphere(1e5, Vec(50,-1e5+81.6,81.6), Vec(),Vec(.75,.75,.75),DIFF) ;// Top
        spheres[6] = Sphere(600, Vec(50,681.6-.27,81.6), Vec(12,12,12),  Vec(), DIFF); // Light
        spheres[7] = Sphere(16.5,Vec(73,16.5,78),        Vec(),Vec(1,1,1)*.999, REFR); // Glas
        spheres[8] = Sphere(16.5,Vec(27,16.5,47),        Vec(),Vec(1,1,1)*.999, SPEC) ;// Mirror
		
        for (int s = 9; s < SPHERES; ++s) {
            // Create weird random spheres
            float radius = 1.0f + 2.0f * erand48(Xi);
            Vec pos = Vec(erand48(Xi) * 100.0 , erand48(Xi) * 100.0 , erand48(Xi) * 100.0 + 50.0);
            Vec c = Vec(erand48(Xi) * 0.8f + 0.1f, erand48(Xi) * 0.8f + 0.1f, erand48(Xi) * 0.8f + 0.1f);
            float reflParam = erand48(Xi);
            Refl_t rt = reflParam < 0.2f ? REFR : (reflParam > 0.8f ? SPEC : DIFF);
        
            spheres[s] = Sphere(radius, pos, Vec(0.05, 0.05, 0.05),  c, rt);
        }
		
		
	//calc bounding box of all objects
	Vec Gbbxmin, Gbbxmax;
		
	//calculate mins and maxs of all objects	
        for (int s = 0; s < SPHERES; ++s) {
            bbxmin[s] = spheres[s].p - spheres[s].rad;
            bbxmax[s] = spheres[s].p + spheres[s].rad;
			
			if(!s)
			{
				Gbbxmin = bbxmin[s];
				Gbbxmax = bbxmax[s];
			}
			else
			{
				if(bbxmin[s].x<Gbbxmin.x) Gbbxmin.x = bbxmin[s].x;
				if(bbxmin[s].y<Gbbxmin.y) Gbbxmin.y = bbxmin[s].y;
				if(bbxmin[s].z<Gbbxmin.z) Gbbxmin.z = bbxmin[s].z;
				
				if(bbxmax[s].x>Gbbxmax.x) Gbbxmax.x = bbxmax[s].x;				
				if(bbxmax[s].y>Gbbxmax.y) Gbbxmax.y = bbxmax[s].y;				
				if(bbxmax[s].z>Gbbxmax.z) Gbbxmax.z = bbxmax[s].z;				
			}
        }
		
		
	//calc number of cells of grid by objects count
	Vec BBlength;
    BBlength.x = ( Gbbxmax.x - Gbbxmin.x );
    BBlength.y = ( Gbbxmax.y - Gbbxmin.y );
    BBlength.z = ( Gbbxmax.z - Gbbxmin.z );
    
    float s = pow( ( ( BBlength.x * BBlength.y * BBlength.z ) / SPHERES ), 1.0f / 3.0f );
    
    int nx = ( int ) ( BBlength.x / s +0.5f);
    int ny = ( int ) ( BBlength.y / s +0.5f);
    int nz = ( int ) ( BBlength.z / s +0.5f);
	
	printf("SPHERES=%d\n", SPHERES);
	printf("GRID nx=%d ny=%d nz=%d\n", nx, ny, nz);
	
	/*
	
	
    
    gridPARAMnx = Aclamp( ( int ) ( BBlength.x / gridPARAMs +0.5f), 1, 128);
    gridPARAMny = Aclamp( ( int ) ( BBlength.y / gridPARAMs +0.5f), 1, 128);
    gridPARAMnz = Aclamp( ( int ) ( BBlength.z / gridPARAMs +0.5f), 1, 128);
	
	
        voxelwx = BBlength.x / gridPARAMnx;
        voxelwy = BBlength.y / gridPARAMny;
        voxelwz = BBlength.z / gridPARAMnz;
        invVoxelwx = 1.0f / voxelwx;
        invVoxelwy = 1.0f / voxelwy;
        invVoxelwz = 1.0f / voxelwz;		
		
		invNx = 1.0f / gridPARAMnx;
        invNy = 1.0f / gridPARAMny;
        invNz = 1.0f / gridPARAMnz;
	
	
	gridARRel gridARRelem;
	
	
	//создать массив grid
	gridARRsz = (sphereCount+2)*6; //приблизительно
	gridARR = (gridARRel*)malloc( gridARRsz * sizeof(gridARRel) );
	gridARRused = 0;
		
	for(i = 0; i < sphereCount; i++) 
	{
		//printf("minX %d lid=%d w=%.5f\n", i, bbx_order[i], bboxes[ bbx_order[i] ].mn.x);		
	
		vassign(abound.mn, bboxes[i].mn);
		vassign(abound.mx, bboxes[i].mx);
			
			vsub(abound.mn, abound.mn, bound.mn);
			vsub(abound.mx, abound.mx, bound.mn);			
			
			int minx = ( int )floor( abound.mn.x * invVoxelwx );
			int maxx = ( int )ceil( abound.mx.x * invVoxelwx );
			int miny = ( int )floor( abound.mn.y * invVoxelwy );
			int maxy = ( int )ceil( abound.mx.y  * invVoxelwy );
			int minz = ( int )floor( abound.mn.z  * invVoxelwz );
			int maxz = ( int )ceil( abound.mx.z  * invVoxelwz );
			
			gridARRelem.objID = i;
			
			for ( int z = minz; z < maxz; z++ )
			{
				for ( int y = miny; y < maxy; y++ )
				{
					for ( int x = minx; x < maxx; x++ )
					{
						gridARRelem.idx = ( z * gridPARAMny * gridPARAMnx ) + ( y * gridPARAMnx ) + x;
						gridARRelem.ordr = gridARRused;
						
						if(gridARRused==gridARRsz)
						{
							//printf("!418\n");
							gridARRsz += gridARRsz;
							gridARR = (gridARRel *)realloc(gridARR, sizeof(gridARRel) * gridARRsz);							
						}
						
						
						gridARR[gridARRused] = gridARRelem;
						gridARRused++;
					}
				}
			}
	}
	
	*/
	
	
	
	
	
  
  for (int y=0; y<h; y++){                       // Loop over image rows
    fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
    for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
      for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
        for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
          for (int s=0; s<samps; s++){
            double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
            double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
            Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                    cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
            r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
        }
  }
  FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}

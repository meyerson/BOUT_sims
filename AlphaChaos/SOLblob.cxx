/*******************************************************************
 * Bare Bones Bout with Python
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <mpi.h>
#include <math.h>

#include <typeinfo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <invert_laplace_gmres.hxx>
#include <boutexception.hxx>

//#include "./StandardAlpha.h"
//#include <python2.6/Python.h>

// Evolving variables 
Field3D u, n; //vorticity, density

//derived variables
Field3D phi,brkt;
int phi_flags;

//other fields
Field3D test1, test2, ReyN, ReyU;
Field2D n0;

//Constrained 
Field3D C_phi;

//other params
BoutReal nu, mu,gam, beta,alpha_c;

Field3D alpha, temp,edgefld,alpha_s;
//solver options
bool use_jacobian, use_precon;

//experimental
bool use_constraint;
bool chaosalpha;

int MZ;

FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner


//BoutReal alphamap(BoutReal x,BoutReal z);
BoutReal alphamap(double x, double Lx, double y,double Ly,
		  double k=1.20,double q0=3.0,double R=100,
		  int max_orbit =4000,double period=1.0,
		  bool count_turn = 0);

BoutReal Ullmann(double x, double Lx, double y,double Ly);
BoutReal Newton_root(double x_in,double y_in,double b = 50.0, double C=.01, double m = 3.0);

const Field3D smooth_xz(const Field3D &f); 

//int alphamapPy();
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int physics_init(bool restarting)
{
  


  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("physics");
  Options *solveropts = globaloptions->getSection("solver");

  OPTION(options, phi_flags, 0);
  //OPTION(options, alpha,3e-5);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options,chaosalpha,false);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  OPTION(options, beta, 6e-4);

  OPTION(globaloptions,MZ,33);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);


  bout_solve(u, "u");
  comms.add(u);
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //Laplacian *lap = Laplacian::create();
  
  bout_solve(n, "n");
  comms.add(n);

  //brute force way to set alpha
  
  if (chaosalpha){
    alpha.allocate();
    BoutReal ***a = alpha.getData();
    BoutReal edge[mesh->ngz];
    BoutReal zoomfactor = 3.0;
    BoutReal lowR = .55;
    BoutReal Lxz = 0;
    //the original 
    for(int jz=0;jz<mesh->ngz;jz++) 
      for(int jx=0;jx<mesh->ngx;jx++){
	Lxz = Ullmann(mesh->GlobalX(jx),1.0,mesh->dz*jz,mesh->zlength);
	for(int jy=0;jy<mesh->ngy;jy++){
	  //output << jy <<endl;
	  //a[jx][jy][jz]=alphamap(lowR+mesh->GlobalX(jx)/zoomfactor,1.0,mesh->dz*jz,mesh->zlength,1.0,4.0,160.0,1000);	
	  a[jx][jy][jz]=Lxz;
	}
      }
    

    alpha = (2.0*.1)/alpha;
    alpha_s = lowPass(alpha,0);
    dump.add(alpha,"alpha",0);
    dump.add(alpha_s,"alpha_smooth",0);
    
    //remap  - NO!
  } else{
    OPTION(options, alpha_c,3e-5);
    alpha = alpha_c;
  }
  
  dump.add(brkt,"brkt",1);
  dump.add(test1,"test1",1);
  dump.add(test2,"test2",1);
  dump.add(ReyN,"ReyN",1);
  dump.add(ReyU,"ReyU",1);
  

  if (use_constraint){
    //solver->setPrecon(precon_phi);
    //solver->setJacobian(jacobian_constrain);
    phi.setBoundary("phi");
    //bout_solve(phi,"phi");
  }else
    dump.add(phi,"phi",1);

  comms.add(phi); //super duper important 

  if (use_jacobian)
    solver->setJacobian(jacobian);

  // if (use_precon)
  //   solver->setPrecon(precon);
    
  output.write("use jacobian %i \n",use_jacobian);
  output.write("use precon %i \n",use_precon);
  output.write("DONE WITH PHYSICS_INIT\n");

  n0 = 1.0;
  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))

int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  
  static Field2D A = 0.0;
  static Field2D C = 1e-24;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);

  phi.applyBoundary("dirichlet");
  // Density
  //f = lowPass(f,1);
  //f = lowPass(g,1);
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
 


  
  ReyU = bracket3D(phi,u)/(nu*LapXZ(u)+1e-5);

 
  ddt(u) -= bracket3D(phi,u);
  ddt(u) += alpha * phi;
  ddt(u) += nu * LapXZ(u);
  ddt(u) += beta* DDZ(n+n0)/(n+n0);

  
  ReyN = bracket3D(phi,n)/(mu * LapXZ(n)+1e-5);
  
  ddt(n)  -= bracket3D(phi,n+n0);
  ddt(n) += mu * LapXZ(n+n0);

 
  return 0;
}


const Field3D mybracket(const Field3D &phi, const Field3D &A)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;

  //output.write("mesh->Bxy = %e\n", (mesh->J*sqrt(mesh->g_22))[2][2]);

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi); 
    
    #pragma omp section
    dpdy = DDY(phi);
    
    #pragma omp section
    dpdz = 0;
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {

    #pragma omp section
    vx = mesh->g_23*dpdz - mesh->g_33*dpdy; 
    //vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_33*dpdx - mesh->g_13*dpdz;
      //vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    //vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
    vz =0;
  }


  // Upwind A using these velocities
  
  Field3D ry, rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);
    
    #pragma omp section
    rz = VDDZ(vz, A);
  }
  //output.write("mesh->g_22: %g  \n" ,vx[4][4]);
  //result = (ry + rz); //usually  convention
  result = (result + ry) / (mesh->J*sqrt(mesh->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}



/* computes Jv, where ddt() is holding v and the system state holds Jv */ 
int jacobian(BoutReal t) {
  mesh->communicate(ddt(u),ddt(n));
  
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 

  mesh->communicate(ddt(phi));

  u=0;
  n=0;

  //u -= mybracket(ddt(phi),ddt(u));
  u += bracket3D(ddt(phi),ddt(u));
  //ddt(u) += alpha * phi;
  u += nu * Delp2(ddt(u));
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  //u -= Grad_par(ddt(n)); 
  //ddt(u).applyBoundary("dirichlet");
  u -= beta* DDZ(n); 
  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 
  n += bracket3D(ddt(phi),ddt(n));
  //n -= mybracket(ddt(phi),ddt(n));
  n += mu * Delp2(ddt(n));
  //n -= alpha* n;
  
  n.applyBoundary();
  u.applyBoundary();
  return 0;
}

/* computes P^-1 r, where ddt() holds (-r) and the system state hold P^-1 r*/

int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  // mesh->communicate(rhscomms);
  mesh->communicate(ddt(n),ddt(u));

  n= 0;
  u =0;

  n += ddt(n);
  // mesh->communicate(n);
  // Ni -= (mesh->Bxy*mesh->Bxy*ddt(Ni) - ddt(rho))/(mesh->Bxy*mesh->Bxy);

  u += ddt(u);
  u -= gamma * Grad_par(ddt(n)); 
 
return 0;
 
}
const Field3D smooth_xz(const Field3D &f){
  Field3D result;
  result.allocate();
  result = f;
  for(int x=2;x<mesh->ngx-2;x++)
    for(int y=0;y<mesh->ngy;y++)
      for(int z=2;z<mesh->ngz-2;z++) {
        result[x][y][z] = 0.5*f[x][y][z] + 0.125*( 0.5*f[x+1][y][z] + 0.125*(f[x+2][y][z] + f[x][y][z] + f[x+1][y][z-1] + f[x+1][y][z+1]) +
                                                   0.5*f[x-1][y][z] + 0.125*(f[x][y][z] + f[x-2][y][z] + f[x-1][y][z-1] + f[x-1][y][z+1]) +
                                                   0.5*f[x][y][z-1] + 0.125*(f[x+1][y][z-1] + f[x-1][y][z-1] + f[x][y][z-2] + f[x][y][z]) +
                                                   0.5*f[x][y][z+1] + 0.125*(f[x+1][y][z+1] + f[x-1][y][z+1] + f[x][y][z] + f[x][y][z+2]));
      }

  mesh->communicate(result);
  return result;
}

BoutReal Ullmann(double x, double Lx, double y,double Ly){
  x = x*(1/Lx);
  y = y*(2.0*M_PI/Ly);
  
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q;

  int max_orbit = 1000;
  
  double L = 0.0;
  double eps = .2;
  double aa = -.05;
  int m = 3;
  double l = 10.0;
  double R = 100;
  double q0 = 3.0;
  double b = 50.0;
  double a= 45.0;
  double x_new;
  double y_new;
  double x_new2;

  double C = ((2*m*l*pow(a,2.0))/(R*q0*pow(b,2.0)))*eps;
     
  while(count < max_orbit and not hit_divert){
    // x_new = x;
    // y_new = y;
    x_new = x/(1-aa*sin(y));
    q = q0*pow((x_new/a),2.0);

    y_new =  (y+ 2*M_PI/q + aa*cos(y));
    y_new = fmod(y_new,2*M_PI);
    
    x_new2 = Newton_root(x_new,y_new,b,C,m);
    
    output<< (-1.0*x_new+ x_new2 +(m*b*C)/(m-1.0) * pow(x_new2/b, m-1.0) * sin(m*y_new))<<endl;

    output<< "old: " <<(-1.0*x_new+ x_new +(m*b*C)/(m-1.0) * pow(x_new/b, m-1.0) * sin(m*y_new))<<endl;

    //chi = (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    //x_new2 = (newton_krylov(func,x_new));
    y_new = (y_new - C*pow(x_new2/b , m-2) * cos(m*y_new));

    //output <<x_new<<endl;
    hit_divert = x_new2>b;
    count++;
    L = L +1.0;
  }
  
  return L;
}

BoutReal alphamap(double x, double Lx,double y,double Ly,
		  double k,double q0 ,double R,int max_orbit,double period,
		  bool count_turns){ 
  

  
  x = x*(2.0*M_PI/Lx)-2.0*M_PI;
  y = y*(2.0*M_PI/Ly);
  //output << "[" << x  << "], "<<endl;
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q;

  double L = 0.0;
  
  //count_turns = true;
  while(count < max_orbit and not hit_divert){
    inSOL = x>0.0; //are we in SOL now?
    
    //update the field line
    x = x+k*sin(period * y);
    y = fmod((x+y),(2.0*M_PI)); //one can argue that x = 2*M_PI*fmod(q(R),1)

    q =q0; //+x/(2*M_PI);
   
    if (count_turns)
      L = L +1.0;
    else
      L = L + q *(R/3.0)*2.0*M_PI; //just assume that the minor radius is 1/3 the major radius, ~ DIII-D
   
    hit_divert = x>0.0;
   
    count++;
    }
  // output << L << endl;
  return L;

}


BoutReal Newton_root(double x_in,double y_in,double b, double C, double m){
  double x_out = x_in;
  double atol = 1.0;
  int iter = 0;
  int max_iter = 100;
  double f;
  double J;

  while (atol > .01 && iter < max_iter)
    {
      f = (-1.0*x_in + x_out +(m*b*C)/(m-1.0) * pow(x_out/b, m-1.0) * sin(m*y_in));
      J = 1.0 + (m*C) * pow(x_out/b,m-2.0) * sin(m*y_in);
      atol = fabs(f);
      x_out = x_out - 1.0*f/J;
      iter++;
      
    }
  //output<< (x_out - x_in)/fabs(x_in)<<endl;
  return x_out;
}

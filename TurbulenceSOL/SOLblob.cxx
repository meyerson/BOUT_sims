/*******************************************************************
 * Bare Bones Bout with Python
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <mpi.h>


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
#include <field_factory.hxx>
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

Field3D alpha, temp,edgefld,alpha_smooth, source;
//solver options
bool use_jacobian, use_precon;

bool withsource;

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

const Field3D smooth_xz(const Field3D &f); 
const Field3D smooth_x2(const Field3D &f,const int sig_pixel); 

const Field3D remap(const Field3D &f, const BoutReal lowR, const BoutReal zoomfactor); 
//int alphamapPy();
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply
#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))
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
  OPTION(options, alpha_c,3e-5);

  OPTION(globaloptions,MZ,33);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);

  OPTION(options,withsource,false);


  bout_solve(u, "u");
  comms.add(u);
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //phi.Appl
  //Laplacian *lap = Laplacian::create();
  
  bout_solve(n, "n");
  comms.add(n);
  
  FieldFactory f(mesh);
  if(withsource){
    //initial_profile("source", v);
    source = f.create3D("gauss(x-0.0,0.02)");
    dump.add(source,"source",0);
    
  }
  //brute force way to set alpha
  if (chaosalpha){
    alpha  = 0;
  }
  else {
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



int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  
  static Field2D A = 0.0;
  static Field2D C = 1e-24;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);

  phi.applyBoundary("neumann");
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
  //ddt(u) += beta* DDZ(n+n0)/(n+n0);
  ddt(u) += beta* DDZ(n);
  
  ReyN = bracket3D(phi,n)/(mu * LapXZ(n)+1e-5);
  
  ddt(n)  -= bracket3D(phi,n);
  ddt(n) += mu * LapXZ(n);
 
  if(withsource)
    ddt(n) += source;

  //apply the boundary    
  n.applyBoundary();
  u.applyBoundary();

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

//const Field3D fftsm


const Field3D smooth_xz(const Field3D &f){
  Field3D result;
  result.allocate();
  result = f;
  for(int x=2;x<mesh->ngx-2;x++)
    for(int y=0;y<mesh->ngy;y++)
      for(int z=0;z<mesh->ngz;z++) {
        result[x][y][z % MZ] = 0.5*f[x][y][z % MZ] + 0.125*( 0.5*f[x+1][y][z % MZ] + 0.125*(f[x+2][y][z % MZ] + f[x][y][z % MZ] + f[x+1][y][(z-1) % MZ] + f[x+1][y][(z+1) % MZ]) +
							     0.5*f[x-1][y][z % MZ] + 0.125*(f[x][y][z % MZ] + f[x-2][y][z % MZ] + f[x-1][y][(z-1) % MZ] + f[x-1][y][(z+1) % MZ]) +
							     0.5*f[x][y][(z-1) % MZ] + 0.125*(f[x+1][y][(z-1) % MZ] + f[x-1][y][(z-1) % MZ] + f[x][y][(z-2) % MZ] + f[x][y][z % MZ]) +
							     0.5*f[x][y][(z+1) % MZ] + 0.125*(f[x+1][y][(z+1) % MZ] + f[x-1][y][(z+1) % MZ] + f[x][y][z % MZ] + f[x][y][(z+2) % MZ]));

	// if (!finite(result[z][y][z % MZ])){
	//   result[x][y][z % MZ] = 0.0;
	  //z--;
      }


  mesh->communicate(result);
  return result;
}

const Field3D smooth_x2(const Field3D &f,const int sig_pixel){
  Field3D result,fs;
  result.allocate();
  fs = f;
  result = 0;
  //int flt_sx = int((mesh->ngx)/4.0);
  int flt_sx = 2.0*sig_pixel;
  //int flt_sz = int((mesh->ngz)/5.0);
  //int flt_sz = 10.0;
  //BoutReal weights[] = {1, 1, 1, 1,1,1,1,.75, .75, .75, .5, .5, .5, .25, .25, .25, .125,.125,.125};
  BoutReal weights[flt_sx];
  BoutReal weights_z[flt_sx];
  
  //BoutReal weights[] = {1.0,.5};
  BoutReal norm =0.0;
  for(int jw=0;jw<flt_sx;jw++){
    norm = norm + exp(-pow(jw,2.0)/(2.0*pow(flt_sx/2.0,2.0)));
    weights[jw] = exp(-pow(jw,2.0)/(2.0*pow(flt_sx/2.0,2.0)));
    weights_z[jw] = exp(-pow(jw,2.0)/(2.0*pow(flt_sx/2.0,2.0)));
    //output<<weights[jw]<<endl;
  }
  BoutReal ahead,behind;

  for(int jy=0;jy<mesh->ngy;jy++)
    for(int jz=0;jz<mesh->ngz;jz++) {
      for(int jx=0;jx<1.0; jx++){
  	result[jx][jy][jz] = fs[jx][jy][jz];
  	result[mesh->ngx-1-jx][jy][jz] = fs[mesh->ngx-1-jx][jy][jz];
      }
    }

  for(int jx=0;jx<mesh->ngx;jx++){
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	for (int jzz=-round(0.0);jzz<round(flt_sx);jzz++)
	  for(int jxx=0;jxx<flt_sx;jxx++){
	    //int jzz = 0;
	//while(jw < flt_sz) 
	    //output<< weights[0];
	    //weights[jxx] = (1.0/norm)*weights[jxx];
	    
	// if ((jx-jw)<0)
	  //   behind = fs[0][jy][jz];
	  // else
	  //   behind = fs[jx-jw][jy][jz];
	    
	    if ((jx-jxx)<0){
	      //behind = fs[0][jy][(jz+jzz)%(MZ)];
	      behind = 0;
	    }
	    else
	      behind = fs[jx-jxx][jy][(jz+jzz)%(MZ)];
	    
	    if ((jx+jxx)>(mesh->ngx-1)){
	      //ahead = fs[mesh->ngx-1][jy][(jz+jzz)%(MZ)];
	      ahead = 0.0;
	    }
	    else
	      ahead = fs[jx+jxx][jy][(jz+jzz)%(MZ)]; 

	    result[jx][jy][jz % MZ] = result[jx][jy][jz % MZ]+ (1.0/norm)*weights[jxx]*(ahead+behind);
	    //jxx++;
	  }
	
	
      }
    // mesh->communicate(result);
  }
  // for(int jy=0;jy<mesh->ngy;jy++)
  //   for(int jz=0;jz<mesh->ngz;jz++) {
  //     for(int jx=0;jx<1.0; jx++){
  // 	result[jx][jy][jz] = fs[jx][jy][jz];
  // 	result[mesh->ngx-1-jx][jy][jz] = fs[mesh->ngx-1-jx][jy][jz];
  //     }
  //   }

  mesh->communicate(result);
  return result;
}

const Field3D remap(const Field3D &f, const BoutReal lowR, const BoutReal zoomfactor){
  Field3D result;
  result.allocate();
  result = f;

 
  //find the edge
  BoutReal edge[mesh->ngz];
  BoutReal sumedge[mesh->ngz];
  //BoutReal localedge[mesh->ngz];
  //int edge[mesh->ngz];
  BoutReal mask[mesh->ngx];
  bool boolmask[mesh->ngx];
  //Field2D mask;
  // mask.allocate();
  
  BoutReal edge_val =0.0;


  int rank, size,i;
  MPI_Comm_rank(BoutComm::get(), &rank);
  MPI_Comm_size(BoutComm::get(), &size);
 
  int current_iter = 0;
  int max_iter = 1; //iterative improvement in the resolvign the border is an improvement
  BoutReal zoom = zoomfactor;
  BoutReal offset = 0;

  for(int z=0;z<mesh->ngz;z++)
    sumedge[z]=0.0; 

  while (current_iter < max_iter){
    offset = .5/zoom; //go to offset  = .25 at fist iter

    result = smooth_xz(result);
    result = smooth_xz(result);
    result = DDX(result)*DDX(result);
    mesh->communicate(result);
  for(int y=1;y<mesh->ngy-1;y++){
    //edge[z]=0.0;
    for(int z=0;z<mesh->ngz;z++){
      //reset for any fixed x,y pair
      edge_val = 0.0;
      edge[z]=0.0; 

      //find the value on edge
      for(int x=2;x<mesh->ngx-2;x++) {
	mask[x] = result[x][y][z];  
        if(result[x][y][z] > edge_val)
	  edge_val = result[x][y][z]; 
      }
      
      //find the value on global edge
      BoutReal localresult = edge_val;
      MPI_Allreduce(&localresult, &edge_val, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get() ); 

      //edge val is the global edge value now

      //get the displacement for the global edge, if not on current
      //proc keep edge[z] = 0, so all BUT the proc holding the edge will
      //have edge[z] = 0
      for(int x=0;x<mesh->ngx;x++)
	if (result[x][y][z] == edge_val)
	  edge[z]= mesh->GlobalX(x);

      
      //scatter the global edge
      BoutReal localedge = edge[z];
      //set the edge to the largest edge displacment value found
      //largest among locals = global

      MPI_Allreduce(&localedge, &edge[z], 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
 
    } // end z loop

    //smooth this edge
    for(int z=0;z<300*mesh->ngz;z++){
      //output<<z %MZ << "  "<<(z + 2) % MZ <<endl;
      edge[z%MZ] = edge[z % MZ] + .5*(edge[(z+1)%MZ] + edge[(z-1)%MZ])+.5*(edge[(z+3)%MZ]+edge[(z-3)%MZ])+ .5*(edge[(z+5)%MZ]+edge[(z-5)%MZ])+ .5*(edge[(z+8)%MZ]+edge[(z-8)%MZ]);
      edge[z%MZ] =  edge[z % MZ]/5.0;
    }
    for(int z=0;z<mesh->ngz;z++){
      if (current_iter == 0 )
	sumedge[z] = edge[z]/(zoom); 
      else
	sumedge[z] += sumedge[z]+edge[z]/(zoom); 
    }
    
    // //   edge[z-1] = .5*edge[z-1] + .25*(edge[z] + edge[z-2]);
    
    for(int z=0;z<mesh->ngz;z++)
      for(int x=0;x<mesh->ngx;x++){
	
    	result[x][y][z]=alphamap((lowR-offset)+sumedge[z] + mesh->GlobalX(x)/(zoom),1.0,mesh->dz*z,mesh->zlength,1.0,3.0,100,20,1);
	if (finite(result[x][y][z]))
	  if(result[x][y][z] <19.0)
	    result[x][y][z] = .00001;
      }
    
 
  } // end y loop
  current_iter++;
  
  //zoom *= 2.0;

  } // end while loop
  
  // for(int y=1;y<mesh->ngy-1;y++)
  //   for(int z=0;z<mesh->ngz;z++)
  //     for(int x=0;x<mesh->ngx;x++)
  //    	result[x][y][z]=alphamap((lowR-1.0/4.0-0.0/8.0)+sumedge[z] + mesh->GlobalX(x)/(zoom/2.0),1.0,mesh->dz*z,mesh->zlength,1.0,3.0,100,100);

  for(int y=1;y<mesh->ngy-1;y++)
    for(int z=0;z<mesh->ngz;z++)
      for(int x=0;x<mesh->ngx;x++)
     	//result[x][y][z]=alphamap((lowR-1.0/32.0-0.0/8.0)+sumedge[z] + mesh->GlobalX(x)/(8*zoom),1.0,mesh->dz*z,mesh->zlength,1.0,3.0,100,1000);
	result[x][y][z]=alphamap(lowR - (lowR)/(8*zoom) + sumedge[z]+mesh->GlobalX(x)/(14*zoom),1.0,mesh->dz*z,mesh->zlength,1.0,3.0,100,1000);

  return result;
}


BoutReal alphamap(double x, double Lx,double y,double Ly,
		  double k,double q0 ,double R,int max_orbit,double period,
		  bool count_turns){ 
  

  
  x = x*(2*M_PI/Lx)-2*M_PI;
  y = y*(2*M_PI/Ly);
  //output << "[" << x  << "], "<<endl;
  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q;

  double L = 0.0;
  
  //count_turns = true;
  while(count < max_orbit and not hit_divert){
    inSOL = x>0; //are we in SOL now?
    
    //update the field line
    x = x+k*sin(period * y);
    y = fmod((x+y),(2*M_PI)); //one can argue that x = 2*M_PI*fmod(q(R),1)

    q =q0 +x/(2*M_PI);
    //q = q0;

    //L = L+q*R*2*M_PI; //one can argue that for x>0, it should be L+q*R*M_PI
    if (count_turns)
      L = L+1.0;
    else
      L = L+q*R*2*M_PI;

    hit_divert = (inSOL and x>0) or (x>M_PI); //did the field line line stay in SOL?
    //cout <<count<<" L:" << L<<endl;
    count++;
    }
  
  return L;

}

/*******************************************************************
 * Bare Bones Bout
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <bout_types.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <invert_laplace_gmres.hxx>
//#include <invert_laplace_XY.hxx>>
#include <boutexception.hxx>

#include <inverter.hxx>
#include <full_gmres.hxx>
#include <utils.hxx>
//lower-level stuff
//extern BoutReal** XYmatrix(int nx, int ny);
//typedef double BoutReal;
// Evolving variables 
Field3D u, n; //vorticity, density

//background density
Field3D n0;

//derived variables
Field3D phi,brkt;
int phi_flags;

//other fields
Field3D test1, test2;


//Constrained 
Field3D C_phi, phibdry;

FieldPerp fu,pho;

//other params
BoutReal alpha, nu, mu,gam, beta;


//inverters
LaplaceGMRES *lapinv;
class Laplacian *lap;
InvertPar *inv; // Parallel inversion class

//solver options
bool use_jacobian, use_precon, use_rootfind;

//experimental
bool use_constraint;

int MZ;

FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner


int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int rootphophi(BoutReal t); //root finding
//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

void ToXZ(const Field3D &in,FieldPerp &out);
void ToXY(const FieldPerp &in,Field3D &out);


int physics_init(bool restarting)
{

  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("physics");
  Options *solveropts = globaloptions->getSection("solver");

  OPTION(options, phi_flags, 0);
  OPTION(options, alpha,0);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  OPTION(options, beta, 6e-4);

  OPTION(globaloptions,MZ,33);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);
  OPTION(solveropts,use_rootfind,false);
 
  bout_solve(u, "u");
  comms.add(u);

  lap  = Laplacian::create(globaloptions->getSection("fullLap"));
  lap->setCoefA(0);
  lap->setCoefC(1e-24);
  lap->setFlags(phi_flags);

  // Initialise parallel inversion class
  inv = InvertPar::Create();
  inv->setCoefA(1e-24);
  inv->setCoefB(1.0);

  //let's copu u to fu in a wierd way
  // BoutReal ***d = u.getData();
  // BoutReal **data  = rmatrix(mesh->ngx, mesh->ngy);

  // for(int jx=0;jx<mesh->ngx;jx++)
  //   for(int jy=0;jy<mesh->ngy;jy++)
  //     data[jx][jy] = d[jx][jy][0];

  // fu.setData(data);

  // ToXZ(u,fu);
  
  phi = lap->solve(u); 
  //phi = inv->solve(phi);
  // pho = lap->solve(fu);

  // ToXY(fu,test1);
  


  //phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //Laplacian *lap = Laplacian::create();
  //gam = full_gmres(u,Laplacian,phi,NULL,0);

  bout_solve(n, "n");
  comms.add(n);
  //comms.add(n0);
  //u.setBoundary("u");
  //n.setBoundary("n");

  //brkt = b0xGrad_dot_Grad(phi, u);

  //dump.add(phi,"phi",1);
  dump.add(brkt,"brkt",1);
  dump.add(test1,"test1",1);
  dump.add(test2,"test2",1);

  if (use_constraint){
    mesh->communicate(phi,u);
    //Add phi equation as a constraint
    if(!bout_constrain(phi, ddt(phi), "phi"))
      throw BoutException("Solver does not support constraints");

    //if (use_precon)
    //  solver->setPrecon(precon_phi);

    if (use_jacobian)
      solver->setJacobian(jacobian_constrain);

    phi.setBoundary("phi");
    phibdry.setBoundary("phi");
    
  }else {
    dump.add(phi,"phi",1);
    //dump.add(pho,"pho",1);
    
    if (use_jacobian)
      solver->setJacobian(jacobian);
    
    if (use_precon)
      solver->setPrecon(precon);

   
    if (use_rootfind){
      output.write("set the rootfindig func\n");
      //bout_rootfind(phi,ddt(phi),"phi"); //almost like a constraint
      //solver->setRootFind(rootphophi);
    }
  }
  
  comms.add(phi); //super duper important 

    
  output.write("use jacobian %i \n",use_jacobian);
  output.write("use precon %i \n",use_precon);
  output.write("DONE WITH PHYSICS_INIT\n");

  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXY(f)(mesh->g11*D2DX2(f) + mesh->g22*D2DY2(f))
 
int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  // output.write("g_33: %g\n",sqrt(mesh->g_33[2][2]));
  //output.write("g_22: %g\n",sqrt(mesh->g_22[2][2]));
  //output.write("J: %g\n",mesh->J[2][2]);

  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-14;
  static Field2D D = 1.0;
  //phi = lap->solve(u);
  //phi = lap->solve(u);
    
  // test1 = Grad2_par2(u);
  // test1 = lap->solve(test1);
  // test1 = lap->solve(test1);
  // test1.applyBoundary("neumann");
  
  // phi = phi - test1;
  //phi.applyBoundary("neumann");
    
  //mesh->communicate(phi,u);
  if (use_constraint){
    phibdry = phi;
    phibdry.applyBoundary();
    phibdry -= phi;
    //ddt(phi) = u - Laplacian(phi);
    //ddt(phi) = ddt(phi)^2;
    ddt(phi) = LapXY(phi)-u;
    //ddt(phi) = Laplacian(phi) - u;
    ddt(phi).setBoundaryTo(phibdry); //removes some boundary errors
  } else {

  //ToXY(pho,phi);
    phi = lap->solve(u);

    phi.applyBoundary();

  }
  

  //phi.applyBoundary("neumann");
  //phi.applyBoundary("dirichlet");
  // Density
  //f = lowPass(f,1);
  //f = lowPass(g,1);
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
 

  test1 = u - LapXY(phi);
  //test2 = mybracket(phi,DDX(n));
  //brkt = mybracket(phi,n);
  
  // output<<"ddt(u)\n";
  ddt(u) -= mybracket(phi,u);
  
  ddt(u) += alpha * phi;
  ddt(u) += nu * LapXY(u);
  //output<<"u2\n";
  //ddt(u).checkData();
  ddt(u) -= beta * DDY(n);
  //ddt(u) -= beta* DDZ(n); 
  //ddt(u) -= Grad_par(n); 
  

  // mesh->communicate(comms); //no don't do this here
  
  ddt(n) -= mybracket(phi,n);
  // output<<"n2\n";
  //ddt(n).checkData();
  //ddt(n)  += bracket3D(phi,n);
  ddt(n) += mu * LapXY(n);
  //ddt(n) += mu * Laplacian(n);
  //output<<"n3\n";
  //ddt(n).checkData();
  //ddt(n) -= alpha* n;
  
  //ddt(n) = lowPass(ddt(n),MZ/8);
  //output<<"n4\n";
  //ddt(n).checkData();
  //ddt(n).applyBoundary("dirichlet");
  //output<<"n5\n";
  //ddt(n).checkData();

  //ddt(n).setBoundaryTo(phibdry);
  //ddt(u).applyBoundary("neumann");

  mesh->communicate(ddt(n),ddt(u));
  //ddt(n) -= VDDZ(n,n) + mu*VDDY(u,n);

  // if (driver){
    
  // }

  

  //ddt(f) -= 10*f;
  //ddt(f) = lowPass(ddt(f),5);

  
 
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
    //vy = -1.0
      //vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    //vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
    vz =0;
  }

  mesh->communicate(vx,vy,vz);
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
  mesh->communicate(result,ry,rz);
  result = (result + ry)/ (mesh->J*sqrt(mesh->g_33));
  
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
  static Field2D C = 1e-4;
  static Field2D D = 1.0;
  
  ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 

  mesh->communicate(ddt(phi));

  u=0;
  n=0;

  u -= mybracket(ddt(phi),ddt(u));
  //u += bracket3D(ddt(phi),ddt(u));
  ddt(u) += alpha * phi;
  u += nu * Laplacian(ddt(u));
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  u -= Grad_par(ddt(n)); 
  //ddt(u).applyBoundary("dirichlet");
  //u -= beta* DDZ(n); 
  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 
  //n += bracket3D(ddt(phi),ddt(n));
  n -= mybracket(ddt(phi),ddt(n));
  n += mu * Laplacian(ddt(n));
  n += alpha* n;
  
  n.applyBoundary();
  u.applyBoundary();
  return 0;
}

int jacobian_constrain(BoutReal t) {
  mesh->communicate(ddt(u),ddt(n),ddt(phi));
  
  // static Field2D A = 0.0;
  // static Field2D C = 1e-12;
  // static Field2D D = 1.0;
  
  // ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 

  //mesh->communicate(ddt(phi));

  u=0;
  n=0;

  phi  = Laplacian(ddt(phi)) - ddt(u);
  
  phibdry = ddt(phi);
  phibdry.applyBoundary();
  phibdry -= ddt(phi); // Contains error in the boundary
  
  phi.setBoundaryTo(phibdry);
  
  u -= mybracket(ddt(phi),ddt(u));
  //u += bracket3D(ddt(phi),ddt(u));
  ddt(u) += alpha * phi;
  u += nu * Laplacian(ddt(u));
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  u -= Grad_par(ddt(n)); 
  //ddt(u).applyBoundary("dirichlet");
  //u -= beta* DDZ(n); 
  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 
  //n += bracket3D(ddt(phi),ddt(n));
  n -= mybracket(ddt(phi),ddt(n));
  n += mu * Laplacian(ddt(n));
  n += alpha* n;
  
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

int precon_phi(BoutReal t, BoutReal gamma, BoutReal delta) {
  // mesh->communicate(rhscomms);
  mesh->communicate(ddt(n),ddt(u),ddt(phi));

  n= 0;
  u =0;

  n += ddt(n);
  // mesh->communicate(n);
  // Ni -= (mesh->Bxy*mesh->Bxy*ddt(Ni) - ddt(rho))/(mesh->Bxy*mesh->Bxy);

  u += ddt(u);

  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  //phi = invert_laplace(ddt(phi)+ ddt(u), phi_flags,&A,&C,&D);
  //phi = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  phi = ddt(phi); 
  return 0;
 
}

// int rootphophi(BoutReal t){
//   output.write("stuff happens here\n");
//   mesh->communicate(ddt(n),ddt(u));
//   //n= 0;
//   //u =0;
//   //u += Laplacian(ddt(phi));
//   //u.applyBoundary();
//   return 0;
 
// }
int rootphophi(BoutReal t) {
  mesh->communicate(ddt(u),ddt(n));
  n= 0;
  u =0;
  //u += Laplacian(ddt(phi));
  
  return 0;
}

void ToXZ(const Field3D &in,FieldPerp &out)
{

  BoutReal **data  = rmatrix(mesh->ngx, mesh->ngy);
  BoutReal ***d = in.getData();

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++){
      data[jx][jy] = d[jx][jy][0];
      //output<<d[jx][jy][0]<<endl;
    }
  out.setData(data);
}

void ToXY(const FieldPerp &in, Field3D &out)
{

  //BoutReal **data  = rmatrix(mesh->ngx, mesh->ngy);
  BoutReal **d = in.getData();
  

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++){
      out.setData(jx,jy,0,&d[jx][jy]);
      //output<<d[jx][jy]<<endl;
    }
      //data[jx][jy][0] = d[jx][jy]; 

  
}

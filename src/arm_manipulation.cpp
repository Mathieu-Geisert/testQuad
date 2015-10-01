/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
  *    \file   examples/ocp/rocket_c.cpp
  *    \author Boris Houska, Hans Joachim Ferreau
  *    \date   2010
  */
#include <iostream>
#include <fstream>
#include <math.h>

#include <acado_optimal_control.hpp>
#include <acado_gnuplot.hpp>
#include <acado/process/process.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <pinocchio/multibody/model.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/kinematics.hpp>

#include <pinocchio/multibody/parser/urdf.hpp>

/* >>> start tutorial code >>> */


// -------------------------------------------------------------------------
//             C-STYLE DEFINITION OF THE PROBLEM DIMENSIONS:
// -------------------------------------------------------------------------


#define  NJ   1    // number of objective functions
#define  NX   13    // number of differential states
#define  NI   12    // number of initial value constraints
#define  NE   2    // number of end-point / terminal constraints
#define  NH   1    // number of inequality path constraints
#define NB_IT 20
#define TPS 10.0


// -----------------------------------------------------------------------------
//   UGLY C-STYLE DEFINITION OF THE OBJECTIVE, MODEL AND CONSTRAINT FUNCTIONS:
// -----------------------------------------------------------------------------

const double Cf = 0.00065;
const double d = 1;
const double c = 0.001;

se3::Model robot = se3::urdf::buildModel("/local/mgeisert/models/ur5_quad.urdf", true);
se3::Data data(robot);

Eigen::VectorXd q(13);
Eigen::VectorXd v(12);
Eigen::VectorXd u(12);
Eigen::VectorXd a(12), a3(3), x3(3), aa(3), gr(4);
Eigen::VectorXd al(3);
Eigen::VectorXd al2(3);
Eigen::VectorXd ei(12);
Eigen::VectorXd Mi;
Eigen::MatrixXd M(12,12), M2(12,12);
Eigen::VectorXd b0;
Eigen::VectorXd b;
Eigen::Matrix3d Mrot;
Eigen::Vector3d MomentRPY;
Eigen::Quaterniond quat;
double r,p,y;
unsigned int i=0,j=0;

void myInequalityPathConstraint( double *x, double *f, void *user_data ){
    for (i=0 ; i<3 ; i++) {
        q[i] = x[i];
    }
    for (i=0 ; i<6 ; i++) {
        q[i+7] = x[i+6];
    }

    //RPY to quaternion
    Mrot = Eigen::AngleAxisd(x[3], Eigen::Vector3d::UnitX())
           *Eigen::AngleAxisd(x[4], Eigen::Vector3d::UnitY())
           *Eigen::AngleAxisd(x[5], Eigen::Vector3d::UnitZ());
    quat = Eigen::Quaterniond(Mrot);

    q[6] = quat.w(); q[3] = quat.x(); q[4] = quat.y(); q[5] = quat.z();

    //set speed vector
    a3[0] = x[12]; a3[1] = x[13]; a3[2] = x[14];
    x3 = (Mrot*a3);
    v[0] = x3[0]; v[1] = x3[1]; v[2] = x3[2];

    a3[0] = x[15]; a3[1] = x[16]; a3[2] = x[17];
    x3 = (Mrot*a3);
    v[3] = x3[0]; v[4] = x3[1]; v[5] = x3[2];

    v[6] = x[18];
    v[7] = x[19];
    v[8] = x[20];
    v[9] = x[21];
    v[10] = x[22];
    v[11] = x[23];

    kinematics(robot, data, q, v);

    //    f[0] = data.oMi[7].translation()[0];
    //    f[1] = data.oMi[7].translation()[1];
    //    f[2] = data.oMi[7].translation()[2];

    gr[0] = 0.; gr[1] = 0.15; gr[2] = 0.; gr[3] = 1.;
    gr = data.oMi[7].toHomogeneousMatrix()*gr;


    Mrot = data.oMi[7].rotation();
    quat = Eigen::Quaterniond(Mrot);
//    std::cout << aa << std::endl;
    f[0] = gr[0];
    f[1] = gr[1];
    f[2] = gr[2];
    f[3] = quat.w();
    f[4] = quat.x();
    f[5] = quat.y();
    f[6] = quat.z();


//    quat = Eigen::Quaterniond(data.oMi[6].rotation());
//    quat = quat*Eigen::Quaterniond(Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX()));
//    r = 2*acos(quat.w());
//    std::cout << quat << std::endl;
//    f[6] = sqrt( (quat.x())*(quat.x()) + (quat.z())*(quat.z()));
    aa = data.v[7].linear();
    aa = Mrot.transpose()*aa;
        f[7] = aa[0];
        f[8] = aa[1];
        f[9] = aa[2];

//        std::cout << data.oMi[6].translation();

}
void myInequalityPathConstraint2( double *x, double *f, void *user_data ){
    for (i=0 ; i<3 ; i++) {
        q[i] = x[i];
    }
    for (i=0 ; i<6 ; i++) {
        q[i+7] = x[i+6];
    }

    //RPY to quaternion
    Mrot = Eigen::AngleAxisd(x[3], Eigen::Vector3d::UnitX())
           *Eigen::AngleAxisd(x[4], Eigen::Vector3d::UnitY())
           *Eigen::AngleAxisd(x[5], Eigen::Vector3d::UnitZ());
    quat = Eigen::Quaterniond(Mrot);

    q[6] = quat.w(); q[3] = quat.x(); q[4] = quat.y(); q[5] = quat.z();

    //set speed vector
    a3[0] = x[12]; a3[1] = x[13]; a3[2] = x[14];
    x3 = (Mrot*a3);
    v[0] = x3[0]; v[1] = x3[1]; v[2] = x3[2];

    a3[0] = x[15]; a3[1] = x[16]; a3[2] = x[17];
    x3 = (Mrot*a3);
    v[3] = x3[0]; v[4] = x3[1]; v[5] = x3[2];

    v[6] = x[18];
    v[7] = x[19];
    v[8] = x[20];
    v[9] = x[21];
    v[10] = x[22];
    v[11] = x[23];

    kinematics(robot, data, q, v);

    //    f[0] = data.oMi[7].translation()[0];
    //    f[1] = data.oMi[7].translation()[1];
    //    f[2] = data.oMi[7].translation()[2];

    gr[0] = 0.; gr[1] = 0.15; gr[2] = 0.; gr[3] = 1.;
    gr = data.oMi[7].toHomogeneousMatrix()*gr;



    f[0] = 4.*(x[0]-3.5)*(x[0]-3.5)+(x[2]+1.)*(x[2]+1.)-4;
//    f[1] = 4*(x[0]-3.5)*(x[0]-3.5)+(x[2]-2)*(x[2]-2)-4;
//    quat = Eigen::Quaterniond(data.oMi[6].rotation());
//    quat = quat*Eigen::Quaterniond(Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX()));
//    r = 2*acos(quat.w());
//    std::cout << quat << std::endl;
//    f[6] = sqrt( (quat.x())*(quat.x()) + (quat.z())*(quat.z()));
//        f[3] = data.v[7].linear()[0];
//        f[4] = data.v[7].linear()[1];
//        f[5] = data.v[7].linear()[2];

//        std::cout << data.oMi[6].translation();
}

void myDifferentialEquation( double *x, double *f, void *user_data )
{

    // --------------------   INITIALIZE INPUT VECTORS -----------------------  //
    //RPY to quaternion
    Mrot = Eigen::AngleAxisd(x[3], Eigen::Vector3d::UnitX())
           *Eigen::AngleAxisd(x[4], Eigen::Vector3d::UnitY())
           *Eigen::AngleAxisd(x[5], Eigen::Vector3d::UnitZ());
    quat = Eigen::Quaterniond(Mrot);

    //set Configuration vector
    q[0] = x[0];
    q[1] = x[1];
    q[2] = x[2];
    q[3] = quat.x();
    q[4] = quat.y();
    q[5] = quat.z();
    q[6] = quat.w();
    q[7] = x[6];
    q[8] = x[7];
    q[9] = x[8];
    q[10] = x[9];
    q[11] = x[10];
    q[12] = x[11];
    //std::cout << "configuration q :   \n" << q << std::endl;

    //set speed vector
    a3[0] = x[12]; a3[1] = x[13]; a3[2] = x[14];
    x3 = (Mrot*a3);
    v[0] = x3[0]; v[1] = x3[1]; v[2] = x3[2];

    a3[0] = x[15]; a3[1] = x[16]; a3[2] = x[17];
    x3 = (Mrot*a3);
    v[3] = x3[0]; v[4] = x3[1]; v[5] = x3[2];

    v[6] = x[18];
    v[7] = x[19];
    v[8] = x[20];
    v[9] = x[21];
    v[10] = x[22];
    v[11] = x[23];

    //set command vector
    u[0] = 0;
    u[1] = 0;
    u[2] = (x[24]+x[25]+x[26]+x[27]);
    u[3] = d*(x[24]-x[25]);;
    u[4] = d*(x[27]-x[26]);
    u[5] = c*(x[24]+x[25]-x[26]-x[27]);
    u[6] = x[28];
    u[7] = x[29];
    u[8] = x[30];
    u[9] = x[31];
    u[10] = x[32];
    u[11] = x[33];

    // --------------------   COMPUTE DYNAMIC -----------------------  //

//    std::cout << " q :    " << q.transpose() << std::endl << std::endl;

//    std::cout << " M rnea :\n" <<  M << std::endl << std::endl;

    M2 = crba(robot, data, q);
    for(i=0 ; i<12 ; i++) {
        for (int j=i ; j<12 ; j++) {
            M(i,j) = M2(i,j);
            M(j,i) = M2(i,j);
        }
    }

    //add mot inertie
    for (i=6 ; i<12 ; i++) {
        M(i,i) = M(i,i) + 0.070756;;
        u[i] = u[i];
    }

//  std::cout << " M crba :\n" <<  M << std::endl << std::endl;

//    std::cout << "v :  " << v.transpose() << std::endl << std::endl;

    a.setZero();

    b = rnea(robot, data, q, v, a);

//    std::cout << " b :    " << b.transpose() << std::endl << std::endl;

    a = M.llt().solve(u-b);
//    std::cout << "F tot moteurs : \n" << x[24]+x[25]+x[26]+x[27] << std::endl << std::endl;
//    std::cout << "Forces : \n" << u << std::endl << std::endl;
//    std::cout << "Acceleration : \n" << a << std::endl << std::endl;

    // --------------------   COMPUTE OUTPUT -----------------------  //
    // --------------------   COMPUTE OUTPUT -----------------------  //

    // vitesses lineaires Freeflyer repere monde
    f[0] = x[12];
    f[1] = x[13];
    f[2] = x[14];

    // vitesses angulaires Freeflyer (transformation RPY)
    f[3] = cos(x[5])/cos(x[4])*x[15] + sin(x[5])/cos(x[4])*x[16];
    f[4] = -sin(x[5])*x[15] + cos(x[5])*x[16];
    f[5] = cos(x[5])*tan(x[4])*x[15] + sin(x[5])*tan(x[4])*x[16] + x[17];

    // vitesses articulations
    for (i=0 ; i<6 ; i++) {
        f[6+i] = v[i+6];
    }

//    std::cout << a << std::endl << std::endl;

    for (i=0 ; i<3 ; i++) {
        al[i] = a[i];
        aa[i] = a[i+3];
    }
    al = Mrot.transpose()*al;
    for (i=0 ; i<3 ;i++) {
        f[i+12] = al[i];
    }

    //accelerations angulaires
    aa = Mrot.transpose()*aa;
    for (i=0 ; i<3 ;i++) {
        f[i+15] = aa[i];
    }

    //acceleration articulation
    for (i=0 ; i<6 ;i++) {
        f[i+18] = a[i+6];
    }

    // cost
    double gripper[10];
    myInequalityPathConstraint(x,gripper,NULL);

    //  Trainer un objet sur le sol
//    f[24]= x[36]*((gripper[0]-3)*(gripper[0]-3)+(gripper[1]+0.)*(gripper[1]+0.)+(gripper[2]+1)*(gripper[2]+1) /*+ 100*(gripper[4]-1.)*(gripper[4]-1.)*/) + x[37]*((gripper[0]-3-(x[35]-5))*(gripper[0]-3-(x[35]-5))+(gripper[1]+0.)*(gripper[1]+0.)+*/(gripper[2]+1)*(gripper[2]+1) /*+ (gripper[4]-1.)*(gripper[4]-1.)*/) + 0.x[38]*((gripper[0]-6)*(gripper[0]-6)+(gripper[1]+0.)*(gripper[1]+0.)+(gripper[2]+1)*(gripper[2]+1) /*+ 100*(gripper[4]-1.)*(gripper[4]-1.)*/) + (u[6]*u[6] + u[7]*u[7] + u[8]*u[8] + u[9]*u[9] + u[10]*u[10] + u[11]*u[11])*0.01;

    // Tourner une manivelle
    f[24]= x[36]*((gripper[0]-3)*(gripper[0]-3)+(gripper[1]+0.)*(gripper[1]+0.)+(gripper[2]+1)*(gripper[2]+1) + 0.*(gripper[4]-1.)*(gripper[4]-1.) )+ x[37]*( (gripper[0]-3-0.2*sin(0.5*(x[35]-5)))*(gripper[0]-3-0.2*sin(0.5*(x[35]-5)))+(gripper[1]+0.)*(gripper[1]+0.)+(gripper[2]+0.8+cos(0.5*(x[35]-5)))*(gripper[2]+0.8+cos(0.5*(x[35]-5))) + 1.*( (gripper[7]-0.2*cos(0.5*(x[35]-5)))*(gripper[7]-0.2*cos(0.5*(x[35]-5)))+(gripper[8])*(gripper[8])+(gripper[9]+0.2*sin(0.5*(x[35]-5)))*(gripper[9])+0.2*sin(0.5*(x[35]-5)) ) + 0.*(gripper[4]-1.)*(gripper[4]-1.) ) + (u[6]*u[6] + u[7]*u[7] + u[8]*u[8] + u[9]*u[9] + u[10]*u[10] + u[11]*u[11])*0.001;
    f[25]=1.;

}






// -------------------------------------------------------------------------
//              USE THE ACADO TOOLKIT TO SOLVE THE PROBLEM:
// -------------------------------------------------------------------------


USING_NAMESPACE_ACADO


int main( ){

    // INTRODUCE THE VARIABLES:
    // --------------------------------------------------
    DifferentialState     xp, yp, zp, phi, theta, psi, spx ,slx, ex, w1x, w2x, w3x, vx, vy, vz, p, q, r, spv, slv, ev, w1v, w2v, w3v, cost, tps;
    Control               u1, u2, u3, u4, spu, slu, eu, w1u, w2u, w3u;
    Parameter             T;
    Disturbance           D1,D2,D3;

    // DEFINE THE DIMENSIONS OF THE C-FUNCTIONS:
    // --------------------------------------------------
    CFunction F( 26, myDifferentialEquation     );
    //CFunction M( NJ, myObjectiveFunction        );
    //CFunction I( NI, myInitialValueConstraint   );
    //CFunction E( NE, myEndPointConstraint       );
    CFunction H( 1, myInequalityPathConstraint2 );


    DifferentialEquation  f;//f2(1.0,2.0)     ;

//    const double h = 0.01;
//    DiscretizedDifferentialEquation  f(h);
//    F.setUserData((void*) &h);


    // DEFINE THE OPTIMIZATION VARIABLES:
    // --------------------------------------------------

    IntermediateState x(39);

     x(0) = xp;  x(1) = yp;  x(2) = zp;  x(3) = phi;  x(4) = theta;  x(5) = psi;
     x(6) = spx;  x(7) = slx;  x(8) = ex;  x(9) = w1x;  x(10) = w2x;  x(11) = w3x;
     x(12) = vx;  x(13) = vy;  x(14) = vz;  x(15) = p; x(16) = q; x(17) = r;
     x(18) = spv;  x(19) = slv;  x(20) = ev;  x(21) = w1v; x(22) = w2v; x(23) = w3v;
     x(24) = u1; x(25) = u2; x(26) = u3; x(27) = u4;
     x(28) = spu; x(29) = slu; x(30) = eu; x(31) = w1u; x(32) = w2u; x(33) = w3u;
     x(34) = cost; x(35) = tps;
     x(36) = D1; x(37) = D2; x(38) = D3;


     // THE EQUATIONS FOR THE JUMP:
    // ---------------------------
//    Transition j;
//    j << xp == xp ;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    Grid timeGrid(0.0,TPS,NB_IT+1);
    OCP ocp( 0.0, TPS, NB_IT );

    VariablesGrid dist1(1,timeGrid), dist2(1,timeGrid), dist3(1,timeGrid);
    dist1.setZero(); dist2.setZero(); dist3.setZero();
    for(i=0 ; i<NB_IT+1 ; i++) {
        if (i>5 && i<10)
            dist1(i,0) = 1.;
        if (i>9 && i<16)
            dist2(i,0) = 1.;
        if (i>15 && i<18)
            dist3(i,0) = 1.;
    }

//    ocp.minimizeMayerTerm(x(34) + 1.*((x(0)-7)*(x(0)-7) + x(1)*x(1) + x(2)*x(2)) );
    ocp.minimizeMayerTerm(x(34));
    ocp.subjectTo( f << F(x) );
    ocp.subjectTo(x(36) == dist1);
    ocp.subjectTo(x(37) == dist2);
    ocp.subjectTo(x(38) == dist3);

    VariablesGrid obs(1,timeGrid);
    obs.setZero();
//    ocp.subjectTo( obs <= H(x) );
//    ocp.subjectTo(4 <= 4*(x(0)-3.5)*(x(0)-3.5)+(x(2)+1)*(x(2)+1) );
//    ocp.subjectTo( 0.0 <= T <= 1);

    // Start constraints
    DVector x0(26);

    for (unsigned int k=0 ; k < 24 ; k++) {
        if (k==7 || k==9) {
            ocp.subjectTo(AT_START, x(k) == -1.57);
            x0(k)=-1.57;
        }
        else {
            ocp.subjectTo(AT_START, x(k) == 0.0);
            x0(k)=0.;
        }
    }
    ocp.subjectTo(AT_START, x(34)==0.0);
    x0(24)=0.;
    ocp.subjectTo(AT_START, x(35)==0.0);
    x0(25)=0.;

    // End constraints
//    ocp.subjectTo( AT_END, x(0) ==  0.0 );
//    ocp.subjectTo( AT_END, x(1) ==  0.0 );
//    ocp.subjectTo(AT_END, x(2) == 0.0 );

    ocp.subjectTo(AT_END, x(3) == 0.0);
    ocp.subjectTo(AT_END, x(4) == 0.0);
    ocp.subjectTo(AT_END, x(5) == 0.0);

    for (unsigned int k=6 ; k < 18 ; k++) {
        if (k==7 || k==9) {
            ocp.subjectTo(AT_END, x(k) == -1.57);
        }
        else
            ocp.subjectTo(AT_END, x(k) == 0.0);
    }

//    VariablesGrid upZ(3,timeGrid), lowZ(3,timeGrid);
//    for (int i = 0 ; i<NB_IT+1 ; i++ ) {
//        if (i == 10 || i==9 || i==11 || i==8 || i==12 || i==13 || i==7) {
//            upZ(i,0) = 1.53;
//            lowZ(i,0) = 1.47;
//            upZ(i,1) = 0.03;
//            lowZ(i,1) = -0.03;
//            upZ(i,2) = -0.97;
//            lowZ(i,2) = -1.03;
//        }
//        else {
//            for (unsigned int k=0; k<3 ; k++) {
//            upZ(i,k) = 10;
//            lowZ(i,k) = -10;
//            }
//        }
//    }
//    ocp.subjectTo(lowZ <= H(x) <= upZ);

    ocp.subjectTo(0 <= x(24) <= 200);
    ocp.subjectTo(0 <= x(25) <= 200);
    ocp.subjectTo(0 <= x(26) <= 200);
    ocp.subjectTo(0 <= x(27) <= 200);
    ocp.subjectTo(-150 <= x(28) <= 150);
    ocp.subjectTo(-150 <= x(29) <= 150);
    ocp.subjectTo(-150 <= x(30) <= 150);
    ocp.subjectTo(-28 <= x(31) <= 28);
    ocp.subjectTo(-28 <= x(32) <= 28);
    ocp.subjectTo(-28 <= x(33) <= 28);

    ocp.subjectTo(-6.28 <= x(6) <= 6.28);
    ocp.subjectTo(-3.14 <= x(7) <= 0);
    for (unsigned int k=2 ; k < 6 ; k++) {
        ocp.subjectTo(-6.28 <= x(k+6) <= 6.28);
    }


    // VISUALIZE THE RESULTS IN A GNUPLOT WINDOW:
    // ------------------------------------------
    GnuplotWindow window1(PLOT_AT_END);
    window1.addSubplot( x(0),"x" );
    window1.addSubplot( x(1) , "y");
    window1.addSubplot( x(2), "z");
    window1.addSubplot( x(3),"phi" );
    window1.addSubplot( x(4) , "theta");
    window1.addSubplot( x(5), "psi");
//    window1.addSubplot( x(12),"xv" );
    window1.addSubplot( x(24),"u1" );
    window1.addSubplot( x(25),"u2" );
    window1.addSubplot( x(26),"u3" );
    window1.addSubplot( x(27),"u4" );
//    window1.addSubplot( x(28),"uArm1" );
//    window1.addSubplot( x(29),"uArm2" );
//    window1.addSubplot( x(29),"uArm2" );
//        window1.addSubplot( x(7),"xArm2" );
//    window1.addSubplot( x(30),"uArm3" );
//    window1.addSubplot( x(31),"uArm4" );
//    window1.addSubplot( x(32),"uArm5" );
//    window1.addSubplot( x(33),"uArm6" );
//    window1.addSubplot( x(17),"w2u" );
//    window1.addSubplot( x(10),"w2v" );
//    window1.addSubplot( x(4),"w2x" );
//    window1.addSubplot( x(17),"w3u" );
//    window1.addSubplot( x(11),"w3v" );
//    window1.addSubplot( x(5),"w3x" );

   // U INITIAL
//    Grid timeGrid(0.0,2.0,NB_IT+1);
    VariablesGrid u_init(10, timeGrid);
    for (int i = 0 ; i<NB_IT+1 ; i++ ) {
      u_init(i,0) = 120;
      u_init(i,1) = 120;
      u_init(i,2) = 120;
      u_init(i,3) = 120;
      u_init(i,4) = 0;
      u_init(i,5) = 0;
      u_init(i,6) = 0;
      u_init(i,7) = 0;
      u_init(i,8) = 0;
      u_init(i,9) = 0;
    }
    VariablesGrid x_init(26, timeGrid);
    x_init.setZero();
    for (int i = 0 ; i<NB_IT+1 ; i++ ) {
        for (unsigned int k=0 ; k < 26 ; k++) {
            if (k==0) {
                x_init(i,k) = 0.;//i* 7./NB_IT;
            }
            else if (k==2)
                x_init(i,k) = 0.;//2.;
            else if (k==24 || k==25) { x_init(i,k) = i*TPS/NB_IT;}
            else if (k==7 || k==9) {x_init(i,k) = -1.57;}
            else
                x_init(i,k) = 0;
        }
    }
    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
    OptimizationAlgorithm algorithm(ocp);
//VariablesGrid u_1, x_1;
//u_1.read("/tmp/controls1.txt");
//u_1.print();
//x_1.read("/tmp/states1.txt");
//x_1.print();
    algorithm.initializeControls(u_init);
    algorithm.initializeDifferentialStates(x_init);

   algorithm.set( INTEGRATOR_TOLERANCE, 1e-3 );
//    algorithm.set( DISCRETIZATION_TYPE, SINGLE_SHOOTING   );
//    algorithm.set( HESSIAN_APPROXIMATION, CONSTANT_HESSIAN );

//            algorithm.set( HESSIAN_APPROXIMATION, BLOCK_BFGS_UPDATE );
    algorithm.set( INTEGRATOR_TYPE, INT_RK45 );
//    algorithm.set( MAX_NUM_INTEGRATOR_STEPS, 100000);
    algorithm.set( MAX_NUM_ITERATIONS, 250 );
    algorithm.set( KKT_TOLERANCE, 1e-12);
//    algorithm.set( SPARSE_QP_SOLUTION, CONDENSING);
    //algorithm.set( ABSOLUTE_TOLERANCE,<double> );
    //algorithm.set( INTEGRATOR_TOLERANCE,<double> );
//    algorithm << window1;
    LogRecord logdiff(LOG_AT_END , PS_DEFAULT);
        logdiff << LOG_DIFFERENTIAL_STATES;
        algorithm << logdiff;
    algorithm.solve();

    algorithm.getLogRecord(logdiff);
    std::ofstream file;
    file.open("/tmp/log_pickup.txt",std::ios::out);
    std::ofstream file2;
    file2.open("/tmp/log_pickup_ocp.txt",std::ios::out);
    logdiff.print(file2);

//    int nb_process= 6 ;
//    Grid timeGrid3(0.0,TPS,nb_process*NB_IT+1);
//    VariablesGrid controli;
//    VariablesGrid state;
//    algorithm.getControls(controli);
//        VariablesGrid control(10,timeGrid3);
//        for (int i=0 ; i<nb_process*NB_IT+1 ; i++) {
//            for (int j=0 ; j<10 ; j++) {
//                control(i,j) = controli(int(double(i)/double(nb_process)),j);
//            }
//        }
//    OutputFcn identity;
//    DynamicSystem dynamicSystem( f,identity );
//    Process process( dynamicSystem,INT_RK78 );
//    process.set(INTEGRATOR_TOLERANCE, 10e-12);
//    process.initializeStartValues(x0);
//    process.init(0.);
//    process.run(control);
//    process.getLast(LOG_DIFFERENTIAL_STATES,state);
//    GnuplotWindow window3(PLOT_AT_END);
//    window3.addSubplot( state(0),"DifferentialState x" );
//    window3.addSubplot( state(1),"DifferentialState y" );
//    window3.addSubplot( state(2),"DifferentialState z" );
//    window3.plot();
//    state.print(file);

//    algorithm.getDifferentialStates("/tmp/states1.txt");
//    algorithm.getControls("/tmp/controls1.txt");
    return 0;
}
/* <<< end tutorial code <<< */



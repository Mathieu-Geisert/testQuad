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


#include <acado_optimal_control.hpp>
//#include <acado_gnuplot.hpp>
//#include <acado_toolkit.hpp>
#include <include/acado_gnuplot/gnuplot_window.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <pinocchio/multibody/model.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>

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


// -----------------------------------------------------------------------------
//   UGLY C-STYLE DEFINITION OF THE OBJECTIVE, MODEL AND CONSTRAINT FUNCTIONS:
// -----------------------------------------------------------------------------

const double Cf = 0.00065;
const double d = 0.250;
const double c = 0.00001;

se3::Model robot = se3::urdf::buildModel("/local/mgeisert/models/ur5_light.urdf", true);
se3::Data data(robot);

Eigen::VectorXd q(13);
Eigen::VectorXd v(12);
Eigen::VectorXd u(12);
Eigen::VectorXd a(12);
Eigen::VectorXd al(3);
Eigen::VectorXd al2(3);
Eigen::VectorXd ei(12);
Eigen::VectorXd Mi;
Eigen::MatrixXd M(12,12);
Eigen::VectorXd b0;
Eigen::VectorXd b;
Eigen::Matrix3d Mrot;
Eigen::Vector3d MomentRPY;
//Eigen::Quaternion quat;
double r,p,y;
unsigned int i,j;

void myDifferentialEquation( double *x, double *f, void *user_data )
{
    //Set vectors
    for (i = 0 ; i<12 ; i++) {
        v[i] = x[i+12];
        a[i] = 0.;
        ei[i] = 0.;
    }

    // Set configuration vector
    for (i=0 ; i<3 ; i++) {
        q[i] = x[i];
    }
    for (i=0 ; i<6 ; i++) {
        q[i+7] = x[i+6];
    }

    //RPY to quaternion
    r = x[3]; p = x[4]; y = x[5];
    Mrot(0,0) = cos(y)*cos(p);
    Mrot(0,1) = cos(y)*sin(p)*sin(r) - sin(y)*cos(r);
    Mrot(0,2) = cos(y)*sin(p)*cos(r) + sin(y)*sin(r);
    Mrot(1,0) = sin(y)*cos(p);
    Mrot(1,1) = sin(y)*sin(p)*sin(r) + cos(y)*cos(r);
    Mrot(1,2) = sin(y)*sin(p)*cos(r) - cos(y)*sin(r);
    Mrot(2,0) = -sin(p);
    Mrot(2,1) = cos(p)*sin(r);
    Mrot(2,2) = cos(p)*cos(r);

    Eigen::Quaternion<double> quat(Mrot);
    q[6] = quat.w(); q[3] = quat.x(); q[4] = quat.y(); q[5] = quat.z();

    //std::cout << "configuration q :   \n" << q << std::endl;

    MomentRPY(0) = x[25];//*0.001;
    MomentRPY(1) = x[26];//*0.001;
    MomentRPY(2) = x[27];//*0.001;

    //MomentRPY = Mrot*MomentRPY;

    u[0] = 0;//- (x[24]+x[25]+x[26]+x[27])*sin(p);
    u[1] = 0;//(x[24]+x[25]+x[26]+x[27])*sin(r)*cos(p);
    u[2] = (x[24]);//*cos(r)*cos(p);
    u[3] = MomentRPY(0);
    u[4] = MomentRPY(1);
    u[5] = MomentRPY(2);

    for (i = 0 ; i<6 ; i++) {
        u[i+6] = x[i+28];//*0.01;
    }


    b0 = rnea(robot, data, q, a, a);

    for (i=0 ; i<12 ; i++) {
        ei(i)=1;
        Mi = rnea(robot, data, q, a, ei);
        for (j=0 ; j<12 ; j++ ) {
            M(i,j) = Mi(j)-b0(j);
        }
        ei(i)=0;
    }

    //std::cout << " M rnea :\n" <<  M << std::endl << std::endl;
    //std::cout << "crba :\n" << crba(robot, data, q) << std::endl << std::endl;
    // self adjoint upper ???

    b = rnea(robot, data, q, v, a);

    a = M.llt().solve(u-b);
//    std::cout << "F tot moteurs : \n" << x[24]+x[25]+x[26]+x[27] << std::endl << std::endl;
//    std::cout << "Forces : \n" << u << std::endl << std::endl;
//    std::cout << "Acceleration : \n" << a << std::endl << std::endl;

    // --------------------   COMPUTE OUTPUT -----------------------  //
    // vitesses lineaires Freeflyer
    for (i=0 ; i<3 ; i++ ) {
        f[i] = v[i];
    }

    // vitesses angulaires Freeflyer (transformation RPY)
    f[3] = v[3] + sin(p)*v[5];
    f[4] = cos(r)*v[4] - sin(r)*cos(p)*v[5];
    f[5] = sin(r)*v[4] + cos(r)*cos(p)*v[5];

    // vitesses articulations
    for (i=0 ; i<6 ; i++) {
        f[6+i] = v[i+6];
    }
    for (i=0 ; i<3 ; i++) {
        al[i] = a[i];
    }
    al2 = Mrot*al;
    for (i=0 ; i<3 ;i++) {
        f[i+12] = al2[i];
    }

    for (i=0 ; i<9 ;i++) {
        f[i+15] = a[i+3];
    }

    // cost
    //f[12]=(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);//*0.01;

}

// -------------------------------------------------------------------------
//              USE THE ACADO TOOLKIT TO SOLVE THE PROBLEM:
// -------------------------------------------------------------------------


USING_NAMESPACE_ACADO


int main( ){

    // INTRODUCE THE VARIABLES:
    // --------------------------------------------------
    DifferentialState     xp, yp, zp, phi, theta, psi, spx ,slx, ex, w1x, w2x, w3x, vx, vy, vz, p, q, r, spv, slv, ev, w1v, w2v, w3v;//  cost;
    Control               u1, u2, u3, u4, spu, slu, eu, w1u, w2u, w3u;
    Parameter             T;

    // DEFINE THE DIMENSIONS OF THE C-FUNCTIONS:
    // --------------------------------------------------
    CFunction F( 24, myDifferentialEquation     );
    //CFunction M( NJ, myObjectiveFunction        );
    //CFunction I( NI, myInitialValueConstraint   );
    //CFunction E( NE, myEndPointConstraint       );
    //CFunction H( NH, myInequalityPathConstraint );


    DifferentialEquation  f      ;

//    const double h = 0.01;
//    DiscretizedDifferentialEquation  f(h);
//    F.setUserData((void*) &h);

    // DEFINE THE OPTIMIZATION VARIABLES:
    // --------------------------------------------------

    IntermediateState x(34);

     x(0) = xp;  x(1) = yp;  x(2) = zp;  x(3) = phi;  x(4) = theta;  x(5) = psi;
     x(6) = spx;  x(7) = slx;  x(8) = ex;  x(9) = w1x;  x(10) = w2x;  x(11) = w3x;
     x(12) = vx;  x(13) = vy;  x(14) = vz;  x(15) = p; x(16) = q; x(17) = r;
     x(18) = spv;  x(19) = slv;  x(20) = ev;  x(21) = w1v; x(22) = w2v; x(23) = w3v;
     x(24) = u1; x(25) = u2; x(26) = u3; x(27) = u4;
     x(28) = spu; x(29) = slu; x(30) = eu; x(31) = w1u; x(32) = w2u; x(33) = w3u;


    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, T, 20 );

    ocp.minimizeMayerTerm(T);

    ocp.subjectTo( f << F(x)*T );

    ocp.subjectTo( 1 <= T <= 20);

    // Start constraints
    for (unsigned int k=0 ; k < 24 ; k++) {
        ocp.subjectTo(AT_START, x(k) == 0.0);
    }

    // End constraints
    ocp.subjectTo( AT_END, x(0) ==  10.0 );
    ocp.subjectTo( AT_END, x(1) ==  0.0 );
    ocp.subjectTo( AT_END, x(2) ==  0.0 );
//    ocp.subjectTo(AT_END, x(6) == 0.0);
//    ocp.subjectTo(AT_END, x(7) == 0.0);
//    ocp.subjectTo(AT_END, x(8) == 0.0);

    for (unsigned int k=3 ; k < 24 ; k++) {
        ocp.subjectTo(AT_END, x(k) == 0.0 );
    }

    ocp.subjectTo(0 <= x(24) <= 120);
    ocp.subjectTo(-5 <= x(25) <= 5);
    ocp.subjectTo(-5 <= x(26) <= 5);
    ocp.subjectTo(-5 <= x(27) <= 5);

    ocp.subjectTo(x(5) == 0.);
    ocp.subjectTo( -1.57 <= x(3) <= 1.57);
    ocp.subjectTo( -1.57 <= x(4) <= 1.57);

    // VISUALIZE THE RESULTS IN A GNUPLOT WINDOW:
    // ------------------------------------------
    GnuplotWindow window1(PLOT_AT_EACH_ITERATION);
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

    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
    OptimizationAlgorithm algorithm(ocp);
    //algorithm.set( INTEGRATOR_TOLERANCE, 1e+2 );
    //algorithm.set( DISCRETIZATION_TYPE, SINGLE_SHOOTING   );
    //algorithm.set( HESSIAN_APPROXIMATION, FULL_BFGS_UPDATE );
    algorithm.set( INTEGRATOR_TYPE, INT_RK78 );
    algorithm.set( MAX_NUM_ITERATIONS, 500 );
    algorithm.set( KKT_TOLERANCE, 1e-12);
    //algorithm.set( ABSOLUTE_TOLERANCE, 10. );
    //algorithm.set( INTEGRATOR_TOLERANCE,<double> );
    algorithm << window1;

//    double testx[34];
//    for (i=0 ; i<34 ; i++) {
//        testx[i] = 0.;
//    }
//    testx[4] = 0;
//    testx[3] = 1.57;
//    testx[13] = 0.;
//    testx[14] = 0.;
//    testx[15] = 0.;
//    testx[16] = 1.;

//    double testf[24];
//    for (i=0 ; i<24 ; i++) {
//        testf[i] = 0;
//    }
//    myDifferentialEquation( testx, testf, testf );

//    for ( unsigned int k=0 ; k<24 ; k++ ) {
//        std::cout << testf[k] << std::endl;
//    }

    algorithm.solve();

    return 0;
}
/* <<< end tutorial code <<< */


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
#include <acado_gnuplot.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

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

se3::Model robot = se3::urdf::buildModel("/local/mgeisert/models/universal_robot/ur_description/urdf/ur5_robot.urdf");
se3::Data data(robot);

void myDifferentialEquation( double *x, double *f, void *user_data )
{
    Eigen::VectorXd q(6);
    q[0] = x[0]; q[1] = x[1]; q[2] = x[2]; q[3] = x[3]; q[4] = x[4]; q[5] = x[5];
    Eigen::VectorXd v(6);
    v[0] = x[6]; v[1] = x[7]; v[2] = x[8]; v[3] = x[9]; v[4] = x[10]; v[5] = x[11];
    Eigen::VectorXd u(6);
    u[0] = x[13]; u[1] = x[14]; u[2] = x[15]; u[3] = x[16]; u[4] = x[17]; u[5] = x[18];
    Eigen::VectorXd a(6);
    a[0] = 0; a[1] = 0; a[2] = 0; a[3] = 0; a[4] = 0; a[5] = 0;

    Eigen::VectorXd ei(6);
    ei[0] = 0.; ei[1] = 0.; ei[2] = 0.; ei[3] = 0.; ei[4] = 0.; ei[5] = 0.;

    Eigen::VectorXd Mi;
    Eigen::MatrixXd M(6,6);
    Eigen::VectorXd b0 = rnea(robot, data, q, a, a);

    for (unsigned int i=0 ; i<6 ; i++) {
        ei(i)=1;
        Mi = rnea(robot, data, q, a, ei);
        for (unsigned int j=0 ; j<6 ; j++ ) {
            M(i,j) = Mi(j)-b0(j);
        }
        ei(i)=0;
    }

    //std::cout << " M rnea :\n" <<  M << std::endl << std::endl;
    //std::cout << "crba :\n" << crba(robot, data, q) << std::endl << std::endl;

    Eigen::VectorXd b = rnea(robot, data, q, v, a);

//    M(1) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;
//    ei[0] = 0; ei[1] = 1;
//    M(1) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;
//    ei[1] = 0; ei[2] = 1;
//    M(2) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;
//    ei[2] = 0; ei[3] = 1;
//    M(3) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;
//    ei[3] = 0; ei[4] = 1;
//    M(4) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;
//    ei[4] = 0; ei[5] = 1;
//    M(5) = rnea(robot, data, q, 0, ei) - rnea(robot, data, q, 0, 0);
//    std::cout << M << std::endl;

    a = M.llt().solve(u-b);
    //std::cout << a << std::endl;

//    for(unsigned int i=0 ; i<6 ; i++) {
//        f[i] = x[i]+x[i+6]*h;
//        f[i+6] = x[i+6] + h*x[12+i];
//    }

    /// DISCRET MODE
//    double h = ((double*) user_data)[0];
//    f[0]=q[0]+v[0]*h;
//    f[1]=q[1]+v[1]*h;
//    f[2]=q[2]+v[2]*h;
//    f[3]=q[3]+v[3]*h;
//    f[4]=q[4]+v[4]*h;
//    f[5]=q[5]+v[5]*h;
//    f[6]=v[0]+a[0]*h;
//    f[7]=v[1]+a[1]*h;
//    f[8]=v[2]+a[2]*h;
//    f[9]=v[3]+a[3]*h;
//    f[10]=v[4]+a[4]*h;
//    f[11]=v[5]+a[5]*h;
//    f[12]=x[12]+(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);//*0.01;

    f[0]=v[0];
    f[1]=v[1];
    f[2]=v[2];
    f[3]=v[3];
    f[4]=v[4];
    f[5]=v[5];
    f[6]=a[0];
    f[7]=a[1];
    f[8]=a[2];
    f[9]=a[3];
    f[10]=a[4];
    f[11]=a[5];
    f[12]=(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);//*0.01;

    std::cout << "ba" << std::endl;

//    for (unsigned int i=0 ; i<21 ; i++ )
//        std::cout << i << "   " << x[i] << std::endl;
}

void myForwardDerive(int number, double* x, double* seed, double* f, double* df, void* userData)
{
    Eigen::VectorXd q(6);
    q[0] = x[0]; q[1] = x[1]; q[2] = x[2]; q[3] = x[3]; q[4] = x[4]; q[5] = x[5];
    Eigen::VectorXd v(6);
    v[0] = x[6]; v[1] = x[7]; v[2] = x[8]; v[3] = x[9]; v[4] = x[10]; v[5] = x[11];
    Eigen::VectorXd u(6);
    u[0] = x[13]; u[1] = x[14]; u[2] = x[15]; u[3] = x[16]; u[4] = x[17]; u[5] = x[18];


}

void myObjectiveFunction( double *x, double *f, void *user_data ){

    f[0] = x[12];
}


void myInitialValueConstraint( double *x, double *f, void *user_data ){

    f[0] =  x[4];
    f[1] =  x[0];
    f[2] =  x[1];
    f[3] =  x[3];
}


void myEndPointConstraint( double *x, double *f, void *user_data ){

    f[0] =  x[4] - 10.0;
    f[1] =  x[0];
}


void myInequalityPathConstraint( double *x, double *f, void *user_data ){

    f[0] =  x[0];
}


// -------------------------------------------------------------------------
//              USE THE ACADO TOOLKIT TO SOLVE THE PROBLEM:
// -------------------------------------------------------------------------


USING_NAMESPACE_ACADO


int main( ){

    // INTRODUCE THE VARIABLES:
    // --------------------------------------------------
    DifferentialState     spx ,slx, ex, w1x, w2x, w3x, spv, slv, ev, w1v, w2v, w3v, cost;
    Control               spu, slu, eu, w1u, w2u, w3u;

    // DEFINE THE DIMENSIONS OF THE C-FUNCTIONS:
    // --------------------------------------------------
    CFunction F( NX, myDifferentialEquation     );
    CFunction M( NJ, myObjectiveFunction        );
    //CFunction I( NI, myInitialValueConstraint   );
    //CFunction E( NE, myEndPointConstraint       );
    //CFunction H( NH, myInequalityPathConstraint );


    DifferentialEquation  f      ;

//    const double h = 0.01;
//    DiscretizedDifferentialEquation  f(h);
//    F.setUserData((void*) &h);

    // DEFINE THE OPTIMIZATION VARIABLES:
    // --------------------------------------------------

    IntermediateState x(19);

     x(0) = spx;  x(1) = slx;  x(2) = ex;  x(3) = w1x;  x(4) = w2x;  x(5) = w3x;
     x(6) = spv;  x(7) = slv;  x(8) = ev;  x(9) = w1v; x(10) = w2v; x(11) = w3v; x(12) = cost;
    x(13) = spu; x(14) = slu; x(15) = eu; x(16) = w1u; x(17) = w2u; x(18) = w3u;


    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, 2.0,10 );

    ocp.minimizeMayerTerm(M(x));

    ocp.subjectTo( f << F(x) );

//    double stats[19];
//    stats[0] = 2; stats[1] = 2; stats[2] = 2; stats[3] = 2; stats[4] = 2; stats[5] = 2; stats[6] = 2; stats[7] = 2; stats[8] = 2;
//    stats[9] = 2; stats[10] = 2; stats[11] = 2; stats[12] = 2; stats[13] = 2; stats[14] = 2; stats[15] = 2; stats[16] = 2;
//    stats[17] = 2; stats[18] = 2;

//    double seed[19];
//    seed[0] = 100; seed[1] = 0; seed[2] = 0; seed[3] = 0; seed[4] = 0; seed[5] = 0; seed[6] = 0; seed[7] = 0; seed[8] = 0;
//    seed[9] = 0; seed[10] = 0; seed[11] = 0; seed[12] = 0; seed[13] = 0; seed[14] = 0; seed[15] = 0; seed[16] = 0;
//    seed[17] = 0; seed[18] = 0;

//    double value[13];
//    double derive[13];

//    F.AD_forward(stats, seed, value, derive);
//    for ( unsigned int i=0 ; i<13 ; i++)
//        std::cout << derive[i] <<std::endl;

    ocp.subjectTo( AT_START, x(0) ==  0.0 );
    ocp.subjectTo( AT_START, x(1) ==  0.0 );
    ocp.subjectTo( AT_START, x(2) ==  1.57 );
    ocp.subjectTo( AT_START, x(3) ==  0.0 );
    ocp.subjectTo( AT_START, x(4) ==  0.0 );
    ocp.subjectTo( AT_START, x(5) ==  0.0 );
    ocp.subjectTo( AT_START, x(6) ==  0.0 );
    ocp.subjectTo( AT_START, x(7) ==  0.0 );
    ocp.subjectTo( AT_START, x(8) ==  0.0 );
    ocp.subjectTo( AT_START, x(9) ==  0.0 );
    ocp.subjectTo( AT_START, x(10) ==  0.0 );
    ocp.subjectTo( AT_START, x(11) ==  0.0 );
    ocp.subjectTo( AT_START, x(12) ==  0.0 );

//    ocp.subjectTo( -2 <= x(18) <= 2 );
//    ocp.subjectTo( -2 <= x(13) <= 2 );
//    ocp.subjectTo( -2 <= x(14) <= 2 );
//    ocp.subjectTo( -2 <= x(15) <= 2 );
//    ocp.subjectTo( -2 <= x(16) <= 2 );
//    ocp.subjectTo( -2 <= x(17) <= 2 );

  //  ocp.subjectTo( AT_START, I(x) ==  0.0 );
  //  ocp.subjectTo( AT_END  , E(x) ==  0.0 );
  //  ocp.subjectTo(  );


    // VISUALIZE THE RESULTS IN A GNUPLOT WINDOW:
    // ------------------------------------------
    GnuplotWindow window1(PLOT_AT_EACH_ITERATION);
    window1.addSubplot( x(13),"spu" );
    window1.addSubplot( x(6),"spv" );
    window1.addSubplot( x(0),"spx" );
    window1.addSubplot( x(14),"slu" );
    window1.addSubplot( x(7),"slv" );
    window1.addSubplot( x(1),"slx" );
    window1.addSubplot( x(15),"eu" );
    window1.addSubplot( x(8),"ev" );
    window1.addSubplot( x(2),"ex" );
    window1.addSubplot( x(16),"w1u" );
    window1.addSubplot( x(9),"w1v" );
    window1.addSubplot( x(3),"w1x" );
//    window1.addSubplot( x(17),"w2u" );
//    window1.addSubplot( x(10),"w2v" );
//    window1.addSubplot( x(4),"w2x" );
//    window1.addSubplot( x(17),"w3u" );
//    window1.addSubplot( x(11),"w3v" );
//    window1.addSubplot( x(5),"w3x" );

    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
    OptimizationAlgorithm algorithm(ocp);
    //algorithm.set( INTEGRATOR_TOLERANCE, 1e-6 );
    //algorithm.set( DISCRETIZATION_TYPE, SINGLE_SHOOTING   );
    algorithm.set( HESSIAN_APPROXIMATION, FULL_BFGS_UPDATE );
    algorithm.set( INTEGRATOR_TYPE, INT_RK45 );
    algorithm.set( MAX_NUM_ITERATIONS, 500 );
    algorithm.set( KKT_TOLERANCE, 1e-16);
    algorithm << window1;
    algorithm.solve();

    return 0;
}
/* <<< end tutorial code <<< */


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
const double d = 0.0001;
const double c = 0.0001;

se3::Model robot = se3::urdf::buildModel("/local/mgeisert/models/freeflyer.urdf", true);
se3::Data data(robot);


Eigen::VectorXd q(7);
Eigen::VectorXd v(6);
Eigen::VectorXd u(6);
Eigen::VectorXd a(6);
Eigen::VectorXd al(3);
Eigen::VectorXd aa(3);
Eigen::VectorXd al2(3);
Eigen::VectorXd aa2(3);
Eigen::VectorXd ei(6);
Eigen::VectorXd Mi;
Eigen::MatrixXd M(6,6);
Eigen::VectorXd b0;
Eigen::VectorXd b;
Eigen::Matrix3d Mrot;
Eigen::Vector3d MomentRPY;
//Eigen::Quaternion quat;
double r,p,y;
unsigned int i,j;

void myDifferentialEquation( double *x, double *f, void *user_data )
{
    for (i = 0 ; i<6 ; i++) {
        v[i] = x[i+6];
        a[i] = 0.;
        ei[i] = 0.;
    }

//    for (i=0 ; i<16 ; i++) {
//        std::cout << " x   : " << x[i] << std::endl;
//    }

    for (i=0 ; i<3 ; i++) {
        q[i] = x[i];
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

//    for (i=0 ; i<12 ; i++) {
//        std::cout << " q   : " << q[i] << std::endl;
//    }

    MomentRPY(0) = d*(x[13]-x[14]);
    MomentRPY(1) = d*(x[16]-x[15]);
    MomentRPY(2) = c*(x[13]+x[14]-x[15]-x[16]);

    //MomentRPY = Mrot*MomentRPY;

    u[0] = 0.;//- (x[13]+x[14]+x[15]+x[16])*sin(p);
    u[1] = 0.;//(x[13]+x[14]+x[15]+x[16])*sin(r)*cos(p);
    u[2] = (x[13]+x[14]+x[15]+x[16]);//(x[13]+x[14]+x[15]+x[16])*cos(r)*cos(p);
    u[3] = MomentRPY(0);
    u[4] = MomentRPY(1);
    u[5] = MomentRPY(2);


    b0 = rnea(robot, data, q, a, a);

  //  std::cout << "b0 : \n" << b0 << std::endl;

    for (i=0 ; i<6 ; i++) {
        ei(i)=1;
        Mi = rnea(robot, data, q, a, ei);
        for (j=0 ; j<6 ; j++ ) {
            M(i,j) = Mi(j)-b0(j);
        }
        ei(i)=0;
    }

//    std::cout << " M rnea :\n" <<  M << std::endl << std::endl;
//    std::cout << "crba :\n" << crba(robot, data, q) << std::endl << std::endl;
//    self adjoint upper ???

    b = rnea(robot, data, q, v, a);

    a = M.llt().solve(u-b);
//    std::cout << "F tot moteurs : \n" << x[12]+x[13]+x[14]+x[15] << std::endl << std::endl;
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

    // accelerations lineaires
    for (i=0 ; i<3 ; i++) {
        al[i] = a[i];
        aa[i] = a[i+3];
    }

    al2 = Mrot*al;
    for (i=0 ; i<3 ;i++) {
        f[i+6] = al2[i];
    }

    //accelerations angulaires
    aa2 = Mrot*aa;
    for (i=0 ; i<3 ;i++) {
        f[i+9] = aa[i];
    }

    // cost
    f[12]=((x[0]-20)*(x[0]-20)+x[2]*x[2]+x[1]*x[1])*0.001;
}

// -------------------------------------------------------------------------
//              USE THE ACADO TOOLKIT TO SOLVE THE PROBLEM:
// -------------------------------------------------------------------------


USING_NAMESPACE_ACADO


int main( ){


    // INTRODUCE THE VARIABLES:
    // --------------------------------------------------
    DifferentialState     xp, yp, zp, phi, theta, psi, vx, vy, vz, p, q, r, cost;
    Control               u1, u2, u3, u4;
    Parameter             T;

    // DEFINE THE DIMENSIONS OF THE C-FUNCTIONS:
    // --------------------------------------------------
    CFunction F( 13, myDifferentialEquation     );
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

    IntermediateState x(17);

     x(0) = xp;  x(1) = yp;  x(2) = zp;  x(3) = phi;  x(4) = theta;  x(5) = psi;
     x(6) = vx;  x(7) = vy;  x(8) = vz;  x(9) = p; x(10) = q; x(11) = r; x(12) = cost;
     x(13) = u1; x(14) = u2; x(15) = u3; x(16) = u4;


    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, T, 20 );

    ocp.minimizeMayerTerm(T);//x(12));

    ocp.subjectTo( f << F(x)*T );

    ocp.subjectTo( 1 <= T <= 20);

    // Start constraints
    for (unsigned int k=0 ; k < 13 ; k++) {
        ocp.subjectTo(AT_START, x(k) == 0.0);
    }

    // End constraints
    ocp.subjectTo( AT_END, x(0) ==  0.0 );
    ocp.subjectTo( AT_END, x(1) ==  0.0 );
    ocp.subjectTo( AT_END, x(2) ==  10.0 );
    for (unsigned int k=3 ; k < 12 ; k++) {
        ocp.subjectTo(AT_END, x(k) == 0.0 );
    }

    // Limits on commands
    double com_min=0;
    double com_max=30;
    ocp.subjectTo(com_min <= x(13) <= com_max);
    ocp.subjectTo(com_min <= x(14) <= com_max);
    ocp.subjectTo(com_min <= x(15) <= com_max);
    ocp.subjectTo(com_min <= x(16) <= com_max);

//    ocp.subjectTo(x(5) == 0.);
//    ocp.subjectTo( -1.57 <= x(4) <= 1.57);
//    ocp.subjectTo( -1.57 <= x(3) <= 1.57);
//    ocp.subjectTo(x(3) == 0.);

    // VISUALIZE THE RESULTS IN A GNUPLOT WINDOW:
    // ------------------------------------------
    GnuplotWindow window1(PLOT_AT_EACH_ITERATION);
    window1.addSubplot( x(0),"x" );
    window1.addSubplot( x(1) , "y");
    window1.addSubplot( x(2), "z");
    window1.addSubplot( x(3),"phi" );
    window1.addSubplot( x(4) , "theta");
    window1.addSubplot( x(5), "psi");
    window1.addSubplot( x(13),"u1" );
    window1.addSubplot( x(14),"u2" );
    window1.addSubplot( x(15),"u3" );
    window1.addSubplot( x(16),"u4" );

    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
    OptimizationAlgorithm algorithm(ocp);
    // algorithm.set( INTEGRATOR_TOLERANCE, 1e-3 );
    //algorithm.set( DISCRETIZATION_TYPE, SINGLE_SHOOTING   );
    //algorithm.set( HESSIAN_APPROXIMATION, FULL_BFGS_UPDATE );
    algorithm.set( INTEGRATOR_TYPE, INT_RK45 );
    algorithm.set( MAX_NUM_ITERATIONS, 500 );
    algorithm.set( KKT_TOLERANCE, 1e-12);
    algorithm << window1;

//    double testx[16];
//    for (i=0 ; i<12 ; i++) {
//        testx[i] = 0.;
//    }
//    testx[4] = 0.;
//    testx[3] = 0.;
//    testx[13] = 12.;
//    testx[14] = 12.;
//    testx[15] = 12.;
//    testx[16] = 12.;

//    double testf[13];
//    for (i=0 ; i<12 ; i++) {
//        testf[i] = 0;
//    }
//    myDifferentialEquation( testx, testf, testf );

//    for ( unsigned int k=0 ; k<13 ; k++ ) {
//        std::cout << testf[k] << std::endl;
//    }

    algorithm.solve();

    return 0;

}
/* <<< end tutorial code <<< */



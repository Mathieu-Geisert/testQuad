
#include <acado_optimal_control.hpp>
#include <acado_gnuplot.hpp>
#include <acado_code_generation.hpp>
#include <acado/process/process.hpp>
#include <acado_toolkit.hpp>

int main( ){

    USING_NAMESPACE_ACADO

    DifferentialState     x,y,z,vx,vy,vz,phi,theta,psi,p,q,r,u1,u2,u3,u4;
    Control               vu1,vu2,vu3,vu4;
    DifferentialEquation  f;


    // Quad constants
    const double c = 0.00001;
    const double Cf = 0.00065;
    const double d = 0.250;
    const double Jx = 0.018;
    const double Jy = 0.018;
    const double Jz = 0.026;
    const double Im = 0.0001;
    const double m = 0.9;
    const double g = 9.81;
    const double Cx = 0.1;


    // Minimization Weights
    double coeffU = 0.000000000000001;
    double coeffX = .00001;
    double coeffX2 = coeffX * 10.;
    double coeffX3 = coeffX * 0.00001;
    double coeffO = -coeffX * 0.1;


    // final values
    double xf = 3., yf = 0., zf = 0.;

    // Temps max
    double tps = 10;
    int nb_it = 20.;
    double tmpc = 0.2;
    int mpc_nb_it = 1;
    int ocp_nb_it = 294;


    OCP ocpMPC( 0.0, tps, nb_it );

    //Cost function
    Function h, hf;
    h << x;
    h << y;
    h << z;

    h << vu1;
    h << vu2;
    h << vu3;
    h << vu4;

    h << p;
    h << q;
    h << r;

    hf << x;
    hf << y;
    hf << z;

    DMatrix Q(10,10), Qf(3,3);
    Q(0,0) = .001;//coeffX;
    Q(1,1) = .001;//coeffX;
    Q(2,2) = .001;//coeffX;

    Q(3,3) = 0.0000001;//coeffU;
    Q(4,4) = 0.0000001;//coeffU;
    Q(5,5) = 0.0000001;//coeffU;
    Q(6,6) = 0.0000001;//coeffU;

    Q(7,7) = 0.010;//coeffX2;
    Q(8,8) = 0.010;//coeffX2;
    Q(9,9) = 0.010;//coeffX2;

    Qf(0,0) = 0.0000001;//coeffX;//*100;
    Qf(1,1) = 0.0000001;//coeffX;//*100;
    Qf(2,2) = 0.0000001;//coeffX;//*100;

    DVector reff(3), ref(10);
    ref(0) = xf;
    ref(1) = yf;
    ref(2) = zf;
    ref(3) = 0.;
    ref(4) = 0.;
    ref(5) = 0.;
    ref(6) = 0.;
    ref(7) = 0.;
    ref(8) = 0.;
    ref(9) = 0.;
//    ref(10) = 0.;
//    ref(11) = 0.;
//    ref(10) = 0.;

    reff(0) = xf;
    reff(1) = yf;
    reff(2) = zf;

    ocpMPC.minimizeLSQ ( Q, h);
    ocpMPC.minimizeLSQEndTerm( Qf, hf);


    // Dynamic equations
    f << dot(x) == vx;
    f << dot(y) == vy;
    f << dot(z) == vz;
    f << dot(vx) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(theta)/m;
    f << dot(vy) == -Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*sin(psi)*cos(theta)/m;
    f << dot(vz) == Cf*(u1*u1+u2*u2+u3*u3+u4*u4)*cos(psi)*cos(theta)/m - g;
    f << dot(phi) == -cos(phi)*tan(theta)*p+sin(phi)*tan(theta)*q+r;
    f << dot(theta) == sin(phi)*p+cos(phi)*q;
    f << dot(psi) == cos(phi)/cos(theta)*p-sin(phi)/cos(theta)*q;
    f << dot(p) == (d*Cf*(u1*u1-u2*u2)+(Jy-Jz)*q*r)/Jx;
    f << dot(q) == (d*Cf*(u4*u4-u3*u3)+(Jz-Jx)*p*r)/Jy;
    f << dot(r) == (c*(u1*u1+u2*u2-u3*u3-u4*u4)+(Jx-Jy)*p*q)/Jz;
    f << dot(u1) == vu1;
    f << dot(u2) == vu2;
    f << dot(u3) == vu3;
    f << dot(u4) == vu4;

    ocpMPC.subjectTo( f );

    ocpMPC.subjectTo( 16 <= u1 <= 95 );
    ocpMPC.subjectTo( 16 <= u2 <= 95 );
    ocpMPC.subjectTo( 16 <= u3 <= 95 );
    ocpMPC.subjectTo( 16 <= u4 <= 95 );
    ocpMPC.subjectTo( -100 <= vu1 <= 100 );
    ocpMPC.subjectTo( -100 <= vu2 <= 100 );
    ocpMPC.subjectTo( -100 <= vu3 <= 100 );
    ocpMPC.subjectTo( -100 <= vu4 <= 100 );

    ocpMPC.subjectTo( -1. <= theta <= 1. );





     OCPexport mpc( ocpMPC );

     mpc.set( HESSIAN_APPROXIMATION,       GAUSS_NEWTON    );
     mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING );
     mpc.set( INTEGRATOR_TYPE,             INT_RK4         );
     mpc.set( NUM_INTEGRATOR_STEPS,        30              );

     mpc.set( QP_SOLVER,                   QP_QPOASES      );
     mpc.set( HOTSTART_QP,                 YES             );
    // 	mpc.set( LEVENBERG_MARQUARDT,         1.0e-4          );
     mpc.set( GENERATE_TEST_FILE,          YES             );
     mpc.set( GENERATE_MAKE_FILE,          YES             );
     mpc.set( CG_USE_ARRIVAL_COST, NO);
    //    mpc.set( GENERATE_MATLAB_INTERFACE,   YES             );
    //    mpc.set( GENERATE_SIMULINK_INTERFACE, YES             );
    // 	mpc.set( USE_SINGLE_PRECISION,        YES             );

        if (mpc.exportCode( "code_generation" ) != SUCCESSFUL_RETURN)
            exit( EXIT_FAILURE );

        mpc.printDimensionsQP( );




}

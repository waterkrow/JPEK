/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7606670772205185248);
void inv_err_fun(double *nom_x, double *true_x, double *out_1293942527211020556);
void H_mod_fun(double *state, double *out_7534678290706100191);
void f_fun(double *state, double dt, double *out_3603120626329206438);
void F_fun(double *state, double dt, double *out_2719231095786819246);
void h_25(double *state, double *unused, double *out_9055706363599892855);
void H_25(double *state, double *unused, double *out_2765573616119153967);
void h_24(double *state, double *unused, double *out_6155136138083892007);
void H_24(double *state, double *unused, double *out_4101441563566573227);
void h_30(double *state, double *unused, double *out_490286340754496395);
void H_30(double *state, double *unused, double *out_7273857388433495179);
void h_26(double *state, double *unused, double *out_7885966461856942546);
void H_26(double *state, double *unused, double *out_8393197622641477367);
void h_27(double *state, double *unused, double *out_4895203673450437546);
void H_27(double *state, double *unused, double *out_240028037382778779);
void h_29(double *state, double *unused, double *out_8550745812536951128);
void H_29(double *state, double *unused, double *out_7637579295250209441);
void h_28(double *state, double *unused, double *out_6129688803721663033);
void H_28(double *state, double *unused, double *out_5707845563936280191);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);

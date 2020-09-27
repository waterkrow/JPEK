
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7606670772205185248) {
   out_7606670772205185248[0] = delta_x[0] + nom_x[0];
   out_7606670772205185248[1] = delta_x[1] + nom_x[1];
   out_7606670772205185248[2] = delta_x[2] + nom_x[2];
   out_7606670772205185248[3] = delta_x[3] + nom_x[3];
   out_7606670772205185248[4] = delta_x[4] + nom_x[4];
   out_7606670772205185248[5] = delta_x[5] + nom_x[5];
   out_7606670772205185248[6] = delta_x[6] + nom_x[6];
   out_7606670772205185248[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1293942527211020556) {
   out_1293942527211020556[0] = -nom_x[0] + true_x[0];
   out_1293942527211020556[1] = -nom_x[1] + true_x[1];
   out_1293942527211020556[2] = -nom_x[2] + true_x[2];
   out_1293942527211020556[3] = -nom_x[3] + true_x[3];
   out_1293942527211020556[4] = -nom_x[4] + true_x[4];
   out_1293942527211020556[5] = -nom_x[5] + true_x[5];
   out_1293942527211020556[6] = -nom_x[6] + true_x[6];
   out_1293942527211020556[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7534678290706100191) {
   out_7534678290706100191[0] = 1.0;
   out_7534678290706100191[1] = 0.0;
   out_7534678290706100191[2] = 0.0;
   out_7534678290706100191[3] = 0.0;
   out_7534678290706100191[4] = 0.0;
   out_7534678290706100191[5] = 0.0;
   out_7534678290706100191[6] = 0.0;
   out_7534678290706100191[7] = 0.0;
   out_7534678290706100191[8] = 0.0;
   out_7534678290706100191[9] = 1.0;
   out_7534678290706100191[10] = 0.0;
   out_7534678290706100191[11] = 0.0;
   out_7534678290706100191[12] = 0.0;
   out_7534678290706100191[13] = 0.0;
   out_7534678290706100191[14] = 0.0;
   out_7534678290706100191[15] = 0.0;
   out_7534678290706100191[16] = 0.0;
   out_7534678290706100191[17] = 0.0;
   out_7534678290706100191[18] = 1.0;
   out_7534678290706100191[19] = 0.0;
   out_7534678290706100191[20] = 0.0;
   out_7534678290706100191[21] = 0.0;
   out_7534678290706100191[22] = 0.0;
   out_7534678290706100191[23] = 0.0;
   out_7534678290706100191[24] = 0.0;
   out_7534678290706100191[25] = 0.0;
   out_7534678290706100191[26] = 0.0;
   out_7534678290706100191[27] = 1.0;
   out_7534678290706100191[28] = 0.0;
   out_7534678290706100191[29] = 0.0;
   out_7534678290706100191[30] = 0.0;
   out_7534678290706100191[31] = 0.0;
   out_7534678290706100191[32] = 0.0;
   out_7534678290706100191[33] = 0.0;
   out_7534678290706100191[34] = 0.0;
   out_7534678290706100191[35] = 0.0;
   out_7534678290706100191[36] = 1.0;
   out_7534678290706100191[37] = 0.0;
   out_7534678290706100191[38] = 0.0;
   out_7534678290706100191[39] = 0.0;
   out_7534678290706100191[40] = 0.0;
   out_7534678290706100191[41] = 0.0;
   out_7534678290706100191[42] = 0.0;
   out_7534678290706100191[43] = 0.0;
   out_7534678290706100191[44] = 0.0;
   out_7534678290706100191[45] = 1.0;
   out_7534678290706100191[46] = 0.0;
   out_7534678290706100191[47] = 0.0;
   out_7534678290706100191[48] = 0.0;
   out_7534678290706100191[49] = 0.0;
   out_7534678290706100191[50] = 0.0;
   out_7534678290706100191[51] = 0.0;
   out_7534678290706100191[52] = 0.0;
   out_7534678290706100191[53] = 0.0;
   out_7534678290706100191[54] = 1.0;
   out_7534678290706100191[55] = 0.0;
   out_7534678290706100191[56] = 0.0;
   out_7534678290706100191[57] = 0.0;
   out_7534678290706100191[58] = 0.0;
   out_7534678290706100191[59] = 0.0;
   out_7534678290706100191[60] = 0.0;
   out_7534678290706100191[61] = 0.0;
   out_7534678290706100191[62] = 0.0;
   out_7534678290706100191[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_3603120626329206438) {
   out_3603120626329206438[0] = state[0];
   out_3603120626329206438[1] = state[1];
   out_3603120626329206438[2] = state[2];
   out_3603120626329206438[3] = state[3];
   out_3603120626329206438[4] = state[4];
   out_3603120626329206438[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3603120626329206438[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3603120626329206438[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2719231095786819246) {
   out_2719231095786819246[0] = 1;
   out_2719231095786819246[1] = 0;
   out_2719231095786819246[2] = 0;
   out_2719231095786819246[3] = 0;
   out_2719231095786819246[4] = 0;
   out_2719231095786819246[5] = 0;
   out_2719231095786819246[6] = 0;
   out_2719231095786819246[7] = 0;
   out_2719231095786819246[8] = 0;
   out_2719231095786819246[9] = 1;
   out_2719231095786819246[10] = 0;
   out_2719231095786819246[11] = 0;
   out_2719231095786819246[12] = 0;
   out_2719231095786819246[13] = 0;
   out_2719231095786819246[14] = 0;
   out_2719231095786819246[15] = 0;
   out_2719231095786819246[16] = 0;
   out_2719231095786819246[17] = 0;
   out_2719231095786819246[18] = 1;
   out_2719231095786819246[19] = 0;
   out_2719231095786819246[20] = 0;
   out_2719231095786819246[21] = 0;
   out_2719231095786819246[22] = 0;
   out_2719231095786819246[23] = 0;
   out_2719231095786819246[24] = 0;
   out_2719231095786819246[25] = 0;
   out_2719231095786819246[26] = 0;
   out_2719231095786819246[27] = 1;
   out_2719231095786819246[28] = 0;
   out_2719231095786819246[29] = 0;
   out_2719231095786819246[30] = 0;
   out_2719231095786819246[31] = 0;
   out_2719231095786819246[32] = 0;
   out_2719231095786819246[33] = 0;
   out_2719231095786819246[34] = 0;
   out_2719231095786819246[35] = 0;
   out_2719231095786819246[36] = 1;
   out_2719231095786819246[37] = 0;
   out_2719231095786819246[38] = 0;
   out_2719231095786819246[39] = 0;
   out_2719231095786819246[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2719231095786819246[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2719231095786819246[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2719231095786819246[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2719231095786819246[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2719231095786819246[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2719231095786819246[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2719231095786819246[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2719231095786819246[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2719231095786819246[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2719231095786819246[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2719231095786819246[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2719231095786819246[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2719231095786819246[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2719231095786819246[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2719231095786819246[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2719231095786819246[56] = 0;
   out_2719231095786819246[57] = 0;
   out_2719231095786819246[58] = 0;
   out_2719231095786819246[59] = 0;
   out_2719231095786819246[60] = 0;
   out_2719231095786819246[61] = 0;
   out_2719231095786819246[62] = 0;
   out_2719231095786819246[63] = 1;
}
void h_25(double *state, double *unused, double *out_9055706363599892855) {
   out_9055706363599892855[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2765573616119153967) {
   out_2765573616119153967[0] = 0;
   out_2765573616119153967[1] = 0;
   out_2765573616119153967[2] = 0;
   out_2765573616119153967[3] = 0;
   out_2765573616119153967[4] = 0;
   out_2765573616119153967[5] = 0;
   out_2765573616119153967[6] = 1;
   out_2765573616119153967[7] = 0;
}
void h_24(double *state, double *unused, double *out_6155136138083892007) {
   out_6155136138083892007[0] = state[4];
   out_6155136138083892007[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4101441563566573227) {
   out_4101441563566573227[0] = 0;
   out_4101441563566573227[1] = 0;
   out_4101441563566573227[2] = 0;
   out_4101441563566573227[3] = 0;
   out_4101441563566573227[4] = 1;
   out_4101441563566573227[5] = 0;
   out_4101441563566573227[6] = 0;
   out_4101441563566573227[7] = 0;
   out_4101441563566573227[8] = 0;
   out_4101441563566573227[9] = 0;
   out_4101441563566573227[10] = 0;
   out_4101441563566573227[11] = 0;
   out_4101441563566573227[12] = 0;
   out_4101441563566573227[13] = 1;
   out_4101441563566573227[14] = 0;
   out_4101441563566573227[15] = 0;
}
void h_30(double *state, double *unused, double *out_490286340754496395) {
   out_490286340754496395[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7273857388433495179) {
   out_7273857388433495179[0] = 0;
   out_7273857388433495179[1] = 0;
   out_7273857388433495179[2] = 0;
   out_7273857388433495179[3] = 0;
   out_7273857388433495179[4] = 1;
   out_7273857388433495179[5] = 0;
   out_7273857388433495179[6] = 0;
   out_7273857388433495179[7] = 0;
}
void h_26(double *state, double *unused, double *out_7885966461856942546) {
   out_7885966461856942546[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8393197622641477367) {
   out_8393197622641477367[0] = 0;
   out_8393197622641477367[1] = 0;
   out_8393197622641477367[2] = 0;
   out_8393197622641477367[3] = 0;
   out_8393197622641477367[4] = 0;
   out_8393197622641477367[5] = 0;
   out_8393197622641477367[6] = 0;
   out_8393197622641477367[7] = 1;
}
void h_27(double *state, double *unused, double *out_4895203673450437546) {
   out_4895203673450437546[0] = state[3];
}
void H_27(double *state, double *unused, double *out_240028037382778779) {
   out_240028037382778779[0] = 0;
   out_240028037382778779[1] = 0;
   out_240028037382778779[2] = 0;
   out_240028037382778779[3] = 1;
   out_240028037382778779[4] = 0;
   out_240028037382778779[5] = 0;
   out_240028037382778779[6] = 0;
   out_240028037382778779[7] = 0;
}
void h_29(double *state, double *unused, double *out_8550745812536951128) {
   out_8550745812536951128[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7637579295250209441) {
   out_7637579295250209441[0] = 0;
   out_7637579295250209441[1] = 1;
   out_7637579295250209441[2] = 0;
   out_7637579295250209441[3] = 0;
   out_7637579295250209441[4] = 0;
   out_7637579295250209441[5] = 0;
   out_7637579295250209441[6] = 0;
   out_7637579295250209441[7] = 0;
}
void h_28(double *state, double *unused, double *out_6129688803721663033) {
   out_6129688803721663033[0] = state[5];
   out_6129688803721663033[1] = state[6];
}
void H_28(double *state, double *unused, double *out_5707845563936280191) {
   out_5707845563936280191[0] = 0;
   out_5707845563936280191[1] = 0;
   out_5707845563936280191[2] = 0;
   out_5707845563936280191[3] = 0;
   out_5707845563936280191[4] = 0;
   out_5707845563936280191[5] = 1;
   out_5707845563936280191[6] = 0;
   out_5707845563936280191[7] = 0;
   out_5707845563936280191[8] = 0;
   out_5707845563936280191[9] = 0;
   out_5707845563936280191[10] = 0;
   out_5707845563936280191[11] = 0;
   out_5707845563936280191[12] = 0;
   out_5707845563936280191[13] = 0;
   out_5707845563936280191[14] = 1;
   out_5707845563936280191[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}

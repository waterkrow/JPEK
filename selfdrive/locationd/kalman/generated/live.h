/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1106718749800000455);
void inv_err_fun(double *nom_x, double *true_x, double *out_6430679954987258979);
void H_mod_fun(double *state, double *out_1497004828671625821);
void f_fun(double *state, double dt, double *out_800714538479493666);
void F_fun(double *state, double dt, double *out_7632238014932592564);
void h_3(double *state, double *unused, double *out_909203479273269773);
void H_3(double *state, double *unused, double *out_3786568584327115003);
void h_4(double *state, double *unused, double *out_6409091132219613490);
void H_4(double *state, double *unused, double *out_3768453196732827335);
void h_9(double *state, double *unused, double *out_111816081692312071);
void H_9(double *state, double *unused, double *out_7751598222144219849);
void h_10(double *state, double *unused, double *out_2097570182416636270);
void H_10(double *state, double *unused, double *out_6861817705678867656);
void h_12(double *state, double *unused, double *out_4847529133592198187);
void H_12(double *state, double *unused, double *out_3693417686661046809);
void h_13(double *state, double *unused, double *out_248610148820804127);
void H_13(double *state, double *unused, double *out_628264136936017662);
void h_14(double *state, double *unused, double *out_111816081692312071);
void H_14(double *state, double *unused, double *out_7751598222144219849);
void h_19(double *state, double *unused, double *out_2933544072579311737);
void H_19(double *state, double *unused, double *out_1687320774814947635);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);
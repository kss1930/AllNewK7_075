/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5302815460919130239);
void inv_err_fun(double *nom_x, double *true_x, double *out_308939119679128541);
void H_mod_fun(double *state, double *out_8205907585195124657);
void f_fun(double *state, double dt, double *out_6664519916989496322);
void F_fun(double *state, double dt, double *out_5244324895196501156);
void h_3(double *state, double *unused, double *out_531919624104925139);
void H_3(double *state, double *unused, double *out_7904903972697666949);
void h_4(double *state, double *unused, double *out_8996109340938877082);
void H_4(double *state, double *unused, double *out_5811366122818198717);
void h_9(double *state, double *unused, double *out_6663746482076330039);
void H_9(double *state, double *unused, double *out_4312338352687437269);
void h_10(double *state, double *unused, double *out_615078502503354550);
void H_10(double *state, double *unused, double *out_2477854920420105880);
void h_12(double *state, double *unused, double *out_4530557800809916787);
void H_12(double *state, double *unused, double *out_1172698024690037189);
void h_13(double *state, double *unused, double *out_6306585579518135669);
void H_13(double *state, double *unused, double *out_208141959279097178);
void h_14(double *state, double *unused, double *out_6663746482076330039);
void H_14(double *state, double *unused, double *out_4312338352687437269);
void h_19(double *state, double *unused, double *out_5153296137457454737);
void H_19(double *state, double *unused, double *out_5962883415731621953);
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
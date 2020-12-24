/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7116643577288957489);
void inv_err_fun(double *nom_x, double *true_x, double *out_6298163448908206639);
void H_mod_fun(double *state, double *out_1465394481245289132);
void f_fun(double *state, double dt, double *out_9040189776974115728);
void F_fun(double *state, double dt, double *out_6410369171251356944);
void h_25(double *state, double *unused, double *out_5868539678396203889);
void H_25(double *state, double *unused, double *out_8320643592696217796);
void h_24(double *state, double *unused, double *out_6701196472411751334);
void H_24(double *state, double *unused, double *out_4467800544752430672);
void h_30(double *state, double *unused, double *out_8837035849074344973);
void H_30(double *state, double *unused, double *out_1863866820764859390);
void h_26(double *state, double *unused, double *out_7692231790682852840);
void H_26(double *state, double *unused, double *out_6608978488032716804);
void h_27(double *state, double *unused, double *out_4534341981671667092);
void H_27(double *state, double *unused, double *out_3047456074598828192);
void h_29(double *state, double *unused, double *out_912993137838891758);
void H_29(double *state, double *unused, double *out_5696739446679887380);
void h_28(double *state, double *unused, double *out_8018911234358471622);
void H_28(double *state, double *unused, double *out_1479682595142529914);
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

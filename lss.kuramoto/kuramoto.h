void
init(double c, int n_grid, double T0,
     int n_chunk, double t_chunk, double dt_max);

double
project_ddt(int i_chunk, double * v);

void
tangent(int i_chunk, double * v0, int inhomo);

void
adjoint(int i_chunk, double * w0, double forcing);

int N_GRID, N_CHUNK, N_STEP;
double C_CONST;
double DT_STEP;
double *** SOLN_U;
double *** SOLN_V;

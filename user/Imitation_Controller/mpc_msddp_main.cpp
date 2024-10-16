#include "MPC_MSDDP_Solver.h"

int main()
{
    MPC_MSDDP_Solver<double> mpc;
    mpc.initialize();
    mpc.run();
}
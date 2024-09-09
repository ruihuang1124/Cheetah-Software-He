/*!
 * @file opt_main.cpp
 * @brief Main Function for the standalone DDP trajectory optimizer
 *
 * The main function initilizes the DDP solver updates upon the request
 * of the main Imitation Controller
 */

#include "MPCSolver.h"
#include "HKDLog.h"

template<typename T>
void MPCSolver<T>::initialize() {
    /* Initialize tunning parameters for DDP solver */
    ddp_options.update_penalty = 7;
    ddp_options.update_relax = 0.1;
    ddp_options.update_ReB = 8;
    ddp_options.update_regularization = 4;
    ddp_options.max_DDP_iter = 10;
    ddp_options.max_AL_iter = 15;
    ddp_options.DDP_thresh = 1e-02;
    ddp_options.AL_active = 1;
    ddp_options.ReB_active = 1;
    ddp_options.pconstr_thresh = .003;
    ddp_options.tconstr_thresh = .003;

    mpc_config.plan_duration = 0.44;
    mpc_config.nsteps_between_mpc = 2;
    mpc_config.timeStep = 0.011;
    dt_mpc = mpc_config.timeStep;
    opt_ref.initialize_referenceData(mpc_config.plan_duration);// first set data in HKDProblem.

    opt_problem.setup(&opt_problem_data, mpc_config);// ref_date in HKD is the same as it in opt_problem_data here.
    opt_problem.initialization();
    opt_problem.print();

    if (!almostEqual_number(dt_mpc, opt_problem_data.ref_data_ptr->dt)) {
        printf(RED);
        printf("problem timestep and reference timestep not match \n");
        printf("problem timestep = %f \n", dt_mpc);
        printf("reference timestep = %f \n", opt_problem_data.ref_data_ptr->dt);
        printf(RESET);
        return;
    }

    // set the initial condition
    xinit.setZero(24);
    body << 0, 0, 0, 0, 0, 0.28, 0, 0, 0, 0, 0, 0;
    qJ << 0, -0.7696, 1.6114, 0, -0.7696, 1.6114, 0, -0.7696, 1.6114, 0, -0.7696, 1.6114;
    pos = body.segment(3, 3);
    eul = body.head(3);

    // get the contact status of the very first phase
    const auto &c = opt_problem_data.ref_data_ptr->contactSeq[0];
    compute_hkd_state(eul, pos, qJ, qdummy, c);

    xinit << body, qdummy;

    // build the solver and solve the TO problem
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (auto phase: opt_problem_data.phase_ptrs) {
        multiple_phases.push_back(phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(xinit);
    solver.solve(ddp_options);
    log_trajectory_sequence(opt_problem_data.trajectory_ptrs);
    mpc_iter = 0;

    update_foot_placement();
    publish_mpc_cmd();
//    publish_rfmpc_cmd();
    publish_debugfoot();
    // opt_problem.lcm_publish();

    // ! for rfmpc!
    horizonLengthMPC_ = 10;
    dt_rfmpc_ = dt_mpc;
    Xt_.resize(30, 1);
    Ut_.resize(12, 1);
    Xd_.resize(30, horizonLengthMPC_);
    Ud_.resize(12, horizonLengthMPC_);
    Xt_.setZero();
    Ut_.setZero();
    Xd_.setZero();
    Ud_.setZero();
}

template<typename T>
void MPCSolver<T>::update() {
    mpc_mutex.lock(); // lock mpc to prevent updating while the previous hasn't finished

    // use less iterations when resolving DDP
    ddp_options.max_AL_iter = 2;
    ddp_options.max_DDP_iter = 2;
    mpc_iter++;

    printf("************************************* \n");
    printf("************Resolving MPC************ \n");
    printf("********MPC Iteration = %d*********** \n", mpc_iter);
    printf("************************************* \n");

    /* update the problem */
    opt_problem.update();
//     opt_problem.print();

//    print_mc_state();
    /* update current state*/
    eul << hkd_data.rpy[2], hkd_data.rpy[1], hkd_data.rpy[0];
    pos << hkd_data.p[0], hkd_data.p[1], hkd_data.p[2];
    vel << hkd_data.vWorld[0], hkd_data.vWorld[1], hkd_data.vWorld[2];
    omega << hkd_data.omegaBody[0], hkd_data.omegaBody[1], hkd_data.omegaBody[2];
    for (int i(0); i < 12; i++) { qJ[i] = hkd_data.qJ[i]; }
    const auto &c = opt_problem_data.ref_data_ptr->contactSeq[0];
    compute_hkd_state(eul, pos, qJ, qdummy, c);
    xinit << eul, pos, omega, vel, qdummy;

    /* build solver and solve the TO problem */
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (auto phase: opt_problem_data.phase_ptrs) {
        multiple_phases.push_back(phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(xinit);
    solver.solve(ddp_options);

    update_foot_placement();
    publish_mpc_cmd();
    publish_debugfoot();
    // opt_problem.lcm_publish();
    mpc_mutex.unlock();
}

template<typename T>
void MPCSolver<T>::mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                                       const hkd_data_lcmt *msg) {
    mpc_mutex.lock();

    printf(GRN);
    printf("Received resolving request\n");
    printf(RESET);


    std::memcpy(&hkd_data, msg, sizeof(hkd_data));
    mpc_time_prev = mpc_time;
    mpc_time = hkd_data.mpctime;    
    
    // get the current foot placements
    const auto& current_pf = hkd_data.foot_placements;
    for (int l = 0; l < 4; l++)
    {
        pf[l] << current_pf[3*l], current_pf[3*l + 1], current_pf[3*l + 2];
    }
    mpc_mutex.unlock();
    std::thread solve_mpc_thread(&MPCSolver::update, this);
    solve_mpc_thread.detach(); // detach the thread from the main thread. The thread would exit once it completes
}

/*
    @brief: Go through the trajectory to find the next foot placement
            Needs to update every time after an MPC update
*/
template<typename T>
void MPCSolver<T>::update_foot_placement()
{
    // mpc_mutex.lock();
    int found_next_ctact[4] = {0};
    const int& n_phases = opt_problem_data.ref_data_ptr->n_phases;
    const auto& ctactSeq = opt_problem_data.ref_data_ptr->contactSeq;

    for (int i(0); i < n_phases - 1; i++)
    {
        const auto &ctact = ctactSeq[i];
        const auto &ctactn = ctactSeq[i+1];
        for (int l(0); l < 4; l++)
        {
            // If we havn't found the next contact
            if (!found_next_ctact[l])
            {
                // Check the search pattern [0, 1]. Iterate through the contact sequence until match the search
                // pattern. qdummy associated with 1 and leg l is then foot placement
                if (ctact[l] == 0 && ctactn[l] == 1)
                {
                    const auto &qdummy =  opt_problem_data.trajectory_ptrs[i+1]->Xbar[0].tail(12);
                    pf[l] = qdummy.segment(3 * l, 3).template cast<float>();
                    found_next_ctact[l] = 1;
                }
            }
        }
        // Break if not found after three phases
        if (i>=3)
        {
            break;
        }
    }
    // mpc_mutex.unlock();
}

template<typename T>
void MPCSolver<T>::publish_mpc_cmd()
{
    int num_controls = mpc_config.nsteps_between_mpc;
    num_controls += 10;     // use 10 more controls than control duration to account for delay

    hkd_cmds.N_mpcsteps = num_controls;
    auto &trajseq = opt_problem_data.trajectory_ptrs;
    auto &ctactSeq = opt_problem_data.ref_data_ptr->contactSeq;
    auto &statusDuration = opt_problem_data.ref_data_ptr->statusDuration;
    int k(0), s(0), i(0);

    while (k < hkd_cmds.N_mpcsteps)
    {
        if (s >= trajseq[i]->horizon)
        {
            s = 0;
            i++;
        }
        for (int j = 0; j < 24; j++)
        {
            hkd_cmds.hkd_controls[k][j] = trajseq[i]->Ubar[s][j];            
        }

        hkd_cmds.mpc_times[k] = mpc_time + (k * dt_mpc);
        for (int l = 0; l < 4; l++)
        {
            hkd_cmds.contacts[k][l] = ctactSeq[i][l];
            hkd_cmds.statusTimes[k][l] = statusDuration(l,i);
        }        
        
        s++;
        k++;
    }
    for (int l = 0; l < 4; l++)
    {
        hkd_cmds.foot_placement[3*l] = pf[l][0];
        hkd_cmds.foot_placement[3*l + 1] = pf[l][1];
        hkd_cmds.foot_placement[3*l + 2] = pf[l][2];
    }
    mpc_lcm.publish("mpc_command", &hkd_cmds);
    printf(GRN);
    printf("published a mpc command message \n");
    printf(RESET);
}

template<typename T>
void MPCSolver<T>::publish_debugfoot()
{
    debug_foot_data.N = 0;
    debug_foot_data.contacts.clear();
    debug_foot_data.qdummy.clear();
    const auto& contacts = opt_problem_data.ref_data_ptr->contactSeq;
    int n_phases = opt_problem_data.ref_data_ptr->n_phases;
    for (int i = 0; i < n_phases; i++)
    {
        auto traj = opt_problem_data.trajectory_ptrs[i];
        auto horizon = traj->horizon;
        for (int k = 0; k < horizon; k++)
        {
            vector<float> qdummy(traj->Xbar[k].data()+12,traj->Xbar[k].data()+23);
            vector<int> ctact(contacts[i].data(), contacts[i].data()+3);
            debug_foot_data.contacts.push_back(ctact);
            debug_foot_data.qdummy.push_back(qdummy);
        }
        debug_foot_data.N += horizon;
    }
    mpc_lcm.publish("debug_foot", &debug_foot_data);
}

template
class MPCSolver<double>;


template<typename T>
void MPCSolver<T>::updateRFMPCSolver() {
    mpc_mutex.lock(); // lock mpc to prevent updating while the previous hasn't finished
    mpc_iter++;

    printf("************************************* \n");
    printf("************Resolving RFMPC************ \n");
    printf("********MPC Iteration = %d*********** \n", mpc_iter);
    printf("************************************* \n");

    /* update rfmpc settings*/
    float alpha;
    float Q[12];
    float Qf[12];
    float R[12];
    float Q_MINI[12] = {1e4, 2e5, 3e6, 5e3, 1e3, 1e3, 1e3, 1e6, 800, 40, 400, 10};  // mini cheetah
    float Q_Mini_f[12] = {1e4, 2e5, 3e6, 5e3, 1e3, 1e3, 1e3, 1e6, 800, 40, 400, 10};
    float R_Mini[12] = {5e2, 1e2, 5e2, 5e2, 1e2, 5e2, 5e2, 1e2, 5e2, 5e2, 1e2, 5e2};
    alpha = 4e-3;                                                             // mini cheetah
    memcpy(Q, Q_MINI, sizeof(Q_MINI));
    memcpy(Qf, Q_Mini_f, sizeof(Q_Mini_f));
    memcpy(R, R_Mini, sizeof(R_Mini));
    float *weights_Q = Q;
    float *weights_Qf = Qf;
    float *weights_R = R;
    mpc_setup_problem(dt_rfmpc_, horizonLengthMPC_, 0.4, 180);
    /* update Xd and Ud from desired trajectory */
    opt_problem.update();

    updateReferenceTrajectoryForVBMPC();

    /* update current state*/
//    eul << hkd_rfmpc_data.rpy[2], hkd_rfmpc_data.rpy[1], hkd_rfmpc_data.rpy[0];
//    pos << hkd_rfmpc_data.p[0], hkd_rfmpc_data.p[1], hkd_rfmpc_data.p[2];
//    vel << hkd_rfmpc_data.vWorld[0], hkd_rfmpc_data.vWorld[1], hkd_rfmpc_data.vWorld[2];
//    omega << hkd_rfmpc_data.omegaBody[0], hkd_rfmpc_data.omegaBody[1], hkd_rfmpc_data.omegaBody[2];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rRobot(i, j) = hkd_rfmpc_data.rRobot[i][j];
        }
    }

    const auto &current_pf = hkd_rfmpc_data.foot_placements;
    for (int l = 0; l < 4; l++) {
        pf[l] << current_pf[3 * l], current_pf[3 * l + 1], current_pf[3 * l + 2];
    }
    for (int i = 0; i < 3; ++i) {
        Xt_(i, 0) = hkd_rfmpc_data.p[i];
        Xt_(i + 3, 0) = hkd_rfmpc_data.vWorld[i];
        for (int j = 0; j < 3; ++j) {
            Xt_(6 + i + j * 3, 0) = rRobot(i, j);
        }
        Xt_(15 + i, 0) = hkd_rfmpc_data.omegaBody[i];
        for (int j = 0; j < 4; ++j) {
            Xt_(18 + i + j * 3, 0) = pf[switch_leg_index_from_left_to_right(j)][i];//TODO, check if right.
            Ut_(i + j * 3, 0) = hkd_rfmpc_data.f_current[switch_leg_index_from_left_to_right(j)][i];
        }
    }

//    print_mc_state();
    /* build solver and solve the TO problem */
    update_vbmpc_problem_data_float_latest(Xt_, Ut_, Xd_, Ud_, weights_Q, weights_Qf, weights_R, alpha);
    update_foot_placement();
    publish_rfmpc_cmd();
    publish_debugfoot();
    // opt_problem.lcm_publish();
    mpc_mutex.unlock();
}

template<typename T>
void MPCSolver<T>::rfmpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                                         const hkd_rfmpc_data_lcmt *msg) {
    mpc_mutex.lock();

    printf(GRN);
    printf("Received resolving request\n");
    printf(RESET);


    std::memcpy(&hkd_rfmpc_data, msg, sizeof(hkd_rfmpc_data));
    mpc_time_prev = mpc_time;
    mpc_time = hkd_rfmpc_data.mpctime;

    // get the current foot placements
    const auto &current_pf = hkd_rfmpc_data.foot_placements;
    for (int l = 0; l < 4; l++) {
        pf[l] << current_pf[3 * l], current_pf[3 * l + 1], current_pf[3 * l + 2];
    }
    mpc_mutex.unlock();
    std::thread solve_mpc_thread(&MPCSolver::updateRFMPCSolver, this);
    solve_mpc_thread.detach(); // detach the thread from the main thread. The thread would exit once it completes
}


template<typename T>
void MPCSolver<T>::publish_rfmpc_cmd() {
    int num_controls = mpc_config.nsteps_between_mpc;
    num_controls += 6;     // use 6 more controls than control duration to account for delay
    hkd_cmds.N_mpcsteps = num_controls;
    auto &XrSeq = opt_problem_data.ref_data_ptr;
    auto &ctactSeq = opt_problem_data.ref_data_ptr->contactSeq;
    auto &statusDuration = opt_problem_data.ref_data_ptr->statusDuration;
    int k(0), s(0), i(0);
    while (k < hkd_cmds.N_mpcsteps) {
        if (s >= XrSeq->horizons[i]) {
            s = 0;
            i++;
        }
        for (int leg = 0; leg < 4; leg++) {
            Vec3<float> f;
            Vec3<float> f_delta;
            Vec3<float> ut;
            for (int axis = 0; axis < 3; axis++) {
                f_delta[axis] = get_vbmpc_solution(leg * 3 + axis);  // delta force from ground to leg in world frame
                ut[axis] = Ut_(leg * 3 + axis, 0);
            }
            f = f_delta + ut; // force from ground to leg in world frame
            for (int j = 0; j < 3; ++j) {
                hkd_cmds.hkd_controls[k][switch_leg_index_from_left_to_right(leg)*3+j] = f[j];
            }
        }
        hkd_cmds.mpc_times[k] = mpc_time + (k * dt_mpc);
        for (int l = 0; l < 4; l++) {
            hkd_cmds.contacts[k][l] = ctactSeq[i][l];
            hkd_cmds.statusTimes[k][l] = statusDuration(l, i);
        }
        std::cout<<"calculated hkdcmd.hkd_controls is"<<hkd_cmds.hkd_controls[k]<<"\n";
        s++;
        k++;
    }
    for (int l = 0; l < 4; l++) {
        hkd_cmds.foot_placement[3 * l] = pf[l][0];
        hkd_cmds.foot_placement[3 * l + 1] = pf[l][1];
        hkd_cmds.foot_placement[3 * l + 2] = pf[l][2];
    }
    mpc_lcm.publish("mpc_command", &hkd_cmds);
    printf(GRN);
    printf("published a rfmpc command message \n");
    printf(RESET);
}


template<typename T>
int MPCSolver<T>::switch_leg_index_from_left_to_right(int vbmpc_leg_index) {// in vbmpc, leg index starts from left while in our robot, leg index starts from right.
    int control_leg_index = 0;
    if (vbmpc_leg_index == 0) {        // 0 FL in vbmpc
        control_leg_index = 1;          // 1 FL in control loop
    } else if (vbmpc_leg_index == 1) {  // 1 FR in vbmpc
        control_leg_index = 0;          // 0 FR in control loop
    } else if (vbmpc_leg_index == 2) {  // 2 HL in vbmpc
        control_leg_index = 3;          // 3 HL in control loop
    } else {                           // 3 HR in vbmpc
        control_leg_index = 2;          // 2 HR in control loop
    }
    return control_leg_index;
}

template<typename T>
void MPCSolver<T>::rpy_to_rotational_matrix_R(Eigen::Matrix<float, 3, 3> &R, float *rpy_in) {
    Eigen::Matrix3f Rz, Ry, Rx;

    Rz.setIdentity();
    Ry.setIdentity();
    Rx.setIdentity();

    Rz(0, 0) = cos(rpy_in[2]);
    Rz(0, 1) = -sin(rpy_in[2]);
    Rz(1, 0) = sin(rpy_in[2]);
    Rz(1, 1) = cos(rpy_in[2]);

    Ry(0, 0) = cos(rpy_in[1]);
    Ry(0, 2) = sin(rpy_in[1]);
    Ry(2, 0) = -sin(rpy_in[1]);
    Ry(2, 2) = cos(rpy_in[1]);

    Rx(1, 1) = cos(rpy_in[0]);
    Rx(1, 2) = -sin(rpy_in[0]);
    Rx(2, 1) = sin(rpy_in[0]);
    Rx(2, 2) = cos(rpy_in[0]);

    R = Rz * Ry * Rx;
}


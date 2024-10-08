#include <memory>
#include "HKDProblem.h"
#include "HKDCost.h"
#include "HSDDP_Utils.h"
#include "HKDConstraints.h"
#include "TrajectoryManagement.h"
#include "SinglePhase.h"
#include <functional> // bind function with fixed arguments

using std::bind;
using std::shared_ptr;
using namespace std::placeholders; // predefined position placeholder for bind

template <typename T>
void HKDProblem<T>::initialization()
{
    deque<int> &horizons = ref_data->horizons;  // the ref_date is the most important parameters we need.!!
    deque<VecM<int, 4>> &ctactSeq = ref_data->contactSeq;
    int n_phases = ref_data->n_phases; // see MPCSolver first set data in HKDProblem. the ref_data n_phases is set
    // check if the size of contactSeq greater
    if (n_phases >= ctactSeq.size())
    {
        printf("Contact sequence size should be greater than n_phase\n");
        return;
    }

    grf_reb_param.delta = 1;
    grf_reb_param.delta_min = 0.001;
    grf_reb_param.eps = 0.2;
    swing_reb_param.delta = 2;
    swing_reb_param.delta_min = 0.01;
    swing_reb_param.eps = 0.02;
    td_al_param.lambda = 0;
    td_al_param.sigma = 100;

    for (int i = 0; i < n_phases; i++)
    {
        shared_ptr<SinglePhase<T,24,24,0>> phase;
        phase = make_shared<SinglePhase<T,24,24,0>>();

        shared_ptr<Trajectory<T,24,24,0>> traj;
        traj = make_shared<Trajectory<T,24,24,0>>(timeStep, horizons[i]);
        
        phase->set_trajectory(traj);      

        create_problem_one_phase(phase,i);

        phase->initialization();

        pdata->trajectory_ptrs.push_back(traj);
        pdata->phase_ptrs.push_back(phase);
    }
}

/*
  @brief:   update every nsteps_between_mpc
 */
template <typename T>
void HKDProblem<T>::update()
{
    // use smaller relaxation parameter
    grf_reb_param.delta = .1;

    /* obtain the previous solutions */
    Xbar_prev.clear();
    Ubar_prev.clear();
    Xbar_prev.resize(ref_data->n_phases);
    Ubar_prev.resize(ref_data->n_phases);
    for (int i = 0; i < ref_data->n_phases; i++)
    {
        std::copy(pdata->trajectory_ptrs[i]->Xbar.begin(), pdata->trajectory_ptrs[i]->Xbar.end(),
                    std::back_inserter(Xbar_prev[i]));
        std::copy(pdata->trajectory_ptrs[i]->Ubar.begin(), pdata->trajectory_ptrs[i]->Ubar.end(),
                    std::back_inserter(Ubar_prev[i]));
    }
    
    for (int j = 0; j < nsteps_between_mpc; j++)
    {
        reference->step();// ref_date is update here too!!!!

        const auto &horizons = ref_data->horizons; // the ref_data in HKD will be update evertime
        const auto &ctactSeq = ref_data->contactSeq;

        auto& trajectories = pdata->trajectory_ptrs;
        auto& phases = pdata->phase_ptrs;

        if (trajectories[0]->size() <= 2)
        {
            trajectories.pop_front();
            phases.pop_front();

            Xbar_prev.pop_front();
            Ubar_prev.pop_front();
        }
        else
        {
            phases[0]->pop_front();
            Xbar_prev[0].pop_front();
            Ubar_prev[0].pop_front();
        }
        // If trajectories is shorter than expected # phases, grow trajectories and phases by one
        if (trajectories.size() < ref_data->n_phases)
        {
            shared_ptr<Trajectory<T, 24, 24, 0>> traj_to_add;
            traj_to_add = make_shared<Trajectory<T, 24, 24, 0>>(timeStep, horizons.back());

            shared_ptr<SinglePhase<T, 24, 24, 0>> phase_to_add;
            phase_to_add = make_shared<SinglePhase<T, 24, 24, 0>>();             

            phase_to_add->set_trajectory(traj_to_add);

            create_problem_one_phase(phase_to_add, ref_data->n_phases-1);

            phase_to_add->initialization();

            trajectories.push_back(traj_to_add);
            phases.push_back(phase_to_add);
        }
        else
        {
            phases.back()->push_back();
        }
    }
    for (auto phase : pdata->phase_ptrs)
    {
        phase->reset_params();
    }
    
}

template<typename T>
void HKDProblem<T>::create_problem_one_phase(shared_ptr<SinglePhase<T,24,24,0>> phase, int idx)

{
    const auto &ctact = ref_data->contactSeq[idx];
    const auto &ctact_n = ref_data->contactSeq[idx+1];
    VecM<int, 4> touchdown_status;
    touchdown_status.setZero();
    for (int leg = 0; leg < 4; leg++)
    {
        if (ctact[leg] == 0 && ctact_n[leg] == 1)
        {
            touchdown_status[leg] = 1;
        }
    }
//    std::cout<<"touchdown_status is: "<< touchdown_status<<"\n";

    /* specialize dynamics and resetmap  */
    auto dynamics_callback = bind(&HKDModel<T>::dynamics, &hkdModel,
                                  _1, _2, _3, _4, ctact, timeStep);                                  
    auto dynamics_partial_callback =
                                bind(&HKDModel<T>::dynamics_partial, &hkdModel,
                                    _1, _2, _3, _4, _5, _6, ctact, timeStep);

    /* set dynamics */
    phase->set_dynamics(dynamics_callback);
    phase->set_dynamics_partial(dynamics_partial_callback);

    /* set resetmap */
    auto resetmap_callback = bind(&HKDReset<T>::resetmap, &hkdReset,
                                  _1, _2, ctact, ctact_n);
    auto resetmap_partial_callback = bind(&HKDReset<T>::resetmap_partial, &hkdReset,
                                          _1, _2, ctact, ctact_n);

    phase->set_resetmap(resetmap_callback);
    phase->set_resetmap_partial(resetmap_partial_callback);

    /* set cost management */
    // tracking cost
    shared_ptr<HKDTrackingCost<T>> track_cost;
    track_cost = make_shared<HKDTrackingCost<T>>(ctact);
    track_cost->set_reference(&(ref_data->Xr[idx]), &(ref_data->Ur[idx]), &(ref_data->Yr[idx]));
    phase->add_cost(track_cost);

    // foot regularization
    shared_ptr<HKDFootPlaceReg<T>> foot_reg;
    foot_reg = make_shared<HKDFootPlaceReg<T>>(ctact);
    foot_reg->set_reference(&(ref_data->Xr[idx]));
    phase->add_cost(foot_reg);

    // previous solution regularization
    if (idx < Xbar_prev.size())
    {
        shared_ptr<PrevSolution_Reg<T>> prev_reg;
        prev_reg = make_shared<PrevSolution_Reg<T>>(ctact);
        prev_reg->set_state_reference(&(Xbar_prev[idx]));
        prev_reg->set_control_reference(&(Ubar_prev[idx]));
        phase->add_cost(prev_reg);
    }
    

    /* add constraints */
    // if there is foot contact, add GRF constraint
    if (find_eigen(ctact, 1).size() > 0)
    {
        shared_ptr<GRFConstraint<T>> grfConstraint;
        grfConstraint = std::make_shared<GRFConstraint<T>>(ctact);
        grfConstraint->update_horizon_len(ref_data->horizons[idx]);
        grfConstraint->create_data();
        grfConstraint->initialize_params(grf_reb_param);
        phase->add_pathConstraint(grfConstraint);
    }
    
    // touchdown constraint
    if (find_eigen(touchdown_status, 1).size() > 0)
    {
        shared_ptr<TouchDownConstraint<T>> tdConstraint;
        tdConstraint = std::make_shared<TouchDownConstraint<T>>(touchdown_status);
        tdConstraint->create_data();
        tdConstraint->initialize_params(td_al_param);
        phase->add_terminalConstraint(tdConstraint);
    }
}

template<typename T>
void HKDProblem<T>::print(){
    printf("***********HKDProblem************\n");
    printf("number of phases = %lu, \n", ref_data->n_phases);
    printf("size of contact sequence = %lu \n", ref_data->contactSeq.size());
    printf("pidx\t dur\t horizon\n");
    for (int i = 0; i < ref_data->n_phases; i++)
    {
        printf("%5i %10.3f %10lu",
               i, ref_data->endTimes[i]-ref_data->startTimes[i], ref_data->horizons[i]);
        printf("\n");
    }
    printf("contact status \n");
    for (int i = 0; i < ref_data->n_phases+1; i++)
    {
        std::cout << ref_data->contactSeq[i].transpose() << std::endl<< std::endl;
    }
    printf("status durations \n");
    std::cout << ref_data->statusDuration.transpose() << std::endl;
    printf("ref trajectory info is \n");
    int xr_point_size = 0;
    int ur_pint_size = 0;
    for (int i = 0; i < ref_data->Xr.size(); ++i) {
        xr_point_size += ref_data->Xr.at(i).size();
    }
    for (int i = 0; i < ref_data->Ur.size(); ++i) {
        ur_pint_size += ref_data->Ur.at(i).size();
    }
    std::cout <<"Xr size is: "<<xr_point_size<<" with total phase: "<<ref_data->Xr.size()<<"\n";
    for (int i = 0; i < ref_data->Xr.size(); ++i) {
        for (int j = 0; j < ref_data->Xr.at(i).size(); ++j) {
            std::cout <<"Xr roll angle value is: "<< ref_data->Xr.at(i).at(j)[0]<<" with phase "<< i+1 <<"\n";
        }

    }
    std::cout <<"Ur size is: "<<ur_pint_size<<" with total phase: "<<ref_data->Ur.size()<<"\n";
    for (int i = 0; i < ref_data->Ur.size(); ++i) {
        for (int j = 0; j < ref_data->Ur.at(i).size(); ++j) {
            std::cout <<"Ur z value for first leg is: "<< ref_data->Ur.at(i).at(j)[2]<<" with phase "<< i+1 <<"\n";
        }
    }
}

template<typename T>
void HKDProblem<T>::lcm_publish()
{
    int k(0);
    float timeStart = ref_data->startTimes[0];
    
    reset_lcm_data();

    for (int i(0); i < ref_data->n_phases ; i++)
    {
        for (int j(0); j < ref_data->horizons[i]; j++)
        {
            lcm_pdata.times.push_back((timeStart + k * (ref_data->dt)) * 1e06);

            for(int l(0); l<4; l++){
                lcm_pdata.contacts[l].push_back(ref_data->contactSeq[i][l]);
                lcm_pdata.qdummy_r[3*l].push_back(ref_data->Xr[i][j][12+3*l]);
                lcm_pdata.qdummy_r[3*l+1].push_back(ref_data->Xr[i][j][12+3*l+1]);
                lcm_pdata.qdummy_r[3*l+2].push_back(ref_data->Xr[i][j][12+3*l+2]);

                lcm_pdata.qdummy[3*l].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3*l]);
                lcm_pdata.qdummy[3*l+1].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3*l+1]);
                lcm_pdata.qdummy[3*l+2].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3*l+2]);
            }

            for (int d = 0; d < 3; d++)
            {
                lcm_pdata.eul_r[d].push_back(ref_data->Xr[i][j][d]);
                lcm_pdata.pos_r[d].push_back(ref_data->Xr[i][j][3+d]);
                lcm_pdata.vel_r[d].push_back(ref_data->Xr[i][j][6+d]);
                lcm_pdata.omega_r[d].push_back(ref_data->Xr[i][j][9+d]);

                lcm_pdata.eul[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][d]);
                lcm_pdata.pos[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3+d]);
                lcm_pdata.vel[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][6+d]);
                lcm_pdata.omega[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][9+d]);
            }            
            k++;
        }      
    }    
    
    lcm_pdata.n_timesteps = k;
    _lcm.publish("DEBUG_HKDMPC", &lcm_pdata);
}

template<typename T>
void HKDProblem<T>::reset_lcm_data()
{
    lcm_pdata.times.clear();
    lcm_pdata.contacts.clear();    
    lcm_pdata.pos_r.clear();
    lcm_pdata.eul_r.clear();
    lcm_pdata.vel_r.clear();
    lcm_pdata.omega_r.clear();
    lcm_pdata.qdummy_r.clear();

    lcm_pdata.pos.clear();
    lcm_pdata.eul.clear();
    lcm_pdata.vel.clear();
    lcm_pdata.omega.clear();
    lcm_pdata.qdummy.clear();

    lcm_pdata.contacts.resize(4);
    lcm_pdata.pos_r.resize(3);
    lcm_pdata.eul_r.resize(3);
    lcm_pdata.vel_r.resize(3);
    lcm_pdata.omega_r.resize(3);
    lcm_pdata.qdummy_r.resize(12);
    lcm_pdata.pos.resize(3);
    lcm_pdata.eul.resize(3);
    lcm_pdata.vel.resize(3);
    lcm_pdata.omega.resize(3);
    lcm_pdata.qdummy.resize(12);
}

template class HKDProblem<double>;
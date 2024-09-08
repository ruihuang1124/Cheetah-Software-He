#ifndef OPT_MAIN_H
#define OPT_MAIN_H

#include <thread>
#include <mutex>
#include <lcm/lcm-cpp.hpp>
#include "hkd_command_lcmt.hpp"
#include "hkd_data_lcmt.hpp"
#include "hkd_rfmpc_data_lcmt.hpp"
#include "opt_sol_lcmt.hpp"

#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"
#include "HKDModel.h"
#include "HKDContactSchedule.h"
#include "HKDProblem.h"
#include "HKDReset.h"
#include "MultiPhaseDDP.h"
#include "cTypes.h"
#include "utilities.h"
#include "Imitation_Reference.h"
#include "HKDReference.h"
#include "MPC_interface.h"
template<typename T>
class MPCSolver
{
public:
    MPCSolver() : mpc_lcm(getLcmUrl(255)), rfmpc_lcm(getLcmUrl(255))
    {
        use_hkd_ = true;
        // Setup reference
        string imitation_path = "../user/Imitation_Controller/PolicyRollout/";
        string contact_fname = imitation_path + "contact_post.csv";
        string state_fname = imitation_path + "state_post.csv";
        imitation_ref.load_contact_data(contact_fname);
        imitation_ref.load_state_data(state_fname);
        imitation_ref.compute_status_duration();

        opt_ref.set_topLevel_reference(imitation_ref.get_data_ptr());

        opt_problem_data.reference_ptr = &opt_ref;
        opt_problem_data.ref_data_ptr = opt_ref.get_referenceData_ptr();// ref_data ptr is set here, which is exactly data(for HKD) in HKDReference reference data used in optimization!!!( See HKD Reference data definition!)

        // Check LCM initialization
        if (!mpc_lcm.good())
        {
            printf(RED);
            printf("Failed to inialize mpc_lcm for hkd command\n");
            return;
        }

        if (use_hkd_){
            mpc_lcm.subscribe("mpc_data", &MPCSolver::mpcdata_lcm_handler, this);
        } else{
            rfmpc_lcm.subscribe("rfmpc_data", &MPCSolver::rfmpcdata_lcm_handler, this);
        }
        
        first_yaw_flip = true;
        yaw_flip_times = 0;
    }
    void mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                             const hkd_data_lcmt *msg);
    void rfmpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                             const hkd_rfmpc_data_lcmt *msg);
    void publish_mpc_cmd();
    void publish_rfmpc_cmd();
    void publish_debugfoot();
    void print_mc_state();
    void initialize();
    void update();
    void update_foot_placement();
    void updateRFMPCSolver();
    void updateReferenceTrajectoryForVBMPC();
    int switch_leg_index_from_left_to_right(int vbmpc_leg_index);
    void run(){
        if (use_hkd_){
            while (mpc_lcm.handle()==0){}
        } else{
            while (rfmpc_lcm.handle()==0){}
        }
    }
    void rpy_to_rotational_matrix_R(Eigen::Matrix<float, 3, 3>& R, float * rpy_in);

public:
    // MPC
    HKDProblem<T> opt_problem;
    HKDProblemData<T> opt_problem_data;
    HKDReference<T> opt_ref;
    Imitation_Reference<T> imitation_ref; 

    HKDPlanConfig<T> mpc_config;
    HSDDP_OPTION ddp_options;

    T dt_mpc;
    T mpc_time;
    T mpc_time_prev;
    int mpc_iter;
    bool first_yaw_flip;
    int yaw_flip_times;
    
    DVec<T> xinit;
    VecM<T, 12> body, qdummy, qJ;
    VecM<T, 3> pos, eul, vel, omega;
    Eigen::Matrix3f rRobot;

    // LCM message
    hkd_data_lcmt hkd_data;
    hkd_rfmpc_data_lcmt hkd_rfmpc_data;
    hkd_command_lcmt hkd_cmds;
    opt_sol_lcmt debug_foot_data;
    lcm::LCM mpc_lcm;
    lcm::LCM rfmpc_lcm;

    // foot placement
    Vec3<float> pf[4];

    // mutex lock
    std::mutex mpc_mutex;

    // use rfmpc solver or HKD solver
    bool use_hkd_ = true;

    Eigen::MatrixXf Xt_, Ut_, Xd_, Ud_;

    int horizonLengthMPC_;
    float dt_rfmpc_;

};


template<typename T>
void MPCSolver<T>::print_mc_state()
{
    printf("***********mc state*************\n");
    printf("eul = %f %f %f\n", hkd_data.rpy[2], hkd_data.rpy[1], hkd_data.rpy[0]);
    printf("pos = %f %f %f\n", hkd_data.p[2], hkd_data.p[1], hkd_data.p[0]);
    printf("omegaB = %f %f %f\n", hkd_data.omegaBody[2], hkd_data.omegaBody[1], hkd_data.omegaBody[0]);
    printf("vWorld = %f %f %f\n\n", hkd_data.vWorld[2], hkd_data.vWorld[1], hkd_data.vWorld[0]);

//    printf("***********More Info about HKDProblem************\n");
//    printf("number of phases = %lu, \n", opt_problem_data.ref_data_ptr->n_phases);
//    printf("size of contact sequence = %lu \n", opt_problem_data.ref_data_ptr->contactSeq.size());
//    printf("pidx\t dur\t horizon\n");
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->n_phases; i++)
//    {
//        printf("%5i %10.3f %10lu",
//               i, opt_problem_data.ref_data_ptr->endTimes[i]-opt_problem_data.ref_data_ptr->startTimes[i], opt_problem_data.ref_data_ptr->horizons[i]);
//        printf("\n");
//    }
//    printf("contact status \n");
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->n_phases+1; i++)
//    {
//        std::cout << opt_problem_data.ref_data_ptr->contactSeq[i].transpose() << std::endl<< std::endl;
//    }
//    printf("status durations \n");
//    std::cout << opt_problem_data.ref_data_ptr->statusDuration.transpose() << std::endl;


//    printf("ref trajectory info is \n");
//    int xr_point_size = 0;
//    int ur_pint_size = 0;
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->Xr.size(); ++i) {
//        xr_point_size += opt_problem_data.ref_data_ptr->Xr.at(i).size();
//    }
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->Ur.size(); ++i) {
//        ur_pint_size += opt_problem_data.ref_data_ptr->Ur.at(i).size();
//    }
//    std::cout <<"Xr size is: "<<xr_point_size<<" with total phase: "<<opt_problem_data.ref_data_ptr->Xr.size()<<"\n";
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->Xr.size(); ++i) {
//        for (int j = 0; j < opt_problem_data.ref_data_ptr->Xr.at(i).size(); ++j) {
//            std::cout <<"Xr roll angle value is: "<< opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[0]<<" with phase "<< i+1 <<"\n";
//        }
//
//    }
//    std::cout <<"Ur size is: "<<ur_pint_size<<" with total phase: "<<opt_problem_data.ref_data_ptr->Ur.size()<<"\n";
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->Ur.size(); ++i) {
//        for (int j = 0; j < opt_problem_data.ref_data_ptr->Ur.at(i).size(); ++j) {
//            std::cout <<"Ur z value for first leg is: "<< opt_problem_data.ref_data_ptr->Ur.at(i).at(j)[2]<<" with phase "<< i+1 <<"\n";
//        }
//    }

//    int point_index = 0;
//    int point_phase_accumulate_index = 0;
//    for (int i = 0; i < opt_problem_data.ref_data_ptr->Xr.size(); ++i) {
//        point_index = 0;
//        if (i >= 1) {
//            point_phase_accumulate_index += opt_problem_data.ref_data_ptr->Xr.at(i - 1).size();
//        }
//        for (int j = 0; j < opt_problem_data.ref_data_ptr->Xr.at(i).size(); ++j) {
//            point_index = point_phase_accumulate_index + j;
//            std::cout << "xr_point_index is: " << point_index << " with phase " << i + 1 << "\n";
//            if (point_index < horizonLengthMPC_) {
//                Eigen::Matrix<float, 3, 3> R_des(3, 3);
//                R_des.setIdentity();
//                float rpy_desired_input[3];
//                for (int k = 0; k < 3; ++k) {
//                    rpy_desired_input[k] = opt_problem_data.ref_data_ptr->Xr.at(i).at(
//                            j)[k];// rpy, xyz, omega_xyz, v_xyz in Xr;
//                }
//                rpy_to_rotational_matrix_R(R_des, rpy_desired_input);
//                for (int k = 0; k < 3; ++k) { // for Xd: xyz, v_xyz,  R_12, omega_xyz, pF_des(xyz1, xyz2, xyz3, xyz4)
//                    Xd_(k, point_index) = opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[k + 3];
//                    Xd_(k + 3, point_index) = opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[k + 9];
//                    for (int l = 0; l < 3; ++l) {
//                        Xd_(6 + k + l * 3, point_index) = R_des(i, j);  //
//                    }
//                    Xd_(15 + k, point_index) = opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[k + 6];
//                }
//
//                Xd_(0, point_index) = opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[0];
//
//                int num_contact_foot = 0;
//                for (int k = 0; k < 4; ++k) {
//                    num_contact_foot += opt_problem_data.ref_data_ptr->contactSeq[i][k];
//                }
//                if (num_contact_foot == 0) {
//                    Ud_.col(point_index).setZero();
//                } else {
//                    float ref_force_z = 80.0 / num_contact_foot;
//                    for (int k = 0; k < 4; ++k) {
//                        Ud_(2 + 3 * k, point_index) =
//                                opt_problem_data.ref_data_ptr->contactSeq[i][switch_leg_index_from_left_to_right(k)] *
//                                ref_force_z;
//                    }
//                }
//                std::cout << "Xd_ vx value is: " << Xd_(3,point_index) << " with phase " << i + 1 << "\n";
//                std::cout <<"Contact status is: "<< opt_problem_data.ref_data_ptr->contactSeq[i].transpose() <<"\n";
//                std::cout<<"Ref Ud_ value is:"<<Ud_.col(point_index).transpose()<<"\n";
//            } else {
//                break;
//            }
//            std::cout <<"Xr roll angle value is: "<< opt_problem_data.ref_data_ptr->Xr.at(i).at(j)[0]<<" with phase "<< i+1 <<"\n";
//        }
//        if (point_index >= horizonLengthMPC_) {
//            break;
//        }
//
//    }

}

#endif
#ifndef IMITATION_CONTROLLER
#define IMITATION_CONTROLLER

#include <RobotController.h>
#include "Controllers/FootSwingTrajectory.h"
#include <deque>
#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"
#include "ImitationUserParameters.h"
#include <mutex>
#include <thread>
#include <lcm/lcm-cpp.hpp>
#include "hkd_command_lcmt.hpp"
#include "hkd_data_lcmt.hpp"
#include "hkd_rfmpc_data_lcmt.hpp"

using std::deque;
class Imitation_Controller : public RobotController
{
public:
    Imitation_Controller();
    virtual void initializeController();
    virtual void runController();
    virtual void updateVisualization() {}
    virtual ControlParameters *getUserControlParameters()
    {
        return &userParameters;
    }

    void update_mpc_if_needed();
    void update_foot_placement();
    void get_a_val_from_solution_bag();

    void standup_ctrl();
    void standup_ctrl_enter();
    void standup_ctrl_run();
    void locomotion_ctrl();
    void getContactStatus();
    void getStatusDuration();

    void handleMPCcommand(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                          const hkd_command_lcmt *msg);
    void handleMPCLCMthread();
    void avoid_leg_collision_CT(int leg);
    void compute_knee_position(Vec3<float>&p, Vec3<float>&q, int leg);

private: //help functions
    void draw_swing();
    void test_mpc_update();
    

protected:
    ImitationUserParameters userParameters;

public:
    // MPC
    deque<VecM<float, 24>> mpc_control_bag;
    double mpc_time;
    int iter;
    int iter_loco;
    int iter_between_mpc_update;
    int iter_between_mpc_node;
    int nsteps_between_mpc_update;
    int nsteps_between_mpc_node;

    VecM<float, 3> f_ff[4];
    bool firstRun;
   
    // Swing Control
    float statusTimes[4];
    float contactStatus[4];
    float swingState[4];
    float swingTimes[4];
    float swingTimesRemain[4];
    Vec4<float> stanceState;
    float stanceTimes[4];
    float stanceTimesRemain[4];
    bool firstSwing[4];
    bool firstStance[4];
    VecM<float, 3> pFoot_des[4];    //desired swing foot positions
    VecM<float, 3> vFoot_des[4];    //desired swing foot velocities
    VecM<float, 3> aFoot_des[4];    //desired swing foot accelerations
    VecM<float, 3> pf[4];           // foothold locations
    VecM<float, 3> pf_filtered[4];           // foothold locations
    VecM<float, 3> pFoot[4];        //actual swing foot positions
    VecM<float, 3> vFoot[4];        //actual swing foot velocities
    FootSwingTrajectory<float> footSwingTrajectories[4];
    Mat3<float> Kp_swing;
    Mat3<float> Kd_swing;
    Mat3<float> Kp_stance;
    Mat3<float> Kd_stance;
    VecM<float, 3> qJ_des[4];

    // foot location filter
    deque<VecM<float, 3>> pf_filter_buffer[4]; // buffer storing re-optimized foot locations
    int filter_window;

    // LCM
    lcm::LCM mpc_data_lcm;
    lcm::LCM rfmpc_data_lcm;
    lcm::LCM mpc_cmds_lcm;
    hkd_command_lcmt mpc_cmds;
    hkd_data_lcmt mpc_data;
    hkd_rfmpc_data_lcmt rfmpc_data;
    std::mutex mpc_cmd_mutex;
    std::mutex mpc_data_mutex;
    std::thread mpcLCMthread;

    // Tracking yaw angle
    int yaw_flip_plus_times;
    int yaw_flip_mins_times;
    float raw_yaw_pre;
    float raw_yaw_cur;
    float yaw;

private:
    VecM<float, 24> mpc_control;
    Vec3<float> init_joint_pos[4];
    bool in_standup;
    int iter_standup;
    int desired_command_mode;
};

#endif

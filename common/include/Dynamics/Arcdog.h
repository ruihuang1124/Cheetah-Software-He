/*! @file Arcdog.h
 *  @brief Utility function to build a Arcdog Quadruped object
 *
 * This file is based on MiniCheetahFullRotorModel_mex.m and builds a model
 * of the Arcdog robot.  The inertia parameters of all bodies are
 * determined from CAD.
 *
 */

#ifndef PROJECT_ARCDOG_H
#define PROJECT_ARCDOG_H

#include "FloatingBaseModel.h"
#include "Quadruped.h"

/*!
 * Generate a Quadruped model of Arcdog
 */
template <typename T>
Quadruped<T> buildArcdog() {
  Quadruped<T> arcdog;
  arcdog._robotType = RobotType::ARCDOG;
  //    mass parameters
  //  Nominal Total Mass = 25.693
  //  Real Total Mass = 22.94 (without Nvidia and outer covering and 2 heads)
  // arcdog._bodyMass = 6.165;//13.777;
  arcdog._bodyMass = 8.598;//with big 24v battery
  //  Other Parts Mass  11.916
  arcdog._abadMass = 0.805;
  arcdog._hipMass = 1.652;
  arcdog._kneeMass = 0.279;
  arcdog._rotorMass = 0.065;

  //    link parameters
  arcdog._bodyLength = 0.548; // distance between the hip forward and backward
  arcdog._bodyWidth = 0.121;  // = 0.237 - 0.058 * 2 (distance between the abad right and left)
  arcdog._bodyHeight = 0.116;

  arcdog._abadLinkLength = 0.0972; // The distance from the center of the abad to the center of the hip
  arcdog._hipLinkLength = 0.2445;  // length of thigh
  arcdog._kneeLinkLength = 0.2635; // length of calf (plus foot, length r=0.02)
  arcdog._maxLegLength = 0.5080; // = _hipLinkLength + _kneeLinkLength

  arcdog._kneeLinkY_offset = 0.000; // bais between the thigh and calf.
  arcdog._abadRotorLocationXOffset = -0.0790; // In X axis, distance from the AbAd joint to Center of AbAd rotor.
  arcdog._hipRotorLocationYOffset = 0.0000;  // In Y axis, distance from the AbAd joint to Center of Hip rotor.
  arcdog._kneeRotorLocationYOffset = 0.0427;// In Y axis, distance from the Hip joint to Center of Knee rotor.

  //    motor parameters
  arcdog._motorTauMax = 18.0f;//TODO, fix the value.
  arcdog._abadGearRatio = 9;
  arcdog._hipGearRatio = 9;
  arcdog._kneeGearRatio = 9;
  arcdog._batteryV = 36;
  arcdog._motorKT = 0.28;  // this is flux linkage(0.0025) * pole pairs(21)*1.5
  arcdog._motorR = 0.43;
  arcdog._jointDamping = 0.02;  // TODO
  arcdog._jointDryFriction = 0.4; // TODO

  // rotor inertia if the rotor is oriented so it spins around the z-axis TODO
  Mat3 <T> rotorRotationalInertiaZ;
  rotorRotationalInertiaZ << 53, 0, 0, 0, 53, 0, 0, 0, 120;
  rotorRotationalInertiaZ = 1e-6 * rotorRotationalInertiaZ;

  Mat3 <T> RY = coordinateRotation<T>(CoordinateAxis::Y, M_PI / 2);
  Mat3 <T> RX = coordinateRotation<T>(CoordinateAxis::X, M_PI / 2);
  Mat3 <T> rotorRotationalInertiaX = RY * rotorRotationalInertiaZ * RY.transpose();
  Mat3 <T> rotorRotationalInertiaY = RX * rotorRotationalInertiaZ * RX.transpose();

  // spatial inertias (base on left leg)
  Mat3 <T> abadRotationalInertia;
  abadRotationalInertia << 572.267, 8.780, 0.046, 8.780,1015.350, -0.015,0.046,-0.015,735.880;
  abadRotationalInertia = abadRotationalInertia * 1e-6;
  Vec3 <T> abadCOM(-0.006155,0.006930,0.000000);  // LEFT
  SpatialInertia<T> abadInertia(arcdog._abadMass, abadCOM, abadRotationalInertia);

  Mat3 <T> hipRotationalInertia;
  hipRotationalInertia << 10621.468, 190.063, 871.080, 190.063, 10208.742, 1376.795, 871.080, 1376.795, 2233.284;
  hipRotationalInertia = hipRotationalInertia * 1e-6;
  Vec3<T> hipCOM(-0.005377, 0.023858, -0.041813);
  SpatialInertia<T> hipInertia(arcdog._hipMass, hipCOM, hipRotationalInertia);

  Mat3<T> kneeRotationalInertia, kneeRotationalInertiaRotated;
  kneeRotationalInertiaRotated << 2485.797,0.000,-93.895,0.000,2512.613,0.000,-93.895,0.000,52.105;
  kneeRotationalInertiaRotated = kneeRotationalInertiaRotated * 1e-6;
  //    kneeRotationalInertia = RY * kneeRotationalInertiaRotated * RY.transpose();
  Vec3 <T> kneeCOM(0.005337,0.000,-0.117807);
  SpatialInertia<T> kneeInertia(arcdog._kneeMass, kneeCOM, kneeRotationalInertiaRotated);

  Vec3 <T> rotorCOM(0, 0, 0);
  SpatialInertia<T> rotorInertiaX(arcdog._rotorMass, rotorCOM, rotorRotationalInertiaX);
  SpatialInertia<T> rotorInertiaY(arcdog._rotorMass, rotorCOM, rotorRotationalInertiaY);

  // Mat3 <T> bodyRotationalInertia;
  // bodyRotationalInertia << 23743.442,-1.324,-165.915,-1.324,127019.601,12.502,-165.915,12.502,140284.756;// inertia from solidworks
  // bodyRotationalInertia = bodyRotationalInertia * 1e-6;
  // //    baseInertia << 0.0996, 0, 0, 0, 0.765, 0, 0, 0, 0.765;
  // Vec3 <T> bodyCOM(-0.000380,-0.000565, -0.007014);
  // Vec3<T> bodyDims(arcdog._bodyLength, arcdog._bodyWidth, arcdog._bodyHeight);
  // SpatialInertia<T> bodyInertia(arcdog._bodyMass, bodyCOM, bodyRotationalInertia); // simplified inertia


  // body with big 24v battery
    Mat3 <T> bodyRotationalInertia;
  bodyRotationalInertia << 31511.770, -2.833, -293.369, -2.833, 142974.387, 11.531, -293.369, 11.531, 158471.459;
  bodyRotationalInertia = bodyRotationalInertia * 1e-6;
  //    baseInertia << 0.0996, 0, 0, 0, 0.765, 0, 0, 0, 0.765;
  Vec3 <T> bodyCOM(-0.000471, -0.0007, -0.008687);
  Vec3<T> bodyDims(arcdog._bodyLength, arcdog._bodyWidth, arcdog._bodyHeight);
  SpatialInertia<T> bodyInertia(arcdog._bodyMass, bodyCOM, bodyRotationalInertia); // simplified inertia

  arcdog._abadInertia = abadInertia;
  arcdog._hipInertia = hipInertia;
  arcdog._kneeInertia = kneeInertia;
  arcdog._abadRotorInertia = rotorInertiaX;
  arcdog._hipRotorInertia = rotorInertiaY;
  arcdog._kneeRotorInertia = rotorInertiaY;
  arcdog._bodyInertia = bodyInertia;

  // locations; base on left leg
  arcdog._abadLocation = Vec3<T>(arcdog._bodyLength, arcdog._bodyWidth, 0) * 0.5;
  arcdog._abadRotorLocation = Vec3<T>(arcdog._bodyLength + arcdog._abadRotorLocationXOffset * 2, arcdog._bodyWidth, 0) * 0.5;
  arcdog._hipLocation = Vec3<T>(0, arcdog._abadLinkLength, 0);
  arcdog._hipRotorLocation = Vec3<T>(0, 0, 0);
  arcdog._kneeLocation = Vec3<T>(0, 0, -arcdog._hipLinkLength);
  arcdog._kneeRotorLocation = Vec3<T>(0, arcdog._kneeRotorLocationYOffset, 0);

  return arcdog;
}

#endif  // PROJECT_ARCDOG_H

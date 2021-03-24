#!/usr/bin/env python3
from cereal import car
from selfdrive.car.tesla.values import CAR
from selfdrive.car import STD_CARGO_KG, gen_empty_fingerprint, scale_rot_inertia, scale_tire_stiffness
from selfdrive.car.interfaces import CarInterfaceBase


class CarInterface(CarInterfaceBase):
  @staticmethod
  def compute_gb(accel, speed):
    # TODO: is this correct?
    return accel

  @staticmethod
  def get_params(candidate, fingerprint=gen_empty_fingerprint(), car_fw=None):
    ret = CarInterfaceBase.get_std_params(candidate, fingerprint)
    ret.carName = "tesla"
    ret.safetyModel = car.CarParams.SafetyModel.tesla
    ret.steerControlType = car.CarParams.SteerControlType.angle
    ret.enableCamera = True
    ret.openpilotLongitudinalControl = False
    ret.communityFeature = True

    ret.steerActuatorDelay = 0.1
    ret.steerRateCost = 0.5

    #safetyParam
    # BIT - MEANING
    # 0   - Has AP?
    # 1   - Has ACC?
    # 2   - OP Long Control?
    # 3   - HUD Integration?
    # 4   - Body Controls?

    if candidate == CAR.AP2_MODELS:
      ret.mass = 2100. + STD_CARGO_KG
      ret.wheelbase = 2.959
      ret.centerToFront = ret.wheelbase * 0.5
      ret.steerRatio = 13.5
      ret.safetyParam = 1 + 2 # has AP, ACC
      ret.openpilotLongitudinalControl = False
    elif candidate == CAR.AP1_MODELS:
      ret.mass = 2100. + STD_CARGO_KG
      ret.wheelbase = 2.959
      ret.centerToFront = ret.wheelbase * 0.5
      ret.steerRatio = 13.5
      ret.safetyParam = 1 + 2 # has AP, ACC
      ret.openpilotLongitudinalControl = True
      ret.communityFeature = True
    else:
      raise ValueError(f"Unsupported car: {candidate}")

    # set safetyParam flag for OP Long Control
    if ret.openpilotLongitudinalControl:
      ret.safetyParam = ret.safetyParam + 4

    ret.rotationalInertia = scale_rot_inertia(ret.mass, ret.wheelbase)
    ret.tireStiffnessFront, ret.tireStiffnessRear = scale_tire_stiffness(ret.mass, ret.wheelbase, ret.centerToFront)

    return ret

  def update(self, c, can_strings):
    self.cp.update_strings(can_strings)
    self.cp_cam.update_strings(can_strings)

    ret = self.CS.update(self.cp, self.cp_cam)
    ret.canValid = self.cp.can_valid and self.cp_cam.can_valid

    events = self.create_common_events(ret)

    ret.events = events.to_msg()
    self.CS.out = ret.as_reader()
    return self.CS.out

  def apply(self, c):
    can_sends = self.CC.update(c.enabled, self.CS, self.frame, c.actuators, c.cruiseControl.cancel)
    self.frame += 1
    return can_sends

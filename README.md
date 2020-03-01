# ComputationallyEfficient3DOTracking
Arduino C/C++ code for calculating 3DO using gyroscope measurements.
All computation is achieved using integer arithemtics only, lookup tables and linear interpolation for square root and in inverse square root, and small-angle approximation for cosine and sine functions estimations.
For each measurement step, the rotation angle and axis are obtained from the three angular velocity measurements using the Simultaneous Orthogonal Rotation Angle (SORA)
Both, the rotation matrix and rotation quaternion implementations are supported.

## Installation and Prerequisites

Runes on an Arduino board

## Usage

StepWise3DOCalculationRM.ino - calculates 3DO for each measurement step using rotation matrices
StepWise3DOCalculationRQ.ino - calculates 3DO for each measurement step using rotation quaternions
OnlyFinal3DOCalculationRQ2RM.ino - calculates only the final 3DO after several consecutive rotation measurements using rotation quaternions


## Author
Sara Stancin

University of Ljubljana
Faculty of electrical engineering
sara.stancin@fe.uni-lj.si

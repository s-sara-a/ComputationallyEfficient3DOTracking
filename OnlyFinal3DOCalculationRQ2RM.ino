// compute only final 3DO (resulting after several consecutive rotations)
// using Rotation quaternions
// author: Sara Stancin; e-mail: sara.stancin@fe.uni-lj.si

#include "configuration.h"

// declaring global variables
int omegaNorm, omegaNormInv, omegaSquare;
int omega[3], v[3];
int S[3][3], Snew[3][3], R[3][3], q[4], qnew[4], qtot[4], qsquare[3], q0q[3], qiq[2], gjq;
int s;

unsigned int timeus;
unsigned int timeustemp;

void setup() {
  while (!Serial);
  Serial.begin(115200);
  delay (250);
  Serial.println("Program started"), Serial.println("");

  fnFillSqrtLookUpTable(N4sqrt, &SqrtLookUpTable[0]);
  fnFillInvSqrtLookUpTable(N4sqrt, &InvSqrtLookUpTable[0]);

  omega[0] = 1000; omega[1] = 1500; omega[2] = 2000;

  // initial orientation
  S[0][0] = 1000000000, S[0][1] = 0, S[0][2] = 0;
  S[1][0] = 0, S[1][1] = 1000000000, S[1][2] = 0;
  S[2][0] = 0, S[2][1] = 0, S[2][2] = 1000000000;

  // initial rotation quaternion
  qtot[0] = fnCos(0);
  qtot[1] = 0;
  qtot[2] = 0;
  qtot[3] = 0;
}

int N = 10000;

void loop(){
  
  timeustemp = micros();
  
  for (int i = 0; i<N; i+=2){
    
    // first measurement step
    
    omegaSquare = fnSquareInts(omega[0]) + fnSquareInts(omega[1]) + fnSquareInts(omega[2]);
    omegaNorm = (fnSqrt4LI(omegaSquare)) >> 1; 
    omegaNormInv = fnInvSqrt4LI(omegaSquare);
    
    v[0] = fnBinaryRShift(fnMultplyInts(omega[0], omegaNormInv), NUMOFIBITSDIFF); 
    v[1] = fnBinaryRShift(fnMultplyInts(omega[1], omegaNormInv), NUMOFIBITSDIFF);
    v[2] = fnBinaryRShift(fnMultplyInts(omega[2], omegaNormInv), NUMOFIBITSDIFF); 

    s = fnSin(omegaNorm); 
    q[0] = fnCos(omegaNorm); 
    q[1] = fnMultiplyAndRShiftInts(v[0], s, NUMOFBITSVXVYVZ);
    q[2] = fnMultiplyAndRShiftInts(v[1], s, NUMOFBITSVXVYVZ);
    q[3] = fnMultiplyAndRShiftInts(v[2], s, NUMOFBITSVXVYVZ);    
  
    qnew[0] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], q[0]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[1], q[1]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[2], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[3], q[3]);
    qnew[1] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[2], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[3], q[2]);
    qnew[2] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[1], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[2], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[3], q[1]);
    qnew[3] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qtot[2], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(qtot[3], q[0]);

    omegaSquare = fnSquareInts(omega[0]) + fnSquareInts(omega[1]) + fnSquareInts(omega[2]);
    omegaNorm = (fnSqrt4LI(omegaSquare)) >> 1; 
    omegaNormInv = fnInvSqrt4LI(omegaSquare);

    // second measurement step
    
    v[0] = fnBinaryRShift(fnMultplyInts(omega[0], omegaNormInv), NUMOFIBITSDIFF); 
    v[1] = fnBinaryRShift(fnMultplyInts(omega[1], omegaNormInv), NUMOFIBITSDIFF);
    v[2] = fnBinaryRShift(fnMultplyInts(omega[2], omegaNormInv), NUMOFIBITSDIFF); 

    s = fnSin(omegaNorm); 
    q[0] = fnCos(omegaNorm); 
    q[1] = fnMultiplyAndRShiftInts(v[0], s, NUMOFBITSVXVYVZ);
    q[2] = fnMultiplyAndRShiftInts(v[1], s, NUMOFBITSVXVYVZ);
    q[3] = fnMultiplyAndRShiftInts(v[2], s, NUMOFBITSVXVYVZ);    
  
    qtot[0] = fnMultiplyAndRShiftIntsWithOverflow(qnew[0], q[0]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[1], q[1]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[2], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[3], q[3]);
    qtot[1] = fnMultiplyAndRShiftIntsWithOverflow(qnew[0], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[2], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[3], q[2]);
    qtot[2] = fnMultiplyAndRShiftIntsWithOverflow(qnew[0], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[1], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[2], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[3], q[1]);
    qtot[3] = fnMultiplyAndRShiftIntsWithOverflow(qnew[0], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(qnew[2], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(qnew[3], q[0]);
 
  }

  qsquare[0] = fnMultiplyAndRShiftIntsWithOverflow(qtot[1], qtot[1]) << 1; qsquare[1] = fnMultiplyAndRShiftIntsWithOverflow(qtot[2], qtot[2]) << 1; qsquare[2] = fnMultiplyAndRShiftIntsWithOverflow(qtot[3], qtot[3]) << 1;
  q0q[0] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], qtot[1]) << 1; q0q[1] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], qtot[2]) << 1; q0q[2] = fnMultiplyAndRShiftIntsWithOverflow(qtot[0], qtot[3]) << 1;
  qiq[0] = fnMultiplyAndRShiftIntsWithOverflow(qtot[1], qtot[2]) << 1; qiq[1] = fnMultiplyAndRShiftIntsWithOverflow(qtot[1], qtot[3]) << 1;
  gjq = fnMultiplyAndRShiftIntsWithOverflow(qtot[2], qtot[3]) << 1;      

  R[0][0] = ONETONUMOFBITSTRIG - qsquare[1] - qsquare[2]; R[0][1] = qiq[0] - q0q[2]; R[0][2] = q0q[1] + qiq[1]; 
  R[1][0] = qiq[0] + q0q[2]; R[1][1] = ONETONUMOFBITSTRIG - qsquare[0] - qsquare[2]; R[1][2] = gjq - q0q[0];
  R[2][0] = qiq[1] - q0q[1]; R[2][1] = q0q[0] + gjq; R[2][2] = ONETONUMOFBITSTRIG - qsquare[0] - qsquare[1];   

  Snew[0][0] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][0]);
  Snew[0][1] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][1]);
  Snew[0][2] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][2]);

  Snew[1][0] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][0]);
  Snew[1][1] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][1]);
  Snew[1][2] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][2]);

  Snew[2][0] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][0]);
  Snew[2][1] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][1]);
  Snew[2][2] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][2]);

 
  timeus = micros() - timeustemp;

  Serial.println();
  Serial.println("Final orientation: ");
  Serial.print("["), Serial.print(Snew[0][0]*pow(10,-9),12), Serial.print(", "), Serial.print(Snew[0][1]*pow(10,-9),12), Serial.print(", "), Serial.print(Snew[0][2]*pow(10,-9),12), Serial.println(";");
  Serial.print(Snew[1][0]*pow(10,-9),12), Serial.print(", "), Serial.print(Snew[1][1]*pow(10,-9),12), Serial.print(", "),Serial.print(Snew[1][2]*pow(10,-9),12), Serial.println(";");
  Serial.print(Snew[2][0]*pow(10,-9),12), Serial.print(", "), Serial.print(Snew[2][1]*pow(10,-9),12), Serial.print(", "),Serial.print(Snew[2][2]*pow(10,-9),12), Serial.println("];");

  Serial.print("Elapsed time: "), Serial.print(timeus); Serial.println(" mikros");
        
  delay(50000);
 
}

// compute stepwise 3DO (for each rotation measurement)
// using Rotation quaternions
// author: Sara Stancin; e-mail: sara.stancin@fe.uni-lj.si

#include "configuration.h"

int omegaNorm, omegaNormInv, omegaSquare;
int omega[3], v[3];
int SRQ[3][3], q[4], temp[3][4];
int s, c;

unsigned int timeus;
unsigned int timeustemp;

void setup() {
  while (!Serial);
  Serial.begin(115200);
  delay (250);
  Serial.println("Program started"), Serial.println("");

  fnFillSqrtLookUpTable(N4sqrt, &SqrtLookUpTable[0]);
  fnFillInvSqrtLookUpTable(N4sqrt, &InvSqrtLookUpTable[0]);

  omega[0] = 500; omega[1] = 0000; omega[2] = 0000;

  // initial orientation
  SRQ[0][0] = 1000000000, SRQ[0][1] = 0, SRQ[0][2] = 0;
  SRQ[1][0] = 0, SRQ[1][1] = 1000000000, SRQ[1][2] = 0;
  SRQ[2][0] = 0, SRQ[2][1] = 0, SRQ[2][2] = 1000000000;
}

int N = 10000;

void loop(){
  
  timeustemp = micros();
  
  for (int i = 0; i<N; i++){

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
  
    temp[0][0] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][0], q[1]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][2], q[3]);
    temp[0][1] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][0], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][1], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][2], q[2]);
    temp[0][2] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][0], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][2], q[1]);
    temp[0][3] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][0], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][1], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[0][2], q[0]);
    
    temp[1][0] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][0], q[1]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][2], q[3]);
    temp[1][1] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][0], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][1], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][2], q[2]);
    temp[1][2] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][0], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][2], q[1]);
    temp[1][3] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][0], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][1], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[1][2], q[0]);
    
    temp[2][0] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][0], q[1]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][2], q[3]);
    temp[2][1] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][0], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][1], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][2], q[2]);
    temp[2][2] = -fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][0], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][2], q[1]);
    temp[2][3] = fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][0], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][1], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(SRQ[2][2], q[0]);
    
    SRQ[0][0] = -fnMultiplyAndRShiftInts(temp[0][0], q[1], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][2], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(temp[0][3], q[2]);
    SRQ[0][1] = -fnMultiplyAndRShiftInts(temp[0][0], q[2], NUMOFBITSTRIG) - fnMultiplyAndRShiftIntsWithOverflow(temp[0][1], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][2], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][3], q[1]);
    SRQ[0][2] = -fnMultiplyAndRShiftInts(temp[0][0], q[3], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(temp[0][2], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(temp[0][3], q[0]);
    
    SRQ[1][0] = -fnMultiplyAndRShiftInts(temp[1][0], q[1], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][2], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(temp[1][3], q[2]);
    SRQ[1][1] = -fnMultiplyAndRShiftInts(temp[1][0], q[2], NUMOFBITSTRIG) - fnMultiplyAndRShiftIntsWithOverflow(temp[1][1], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][2], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][3], q[1]);
    SRQ[1][2] = -fnMultiplyAndRShiftInts(temp[1][0], q[3], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(temp[1][2], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(temp[1][3], q[0]);
    
    SRQ[2][0] = -fnMultiplyAndRShiftInts(temp[2][0], q[1], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][1], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][2], q[3]) - fnMultiplyAndRShiftIntsWithOverflow(temp[2][3], q[2]);
    SRQ[2][1] = -fnMultiplyAndRShiftInts(temp[2][0], q[2], NUMOFBITSTRIG) - fnMultiplyAndRShiftIntsWithOverflow(temp[2][1], q[3]) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][2], q[0]) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][3], q[1]);
    SRQ[2][2] = -fnMultiplyAndRShiftInts(temp[2][0], q[3], NUMOFBITSTRIG) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][1], q[2]) - fnMultiplyAndRShiftIntsWithOverflow(temp[2][2], q[1]) + fnMultiplyAndRShiftIntsWithOverflow(temp[2][3], q[0]);

  }
  
  timeus = micros() - timeustemp;

  Serial.println();
  Serial.println("Final orientation: ");
  Serial.print(SRQ[0][0]), Serial.print(", "), Serial.print(SRQ[0][1]), Serial.print(", "),Serial.println(SRQ[0][2]);
  Serial.print(SRQ[1][0]), Serial.print(", "), Serial.print(SRQ[1][1]), Serial.print(", "),Serial.println(SRQ[1][2]);
  Serial.print(SRQ[2][0]), Serial.print(", "), Serial.print(SRQ[2][1]), Serial.print(", "),Serial.println(SRQ[2][2]);
  
  Serial.print("Elapsed time: "), Serial.print(timeus); Serial.println(" mikros");
        
  delay(50000);
 
}

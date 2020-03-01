// compute stepwise 3DO (for each rotation measurement)
// using Rotation matrices
// author: Sara Stancin; e-mail: sara.stancin@fe.uni-lj.si

#include "configuration.h"

// declaring global variables
int omegaNorm, omegaNormInv, omegaSquare;
int v[3], vs[3], omega[3], vminc[3], vminc_m[3];
int S[3][3], Snew[3][3], R[3][3], Rtot[3][3];
int s, c, minc;

unsigned int timeus;
unsigned int timeustemp;

void setup() {
  while (!Serial);
  Serial.begin(115200);
  delay (250);
  Serial.println("Program started"), Serial.println("");

  fnFillSqrtLookUpTable(N4sqrt, &SqrtLookUpTable[0]);
  fnFillInvSqrtLookUpTable(N4sqrt, &InvSqrtLookUpTable[0]);

  omega[0] = 2000; omega[1] = 2000; omega[2] = 2000;
 
  // initial orientation
  S[0][0] = 1000000000, S[0][1] = 0, S[0][2] = 0;
  S[1][0] = 0, S[1][1] = 1000000000, S[1][2] = 0;
  S[2][0] = 0, S[2][1] = 0, S[2][2] = 1000000000;
}

int N = 10000;

void loop(){

  timeustemp = micros();
  
  for (int i = 0; i<N; i+=2){
    
    // first measurement step
    
    omegaSquare = fnSquareInts(omega[0]) + fnSquareInts(omega[1]) + fnSquareInts(omega[2]);
    omegaNorm = fnSqrt4LI(omegaSquare); 
    omegaNormInv = fnInvSqrt4LI(omegaSquare);
    
    v[0] = fnBinaryRShift(fnMultplyInts(omega[0], omegaNormInv), NUMOFIBITSDIFF);  
    v[1] = fnBinaryRShift(fnMultplyInts(omega[1], omegaNormInv), NUMOFIBITSDIFF);
    v[2] = fnBinaryRShift(fnMultplyInts(omega[2], omegaNormInv), NUMOFIBITSDIFF); 

    s = fnSin(omegaNorm), c = fnCos(omegaNorm), minc = (ONETONUMOFBITSTRIG - c); 

    vminc[0] = fnMultiplyAndRShiftInts(minc, v[0], NUMOFBITSVXVYVZ), vminc[1] = fnMultiplyAndRShiftInts(minc, v[1], NUMOFBITSVXVYVZ), vminc[2] = fnMultiplyAndRShiftInts(minc, v[2], NUMOFBITSVXVYVZ);   
  
    vs[0]  = fnMultiplyAndRShiftInts(s, v[0], NUMOFBITSVXVYVZ), vs[1]  = fnMultiplyAndRShiftInts(s, v[1], NUMOFBITSVXVYVZ), vs[2]  = fnMultiplyAndRShiftInts(s, v[2], NUMOFBITSVXVYVZ);   
    
    vminc_m[0] = fnMultiplyAndRShiftInts(vminc[0], v[1], NUMOFBITSVXVYVZ), vminc_m[1] = fnMultiplyAndRShiftInts(vminc[0], v[2], NUMOFBITSVXVYVZ), vminc_m[2] = fnMultiplyAndRShiftInts(vminc[1], v[2], NUMOFBITSVXVYVZ); 
 
    R[0][0] = (c + fnMultiplyAndRShiftInts(vminc[0], v[0], NUMOFBITSVXVYVZ)); R[0][1] = (vminc_m[0] - vs[2]); R[0][2] = (vminc_m[1] + vs[1]);
    R[1][0] = (vminc_m[0] + vs[2]); R[1][1] = (c + fnMultiplyAndRShiftInts(vminc[1], v[1], NUMOFBITSVXVYVZ)); R[1][2] = (vminc_m[2] - vs[0]);
    R[2][0] = (vminc_m[1] - vs[1]); R[2][1] = (vminc_m[2] + vs[0]); R[2][2] = (c + fnMultiplyAndRShiftInts(vminc[2], v[2], NUMOFBITSVXVYVZ));

    Snew[0][0] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][0]);
    Snew[0][1] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][1]);
    Snew[0][2] = fnMultiplyAndRShiftIntsWithOverflow(S[0][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[0][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[0][2], R[2][2]);

    Snew[1][0] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][0]);
    Snew[1][1] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][1]);
    Snew[1][2] = fnMultiplyAndRShiftIntsWithOverflow(S[1][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[1][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[1][2], R[2][2]);

    Snew[2][0] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][0]);
    Snew[2][1] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][1]);
    Snew[2][2] = fnMultiplyAndRShiftIntsWithOverflow(S[2][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(S[2][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(S[2][2], R[2][2]);

    // second measurement step
    
    omegaSquare = fnSquareInts(omega[0]) + fnSquareInts(omega[1]) + fnSquareInts(omega[2]);
    omegaNorm = fnSqrt4LI(omegaSquare); 
    omegaNormInv = fnInvSqrt4LI(omegaSquare);
    
    v[0] = fnBinaryRShift(fnMultplyInts(omega[0], omegaNormInv), NUMOFIBITSDIFF); 
    v[1] = fnBinaryRShift(fnMultplyInts(omega[1], omegaNormInv), NUMOFIBITSDIFF);
    v[2] = fnBinaryRShift(fnMultplyInts(omega[2], omegaNormInv), NUMOFIBITSDIFF);
  
    s = fnSin(omegaNorm), c = fnCos(omegaNorm), minc = (ONETONUMOFBITSTRIG - c); 

    vminc[0] = fnMultiplyAndRShiftInts(minc, v[0], NUMOFBITSVXVYVZ), vminc[1] = fnMultiplyAndRShiftInts(minc, v[1], NUMOFBITSVXVYVZ), vminc[2] = fnMultiplyAndRShiftInts(minc, v[2], NUMOFBITSVXVYVZ);   // prave vrednosti so enake izracunanim * 2^(-(numOfBitsTrig)) = 2^(-23)
    vs[0]  = fnMultiplyAndRShiftInts(s, v[0], NUMOFBITSVXVYVZ), vs[1]  = fnMultiplyAndRShiftInts(s, v[1], NUMOFBITSVXVYVZ), vs[2]  = fnMultiplyAndRShiftInts(s, v[2], NUMOFBITSVXVYVZ);  // prave vrednosti so enake izracunanim * 2^(-(numOfBitsTrig)) = 2^(-23)  
    vminc_m[0] = fnMultiplyAndRShiftInts(vminc[0], v[1], NUMOFBITSVXVYVZ), vminc_m[1] = fnMultiplyAndRShiftInts(vminc[0], v[2], NUMOFBITSVXVYVZ), vminc_m[2] = fnMultiplyAndRShiftInts(vminc[1], v[2], NUMOFBITSVXVYVZ); // prave vrednosti so enake izracunanim * 2^(-(numOfBitsTrig)) = 2^(-23)

    R[0][0] = (c + fnMultiplyAndRShiftInts(vminc[0], v[0], NUMOFBITSVXVYVZ)); R[0][1] = (vminc_m[0] - vs[2]); R[0][2] = (vminc_m[1] + vs[1]);
    R[1][0] = (vminc_m[0] + vs[2]); R[1][1] = (c + fnMultiplyAndRShiftInts(vminc[1], v[1], NUMOFBITSVXVYVZ)); R[1][2] = (vminc_m[2] - vs[0]);
    R[2][0] = (vminc_m[1] - vs[1]); R[2][1] = (vminc_m[2] + vs[0]); R[2][2] = (c + fnMultiplyAndRShiftInts(vminc[2], v[2], NUMOFBITSVXVYVZ));  

    
    S[0][0] = fnMultiplyAndRShiftIntsWithOverflow(Snew[0][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][2], R[2][0]);
    S[0][1] = fnMultiplyAndRShiftIntsWithOverflow(Snew[0][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][2], R[2][1]);
    S[0][2] = fnMultiplyAndRShiftIntsWithOverflow(Snew[0][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[0][2], R[2][2]);

    S[1][0] = fnMultiplyAndRShiftIntsWithOverflow(Snew[1][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][2], R[2][0]);
    S[1][1] = fnMultiplyAndRShiftIntsWithOverflow(Snew[1][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][2], R[2][1]);
    S[1][2] = fnMultiplyAndRShiftIntsWithOverflow(Snew[1][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[1][2], R[2][2]);

    S[2][0] = fnMultiplyAndRShiftIntsWithOverflow(Snew[2][0], R[0][0]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][1], R[1][0])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][2], R[2][0]);
    S[2][1] = fnMultiplyAndRShiftIntsWithOverflow(Snew[2][0], R[0][1]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][1], R[1][1])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][2], R[2][1]);
    S[2][2] = fnMultiplyAndRShiftIntsWithOverflow(Snew[2][0], R[0][2]) + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][1], R[1][2])  + fnMultiplyAndRShiftIntsWithOverflow(Snew[2][2], R[2][2]);

  }
  
  timeus = micros() - timeustemp;

  Serial.println();
  Serial.println("Final orientation: ");
  Serial.print("["), Serial.print(S[0][0]*pow(10,-9),12), Serial.print(", "), Serial.print(S[0][1]*pow(10,-9),12), Serial.print(", "), Serial.print(S[0][2]*pow(10,-9),12), Serial.println(";");
  Serial.print(S[1][0]*pow(10,-9),12), Serial.print(", "), Serial.print(S[1][1]*pow(10,-9),12), Serial.print(", "),Serial.print(S[1][2]*pow(10,-9),12), Serial.println(";");
  Serial.print(S[2][0]*pow(10,-9),12), Serial.print(", "), Serial.print(S[2][1]*pow(10,-9),12), Serial.print(", "),Serial.print(S[2][2]*pow(10,-9),12), Serial.println("];");
  
  Serial.print("Elapsed time: "), Serial.print(timeus); Serial.println(" mikros");
  Serial.print("SqrtLookUpTable[0]: "), Serial.println(SqrtLookUpTable[0]); 
  Serial.print("SqrtLookUpTable[100]: "), Serial.println(SqrtLookUpTable[100]);
  Serial.print("SqrtLookUpTable[1000]: "), Serial.println(SqrtLookUpTable[1000]); 
  Serial.print("SqrtLookUpTable[3434]: "), Serial.println(SqrtLookUpTable[3433]); 
  
        
  delay(500000);
 
}

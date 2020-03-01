#include <math.h>
#include <stdlib.h>

int fnMultplyInts(int number1, int number2){return number1*number2;}; 
int fnSquareInts(int number){return fnMultplyInts(number, number);};
int fnCalculateOmNorm(int omegax, int omegay, int omegaz){return fnSqrt4LI(fnSquareInts(omegax) + fnSquareInts(omegay) + fnSquareInts(omegaz));};
inline int fnBinaryRShift(int number, int bits){return (number >> bits) + ((number >> (bits - 1)) & 1);};   
int fnSqrt4LI(int number){
    int numberShift = number >> BITS4SQRT;
    int sqrt1 = SqrtLookUpTable[numberShift];
    int sqrt2 = SqrtLookUpTable[numberShift+1];
    return sqrt1 + ((sqrt2-sqrt1)*(number & TWOTOBITS4SQRTMIN1) >> BITS4SQRT); 
}
int fnInvSqrt4LI(int number){
    int numberShift = number >> BITS4SQRT;
    int sqrt1 = InvSqrtLookUpTable[numberShift];
    int sqrt2 = InvSqrtLookUpTable[numberShift+1];
    return sqrt1 + ((sqrt2-sqrt1)*(number & TWOTOBITS4SQRTMIN1) >> BITS4SQRT); 
}
int fnMultiplyAndRShiftInts(int number1, int number2, int bits){return fnBinaryRShift(fnMultplyInts(number1,number2),bits);};  
int fnMultiplyAndRShiftIntsWithOverflow(int number1, int number2){
  int number1a = number1 >> NUMOFINTBITSHALF;
  int number1b = number1 & TWOTONUMOFBITSHALFMIN1;  
  int number2a = number2  >> NUMOFBITSREMAINDER;
  int number2b = number2 & TWOTONUMBEROFBITSREMAINDERMIN1;  
  return (number1a*number2a) + ((number1a*number2b) >> NUMOFBITSREMAINDER) + ((number2a*number1b) >> NUMOFINTBITSHALF);
}
void fnFillSqrtLookUpTable(int N, int* resultTable){
  double resultInDouble;
  for (int i = 0; i < N; i++){
    resultInDouble = sqrt((i*pow(2,BITS4SQRT))*pow(10,-2*DECEX))*pow(10,DECEX);
    *(resultTable+i) = resultInDouble;
    if (resultInDouble - *(resultTable+i) > 0.5)
      (*(resultTable+i))++;
  }
}
void fnFillInvSqrtLookUpTable(int N, int* resultTable){
  double resultInDouble;
  for (int i = 0; i < N; i++){
    resultInDouble = pow(2,NUMOFBITSINVSQRT)/sqrt(i*pow(2,BITS4SQRT));
    *(resultTable+i) = pow(2,NUMOFBITSINVSQRT)/sqrt(i*pow(2,BITS4SQRT));
    if (resultInDouble - *(resultTable+i) > 0.5)
      (*(resultTable+i))++;
  }
}
int fnSin(int number){
    return fnBinaryRShift(number*phiDegRadTxNUMOFBITSCONSTSin, NUMOFBITSTRIGMINCONSTSin);
}
int fnCos(int number){
    int vmesniRezultat = fnBinaryRShift(number*phiDegRadTxNUMOFBITSCONSTSin, NUMOFBITSTRIGMINCONSTSin); 
    return ONETONUMOFBITSTRIG - fnBinaryRShift(fnBinaryRShift(vmesniRezultat, 4)*vmesniRezultat, 17); 
}


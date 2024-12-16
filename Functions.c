/*
 * Functions.c
 *
 *  Created on: 30 jun. 2023

 */

#include <stdint.h>
#include <math.h>
#include <intrin.h>
#include <stdbool.h>
#include <time.h>

/* Global variables */
uint32_t NumBits;

/*----------------------------------------------------------------------------*/
static inline uint32_t not(uint32_t a){
	return (a^1);
}
/*----------------------------------------------------------------------------*/
static inline uint32_t isNotZero(uint32_t a){
	return((a|-a)>>31);
}
/*----------------------------------------------------------------------------*/
static inline uint32_t isZero(uint32_t a){
	return(not((a|-a)>>31));
}
/*----------------------------------------------------------------------------*/
static inline uint32_t MaX (int32_t a, int32_t b) {

	return( (b&((a - b)>>31)) | (a&(((a - b)>>31)^0xFFFFFFFF)));
}
/*----------------------------------------------------------------------------*/
static inline uint32_t bitMux(uint32_t c, uint32_t a, uint32_t b) {
	return b ^ (c & (a ^ b));
}
/*----------------------------------------------------------------------------*/
static inline uint32_t isLessThan (int32_t a, int32_t b) {
uint32_t c = a - b;
	return bitMux(a^b, a, c) >> 31;
}
/*----------------------------------------------------------------------------*/

/* Returns the maximum number of bits needed by LongDivision algorithm */
uint32_t NumMaxBits(int Lmin, int T){

	uint32_t Nnum=32;
	uint32_t Value=(Lmin+T-1);
	for (int i=31; i>=0; i--){
		if ((Value & 0x80000000)==0x80000000)
			return(Nnum);
		Nnum--;
		Value=Value<<1;
	}
	return(Nnum);
}

/*----------------------------------------------------------------------------*/
/* Initial values for NumBits and "seed" used in random generation */
void Initialization (int Lmin, int T){

	NumBits=NumMaxBits(Lmin, T);
	srand(time(NULL));
	return;
}

/*----------------------------------------------------------------------------*/
/* Generates a random BitString with a number of bits equal to NumberOfBits*/
void GenerateRandomBitString(int NumberOfBits, uint32_t *BitString){

	int NRounds=ceil(NumberOfBits/32.0f);

	for (int j=0; j<NRounds; j++){
		BitString[j]=0;
		for (int k=0; k<4; k++){
			BitString[j]=BitString[j] ^ ((rand() % (256))<<(k<<3));
		}
	}
	return;
}

/*----------------------------------------------------------------------------*/
/* Read one bit from BitString at the position indicated by "Pointer" */
int Read (int Pointer, uint32_t *BitString){

	uint32_t Word;
	uint32_t Element;
	uint32_t Value;
	uint32_t Mask;
	uint32_t test;

	Word=Pointer>>5;
	Element=Pointer-(Word<<5);
	Value=BitString[Word];
	Mask=1<<Element;
	test=Mask&Value;

	return(isNotZero(test));
}

/*----------------------------------------------------------------------------*/
/* Write the integer "Value" into string "BitString" at the position indicated by "Pointer" */
void Write (int Pointer, uint32_t *BitString, int Value, int NumberOfbits){

	uint32_t Word;
	uint32_t Element;
	uint32_t MyValue;

	MyValue=Value;
	for (int j=NumberOfbits-1; j>=0; j--){
		Word=(Pointer+j)>>5;
		Element=(Pointer-(Word<<5))+j;
		BitString[Word]=(BitString[Word] & (0xffffffff ^ (1<<Element)) )*isZero(((MyValue)& 1))+(BitString[Word]| (1<<Element))*isNotZero(((MyValue)& 1));
		MyValue=MyValue>>1;
	}

	return;
}

/*----------------------------------------------------------------------------*/
/* Make the integer division of two integers a/b where a = quotient * b + rest*/
void LongDivision(int a, int b, int* rest, int* quotient){

	uint32_t LessThan;
	*rest=0; *quotient=0;

	for (int i=NumBits; i>=0; i--){
		*rest=*rest<<1;
		*rest=*rest|((a>>i) & 0x00000001);
		LessThan=isLessThan (*rest,b);
		*rest-=(b & (LessThan-1));
		*quotient=(*quotient | (not(LessThan)<<i));
	}

	return;
}

/*----------------------------------------------------------------------------*/
/* Proposed algorithm to find the value of d and log(d) */
void best_d(int t, int L, uint32_t *logb2, uint32_t *d){

	int rest,ceil_div;
	uint32_t Maxim;

	LongDivision(L+t-1,t,&rest,&ceil_div);
	Maxim=MaX(1,ceil_div)-1;
	*d=1<<Maxim;
	*logb2=Maxim;

	return ;
}

/*----------------------------------------------------------------------------*/
/* Converting a Constant Weight Word into a Binary string */
void CW2Bin(int N, int T, uint32_t *BitString, uint32_t *Tupple, int Lmin){

	int index,nn,tt;
	int Pointer;
	uint32_t u,d,Deltaindex;
	uint32_t *TuppleCopy;
	TuppleCopy=malloc(T*sizeof(uint32_t));

	nn=N;tt=T;
	index=0;
	Pointer=0;

	for (int j=0; j<T; j++) {
		TuppleCopy[j]=Tupple[j];
	}
	while(tt!=0 && nn>tt){
		best_d(tt,Lmin-Pointer,&u,&d);
		if(TuppleCopy[index]>=d){
			nn-=d; TuppleCopy[index]-=d;
			Write(Pointer,BitString,1,1);
			Pointer++;
		}
		else {
			Write(Pointer,BitString,0,1);
			Pointer++;
			/* encodefd */
			Deltaindex=TuppleCopy[index];
			if (Deltaindex<((1<<u)-d))
				u--;
			else
				Deltaindex=Deltaindex+(1<<u)-d;

			if (d>1){  // Write if not unary
				Write(Pointer,BitString,Deltaindex,u);
				Pointer+=u;
			}
			nn-=TuppleCopy[index]+1; tt--;index++;
		}
	}
	free(TuppleCopy);

	return;
}

/*----------------------------------------------------------------------------*/
/* Converting binary string into constant weight word*/
int Bin2CW (int N, int T, uint32_t *BitString, uint32_t *Tupple, int *Lambda,int Lmin){

	uint32_t delta,index,d;
	uint32_t NumberOfBitsRead,ReadValue,i;
	uint32_t n, t;
	uint32_t u,ud;
	uint32_t pt;
	uint32_t bitstop;
	uint32_t IsUZero;

	n=N;t=T;
	delta=0;index=0;NumberOfBitsRead=0;pt=0;bitstop=0,u=2;
	for (int j=0; j<T; j++){
		Tupple[j]=0;
	}
	for (int j=0; j<Lmin; j++){
		best_d(t,Lmin-NumberOfBitsRead,&ud,&d);
		ReadValue=Read(NumberOfBitsRead, BitString);
		NumberOfBitsRead++;
		pt++;
		u=(ud*((isZero(bitstop))))+(u-1)*isNotZero(bitstop);
		IsUZero=isZero(u);
		delta+=((1<<u)*ReadValue*isNotZero(bitstop))+d*ReadValue*isZero(bitstop);
		Tupple[index]=delta*((IsUZero)*(ReadValue*isNotZero(bitstop)+(isZero(ReadValue) ))+ReadValue*isZero(bitstop));
		delta-=delta*(IsUZero)*(ReadValue*isNotZero(bitstop) +(isZero(ReadValue)));
		index+=(IsUZero)*(ReadValue*isNotZero(bitstop)+isZero(ReadValue));
		t-=IsUZero*(ReadValue*isNotZero(bitstop)+isZero(ReadValue));
		bitstop=(isNotZero(ReadValue)*isZero(IsUZero)*isNotZero(bitstop)) + (isZero(ReadValue)*isZero(IsUZero));
		n--;
	}
	/*Finally, update Lambda (positions of 1's). If Lambda is not needed, this part of the code can be removed */
	for (int j=0; j<N; j++){
		Lambda[j]=0;
	}
	i=0;
	for (int j=0; j<T;j++){
		Lambda[Tupple[j]+i]=1;
		i+=(Tupple[j]+1);
	}
	return NumberOfBitsRead;
}

















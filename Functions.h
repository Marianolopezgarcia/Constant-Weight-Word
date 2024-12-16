/*
 * Functions.h
 *
 * Created on: 30 jun. 2023
 *
 */

#ifndef FUNCTIONS_H_
	#define FUNCTIONS_H_
	void Initialization (int Lmin, int T);
	void GenerateRandomBitString(int NumBits, uint32_t *BitString);
	int Bin2CW(int N, int T, uint32_t *BitString, uint32_t *Tupple, int *Lambda, int Lmin);
	void CW2Bin(int N, int T, uint32_t *BitString, uint32_t *Tupple, int Lmin);
#endif /* FUNCTIONS_H_ */

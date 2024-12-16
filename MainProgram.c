/*
 * Software related to the paper "Converting Fixed-Length Binary Strings
 * into Constant Weight Words: Application on Post-Quantum Cryptography"
 * Presented at IEEE Transactions on Dependable and Secure Computing, 2024
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Created on: 30 jun. 2023
 */

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <intrin.h>
#include <stdbool.h>
#include "Functions.h"

#define T 64  // Hamming weight
#define Lmin 365 // Length of the input binary string
#define N ((Lmin+2*T-T*ceil((float)Lmin/T))*pow(2,ceil((float)Lmin/T)-1)) // Length of the output binary string.

int main(void){

	uint32_t InputString[(int)(ceil(N/32.0f))]; // Input binary string represented by 32-bit words
	uint32_t OutputString[(int)(ceil(N/32.0f))]; // Used only for checking purpose
	uint32_t Tupple[T]; // t-Tupple delta
	int Lambda[(int)N]; // N-bit array representing the constant weight word. Used only for checking purpose: each bit is represented in one position.
	uint32_t ReadBits;  //Number of bits of InputString: it should be equal to Lmin.
	uint32_t Mask, NWords; //General variables

	if(Lmin<T){
		printf("Wrong Values for L, N or T (Lmin<T)");
		printf("\nValues: N=%d, Lmin=%d, T=%d\n",(int)N,(int)Lmin,(int)T);
		return 0;
	}
	printf("***********************************\n");
	printf("Values: N=%d, Lmin=%d, t=%d\n",(int)N,(int)Lmin,(int)T);
	printf("***********************************\n");

	/* Initialization function */
	Initialization (Lmin, T);

	for (int j=0; j<100; j++){ // Test 100 random InputStrings

		/* Generating a N-bit random string*/
		GenerateRandomBitString((int)N, InputString);

		/* Converting a binary input string into a constant weight word represented by Tupple or Lambda*/
		ReadBits=Bin2CW(N,T,InputString,Tupple,Lambda,Lmin);

		/* Inverse process: converting a constant weight word (Tupple) into a binary string */
		CW2Bin(N,T,OutputString,Tupple,Lmin);

		/* Checking results: if InputString[k]==OutputString[k] the result is correct */
		NWords=ceil(ReadBits/32.0f);
		printf("Round no.%d: ", j);
		for(int k=0; k<NWords;k++){
			Mask=0xffffffff;
			if (k==NWords-1 && (ReadBits%32)!=0){
				Mask=(1<<(ReadBits-(k<<5)))-1;
			}
			if((InputString[k] & Mask)!=(OutputString[k]& Mask) || ReadBits<Lmin){
				printf(" Error when checking if input and output data are identical");
				return 0;
			}
			else {
				printf ("InputString[%d]=%d-OutputString[%d]=%d; ", k,InputString[k]& Mask, k, OutputString[k]& Mask);
			}
		}
		printf("\nRound %d Successful checking\n", j);
	}
	return 0;
}







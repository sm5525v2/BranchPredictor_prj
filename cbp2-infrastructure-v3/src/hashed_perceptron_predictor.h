// my_predictor.h

#include <cmath>
#include <bitset>
#include <iostream>
#include <unistd.h>
// total Bits in history register 127
// total number of weight per perceptron 128
// total bits to store weight: 9 bits
//total numbe rof peceptrons tables size = 256
// total perceptrons 256 * 256 * 9 = 65 KiloBytes
// total memory budgets will be 73.7 KBytes + 256 Bytes + 256 Bits = 74.2KBytes
class perceptron_update : public branch_update {
public:
	unsigned int index;
};
using namespace std;
class perceptron : public branch_predictor {
public:
#define HISTORY_LENGTH	256 // the number of history 
#define TABLE_SIZE 256 // for each of table, the # of weights
#define SEGMENT_LENGTH 8
	perceptron_update u;
	branch_info bi;

	std::bitset<HISTORY_LENGTH> history; // the global history table which stores the last result bits
	unsigned int index[HISTORY_LENGTH]; // the index table which stores 
	signed long W[HISTORY_LENGTH][TABLE_SIZE]; //2D-List which hold each weight of each perceptron
	int lastout; // 
	unsigned int theta;

	bool confident;

	perceptron (void) {  
		memset (index, 0, sizeof (index));
		lastout = 0;
		theta = floor(1.93 * HISTORY_LENGTH + (HISTORY_LENGTH * 1.0) / 2);
		for(int a = 0; a< HISTORY_LENGTH;a++){ //initialize the Weight table with all value starting with one, gives a taken intuition
			for(int b = 0; b <=TABLE_SIZE; b++){
				W[a][b] = 1;
			}
		}
		history.set(); //set all the bits in the history table to be 1, indicating all previous history from start is taken 

	}

	branch_update *predict (branch_info & b) { 
		
		if(b.br_flags & BR_CONDITIONAL){
			
			index[0] = b.address % TABLE_SIZE;// calulating the index of the first table using only branch address history. 
			int out = W[0][index[0]];
			std::bitset<SEGMENT_LENGTH> History_segment; // craete a hisotry_segment bitset to store log(# weight) bits of history
		
			for (int j = 1; j < HISTORY_LENGTH; j++) {
				
				//assgin the value based on the whole history tables.
				for(int i = 0; i< SEGMENT_LENGTH;i++){
					if(history[i+j] == 1){ 
						History_segment.set(i,1);
					}
					else{
						History_segment.set(i,0);
					}
				}
				// change it from bitset to unsigned long for XOR calculation
				unsigned long history_seg = History_segment.to_ulong();

				// for each of the table, calculate the location of the weight from the table to use. 
				index[j] = (history_seg ^ b.address) % TABLE_SIZE;
				
				//weighted sum of the all the weights from the table
				out = out + W[j][index[j]];
				
			}
			//record the last time out data
			lastout = out;

			//check confidence
			confident = true;
			if(abs(lastout) <= theta) confident = false;
			
			//the conditional statement based on the literature to determine take or not
			//If the out is less than 0, it will be predicted as not taken.
			bool taken = 1;
			if (out < 0) {
				taken = 0;
			}
			u.direction_prediction(taken);
			
		}
		else {
			u.direction_prediction(1);
		}
		return &u;
	}

	void update (branch_update *u, bool taken, unsigned int target) {
		// update the history, if the prediction is not equal to result or the predicition is not reach the theta
		if (u->direction_prediction() != taken || abs(lastout) <= theta) {
			// update each of the weight table 
			for (int j = 0; j < HISTORY_LENGTH; j++) {
				//if the result is taken, we will increase the positive confidence of that weight.
				if (taken) {
					//the size of the weight is limited to be from -255 to 256, which is 9 bits 
					if(W[j][index[j]] < 256){
						W[j][index[j]] += 1;
					}
					else {
						W[j][index[j]]  = 256;
					}
				//
				} else {
					if(W[j][index[j]] > -255){
						W[j][index[j]]  -= 1;
					} else {
						W[j][index[j]]  = -255;
					}
						
				}
			}
		}
		//update the global history tables and store new result.
		history = history << 1;

		if(taken){
			history.set(0,1);

		}else {
			history.set(0,0);
		}
	}
};
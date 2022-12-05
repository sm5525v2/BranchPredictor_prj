// my_predictor.h
// This file contains a sample my_predictor class.
// It has a simple 32,768-entry gshare with a history length of 15 and a
// simple direct-mapped branch target buffer for indirect branch prediction. 
#include <cmath>
#include <bitset>

class perceptron_update : public branch_update {
public:
	unsigned int index;
};

class perceptron : public branch_predictor {
public:
#define HISTORY_LENGTH	128
#define TABLE_SIZE 512
	perceptron_update u;
	branch_info bi;

	std::bitset<HISTORY_LENGTH> history;
	unsigned int index[TABLE_SIZE];
	unsigned int W[TABLE_SIZE][HISTORY_LENGTH];
	int lastout;
	unsigned int theta;

	perceptron (void) {  
		memset (index, 0, sizeof (index));
		lastout = 0;
		theta = floor(1.93 * HISTORY_LENGTH + (HISTORY_LENGTH*1.0) / 2);
	}

	branch_update *predict (branch_info & b) { 

		if(b.br_flags & BR_CONDITIONAL){
			index[0] = b.address % TABLE_SIZE;
			int out = W[0][index[0]];
			for (int j = 1; j <= HISTORY_LENGTH; j++) {

				index[j] = (history[j - 1] ^ b.address) % TABLE_SIZE;
				out = out + W[j][index[j]];
			}
			//printf("%d\n", out);
			lastout = out;
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
		// update the history 
		if (u->direction_prediction() != taken || abs(lastout) <= theta) {
			for (int j = 0; j <= HISTORY_LENGTH; j++) {
				if (taken) {
					W[j][index[j]] = W[j][index[j]] + 1;
				} else {
					W[j][index[j]] = W[j][index[j]] - 1;
				}
			}
		}
		history = history << 1;
		if(taken){
			history.set(0,1);

		}else {
			history.set(0,0);
		}
	}
};

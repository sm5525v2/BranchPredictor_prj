// my_predictor.h
// This file contains a sample my_predictor class.
// It has a simple 32,768-entry gshare with a history length of 15 and a
// simple direct-mapped branch target buffer for indirect branch prediction.


#include "batage_predictor.h"
// #include "hashed_perceptron_predictor.h"
#include "chooser.h"

class my_update : public branch_update {
public:
    //add variable if need
};

class my_predictor : public branch_predictor{

 private:

  //batage
  batage batage_pred;
  histories hist;

  //hashed_perceptron
  //perceptron perceptron_pred

  //chooser
  chooser chsr;

  my_update u; //predict.cc use branch_update to get result
  long long int PC; //PC is need to predict and update batage, so store it in member variable and pass to pred.predict method

 public:
  my_predictor(void)
  {
	//initialize batage and perectron
    hist.printconfig();
    printf("total bits = %d\n",batage_pred.size()+hist.size());
  }

  branch_update* predict (branch_info & b) {
      PC = b.address;
    //   u.direction_prediction(batage_pred.predict(PC, hist));

	  u.direction_prediction(chsr.predict(batage_pred, PC, hist));

	  //bool pred_taken = chsr.predict(batage_pred, perceptron_pred, PC, hist);
	  //u.direction_prediction(pred_taken);
      return &u;
  }

  void update (branch_update *u, bool taken, unsigned int target) {
      //update batage
	  batage_pred.update(PC,taken,hist, false);

	  //update hashed_perceptron
	  //perceptron perceptron_pred.update()

      hist.update(target,taken);
  }
};
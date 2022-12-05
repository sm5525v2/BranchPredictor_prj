#include "batage_predictor.h"
// #include hashed_perceptron_predictor.h

class chooser {
public:
    bool choose(batage& batage_pred, perceptron& perceptron_pred, uint32_t pc, histories & p, branch_info& b) {
        // return batage_pred.predict(pc, p);
        bool taken = batage_pred.predict(pc, p);
        if(batage_pred.curConfLevel == 0) return taken;
        else return perceptron_pred.predict(b)->direction_prediction();
    }
};
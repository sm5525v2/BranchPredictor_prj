class chooser {
public:
    bool choose(batage& batage_pred, perceptron& perceptron_pred, uint32_t pc, histories & p, branch_info& b) {
        bool batageTaken = batage_pred.predict(pc, p);
        bool perceptronTaken = perceptron_pred.predict(b);
        if(batage_pred.curConfLevel != 0 && perceptron_pred.confident) return perceptronTaken;
        else return batageTaken;
    }
};
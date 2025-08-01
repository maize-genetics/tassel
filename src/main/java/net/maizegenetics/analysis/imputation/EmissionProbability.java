package net.maizegenetics.analysis.imputation;

public class EmissionProbability {
	double[][] probObsGivenState; //state in rows, observation in columns
	
	public EmissionProbability() {
		
	}
	
	public double getProbObsGivenState(int state, int obs) {
		return probObsGivenState[state][obs];
	}
	
	public double getLnProbObsGivenState(int state, int obs) {
	    double prob = getProbObsGivenState(state, obs);
	    if (prob == 0.0) prob = Double.MIN_VALUE;
		return Math.log(prob);
	}

	public double getProbObsGivenState(int state, int obs, int node) {
		return probObsGivenState[state][obs];
	}
	
	public double getLnProbObsGivenState(int state, int obs, int node) {
        double prob = getProbObsGivenState(state, obs, node);
        if (prob == 0.0) prob = Double.MIN_VALUE;
        return Math.log(prob);
	}

	public void setEmissionProbability(double[][] probabilityMatrix) {
		probObsGivenState = probabilityMatrix;
	}
}

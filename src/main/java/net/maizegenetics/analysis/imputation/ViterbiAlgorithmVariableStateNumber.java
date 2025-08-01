package net.maizegenetics.analysis.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class ViterbiAlgorithmVariableStateNumber {

	//adapted from Rabiner Proceedings of the IEEE 77(2):257-286
	//initialize
	//d(0, i) = p(Obs-0|S(i)*p(S(i))
	//where d0 = path length, Obs-0 = observation 0, S0 = true state 0
	//iterate:
	//for t = 1 to n
	//d(t,S(j)) = max(j){d(t-1,S(i)) * p[S(j)|S(i)] * p[Obs(t)|S(j)]}
	//h(t, j) = the value of i that maximizes distance
	//where h() = path history
	//termination
	//choose state that maximizes path length
	//back tracking
	//S(t) = h(t+1, S(t+1)), to decode best sequence
	
	TransitionProbability myTransitionMatrix;
	EmissionProbability probObservationGivenState;
	byte[] obs;
	ArrayList<int[]> history;
	double[] distance;
	double[] probTrueStates; //ln of probabilities
	int numberOfObs;
	byte[] finalState;
	int numberOfCurrentNodeStates;
	int numberOfPreviousNodeStates;
	
	public ViterbiAlgorithmVariableStateNumber(byte[] observations, TransitionProbability transitionMatrix, EmissionProbability obsGivenTrue, double[] pTrue) {
		obs = observations;
		numberOfObs = obs.length;
//		numberOfStates = transitionMatrix.getNumberOfStates();
		
		myTransitionMatrix = transitionMatrix;
		probObservationGivenState = obsGivenTrue;
		probTrueStates = new double[pTrue.length];
		for (int i = 0; i < pTrue.length; i++) {
			probTrueStates[i] = Math.log(pTrue[i]);
		}
		
		history = new ArrayList<>(numberOfObs);
//		distance = new double[numberOfStates];
	}
	
	public void calculate() {
		initialize();
		for (int i = 1; i < numberOfObs; i++) {
			updateDistanceAndHistory(i);
		}
	}
	
	public void initialize() {
		int n = probTrueStates.length;
		numberOfCurrentNodeStates = numberOfPreviousNodeStates = n;
		distance = new double[n];
		for (int i = 0; i < n; i++) {
			try{
				distance[i] = probObservationGivenState.getLnProbObsGivenState(i, obs[0], 0) + probTrueStates[i];
			} catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
	
	public void updateDistanceAndHistory(int node) {
		
		myTransitionMatrix.setNode(node);
		numberOfPreviousNodeStates = numberOfCurrentNodeStates;
		numberOfCurrentNodeStates = myTransitionMatrix.getNumberOfStates();
		double[][] candidateDistance = new double[numberOfPreviousNodeStates][numberOfCurrentNodeStates];
		
		int distanceLength = distance.length;
		
		try { //this try block for debugging
		for (int i = 0; i < numberOfPreviousNodeStates; i++) { //this is the number of nodes for the previous anchor
			for (int j = 0; j < numberOfCurrentNodeStates; j++) {
				candidateDistance[i][j] = distance[i] + myTransitionMatrix.getLnTransitionProbability(i, j) + probObservationGivenState.getLnProbObsGivenState(j, obs[node], node);
			}
		}
		} catch(Exception e) {
			System.out.println(String.format("at node %d, distance.length = %d, numberOfPreviousNodeStates = %d, and numberOfCurrentNodeStates = %d",
					node, distance.length, numberOfPreviousNodeStates, numberOfCurrentNodeStates));
			throw new RuntimeException(e);
		}
		
		//find the maxima
		int[] max = new int[numberOfCurrentNodeStates];
		for (int i = 0; i < numberOfPreviousNodeStates; i++) {
			for (int j = 0; j < numberOfCurrentNodeStates; j++) {
				if (candidateDistance[i][j] > candidateDistance[max[j]][j]) max[j] = i;
			}
		}

		//update distance and history
		distance = new double[numberOfCurrentNodeStates];
		int[] nodeHistory = new int[numberOfCurrentNodeStates];
		history.add(nodeHistory);
		for (int j = 0; j < numberOfCurrentNodeStates; j++) {
			distance[j] = candidateDistance[max[j]][j];
			nodeHistory[j] = (byte) max[j];
		}
		
		//debug
//		System.out.print("distance: ");
//		for (double dbl : distance) System.out.print(String.format("  %1.2e", dbl));
//		System.out.println();
//		System.out.println("candidate distance:");
//		for (int i = 0; i < numberOfPreviousNodeStates; i++) {
//			for (int j = 0; j < numberOfCurrentNodeStates; j++) {
//				System.out.printf("%1.1e  ", candidateDistance[i][j]);
//			}
//			System.out.println();
//		}
//		System.out.print("node history: ");
//		for (byte b : nodeHistory) System.out.print(String.format("  %d", b));
//		System.out.println();
//		if (node > 5) System.exit(0);
		//end debug
		
		//if the min distance is less than -1e100, subtract the min distance;
		double maxd = distance[0];
		double mind = 0;
		for (int i = 0; i < numberOfCurrentNodeStates; i++) {
			if (distance[i] > maxd) maxd = distance[i];
			if (distance[i] != Double.NEGATIVE_INFINITY && distance[i] < mind) mind = distance[i];
		}
		if (mind < -1e100) {
			for (int i = 0; i < numberOfCurrentNodeStates; i++) {
				distance[i] -= maxd;
			}
		}
	}

	/**
	 *
	 * @return	an int array representing the sequence of most likely states given the data
	 */
	public int[] getMostProbableIntegerStateSequence() {
		int[] seq = new int[numberOfObs];
		byte finalState = 0;
		
		for (int i = 1; i < distance.length; i++) {
			if (distance[i] > distance[finalState]) finalState = (byte) i;
		}
		
		//S(t) = h(t+1, S(t+1)), to decode best sequence
		seq[numberOfObs - 1] = finalState;
		for (int i = numberOfObs - 2; i >= 0; i--) {
			 seq[i] = history.get(i)[seq[i + 1]];
		}
		return seq;
	}

	/**
	 *
	 * @return	a byte array representing the sequence of most likely states given the data
	 *
	 * @deprecated Using this method when the number of states at a position > 127 will yield
	 * in unpredictable results. Replace with {@link #getMostProbableIntegerStateSequence()}
	 */
	@Deprecated
	public byte[] getMostProbableStateSequence() {
		byte[] seq = new byte[numberOfObs];
		byte finalState = 0;

		for (int i = 1; i < distance.length; i++) {
			if (distance[i] > distance[finalState]) finalState = (byte) i;
		}

		//S(t) = h(t+1, S(t+1)), to decode best sequence
		seq[numberOfObs - 1] = finalState;
		for (int i = numberOfObs - 2; i >= 0; i--) {
			seq[i] = (byte) history.get(i)[seq[i + 1]];
		}
		return seq;
	}

	public void setStateProbability(double[] probTrueState) {
		int n = probTrueState.length;
		probTrueStates = new double[n];
		for (int i = 0; i < n; i++) {
			probTrueStates[i] = Math.log(probTrueState[i]);
		}
	}
}

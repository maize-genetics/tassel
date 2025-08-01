package net.maizegenetics.analysis.imputation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BackwardForwardVariableStateNumber {
	/**
	 * @author pbradbury
	 * The BackwardForward algorithm is an HMM that is used to estimate the probability of every possible state at each position in a sequence.
	 * This implementation is based on the description in LR Rabiner (1986) A tutorial on hidden Markov models and select applications in speech recognition. Proceedings of the IEEE 77(2):257-286.
	 * To use this class, first supply observations, positions, a TransitionProbability object, an EmissionProbability object, and the initial state probabilities. The calculate alpha and beta. Typical usage will be:
	 *    BackwardForwardAlgorithm myBackwardForward = new BackwardForwardAlgorithm()
	 *       .observations(obs)
	 *       .positions(pos)
	 *       .transition(transition)
	 *       .emission(emission)
	 *       .initialStateProbability(probs)
	 *       .calculateAlpha()
	 *       .calculateBeta()
	 *       
	 * Either use gamma() to retrieve a List<double[]> of values or use writeGamma(String outputFile) to save the values to an output file.
	 */
		private static final Logger myLogger = LogManager.getLogger(BackwardForwardAlgorithm.class);
		
		private int[] myObservations;
		private int[] myPositions;
		private TransitionProbability myTransitions;
		private EmissionProbability myEmissions;
		private double[] initialStateProbability;
		private List<double[]> alpha;
		private List<double[]> beta;
		
		/**
		 * The BackwardForward algorithm is used to calculate the probability of each state at each position.
		 * 
		 */
		public BackwardForwardVariableStateNumber() {
			
		}
		
		public BackwardForwardVariableStateNumber calculateAlpha() {
			myTransitions.setNode(0);
			int nStates = myTransitions.getNumberOfStates();
			int nObs = myObservations.length;
			alpha = new LinkedList<>();
			
			//1. initialize: alpha[1](i) = p[i]b[i](O[1]), i = state i
			double[] aPrior = new double[nStates];
			for (int s = 0; s < nStates; s++) 
				aPrior[s] = initialStateProbability[s] * myEmissions.getProbObsGivenState(s, myObservations[0], 0);
			alpha.add(aPrior);

			//2. induction: alpha[t+1](j) = {sum[i=1 to N] alpha[t](i)a[ij]} b[j](O[t+1])
			//alpha[t+1](j) is the value of alpha at node t+1, state j
			//alpha[t](i) is the value of alpha at the previous node for state i
			//a[ij] is the transition probability from node t, state i to node t+1, state j
			//b[j](O[t+1]) is the emission probability for node t+1, state j given the set of observations for node t+1
			for (int t = 1; t < nObs; t++) { //this t is the t+1 in the formula, aPrior = alpha[t]
				myTransitions.setNode(t);
				int nPreviousStates = nStates;
				nStates = myTransitions.getNumberOfStates(); //number of states for current node
				double[] aT = new double[nStates]; //this is alpha[t+1]
				for (int j = 0; j < nStates; j++) {
					double sumTrans = 0;
					for (int i = 0; i < nPreviousStates; i++) sumTrans += aPrior[i] * myTransitions.getTransitionProbability(i, j);
					aT[j] = sumTrans * myEmissions.getProbObsGivenState(j, myObservations[t], t);
				}
				
				aT = multiplyArrayByConstantIfSmall(aT);
				alpha.add(aT);
				aPrior = aT;
			}
			
			return this;
		}
		
		private double[] multiplyArrayByConstantIfSmall(double[] dblArray) {
			double minval = Arrays.stream(dblArray).min().getAsDouble();
			if (minval < 1e-50) {
				double maxval = Arrays.stream(dblArray).max().getAsDouble();
				if (maxval < 1e-25) return Arrays.stream(dblArray).map(d -> d*1e25).toArray();
			}
			
			return dblArray;
		}
		
		public BackwardForwardVariableStateNumber calculateBeta() {
			int nObs = myObservations.length;
			LinkedList<double[]> betaTemp = new LinkedList<>();
			myTransitions.setNode(nObs - 1);
			int nStates = myTransitions.getNumberOfStates();
			
			//initialization: beta[T](i) = 1
			double[] bNext = new double[nStates];
			Arrays.fill(bNext, 1.0);
			betaTemp.add(bNext);
			
			//induction: beta[t](i) = sum(j=1 to N): a[i][j]*b[j](O[t+1])*beta[t+1](j)
			for (int t = nObs - 2; t >= 0; t--) {
				myTransitions.setNode(t); //in order to get the number of states
				int nNextStates = nStates;
				nStates = myTransitions.getNumberOfStates();
				double[] bT = new double[nStates];
				myTransitions.setNode(t+1); //because want transitions from current State to next State not previous to current
				for (int i = 0; i < nStates; i++) {
					double sumStates = 0;
					for (int j = 0; j < nNextStates; j++) {
						sumStates += myTransitions.getTransitionProbability(i, j) * myEmissions.getProbObsGivenState(j, myObservations[t + 1], t + 1) * bNext[j];
					}
						
					bT[i] = sumStates;
				}
				bT = multiplyArrayByConstantIfSmall(bT);
				betaTemp.addFirst(bT);
				bNext = bT;
			}
			beta = betaTemp;
			
			return this;
		}
		
		public List<double[]> gamma() {
			List<double[]> gamma = new ArrayList<>();
			Iterator<double[]> itAlpha = alpha.iterator();
			Iterator<double[]> itBeta = beta.iterator();
			
			//gamma[t](i) = P(q[t] = S[i] | O,model)
			//gamma[t](i) = alpha[t](i)*beta[t](i) / {sum(j=1 to N): alpha[t](j)*beta[t](j)}
			while(itAlpha.hasNext()) {
				double[] alphaArray = itAlpha.next();
				double[] betaArray = itBeta.next();
				int n = alphaArray.length;
				double[] gammaArray = new double[n];
				for (int i = 0; i < n; i++) gammaArray[i] = alphaArray[i] * betaArray[i];
				
				double divisor = Arrays.stream(gammaArray).sum();
				for (int i = 0; i < n; i++) gammaArray[i] /= divisor;
				gamma.add(gammaArray);
			}
			
			return gamma;
		}
		
		public BackwardForwardVariableStateNumber emission(EmissionProbability emission) {
			myEmissions = emission;
			return this;
		}
		
		public BackwardForwardVariableStateNumber transition(TransitionProbability transition) {
			myTransitions = transition;
			return this;
		}

		public  BackwardForwardVariableStateNumber observations(int[] observations) {
			myObservations = observations;
			return this;
		}
		
		public  BackwardForwardVariableStateNumber positions(int[] positions) {
			myPositions = positions;
			return this;
		}
		
		public BackwardForwardVariableStateNumber initialStateProbability(double[] probs) {
			initialStateProbability = probs;
			return this;
		}
		
		public List<double[]> alpha() {return alpha;}
		
		public List<double[]> beta() {return beta;}

}

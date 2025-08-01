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

import net.maizegenetics.analysis.imputation.EmissionProbability;
import net.maizegenetics.analysis.imputation.TransitionProbability;

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
 * The either use gamma() to retrieve a List<double[]> of values or use writeGamma(String outputFile) to save the values to an output file.
 */
public class BackwardForwardAlgorithm {
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
	public BackwardForwardAlgorithm() {
		
	}
	
	public BackwardForwardAlgorithm calculateAlpha() {
		int nStates = myTransitions.getNumberOfStates();
		int nObs = myObservations.length;
		alpha = new LinkedList<>();
		
		//1. initialize: alpha[1](i) = p[i]b[i](O[1]), i = state i
		double[] aPrior = new double[nStates];
		for (int s = 0; s < nStates; s++) 
			aPrior[s] = initialStateProbability[s] * myEmissions.getProbObsGivenState(s, myObservations[0], 0);
		alpha.add(aPrior);

		//2. induction: alpha[t+1](j) = {sum[i=1 to N] alpha[t](i)a[ij]} b[j](O[t+1])
		for (int t = 1; t < nObs; t++) { //this t is the t+1 in the formula, aPrior = alpha[t]
			double[] aT = new double[nStates]; //this is alpha[t+1]
			myTransitions.setNode(t);
			for (int j = 0; j < nStates; j++) {
				double sumTrans = 0;
				for (int i = 0; i < nStates; i++) sumTrans += aPrior[i] * myTransitions.getTransitionProbability(i, j);
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
	
	public BackwardForwardAlgorithm calculateBeta() {
		LinkedList<double[]> betaTemp = new LinkedList<>();
		int nStates = myTransitions.getNumberOfStates();
		int nObs = myObservations.length;
		
		//initialization: beta[T](i) = 1
		double[] bNext = new double[nStates];
		Arrays.fill(bNext, 1.0);
		betaTemp.add(bNext);
		
		//induction: beta[t](i) = sum(j=1 to N): a[i][j]*b[j](O[t+1])*beta[t+1](j)
		for (int t = nObs - 2; t >= 0; t--) {
			double[] bT = new double[nStates];
			myTransitions.setNode(t+1);
			for (int i = 0; i < nStates; i++) {
				double sumStates = 0;
				for (int j = 0; j < nStates; j++) {
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
	
	private void printSite(int pos, double[] values) {
		System.out.print(pos + ": ");
		Arrays.stream(values).mapToObj(d -> String.format("%1.4f ", d)).forEach(System.out::print);
		System.out.println();
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
	
	public void writeGamma(String outputFile, String formatString) {
		
		Iterator<double[]> itAlpha = alpha.iterator();
		Iterator<double[]> itBeta = beta.iterator();
		int counter = 0;
		
		try(BufferedWriter bw = Files.newBufferedWriter(Paths.get(outputFile))) {
			while(itAlpha.hasNext()) {
				double[] alphaArray = itAlpha.next();
				double[] betaArray = itBeta.next();
				int n = alphaArray.length;
				double[] gammaArray = new double[n];
				for (int i = 0; i < n; i++) gammaArray[i] = alphaArray[i] * betaArray[i];
				
				double divisor = Arrays.stream(gammaArray).sum();
				double[] normalizedGamma = Arrays.stream(gammaArray).map(g -> g/divisor).toArray();
				
				bw.write(myPositions[counter] + "\t");
				bw.write(Arrays.stream(normalizedGamma)
						.mapToObj(dbl -> String.format(formatString, dbl))
						.collect(Collectors.joining("\t", "", "\n")));
				
				counter++;
			}

		} catch(IOException ioe) {
			throw new RuntimeException("Unable to write " + outputFile, ioe);
		}
	}
	
	public void writeGamma(String outputFile) {
		writeGamma(outputFile, "%1.2e");
	}
	
	public BackwardForwardAlgorithm emission(EmissionProbability emission) {
		myEmissions = emission;
		return this;
	}
	
	public BackwardForwardAlgorithm transition(TransitionProbability transition) {
		myTransitions = transition;
		return this;
	}

	public  BackwardForwardAlgorithm observations(int[] observations) {
		myObservations = observations;
		return this;
	}
	
	public  BackwardForwardAlgorithm positions(int[] positions) {
		myPositions = positions;
		return this;
	}
	
	public BackwardForwardAlgorithm initialStateProbability(double[] probs) {
		initialStateProbability = probs;
		return this;
	}
	
	public List<double[]> alpha() {return alpha;}
	
	public List<double[]> beta() {return beta;}
}

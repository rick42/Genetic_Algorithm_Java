import java.util.*;

public class GeneticAlgorithm {
	
	private Map<Double[],Double> population;
	private int chromosomeLength;
	private int populationSize;
    private double percentChild = 0.6;
    private double percentMutate = 0.1;
    
    public GeneticAlgorithm(int popSize, int chrmLen)
    {
    	// TODO: Complete this constructor

    	populationSize = popSize;
    	chromosomeLength = chrmLen;
    }

    public void initializePopulation()
    {
    	// TODO: Write this method
    }

    public void nextGeneration()
    {
    	// TODO: Write this method
    }

    public void train(int maxGenerations)
    {
    	// TODO: Write this method
    }

    private double[] crossOver(double[] parent1, double[] parent2)
    {
    	// TODO: Write this method

    	return new double[0];
    }

    private double[] mutate(double[] chromosome)
    {
    	// TODO: Write this method

    	return new double[0];
    }
}
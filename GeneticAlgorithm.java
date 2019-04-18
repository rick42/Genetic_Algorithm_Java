import java.util.Map;

public class GeneticAlgorithm {

	///// ATTRIBUTES /////

    private Map<Double[],Double> population;
    private int chromosomeLength;
    private int populationSize;
    private double percentChild = 0.6;
    private double percentMutate = 0.1;
    
    ///// CONSTRUCTOR /////

    public GeneticAlgorithm(int popSize, int chrmLen)
    {
    	// TODO: Complete this constructor

    	populationSize = popSize;
    	chromosomeLength = chrmLen;
    }

    ///// METHODS /////

    /** 
    * initializes the private variable population with random values
    */
    public void initializePopulation()
    {
    	// TODO: Write this method
    }

    /**
    * nextGeneration method updates the population table by a single generation 
    */
    private void nextGeneration()
    {
    	// TODO: Write this method
    }

    /*
    * train method performs a number of generations specified by parameter 
    * maxGenerations
    */
    public void train(int maxGenerations)
    {
    	// TODO: Write this method
    }

    /*
    * crossOver method performs cross over on parameters parent1 and parent2 and 
    * producing a single child chromosome.  There is a random chance that mutation
    * will occur in the child specified by attribute percentMutate.  If the child
    * by crossOver already exists in the population, perform crossover again,
    * else return the child chromosome.
    */
    private double[] crossOver(double[] parent1, double[] parent2)
    {
    	// TODO: Write this method

    	return new double[0];
    }

    /**
    * mutate randomly chooses a gene in chromosome and chooses a new random value
    * for that gene.  The resulting chromosome is returned.
    */
    private double[] mutate(double[] chromosome)
    {
    	// TODO: Write this method

    	return new double[0];
    }
}
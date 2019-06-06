import java.util.ArrayList;
import java.util.Random;
import java.lang.Math;

public class GeneticAlgorithm {

    ///// ATTRIBUTES /////

    private ArrayList<Individual> population;
    private double[] fitnesses;
    private int chromosomeLength;
    private int populationSize;
    private double percentChild = 0.6;
    private double percentMutate = 0.1;
    private double geneRange = 200.0;
    private double geneOffset = -100.0;
    
    ///// CONSTRUCTOR /////

    public GeneticAlgorithm(int popSize, int chrmLen)
    {
    	// TODO: Complete this constructor

    	populationSize = popSize;
    	chromosomeLength = chrmLen;
        fitnesses = new double[populationSize];
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

    /**
    * train method performs a number of generations specified by parameter 
    * maxGenerations
    */
    public void train(int maxGenerations)
    {
    	// TODO: Write this method
    	
    }

    /**
    * crossOver method performs cross over on parameters parent1 and parent2 and 
    * producing a two child chromosomes.  There is a random chance that mutation
    * will occur in the children specified by attribute percentMutate. Returns 
    * two children.
    */
    public double[][] crossOver(double[] parent1, double[] parent2)
    {
        Random rand = new Random();
        
        int crossOverPoint;   //point where child will split parent 1 and parent 2
        double[][] children = new double [2][chromosomeLength]; // two dimentional array for two children
    	if (chromosomeLength >2)
    	{
            crossOverPoint = rand.nextInt(chromosomeLength-2)+1;
        }
        else
        {
            crossOverPoint = 1;
        }
    	
    	for (int i=0; i < chromosomeLength; i++)
    	{
            if (i < crossOverPoint)          // if coin lands on 0, choose parent 1
            {
                children[0][i] = parent1[i];
                children[1][i] = parent2[i];
            }
            else if (i >= crossOverPoint)     // if coin lands on 1, choose parent 2
            {
                children[0][i] = parent2[i];
                children[1][i] = parent1[i];
            }
            else{}
        }
        
        
        if (rand.nextDouble() <= percentMutate)
        {
            mutate(children[0]);
        }
        else{}
        if (rand.nextDouble() <= percentMutate)
        {
            mutate(children[1]);
        }
        else{}

    	return children;
    }

    /**
    * mutate randomly chooses a gene in chromosome and chooses a new random value
    * for that gene.  The resulting chromosome is returned.
    */
    private double[] mutate(double[] chromosome)
    {
    	Random rand = new Random();
    	int randomGene = rand.nextInt(chromosomeLength);
    	chromosome[randomGene] = (rand.nextDouble() * geneRange) + geneOffset;

    	return chromosome;
    }
}

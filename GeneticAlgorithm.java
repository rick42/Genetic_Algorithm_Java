import java.util.ArrayList;
import java.lang.Math;
import java.lang.Comparable;
import java.util.Random;

public class GeneticAlgorithm {

    ///// ATTRIBUTES /////

    public ArrayList<Individual> population;
    private String[] geneDataType;
    private int chromosomeLength;
    private int populationSize;
    private double percentChild;
    private double [] geneRange;
    private double [] geneOffset;
    private double percentMutate;
    private FitnessFunction fitnessFunction;
    private String fitnessPriority;
    Random rand = new Random(System.nanoTime());
    private long seed = System.currentTimeMillis();


    ///// CONSTRUCTOR /////

    public GeneticAlgorithm(int popSize, int chrmLen, String [] chrmDT, FitnessFunction fitFunc, double [] minGeneValue, double [] maxGeneValue, double chldPrcnt, double mttPrcnt, String fitPrio)
    {
    	// TODO: Complete this constructor
    	
    	population = new ArrayList<Individual>();
    	populationSize = popSize;
    	chromosomeLength = chrmLen;
        fitnessFunction = fitFunc;
        percentMutate = mttPrcnt;
        percentChild = chldPrcnt;
        geneDataType = new String [chromosomeLength];
        geneRange = new double [chromosomeLength];
        geneOffset = new double [chromosomeLength];
        
        for (int i = 0; i < chromosomeLength; i++)
        {
            geneDataType[i] = chrmDT[i];
            geneRange[i] = (maxGeneValue[i] - minGeneValue[i]);
            geneOffset[i] = minGeneValue[i];
        
        }
       // geneRange = (maxGeneValue - minGeneValue);
       // geneOffset = (minGeneValue);
        fitnessPriority = fitPrio;
        initializePopulation();

    }

    ///// METHODS /////

    /** 
    * initializes the private variable population with random values
    */
    private void initializePopulation()
    {
    	// TODO: Write this method
    	double[] tempChrome;
    	double[] tempChrome2 = new double [chromosomeLength];
    	int check= 1;
    	System.out.println("\nFitness:\n");
    	for (int i = 0; i < populationSize; i++)
    	{
            tempChrome = new double [chromosomeLength];
            for (int j = 0; j < chromosomeLength; j++)
            {   
                seed=seed+getInt(1001);
                tempChrome[j]= getGene(j);
                seed=seed+getInt(1002);
                //System.out.printf("%2.10f\n ",tempChrome[j]);

            }
                
            if (i == 0)
            {
                population.add(new Individual(tempChrome,fitnessFunction));
                //System.out.printf("%2.10f\n ",tempChrome[0]);

                //System.out.println(population.get(i).getFitness());
            }
            else
            {
                check = (findChromosome(population, tempChrome));
                //System.out.printf("%d <--------\n",check);
            }
            
            
            if ((i > 0) && (check == 0)) 
            {
                //tempChrome2 = population.get(i-1).getChromosome();
                //System.out.printf("\nthis is i-1 %2.10f\n ",tempChrome2[0]);
                population.add(new Individual(tempChrome,fitnessFunction));
                //tempChrome = population.get(i).getChromosome();
                //System.out.printf("this is i %2.10f\n ",tempChrome[0]);
                //System.out.println(population.get(i).getFitness()[0]);


            }
            else if (i != 0)    
            {   
                System.out.println("bad");
                i--;
            }
    	}
    	for (int k = 0; k < populationSize; k++)
    	{
            population.get(k).printChromosome();
    	
    	}
    	sortPopulation();
    }

    /**
    * nextGeneration method updates the population table by a single generation 
    */
    private void nextGeneration()
    {
    	// TODO: Write this method
        double [][] parents;
        double [][] children;
        int newGeneration = (int)(percentChild * populationSize);
        int maxPop = newGeneration + populationSize;
        
        while (population.size() < maxPop)
        {
            parents = new double [2][chromosomeLength];
            children = new double [2][chromosomeLength];
            parents = createMatingPool();
            children = crossOver(parents[0] ,parents[1]);
            
            if ( (population.size() < maxPop) && (findChromosome(population,children[0]) == 0 ) )
            {
                population.add(new Individual(children[0], fitnessFunction));
            }
            if ( (population.size() < maxPop) && (findChromosome(population,children[1]) == 0 ) )
            {
                population.add(new Individual(children[1], fitnessFunction));
            }
        }
        sortPopulation();
        
        for (int i = maxPop-1; i >= populationSize; i--)
        {
            population.remove(i);
        }
    }

    /**
    * train method performs a number of generations specified by parameter 
    * maxGenerations
    */
    public void train(int maxGenerations)
    {
    	// TODO: Write this method
    	for ( int i = 0; i < maxGenerations; i++)
    	{
            nextGeneration();
            //System.out.printf("this is Generation %d\n ",i+1);
        }
    	
    }
    
    /**
    * selects chromosomes to create the next generation
    *
    */
    
    private double[][] createMatingPool()
    {
    	
        final int poolSize = population.size();
        int tickets = (poolSize * (poolSize + 1)) / 2;
        int [] indexPool = new int [tickets];
        int k = 0;
        for (int i = 0; i < poolSize; i++)
        {   
            for (int j = i; j < poolSize; j++)
            {
                indexPool[k] = i;
                //System.out.println(indexPool[k]);
                k++;
            }
            
        }
        
        double [][] parents = new double [2][chromosomeLength];
        int lottery1,lottery2;
        lottery1 = getInt(tickets);
        seed=seed+1003;
        
        do{
            lottery2 = getInt(tickets);
            seed=seed+1004;
        }while (indexPool[lottery1] == indexPool[lottery2] );
        
        parents[0] = population.get(indexPool[lottery1]).getChromosome();
        parents[1] = population.get(indexPool[lottery2]).getChromosome();
        
        
        return parents;
    }
    /**
    * crossOver method performs cross over on parameters parent1 and parent2 and 
    * producing a two child chromosomes.  There is a random chance that mutation
    * will occur in the children specified by attribute percentMutate. Returns 
    * two children.
    */
    private double[][] crossOver(double[] parent1, double[] parent2)
    {
        
        
        int crossOverPoint;   //point where child will split parent 1 and parent 2
        double[][] children = new double [2][chromosomeLength]; // two dimentional array for two children
        double mutationBonusMultiplier = 5; 
        double mutationCheck = percentMutate;

    	if (chromosomeLength >2)
    	{
             crossOverPoint = getInt(chromosomeLength-2)+1;
             seed=seed+getInt(1005);
        }
        else
        {
            crossOverPoint = 1;
            mutationCheck = (mutationBonusMultiplier * percentMutate);
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
        
        
        if (getDouble() <= mutationCheck)
        {
            mutate(children[0]);
            seed+=getInt(1009);
        }
        else{}
        if (getDouble() <= mutationCheck)
        {
            mutate(children[1]);
            seed+=getInt(1010);
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
        rand.setSeed(seed);
        int randomGene = rand.nextInt(chromosomeLength);
    	chromosome[randomGene] = (rand.nextDouble() * geneRange[randomGene]) + geneOffset[randomGene];

    	return chromosome;
    }
    
    /**
    * find if chromosome already exist in population, true if chromosome exsist, false if it does not.
    */
    private int findChromosome( ArrayList<Individual> arr, double [] chromosome)
    {
        //TODO: Write this method
        double [] popChromosome = new double [chromosomeLength];
        int currentPopulationSize = arr.size();
        int check =  0;

        for (int i = 0; i < currentPopulationSize; i++)
        {
            popChromosome = (arr.get(i).getChromosome());
            
            for (int j = 0; j < chromosomeLength; j++)
            {
                if (popChromosome[j] != chromosome[j])
                {
                    j = chromosomeLength;       //breaks out of j loop
                }
                else if ((j == chromosomeLength-1) && (popChromosome[j] == chromosome[j]) && (j != 0))    
                {
                    check = 1;
                    //.out.println("bad1");
                    ///System.out.printf("%2.10f %2.10f",popChromosome[j],chromosome[j]);

                    i = currentPopulationSize;  //breaks out of i loop
                    
                }
            }
        
        }
        
        return check;
    }
    
    /** 
    * sort population according to fitness priority to maximize, or minimize
    * fitness value.
    */
    private void sortPopulation()
    {

        if (fitnessPriority == "minimize")
        {
            sortHeap(population);
        }
        else if (fitnessPriority == "maximize")
        {
            reverseSortHeap(population);
        }
        else 
        {
            System.out.println("\nThis is bad");
        }
    }
    
    /**
    *the following are functions to sort by heap.
    *
    */
    private void heapify(ArrayList<Individual> arr, int i, int end){
        
        int listEnd = arr.size();
        int r= (2*i)+2;
        int l= (2*i)+1;
        double leftChild, rightChild, parent;
        if (l < listEnd)
        {
            leftChild = arr.get(l).getFitness();
        }
        else
        {
            leftChild = -3000;
        }
        if (r < listEnd)
        {
            rightChild = arr.get(r).getFitness();
        }
        else 
        {
            rightChild = -3000;
        }
        parent = arr.get(i).getFitness();

        if ( (l <= end) && (r <= end) && (parent < leftChild)  && (leftChild >= rightChild))
        {

            swap(arr, i, l);
        }
        else if ( (r <= end) && parent < rightChild)
        {
            swap(arr, i, r);
        }
        else if ( (r> end) && (parent < leftChild) )
        {

        //System.out.println("\nno right child\n");

            swap(arr, i, l);

        }
        //System.out.println(parent);
        //System.out.printf("\n%2.2f %2.2f\n",leftChild,rightChild);

    }

    private void swap(ArrayList<Individual> arr, int i, int j)
    {
        //final int chromosomeEnd = arr.get(i).getChromosome().length;
        //double [] tempChrome = new double[chromosomeEnd];
        //double tempFitness;
        
        //tempChrome = arr.get(i).getChromosome();
        ArrayList<Individual> temp = new ArrayList<Individual>();
        temp.add(arr.get(i));
        arr.set(i, arr.get(j));
        arr.set(j, temp.get(0));

    }
    
    private void buildHeap(ArrayList<Individual> arr)
    {
        final int listEnd = arr.size();
        int i;

        for (int e = 0; e < listEnd; e++)
        {
            i = e;

            while ( i>=0 && (arr.get(i).getFitness() > arr.get( ( ( (i-1) /2) ) ).getFitness() ) )
            {
                i = ((i-1)/2);
                heapify(arr,i,e);
            }
        }
        
    }
    
    private void sortHeap(ArrayList<Individual> arr ) 
    {   
        buildHeap(arr);
        final int listEnd = arr.size();
        int i;
        int r;
        int l;
        
        
        for (int e = listEnd-1; e >0 ;)
        {
            swap (arr, e, 0);
            i=0;
            r=2;
            l=1;
            e--;
            while ( (i <= e) && ( (l <= e) || (r <= e) ) && ((arr.get(i).getFitness() < arr.get(l).getFitness()) || (arr.get(i).getFitness() < arr.get(r).getFitness()) ))
            {
                if ( (r <= e) && ( arr.get(l).getFitness() < arr.get(r).getFitness() ) )
                {
                    heapify(arr,i,e);
                    i=r;
                    r= (2*i)+2;
                    l= (2*i)+1;
                }
                else
                {
                    heapify(arr,i,e);
                    i=l;
                    r= (2*i)+2;
                    l= (2*i)+1;
                
                }
            
            }
            
        
        }
        
    }
    private void reverseSortHeap(ArrayList<Individual> arr)
    {   
        buildHeap(arr);
        sortHeap(arr);
        int i = 0;
        int j = arr.size()-1;
        
        while (i < j)
        {   
            swap(arr, i, j);
            i++;
            j--;
            
        }
    }
    
    private double getDouble()
    {
        double duece;
        rand.setSeed(System.currentTimeMillis());

        duece = rand.nextDouble();
        
        return duece;
    }
    
    private double getGene(int i)
    {
        double gene;
        rand.setSeed(seed);

        gene = (rand.nextDouble() * geneRange[i]) + geneOffset[i];
        
        return gene;
    }
    private int getInt(int range)
    {
        rand.setSeed(seed);

        int rando = rand.nextInt(range);
        
        return rando;
    }

}

public class Individual{
    private double[] chromosome;
    private double fitness;

    public Individual(double[] chrom, double fit)
    {
        chromosome = chrom;
        fitness = fit;
    }

    public double[] getChromosome(){return chromosome;}
    public double getFitness(){return fitness;}

    public static int cmpFitness(Individual ind1, Individual ind2)
    {
        if(ind1.getFitness() < ind2.getFitness())
        {
            return -1;
        }
        else if (ind1.getFitness() > ind2.getFitness())
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}
// Use the following interface to match with a Lambda expression.
// That Lambda expression should take in a chromosome and returns
// the chromosome's fitness.

@FunctionalInterface
public interface FitnessFunction {

	public double calculateFitness(double[] chromosome);

}

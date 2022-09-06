package toShare;
public interface PairPot {

	/**computes the pair energy
	 * @param r double, the pair distance
	 * @return double, the pair energy */
	public double pairE(double r);

	public double getRange();

	/**computes the next guess
 * 	@param boxArr Box[]
 * 	@param t double, tempearature
 * 	@param skTarget double[][], the target structure factor */
	public void nextGuess(Box[] boxArr, double t, double[][] skTarget);

	/**to String
 * 	@return a String describing the potential */
	public String toString();
}

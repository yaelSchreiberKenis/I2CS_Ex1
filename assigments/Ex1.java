package assigments;

import static java.lang.Math.*;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}

    /**
     * This function returns the reduces polynomial - it removes
     * all the redundant zeroes that are located in the highest degrees
     * @param p The polynomial
     * @return Reduces polynomial
     */
    public static double [] ReducePolynomial(double [] p) {
        int zeroes = 0;
        for (int i = p.length - 1; i > 0; i--){
            if (p[i] == 0) {
                zeroes++;
            }
            else {
                break;
            }
        }

        if (zeroes == 0) {
            return p;
        }

        double [] reduced = new double[p.length - zeroes];
        for (int i = 0; i < p.length - zeroes; i++){
            reduced[i] = p[i];
        }

        return reduced;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p1(x) -p2(x)| < eps,
     * assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     *
     * @param p1 - first polynomial function
     * @param p2 - second polynomial function
     * @param x1 - minimal value of the range
     * @param x2 - maximal value of the range
     * @param eps - epsilon (positive small value, like 10^-3 or 10^-6)
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans;
        // Decide which polynomial has fewer coefficients
        double[] minarr = p2;
        double[] maxarr = p1;

        // Array to store differences between polynomials
        double[] diffArr = new double[max(p1.length, p2.length)];

        // If p1 is smaller, switch minarr and maxarr
        if (p1.length < p2.length) {
            minarr = p1;
            maxarr = p2;
        }

        // Subtract coefficients that exist in both polynomials
        for (int i = 0; i < minarr.length; i++) {
            diffArr[i] = maxarr[i] - minarr[i];
        }

        // Copy remaining coefficients from longer polynomial
        for (int i = minarr.length; i < maxarr.length; i++) {
            diffArr[i] = maxarr[i];
        }

        // Find root of the difference polynomial in range [x1, x2]
        ans = root_rec(diffArr, x1, x2, eps);

        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n inner sample points and computing the segment-path between them.
     * assuming x1 < x2.
     * This function should be implemented iteratively (none recursive).
     *
     * @param p - the polynomial function
     * @param x1 - minimal value of the range
     * @param x2 - maximal value of the range
     * @param numberOfSegments - (A positive integer value (1,2,...))
     * @return the length approximation of the function between f(x1) and f(x2)
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        // Compute width of each segment along x-axis
        double segmentLength = (x2 - x1) / numberOfSegments;

        // Variable to store total length
        double sum = 0;

        // Loop over each segment
        for (double start = x1; start < x2; start += segmentLength) {
            // Compute y-value at start and end of segment
            double fp1 = f(p, start);
            double fp2 = f(p, start + segmentLength);

            // Compute distance using Pythagoras: sqrt(dx^2 + dy^2)
            double x2y2 = (segmentLength * segmentLength) + (fp1 - fp2) * (fp1 - fp2);

            // Add segment length to total sum
            sum += sqrt(x2y2);
        }

        return sum;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids
     * between the functions (number of samples on each polynomial).
     * This function computes an approximation of the area between the polynomial functions within the x-range.
     * The area is computed using Riemann's like integral.
     *
     * @param p1 - first polynomial function
     * @param p2 - second polynomial function
     * @param x1 - minimal value of the range
     * @param x2 - maximal value of the range
     * @param numberOfTrapezoid - number of trapezoids to approximate area
     * @return the approximated area between the two polynomial functions within the [x1,x2] range
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double segmentLength = (x2 - x1) / numberOfTrapezoid;
        double sum = 0;

        // Loop over each segment in x
        for (double left = x1; left < x2; left += segmentLength) {
            double right = left + segmentLength;

            // Compute height difference at left and right
            double h1 = f(p1, left) - f(p2, left);
            double h2 = f(p1, right) - f(p2, right);

            if (h1 * h2 >= 0) {
                // If heights have same sign, compute trapezoid area normally
                double trapeze = ((abs(h1) + abs(h2)) / 2.0) * segmentLength;
                sum += abs(trapeze);
            } else {
                // If polynomials cross, find x where they are equal
                double equalityPoint = sameValue(p1, p2, left, right, EPS);

                // Compute area of right triangle
                double rightTriangle = (right - equalityPoint) * abs(h2) / 2;
                sum += rightTriangle;

                // Compute area of left triangle
                double leftTriangle = (equalityPoint - left) * abs(h1) / 2;
                sum += leftTriangle;
            }
        }

        return sum;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Given a polynomial as a double array, getPolynomFromString(poly(p)) should
     * return an array equal to p.
     *
     * @param p - a String representing polynomial function.
     * @return array of doubles representing polynomial coefficients
     */
    public static double[] getPolynomFromString(String p) {
        // Replace '-' with '+-' for future splitting by '+'
        p = p.replace("-", "+-");

        // Remove all spaces
        p = p.replace(" ", "");

        // In this case for helping us to interpret -x as -1*x or x as 1*x
        p = p.replace("-x", "-1x");
        p = p.replace("+x", "1x");
        if (p.charAt(0) == 'x') {
            p = "1" + p;
        }

        // Remove leading '+' if exists
        if (p.charAt(0) == '+') {
            p = p.substring(1);
        }

        // Split the string into individual terms
        String[] parts = p.split("\\+");

        // Array size is the highest degree + 1
        int arrSize = getDegreeOfCoefficient(parts[0]) + 1;
        double[] poly = new double[arrSize];

        // Loop over terms
        for (int i = 0; i < parts.length; i++) {
            int degree = getDegreeOfCoefficient(parts[i]);
            double coefficient;

            // If term contains 'x', extract coefficient before 'x'
            if (parts[i].contains("x")) {
                coefficient = Double.parseDouble(parts[i].substring(0, parts[i].indexOf("x")));
            } else {
                coefficient = Double.parseDouble(parts[i]);
            }

            poly[degree] = coefficient;
        }

        return poly;
    }

    /**
     * Computes degree of a single term in a polynomial string
     * Examples: "3x^2" -> 2, "4x" -> 1, "5" -> 0
     */
    public static int getDegreeOfCoefficient(String p) {
        if (!p.contains("x")) return 0;   // No x means degree 0
        if (!p.contains("^")) return 1;   // x without ^ means degree 1
        String degreeString = p.substring(p.indexOf("^") + 1);
        return Integer.parseInt(degreeString);
    }

    /**
     * Adds two polynomials p1 and p2
     */
    public static double[] add(double[] p1, double[] p2) {
        // Determine which polynomial is longer
        double[] maxarr = p1;
        double[] minarr = p2;
        double[] result = new double[max(p1.length, p2.length)];

        if (p1.length < p2.length) {
            minarr = p1;
            maxarr = p2;
        }

        // Add coefficients where both polynomials have terms
        for (int i = 0; i < minarr.length; i++) {
            result[i] = maxarr[i] + minarr[i];
        }

        // Copy remaining coefficients from longer polynomial
        for (int i = minarr.length; i < maxarr.length; i++) {
            result[i] = maxarr[i];
        }

        return result;
    }

    /**
     * Multiplies two polynomials p1 and p2
     */
    public static double[] mul(double[] p1, double[] p2) {
        double[] ans = new double[p1.length + p2.length - 1];

        // Multiply each term of p1 by each term of p2
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[j + i] += p1[i] * p2[j];
            }
        }

        return ReducePolynomial(ans);
    }

    /**
     * Computes derivative of polynomial po
     */
    public static double[] derivative(double[] po) {
        // Degree 0 polynomial derivative is 0
        if (po.length < 2) return new double[]{0};

        double[] ans = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = po[i] * i;
        }

        return ans;
    }

    /**
     * Computes a polynomial from 2 or 3 points
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        // Check that points are valid (2 or 3 points)
        int lengthOfx = xx.length;
        int lengthOfy = yy.length;
        if (!(xx != null && yy != null && lengthOfx == lengthOfy && lengthOfx > 1 && lengthOfx < 4)) {
            return null;
        }

        if (lengthOfx == 2) {
            // Avoid deviding in null
            if (xx[1] == xx[0]) {
                return null;
            }
            // Linear polynomial: y = a + bx
            double b = (yy[1] - yy[0]) / (xx[1] - xx[0]);
            double a = yy[0] - b * xx[0];
            return new double[]{a, b};
        }

        // Quadratic polynomial for 3 points
        double x1 = xx[0], x2 = xx[1], x3 = xx[2];
        double y1 = yy[0], y2 = yy[1], y3 = yy[2];

        double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
        // Avoid deviding in null
        if (denom == 0) {
            return null;
        }

        double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
        double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
        double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

        return new double[]{c, b, a};
    }

    /**
     * Two polynomial functions are equal if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2), up to EPS.
     */
    public static boolean equals(double[] p1, double[] p2) {
        int maxLength = max(p1.length, p2.length);

        for (double i = 0; i < maxLength; i++) {
            // If values differ more than EPS, polynomials are not equal
            if (abs(f(p1, i) - f(p2, i)) > EPS) return false;
        }

        return true;
    }

    /**
     * Converts polynomial array into readable string
     * Example: {2,0,3.1,-1.2} -> "-1.2x^3 +3.1x^2 +2.0"
     */
    public static String poly(double[] poly) {
        StringBuilder ans = new StringBuilder();
        if (poly.length == 0) return "0";

        for (int i = poly.length - 1; i >= 0; i--) {
            if (poly[i] != 0) {
                if (poly[i] > 0 && i != poly.length - 1) {
                    ans.append(" +");
                }
                else {
                    ans.append(" ");
                }
                ans.append(poly[i]);
                if (i == 1) {
                    ans.append("x");
                }
                else if (i > 1) {
                    ans.append("x^").append(i);
                }
            }
        }

        if (ans.length() == 0) return "0";

        // Remove first space
        return ans.toString().substring(1);
    }
}

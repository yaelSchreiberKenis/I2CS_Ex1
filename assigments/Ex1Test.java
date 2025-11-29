package assigments;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
    static double[] xPoints0 = {2, 4};
    static double[] xPoints1 = {1, 7, -0.3};
    static double[] yPoints0 = {-1, 5.5};
    static double[] yPoints1 = {6, 3, 0};

    // -------------------- Additional tests --------------------
    @Test
    void testSameValueAdditional() {
        double[] p1 = {1, -2, 1}; // (x-1)^2
        double[] p2 = {0, 1};     // x
        double x = Ex1.sameValue(p1, p2, 1, 3, Ex1.EPS);
        assertEquals(2.618, x, 0.01); // approximate root of (x-1)^2 - x = 0
        x = Ex1.sameValue(p1, p2, 0, 1, Ex1.EPS);
        assertEquals(0.382, x, 0.01); // another approximate root of (x-1)^2 - x = 0

        double[] p3 = {0, 1};     // x
        double[] p4 = {0, -1};    // -x
        double y = Ex1.sameValue(p3, p4, -1, 1, Ex1.EPS);
        assertEquals(0, y, Ex1.EPS);

        double[] p5 = {0, 0, 1};  // x^2
        double[] p6 = {1};        // 1
        double z = Ex1.sameValue(p5, p6, 0, 2, Ex1.EPS);
        assertEquals(1, z, 0.01);
    }

    @Test
    void testLengthAdditional() {
        double[] p1 = {0, 1}; // y = x
        double len1 = Ex1.length(p1, 0, 3, 3);
        assertEquals(3 * Math.sqrt(2), len1, 0.001);

        double[] p2 = {0, 0, 1}; // y = x^2
        double len2 = Ex1.length(p2, 0, 1, 10);
        assertTrue(len2 > 1);

        double[] p3 = {5}; // y = 5
        double len3 = Ex1.length(p3, 0, 2, 4);
        assertEquals(2, len3, Ex1.EPS);
    }

    @Test
    void testAreaAdditional() {
        double[] p1 = {0, 1}; // y = x
        double[] p2 = {0};    // y = 0
        double a1 = Ex1.area(p1, p2, 0, 2, 100);
        assertEquals(2.0, a1, 0.01);

        double[] p3 = {0, 0, 1}; // y = x^2
        double[] p4 = {0, 0, 0.5}; // y = 0.5 x^2
        double a2 = Ex1.area(p3, p4, 0, 2, 50);
        assertEquals(4.0 / 3.0, a2, 0.05);

        double[] p5 = {0}; // y = 0
        double[] p6 = {1}; // y = 1
        double a3 = Ex1.area(p5, p6, 0, 5, 10);
        assertEquals(5, a3, 0.01);
    }

    @Test
    void testGetPolynomFromStringAdditional() {
        String s1 = "2x^3 -3x^2 +4x -5";
        double[] p1 = Ex1.getPolynomFromString(s1);
        assertArrayEquals(new double[]{-5, 4, -3, 2}, p1, 0);

        String s2 = "-x^2 + 2x - 3";
        double[] p2 = Ex1.getPolynomFromString(s2);
        assertArrayEquals(new double[]{-3, 2, -1}, p2, 0);

        String s3 = "5";
        double[] p3 = Ex1.getPolynomFromString(s3);
        assertArrayEquals(new double[]{5}, p3, 0);
    }

    @Test
    void testAddAdditional() {
        double[] p1 = {1, 2};
        double[] p2 = {3, 4};
        double[] sum1 = Ex1.add(p1, p2);
        assertArrayEquals(new double[]{4, 6}, sum1, 0);

        double[] p3 = {0, 1, 2};
        double[] p4 = {1};
        double[] sum2 = Ex1.add(p3, p4);
        assertArrayEquals(new double[]{1, 1, 2}, sum2, 0);

        double[] p5 = Ex1.ZERO;
        double[] sum3 = Ex1.add(p5, p4);
        assertArrayEquals(p4, sum3, 0);
    }

    @Test
    void testMulAdditional() {
        double[] p1 = {1, 1}; // 1 + x
        double[] p2 = {1, -1}; // 1 - x
        double[] m1 = Ex1.mul(p1, p2);
        assertArrayEquals(new double[]{1, 0, -1}, m1, 0);

        double[] p3 = {0, 1}; // x
        double[] p4 = {0, 0, 1}; // x^2
        double[] m2 = Ex1.mul(p3, p4);
        assertArrayEquals(new double[]{0, 0, 0, 1}, m2, 0);
    }

    @Test
    void testDerivativeAdditional() {
        double[] p1 = {1, -3, 2}; // 2x^2 -3x +1
        double[] d1 = Ex1.derivative(p1);
        assertArrayEquals(new double[]{-3, 4}, d1, 0);

        double[] p2 = {5};
        double[] d2 = Ex1.derivative(p2);
        assertArrayEquals(new double[]{0}, d2, 0);

        double[] p3 = {0, 0, 1}; // x^2
        double[] d3 = Ex1.derivative(p3);
        assertArrayEquals(new double[]{0, 2}, d3, 0);
    }

    @Test
    void testEqualsAdditional() {
        double[] p1 = {1, 2};
        double[] p2 = {1, 2, 0};
        assertTrue(Ex1.equals(p1, p2));

        double[] p3 = {1, 2};
        double[] p4 = {1, 2.002};
        assertFalse(Ex1.equals(p3, p4));
    }

    @Test
    void testPolyAdditional() {
        double[] p1 = {1, 0, -2}; // -2x^2 +1
        String s1 = Ex1.poly(p1);
        assertEquals("-2.0x^2 +1.0", s1);

        double[] p2 = {0};
        String s2 = Ex1.poly(p2);
        assertEquals("0", s2);

        double[] p3 = {3, -1, 0, 2};
        String s3 = Ex1.poly(p3);
        assertEquals("2.0x^3 -1.0x +3.0", s3);
    }


    @Test
    void testPolynomFromPoints() {
        double []poly1 = Ex1.PolynomFromPoints(xPoints0, yPoints0);
        double []poly2 = Ex1.PolynomFromPoints(xPoints1, yPoints1);
        assertEquals(-7.5, poly1[0], 0);
        assertEquals(3.25, poly1[1], 0);

        double[] xs1 = {0, 1};
        double[] ys1 = {1, 3};
        double[] poly3 = Ex1.PolynomFromPoints(xs1, ys1);
        assertArrayEquals(new double[]{1, 2}, poly3, 0);

        double[] xs2 = {0, 1, 2};
        double[] ys2 = {1, 3, 7};
        double[] poly4 = Ex1.PolynomFromPoints(xs2, ys2);
        assertArrayEquals(new double[]{1, 1, 1}, poly4, 0);
    }

    // -------------------- Edge case / additional tests --------------------

    @Test
    void testSameValueEdgeCases() {
        // Same polynomial, should return any x in range
        double[] p1 = {0, 1}; // x
        double[] p2 = {0, 1}; // x
        double x = Ex1.sameValue(p1, p2, -10, 10, Ex1.EPS);
        assertTrue(x >= -10 && x <= 10);

        // Constant polynomials
        double[] p3 = {5};
        double[] p4 = {5};
        assertEquals(0, Ex1.sameValue(p3, p4, -1, 1, Ex1.EPS), Ex1.EPS);
    }

    @Test
    void testLengthEdgeCases() {
        // Single segment
        double[] p1 = {0, 1}; // x
        double len = Ex1.length(p1, 0, 2, 1);
        assertEquals(2.828, len, 0.01);

        // Zero polynomial
        double[] p2 = {0};
        double len2 = Ex1.length(p2, 0, 10, 5);
        assertEquals(10, len2, 0.01); // horizontal distance is sum of dx, dy=0
    }

    @Test
    void testAreaEdgeCases() {
        // Overlapping polynomials
        double[] p1 = {0};
        double[] p2 = {0};
        double area = Ex1.area(p1, p2, -5, 5, 10);
        assertEquals(0, area, 0);

        // One polynomial zero
        double[] p3 = {0};
        double[] p4 = {2};
        double area2 = Ex1.area(p3, p4, 0, 3, 3);
        assertEquals(6, area2, 0.01);

        // Negative values
        double[] p5 = {-1, 0, 1}; // x^2 -1
        double[] p6 = {0};
        double area3 = Ex1.area(p5, p6, -1, 1, 100);
        assertTrue(area3 > 0); // area is always positive
    }

    @Test
    void testGetPolynomFromStringEdgeCases() {
        // Single term
        String s1 = "-5x^3";
        double[] p1 = Ex1.getPolynomFromString(s1);
        assertArrayEquals(new double[]{0,0,0,-5}, p1, 0);

        // Term with implicit coefficient 1
        String s2 = "x^2 - x + 3";
        double[] p2 = Ex1.getPolynomFromString(s2);
        assertArrayEquals(new double[]{3, -1, 1}, p2, 0);

        // Leading '+'
        String s3 = "+3x^2 + 2x +1";
        double[] p3 = Ex1.getPolynomFromString(s3);
        assertArrayEquals(new double[]{1, 2, 3}, p3, 0);

        // Empty polynomial string
        String s4 = "0";
        double[] p4 = Ex1.getPolynomFromString(s4);
        assertArrayEquals(new double[]{0}, p4, 0);
    }

    @Test
    void testAddEdgeCases() {
        // Adding zero polynomial
        double[] p1 = {1,2,3};
        double[] sum = Ex1.add(p1, Ex1.ZERO);
        assertArrayEquals(p1, sum, 0);

        // Different lengths
        double[] p2 = {1};
        double[] p3 = {1,2,3,4};
        double[] result = Ex1.add(p2, p3);
        assertArrayEquals(new double[]{2,2,3,4}, result, 0);

        // Both zero
        double[] sum2 = Ex1.add(Ex1.ZERO, Ex1.ZERO);
        assertArrayEquals(Ex1.ZERO, sum2, 0);
    }

    @Test
    void testMulEdgeCases() {
        // Multiply by zero
        double[] p1 = {1, 2, 3};
        double[] m1 = Ex1.mul(p1, Ex1.ZERO);
        assertArrayEquals(Ex1.ZERO, m1, 0);

        // Multiply by one-term polynomial
        double[] p2 = {2};
        double[] m2 = Ex1.mul(p1, p2);
        assertArrayEquals(new double[]{2,4,6}, m2, 0);

        // Multiply two zero polynomials
        double[] m3 = Ex1.mul(Ex1.ZERO, Ex1.ZERO);
        assertArrayEquals(Ex1.ZERO, m3, 0);
    }

    @Test
    void testDerivativeEdgeCases() {
        // Constant polynomial
        double[] p1 = {5};
        double[] d1 = Ex1.derivative(p1);
        assertArrayEquals(new double[]{0}, d1, 0);

        // Single-term x
        double[] p2 = {0,1};
        double[] d2 = Ex1.derivative(p2);
        assertArrayEquals(new double[]{1}, d2, 0);

        // Higher degree polynomial
        double[] p3 = {0,0,3}; // 3x^2
        double[] d3 = Ex1.derivative(p3);
        assertArrayEquals(new double[]{0,6}, d3, 0);
    }

    @Test
    void testPolynomFromPointsEdgeCases() {
        // Linear polynomial
        double[] xs1 = {1,1};
        double[] ys1 = {2,2};
        assertNull(Ex1.PolynomFromPoints(xs1, ys1)); // same x points invalid

        // Three points giving quadratic
        double[] xs2 = {0,1,2};
        double[] ys2 = {0,1,4};
        double[] poly2 = Ex1.PolynomFromPoints(xs2, ys2);
        assertArrayEquals(new double[]{0,0,1}, poly2, 0);
    }

    @Test
    void testEqualsEdgeCases() {
        // Polynomials with trailing zeros
        double[] p1 = {1,2};
        double[] p2 = {1,2,0,0};
        assertTrue(Ex1.equals(p1, p2));

        // Polynomials differing by EPS
        double[] p3 = {1,2};
        double[] p4 = {1,2+Ex1.EPS/2};
        assertTrue(Ex1.equals(p3, p4));

        // Polynomials differing more than EPS
        double[] p5 = {1,2};
        double[] p6 = {1,2+2*Ex1.EPS};
        assertFalse(Ex1.equals(p5, p6));
    }

    @Test
    void testPolyEdgeCases() {
        // Zero polynomial
        double[] p1 = {0};
        assertEquals("0", Ex1.poly(p1));

        // Single term negative
        double[] p2 = {-3};
        assertEquals("-3.0", Ex1.poly(p2));

        // Coefficients with zeros in middle
        double[] p3 = {1,0,2}; // 2x^2 +1
        assertEquals("2.0x^2 +1.0", Ex1.poly(p3));

        // Large degree
        double[] p4 = new double[10];
        p4[9] = 1;
        assertEquals("1.0x^9", Ex1.poly(p4));
    }


    @Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
    public void testArea() {
        double x1 = -4, x2 = 0;
        double a1 = Ex1.area(po1, po2, x1, x2, 100);
        double a2 = Ex1.area(po2, po1, x1, x2, 100);

        assertEquals(a1, a2, Ex1.EPS);
    }

    @Test
    /**
     * Test the area when f1(x)=0, f2(x)=x.
     * The area between x=−1 and x=2 is:
     * ∫ |x| dx = 0.5 + 2 = 2.5
     */
    public void testArea2() {
        double[] po_a = Ex1.ZERO;  // f(x)=0
        double[] po_b = {0,1};     // f(x)=x

        double x1 = -1;
        double x2 = 2;
        double area = 2.5;

        double a1 = Ex1.area(po_a, po_b, x1, x2, 1);
        double a2 = Ex1.area(po_a, po_b, x1, x2, 2);
        double a3 = Ex1.area(po_a, po_b, x1, x2, 3);
        double a100 = Ex1.area(po_a, po_b, x1, x2, 100);

        assertEquals(area, a1, Ex1.EPS);
        assertEquals(area, a2, Ex1.EPS);
        assertEquals(area, a3, Ex1.EPS);
        assertEquals(area, a100, Ex1.EPS);
    }

    @Test
    /**
     * Test the area between two polynomials on a computed interval.
     */
    public void testArea3() {
        double[] po_a = {2,1,-0.7, -0.02,0.02};
        double[] po_b = {6, 0.1, -0.2};

        // find intersection
        double x1 = Ex1.sameValue(po_a, po_b, -10, -5, Ex1.EPS);

        double a1 = Ex1.area(po_a, po_b, x1, 6, 8);
        double expected = 58.5658;

        assertEquals(expected, a1, Ex1.EPS);
    }
}

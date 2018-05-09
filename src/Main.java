
import java.math.BigDecimal;

import core.GeneratorCreator;
import kurs.BigDecimalMatrix;
import kurs.Matrix;

public class Main {

    static BigDecimal sqrt(BigDecimal value) {
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        return x.add(new BigDecimal(value.subtract(x.multiply(x)).doubleValue() / (x.doubleValue() * 2.0)));
    }

    public static void main(String[] args) {
            BigDecimal gamma = new BigDecimal(0.4);
            BigDecimal lambda = new BigDecimal(0.3);
            int scale = 20;
            BigDecimal T = new BigDecimal(2);
            BigDecimal accuracy = new BigDecimal(0.00001);
            BigDecimal K = new BigDecimal(5);

            BigDecimal[][] d0 = {
                    {new BigDecimal(-0.81156/2), new BigDecimal(0)},
                    {new BigDecimal(0), new BigDecimal(-0.026346/2)}
            };
            BigDecimal[][] d1 = {
                    {new BigDecimal(0.80616/2), new BigDecimal(0.0054/2)},
                    {new BigDecimal(0.014676/2), new BigDecimal(0.01167/2)}
            };
            BigDecimalMatrix one = new BigDecimalMatrix(d0, scale);
            BigDecimalMatrix two = new BigDecimalMatrix(d1, scale);
            GeneratorCreator gC = new GeneratorCreator(gamma, lambda, K.intValue(), one, two, T, accuracy);
            BigDecimalMatrix mn = gC.funcM(2);
            BigDecimalMatrix mWithHat = gC.funcMWithHat(2);


            System.out.println(gC.funcPhi(new BigDecimal(3), new BigDecimal(2)));
            //System.out.println(gC.funcPHi(0));
            //System.out.println(gC.funcP(one, two, 0));
            System.out.println("Matrix Mn:");
            System.out.println(mn);
            System.out.println("Matrix mWithHat:");
            System.out.println(mWithHat);
    }
}

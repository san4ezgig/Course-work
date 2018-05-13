
import java.math.BigDecimal;
import java.math.BigInteger;

import core.GeneratorCreator;
import kurs.BigDecimalMatrix;
import core.*;
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
            BigDecimal T = new BigDecimal(1);
            BigDecimal accuracy = new BigDecimal(0.00001);
            BigDecimal K = new BigDecimal(20);

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

            BigDecimalMatrix g0 = BigDecimalMatrix.identity(2*K.intValue() + 2);
            GeneratorCreator generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), one, two, T, accuracy);
            GMatrixCreator gMatrixCreator = new GMatrixCreator(generatorCreator);
            PSlashMatrixCreator pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
            PhiMatrixCreator phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
            StationaryDistributionCreator sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
            PerformanceParameters pParameters = new PerformanceParameters(sdCreator.getPiVectors());

            System.out.println(pParameters.getAverageNumberOfRequests().toString());
            System.out.println(pParameters.getAverageNumberOfEnergyUnits().toString());
            //System.out.println(sdCreator.getPiVectors());
            //System.out.println(sdCreator.getPiVectors());
            //System.out.println(generatorCreator.funcPhi(new BigDecimal(5), new BigDecimal(3)));
            //System.out.println(gMatrixCreator.create(g0));
    }
}

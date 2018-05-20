
import java.math.BigDecimal;
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
            BigDecimal T = new BigDecimal(1.5);
            BigDecimal accuracy = new BigDecimal(0.001);
            BigDecimal K = new BigDecimal(4);

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
            BigDecimal HALF = new BigDecimal("0.5");

           /* BigDecimalMatrix g0 = BigDecimalMatrix.identity(2*K.intValue() + 2);
            GeneratorCreator generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), one, two, T, accuracy);
            GMatrixCreator gMatrixCreator = new GMatrixCreator(generatorCreator);
            PSlashMatrixCreator pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
            PhiMatrixCreator phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
            StationaryDistributionCreator sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
            PerformanceParameters pParameters = new PerformanceParameters(sdCreator.getPiVectors());

            System.out.println(pParameters.getAverageNumberOfEnergyUnits());
            MatrixContainer.reInit();*/
            for(BigDecimal t = HALF; t.compareTo(T) <= 0; t = t.add(HALF)) {
                System.out.println("T = " + t);
                for(int k = 1; k < 20; k++) {
                    K = new BigDecimal(k);
                    BigDecimalMatrix g0 = BigDecimalMatrix.identity(2*K.intValue() + 2);
                    GeneratorCreator generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), one, two, t, accuracy);
                    GMatrixCreator gMatrixCreator = new GMatrixCreator(generatorCreator);
                    PSlashMatrixCreator pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                    PhiMatrixCreator phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
                    StationaryDistributionCreator sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
                    PerformanceParameters pParameters = new PerformanceParameters(sdCreator.getPiVectors());
                    System.out.println(pParameters.getEnergyUnitLooseProbability().toString().substring(0, 11));
                    MatrixContainer.reInit();
                }
            }
    }
}

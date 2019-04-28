package core;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Objects;
import kurs.BigDecimalMatrix;

public class StationaryDistributionCreator {
    private PSlashMatrixCreator pSlashCreator;
    private ArrayList<BigDecimalMatrix> phiMatrices;
    private Integer K;
    private ArrayList<BigDecimalMatrix> piVectors;

    public StationaryDistributionCreator(PSlashMatrixCreator pSlashCreator, ArrayList<BigDecimalMatrix> phiMatrices, Integer K) {
        this.pSlashCreator = pSlashCreator;
        this.phiMatrices = phiMatrices;
        this.K = K;
        this.piVectors = new ArrayList<>(phiMatrices.size());
    }

    public ArrayList<BigDecimalMatrix> getPiVectors() {
        int size = 2 * (K + 1);
        BigDecimalMatrix A = BigDecimalMatrix.identity(size).subtract(pSlashCreator.create(0, 0));
        BigDecimalMatrix IWaved = BigDecimalMatrix.identity(size);
        IWaved.setElement(0, 0, BigDecimal.ZERO);

        BigDecimalMatrix eWaved = BigDecimalMatrix.eRow(size, BigDecimal.ZERO);
        eWaved.setElement(0, 0, BigDecimal.ONE);

        BigDecimalMatrix e = BigDecimalMatrix.eCol(size, BigDecimal.ONE);
        BigDecimalMatrix sum = BigDecimalMatrix.zeroMatrix(size);

        for (BigDecimalMatrix phi : phiMatrices) {
            sum = sum.add(phi.multiply(e).multiply(eWaved));
        }

        BigDecimalMatrix AWaved = A.multiply(IWaved).add(sum);

        BigDecimalMatrix pi0 = eWaved.multiply(AWaved.inverse());
        piVectors.add(0, pi0);
        MatrixContainer.getPiVectors().put(Objects.hash(0), pi0);

        for (int i = 1; i < phiMatrices.size(); i++) {
            MatrixContainer.getPiVectors().put(Objects.hash(i), pi0.multiply(phiMatrices.get(i)));
            piVectors.add(i, pi0.multiply(phiMatrices.get(i)));
        }

        // System.out.println(piVectors.size());

        return piVectors;
    }

    public void checkPIMatrix(BigDecimal T, BigDecimal gamma, BigDecimalMatrix d0, BigDecimal lambda) {
        ArrayList<BigDecimalMatrix> vectors = this.getPiVectors();
        BigDecimal result = BigDecimal.ZERO;
        BigDecimal sum = BigDecimal.ZERO;
        BigDecimalMatrix eCol = BigDecimalMatrix.eCol(2, BigDecimal.ONE);
        BigDecimalMatrix negativeD0 = d0.multiply(new BigDecimal(-1));
        BigDecimalMatrix piRow;

        result = result.add(T);

        for(int i = 0; i < phiMatrices.size(); i++) {
            piRow = BigDecimalMatrix.eRow(2, vectors.get(i).getElement(0, 0));
            piRow.setElement(0, 1, vectors.get(i).getElement(0, 1));
            sum = sum.add(piRow.multiply(new BigDecimal(1 / gamma.doubleValue())).multiply(eCol).getElement(0, 0));
        }
        result = result.add(sum);

        sum = BigDecimal.ZERO;
        for(int k = 0; k <= K; k++) {
            piRow = BigDecimalMatrix.eRow(2, vectors.get(0).getElement(0, k * 2));
            piRow.setElement(0, 1, vectors.get(0).getElement(0, 2 * k + 1));

            sum = sum.add(piRow.multiply(negativeD0.inverse().multiply(eCol)).getElement(0, 0));
        }
        result = result.add(sum);

        piRow = BigDecimalMatrix.eRow(2, vectors.get(0).getElement(0, 0));
        piRow.setElement(0, 1, vectors.get(0).getElement(0, 1));
        sum = piRow
                .multiply(new BigDecimal(2))
                .multiply(negativeD0.add(BigDecimalMatrix.identity(2).multiply(gamma)).inverse())
                .multiply(eCol)
                .getElement(0 ,0);
        result = result.subtract(sum);

        System.out.println(result);
        System.out.println(1 / lambda.doubleValue());
    }
}


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
        BigDecimalMatrix A = BigDecimalMatrix.identity(2*K + 2).subtract(pSlashCreator.create(0, 0));
        BigDecimalMatrix IWaved = BigDecimalMatrix.identity(2*K + 2);
        IWaved.setElement(0, 0, BigDecimal.ZERO);

        BigDecimalMatrix eWaved = BigDecimalMatrix.eRow(2*K + 2, BigDecimal.ZERO);
        eWaved.setElement(0, 0, BigDecimal.ONE);

        BigDecimalMatrix e = BigDecimalMatrix.eCol(2*K + 2, BigDecimal.ONE);
        BigDecimalMatrix sum = BigDecimalMatrix.zeroMatrix(2*K + 2);

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

        return piVectors;
    }
}


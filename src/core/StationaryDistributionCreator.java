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
        //System.out.println(pSlashCreator.create(0, 0).getWidth());
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
        BigDecimalMatrix a = new BigDecimalMatrix(2, K + 1, 30);
        for(int j = 0; j < pi0.getWidth(); j += 2) {
            a.setElement(0, j/2, pi0.getElement(0, j));
            a.setElement(1, j/2, pi0.getElement(0, j + 1));
        }
        piVectors.add(0, a);
        MatrixContainer.getPiVectors().put(Objects.hash(0), pi0);

        for (int i = 1; i < phiMatrices.size(); i++) {
            MatrixContainer.getPiVectors().put(Objects.hash(i), pi0.multiply(phiMatrices.get(i)));
            BigDecimalMatrix b = pi0.multiply(phiMatrices.get(i));
            BigDecimalMatrix val = new BigDecimalMatrix(2, K + 1, 30);
            for(int j = 0; j < b.getWidth(); j += 2) {
                val.setElement(0, j/2, b.getElement(0, j));
                val.setElement(1, j/2, b.getElement(0, j + 1));
            }

            piVectors.add(i, val);
        }
        return piVectors;
    }
}


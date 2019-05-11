package core;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Objects;

import core.MatrixContainer;
import kurs.BigDecimalMatrix;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Objects;

public class PhiMatrixCreator {
    private ArrayList<BigDecimalMatrix> phiMatrices;
    private PSlashMatrixCreator pSlashCreator;
    private BigDecimal accuracy;
    private Integer K;
    private int systemSize;

    //Approved
    public PhiMatrixCreator(PSlashMatrixCreator pSlashCreator, Integer K, int systemSize) {
        phiMatrices = new ArrayList<>();
        this.pSlashCreator = pSlashCreator;
        this.K = K;
        this.accuracy = pSlashCreator.getAccuracy();
        this.systemSize = systemSize;
    }

    //Approved
    public ArrayList<BigDecimalMatrix> getPhiMatrices() {
        initPsiMatrices();
        return phiMatrices;
    }

    //Approved
    private void initPsiMatrices() {
        final BigDecimalMatrix I = BigDecimalMatrix.identity(systemSize*(K + 1));
        phiMatrices.add(0, I);
        int hash = Objects.hash(0);
        try {
            if (!MatrixContainer.getPhiMatrices().containsKey(hash)) {
                MatrixContainer.getPhiMatrices().put(hash, I.clone());
            }
            int index = 1;
            BigDecimalMatrix leftTerm;
            BigDecimalMatrix rightTerm;
            BigDecimalMatrix result;
            do {
                hash = Objects.hash(index);
                if (MatrixContainer.getPhiMatrices().containsKey(hash)) {
                    phiMatrices.add(index, MatrixContainer.getPhiMatrices().get(hash).clone());
                } else {
                    rightTerm = I.subtract(pSlashCreator.create(index, index)).inverse();

                    leftTerm = BigDecimalMatrix.zeroMatrix(systemSize *(K + 1));

                    for(int i = 1; i <= index - 1; i++) {
                        leftTerm = leftTerm.add(phiMatrices.get(i).multiply(pSlashCreator.create(i, index)));
                    }

                    leftTerm = leftTerm.add(pSlashCreator.create(0, index));
                    result = leftTerm.multiply(rightTerm);
                    phiMatrices.add(index, result);
                    MatrixContainer.getPhiMatrices().put(hash, result.clone());
                }
                index++;
            } while (phiMatrices.get(index - 1).subtract(phiMatrices.get(index - 2)).squaredEuclidianNorm().doubleValue()
            > accuracy.doubleValue());
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
    }
}

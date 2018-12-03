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

    //Approved
    public PhiMatrixCreator(PSlashMatrixCreator pSlashCreator, Integer K) {
        phiMatrices = new ArrayList<>();
        this.pSlashCreator = pSlashCreator;
        this.K = K;
        this.accuracy = pSlashCreator.getAccuracy();
    }

    //Approved
    public ArrayList<BigDecimalMatrix> getPhiMatrices() {
        initPsiMatrices();
        return phiMatrices;
    }

    //Approved
    private void initPsiMatrices() {
        final BigDecimalMatrix I = BigDecimalMatrix.identity(2*K + 2);
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
                    //.inverse()
                    rightTerm = I.subtract(pSlashCreator.create(index, index)).inverse();


                    int i = 1;
                    leftTerm = BigDecimalMatrix.zeroMatrix(K*2 + 2);
                    //System.out.println(1);
                    while (i < index) {
                        leftTerm = leftTerm.add(phiMatrices.get(i).multiply(pSlashCreator.create(i, index)));
                        i++;
                    }
                    leftTerm = leftTerm.add(pSlashCreator.create(0, index));
                    result = leftTerm.multiply(rightTerm);
                    //System.out.println(result.squaredEuclidianNorm());
                    phiMatrices.add(index, result);
                    MatrixContainer.getPhiMatrices().put(hash, result.clone());
                }
                index++;
            } while (phiMatrices.get(index - 1).subtract(phiMatrices.get(index - 2)).squaredEuclidianNorm().compareTo(accuracy) > 0);
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
    }
}

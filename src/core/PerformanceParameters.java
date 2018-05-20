package core;


import java.math.BigDecimal;
import java.util.ArrayList;
import kurs.BigDecimalMatrix;
import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;

/**
 * Created by user on 04.06.2017.
 */
public class PerformanceParameters {
    private ArrayList<BigDecimalMatrix> piVector;
    private int piVectorSize;
    private int piSize;

    public PerformanceParameters(ArrayList<BigDecimalMatrix> piVector) {
        this.piVector = piVector;
        this.piVectorSize = piVector.size();
        this.piSize = piVector.get(0).getWidth();
    }

    public BigDecimal getAverageNumberOfRequests() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix e = BigDecimalMatrix.eCol(this.piSize, ONE);
        for (int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(valueOf(i).multiply(piVector.get(i).multiply(e).getElement(0, 0).add(piVector.get(i).multiply(e).getElement(1, 0))));
        }
        return sum.add(new BigDecimal(2));
    }

    public BigDecimal getAverageNumberOfEnergyUnits() {
        BigDecimal sum = ZERO;
        for(int i = 0; i < this.piVectorSize; i++) {
            for (int k = 0; k < this.piSize; k++) {
                BigDecimalMatrix val = new BigDecimalMatrix(1, 2, 30);
                val.setElement(0, 0, piVector.get(i).getElement(0, k));
                val.setElement(0, 1, piVector.get(i).getElement(1, k));
                sum = sum.add(valueOf(k).multiply(val.multiply(BigDecimalMatrix.eCol(2, ONE)).getElement(0, 0)));
            }
        }
        return sum.abs().multiply(new BigDecimal(10));
    }

    public BigDecimal getEnergyUnitLooseProbability() {
        BigDecimal sum = ZERO;
        for(int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).getElement(0, this.piSize - 1));
        }
        return sum;
    }

    public BigDecimal getNoRequestsProbability() {
        return piVector.get(0).multiply(BigDecimalMatrix.eCol(this.piSize, ONE)).getElement(0, 0).add(piVector.get(0).multiply(BigDecimalMatrix.eCol(this.piSize, ONE)).getElement(1, 0));
    }

    public BigDecimal getNoEnergyUnitsProbability() {
        BigDecimalMatrix e = BigDecimalMatrix.eCol(2, ONE);
        BigDecimal sum = ZERO;
        for(int i = 0; i < this.piVectorSize; i++) {
            BigDecimalMatrix val = new BigDecimalMatrix(1, 2, 30);
            val.setElement(0, 0, piVector.get(i).getElement(0, 0));
            val.setElement(0, 1, piVector.get(i).getElement(1, 0));
            sum = sum.add(val.multiply(e).getElement(0, 0));
        }
        return sum;
    }

    public BigDecimal getNoEnergySystemIdleProbability() {
        BigDecimalMatrix e = BigDecimalMatrix.eCol(2, ONE);
        BigDecimal sum = ZERO;
        for(int i = 1; i < this.piVectorSize; i++) {
            BigDecimalMatrix val = new BigDecimalMatrix(1, 2, 30);
            val.setElement(0, 0, piVector.get(i).getElement(0, 0));
            val.setElement(0, 1, piVector.get(i).getElement(1, 0));
            sum = sum.add(val.multiply(e).getElement(0, 0));
        }
        return sum;
    }
}

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
    private BigDecimal lambda;
    private BigDecimalMatrix d1;

    public PerformanceParameters(ArrayList<BigDecimalMatrix> piVector, BigDecimal lambda, BigDecimalMatrix d1) {
        this.piVector = piVector;
        this.piVectorSize = piVector.size();
        this.d1 = d1;
        this.piSize = piVector.get(0).getWidth();
        this.lambda = lambda;
    }

    public BigDecimal getAverageNumberOfRequests() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix e = BigDecimalMatrix.eCol(this.piSize, ONE);
        for (int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(valueOf(i).multiply(piVector.get(i).multiply(e).getElement(0, 0)));
        }
        return sum;
    }

    public BigDecimal getAverageNumberOfEnergyUnits() {
        BigDecimal sum = ZERO;
        for (int i = 0; i < this.piVectorSize; i++) {
            for (int k = 0; k < this.piSize; k++) {
                sum = sum.add(valueOf(k).multiply(piVector.get(i).getElement(0, k)));
            }
        }
        return sum;
    }

    // It doesnt work
    /*public BigDecimal getEnergyUnitLossProbability() {
        BigDecimal sum = ZERO;
        for(int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).getElement(0, this.piSize - 1));
        }
        return sum;
    }*/

    public BigDecimal getNoRequestsProbability() {
        return piVector.get(0).multiply(BigDecimalMatrix.eCol(this.piSize, ONE)).getElement(0, 0);
    }

    public BigDecimal getNoEnergyUnitsProbability() {
        BigDecimal sum = ZERO;
        for (int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).getElement(0, 0));
        }
        return sum;
    }

    public BigDecimal getNoEnergySystemIdleProbability() {
        BigDecimal sum = ZERO;
        for (int i = 1; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).getElement(0, 0));
        }
        return sum;
    }

    //Запрос начнет обслуживаться в момент прихода
    public BigDecimal getRequestStartWorkedWhenHisWent() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix eCol = BigDecimalMatrix.eCol(2, ONE);
        BigDecimalMatrix eRow = BigDecimalMatrix.eRow(2, ONE);

        for (int k = 0; k < this.piSize; k++) {
            sum = eRow.multiply(d1).multiply(piVector.get(0).getElement(0, k)).multiply(eCol).getElement(0, 0);
        }

        return new BigDecimal(1 / lambda.doubleValue() * sum.doubleValue());
    }

    public void check(BigDecimalMatrix tetta) {
        BigDecimal sum = ZERO;
        for (int i = 0; i < this.piVectorSize; i++) {
            for (int k = 0; k < this.piSize; k++) {
                sum = sum.add(valueOf(k).multiply(piVector.get(i).getElement(0, k)));
            }
        }
        System.out.println(tetta);
        System.out.println(sum);
    }
}

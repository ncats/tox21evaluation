
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * This is the AUC ROC code used to evaluate the Tox21 challenge. Comment
 * and questions about this code should be addressed to Dr. Ruili Huang at
 * huangru@mail.nih.gov.
 */
public class ROC {
        
    //vals-->real values, refVals-->submitted values, activeVal=1, prec=0.001
    public static Map<String,Object> getROC(double vals[], double refVals[], int len, double activeVal, double prec)  {
                 
        /*System.out.println("*");
                
          System.out.println("*vals[]="+Arrays.toString(vals));
          System.out.println("*refVals[]="+Arrays.toString(refVals));*/
                
        double roc, val, pVal;
        int pValCnt;
        double specVal, sensVal, specPre, sensPre;
        double pVals[],tmp[];
            
        double inactiveVal = 0.1;
            
        double[] x = new double[len], y=new double[len];

        tmp = new double[len];
                
        pValCnt = 0;
        for(int i=0; i<len; i++) {
                 
            if(!Double.isNaN(refVals[i])) {
                tmp[pValCnt] = refVals[i];
                pValCnt++;
            }
        } 
        //System.out.println("*pValCnt="+pValCnt);
                
        pVals = sortHeap(tmp, pValCnt, prec);
            
        /*
         * 
         */
        Double[] doubleArray = ArrayUtils.toObject(pVals);
        List<Double> list = Arrays.asList(doubleArray); 
        Collections.sort(list); 
        pVals= ArrayUtils.toPrimitive(list.toArray(new Double[0]));
            
        /*for(int m=0;m<pVals.length;m++){
          System.out.println((m+1)+"\t"+pVals[m]);
          }*/

        int TP[] = new int[pValCnt];
        int FP[] = new int[pValCnt];
        int TN[] = new int[pValCnt];
        int FN[] = new int[pValCnt];

        for(int i=0; i<pValCnt; i++) {
            TP[i] = 0;
            FP[i] = 0;
            TN[i] = 0;
            FN[i] = 0;          
        }

        for (int j = 0; j < len; j++) {
            pVal = refVals[j];
            val = vals[j];

            if (!Double.isNaN(val)) {
                for (int i = 0; i < pValCnt; i++) {

                    if (val > activeVal) {
                        if (Double.isNaN(pVal))
                            FN[i]++;
                        else if (pVal >= pVals[i])
                            TP[i]++;
                        else
                            FN[i]++;
                    } else {
                        if (val < inactiveVal) {
                            if (Double.isNaN(pVal))
                                FP[i]++;
                            else if (pVal >= pVals[i])
                                FP[i]++;
                            else
                                TN[i]++;
                        }
                    }
                }
            }
        }

        roc = 0;

        specPre = (double)(FP[0])/(double)(FP[0]+TN[0]);
        sensPre = (double)(TP[0])/(double)(TP[0]+FN[0]);

        for(int p=0; p<pValCnt; p++) {
                
            if (p > 0) {
                specVal = (double) (FP[p]) / (double) (FP[p] + TN[p]);// x
                x[p] = specVal;

                sensVal = (double) (TP[p]) / (double) (TP[p] + FN[p]);// y
                y[p] = sensVal;

                roc += (sensVal + sensPre) * (specPre - specVal) / 2.0;
                specPre = specVal;
                sensPre = sensVal;
                // System.out.println(p+"\t"+x[p]+"\t"+y[p]);
            } 
        }

        Map<String,Object> map = new HashMap<String, Object>();
            
        map.put("x", x);
        map.put("y", y); 
        map.put("roc", roc);
            
        return map; 
    }

    private static double[] sortHeap(double vals[], int len, double prec) {
        int HeapSize = len;
        double tmpVal, tmpVals[];

        tmpVals = new double[len];

        for (int i = 0; i < len; i++) {
            tmpVals[i] = prec * Math.round(vals[i] / prec);
        }

        for (int i = HeapSize / 2 - 1; i >= 0; i--)
            sortHeapify(tmpVals, i, HeapSize);

        for (int i = len - 1; i >= 1; i--) {

            tmpVal = tmpVals[i];
            tmpVals[i] = tmpVals[0];
            tmpVals[0] = tmpVal;

            HeapSize--;
            sortHeapify(tmpVals, 0, HeapSize);
        }

        return tmpVals;
    }

    private static void sortHeapify(double A[], int i, int HeapSize) {

        int left=2*i, right=2*i+1;
        int largest;
        double tmp;


        if ((left<HeapSize) && (A[left]>A[i])) // sort asc
            largest = left;
        else
            largest = i;

        if ((right<HeapSize) && (A[right]>A[largest])) // sort asc
            largest = right;

        if (largest != i) {

            tmp = A[i];
            A[i] = A[largest];
            A[largest] = tmp;

            sortHeapify(A, largest, HeapSize);
        }
    }
}

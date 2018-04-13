/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class NavierStokesBodyForced extends NavierStokes {
    double[] acceleration;
    
    @Override
    public void init(FlowProProperties props) throws IOException {
        super.init(props);
        
        // acceleration
        if(props.containsKey("acceleration")){
            acceleration = props.getDoubleArray("acceleration");
        } else {
            acceleration = new double[dim];
        }
    }
    
    @Override
    public boolean isSourcePresent() {
        return true;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        double[] source = new double[nEqs];
        double[] avg = elem.integralMonitor;
        double uMean = 0.2;
        source[1] = -W[0]*(avg[0] - uMean);
        return source;
    }
}

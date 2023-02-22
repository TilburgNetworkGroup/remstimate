# testing tie-oriented modeling

    Code
      tie_mle
    Output
      Relational Event Model (tie oriented) 
      
      Coefficients:
      
            baseline indegreeSender        inertia    reciprocity 
         -4.99153889     0.01349127    -0.16034529     0.02923781 
      
      Null deviance: 11634.18 
      Residual deviance: 11593.73 
      AIC: 11601.73 AICC: 11601.77 BIC: 11621.36 
      

---

    Code
      summary(tie_mle)
    Output
      Relational Event Model (tie oriented) 
      
      Call:
      ~baseline + indegreeSender + inertia + reciprocity
      
      
      Coefficients (MLE with interval likelihood):
      
                        Estimate    Std. Err     z value Pr(>|z|)    Pr(=0)
      baseline        -4.9915389   0.0617953 -80.7753244   0.0000 < 2.2e-16
      indegreeSender   0.0134913   0.0029837   4.5217287   0.0000 0.0011474
      inertia         -0.1603453   0.0320125  -5.0088263   0.0000 0.0001127
      reciprocity      0.0292378   0.0316181   0.9247164   0.3551 0.9537494
      Null deviance: 11634.18 on 1000 degrees of freedom
      Residual deviance: 11593.73 on 996 degrees of freedom
      Chi-square: 40.45225 on 4 degrees of freedom, asymptotic p-value 3.489608e-08 
      AIC: 11601.73 AICC: 11601.77 BIC: 11621.36 


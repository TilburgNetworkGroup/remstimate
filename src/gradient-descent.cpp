
//if loss function oscillates in later epochs, then reduce learning rate for better convergence
// [[Rcpp::export]]
Rcpp::List GD(const arma::vec & pars,const arma::cube & stats, const arma::mat &event_binary,const arma::vec & interevent_time,int epochs, double learning_rate){
    arma::vec pars_prev = pars;
    arma::vec pars_next;
    arma::uword P = pars.n_elem;

    //for output
    arma::mat beta(P,epochs,arma::fill::zeros);
    arma::vec loss(epochs,arma::fill::zeros);

    for(int i =0; i<epochs;i++){
        Rcpp::List derv = remDerivatives(pars_prev,stats,event_binary,interevent_time,true,false);
        double loglik = derv["value"];
        arma::vec grad = derv["gradient"];
        pars_next = pars_prev - (learning_rate * grad);
        beta.col(i) = pars_next;
        loss(i) = loglik;
        pars_prev = pars_next;
    }
    return Rcpp::List::create(Rcpp::Named("loss") = loss, Rcpp::Named("coef") = pars_next,Rcpp::Named("betas")=beta);
}

// default values for hyper-parameters recommended in literature
// β1 =  0.9 for β2 =  0.999
// eta = 10^-8
// [[Rcpp::export]]
Rcpp::List GD_ADAM(const arma::vec & pars,const arma::cube & stats, const arma::mat &event_binary,const arma::vec & interevent_time,int epochs=50, double learning_rate=0.002,double beta1=0.9, double beta2=0.999,double eta = 0.00000001){

    arma::uword P = pars.n_elem;
    arma::vec pars_prev = pars;
    arma::vec pars_next;


    arma::vec moment1(P,arma::fill::zeros);
    arma::vec moment2(P,arma::fill::zeros);

    arma::vec moment1_est(P,arma::fill::zeros);
    arma::vec moment2_est(P,arma::fill::zeros);

    //for output
    arma::mat beta(P,epochs,arma::fill::zeros);
    arma::vec loss(epochs,arma::fill::zeros);

    for(int i =0; i<epochs;i++){
        Rcpp::List derv = remDerivatives(pars_prev,stats,event_binary,interevent_time,true,false);
        double loglik = derv["value"];
        
        arma::vec grad = derv["gradient"];
        moment1 = (beta1*moment1) + (1-beta1)*grad;
        moment2 = (beta2 * moment2) + (1-beta2) * (grad%grad);
        // bias-corrected first and second moment estimates:
        moment1_est = moment1 / (1-pow(beta1,i+1));
        moment2_est = moment2/  (1-pow(beta2,i+1));
        pars_next = pars_prev - (learning_rate * moment1_est)/(sqrt(moment2_est)+eta);
        
        beta.col(i) = pars_next;
        loss(i) = loglik;
        
        pars_prev = pars_next;
    }
    return Rcpp::List::create(Rcpp::Named("loss") = loss, Rcpp::Named("coef") = pars_next,Rcpp::Named("betas")=beta);
}


#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int sample_int(int n) {
  Rcpp::IntegerVector pool = Rcpp::seq(1, n);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[0];
}

double lfactorial_single(int n){
  
  NumericVector xx(1); xx(0) = n;
  return(lfactorial(xx)[0]);
}

void quicksort(IntegerVector order, NumericVector x, int l, int r)
{
  // Base case: No need to sort arrays of length <= 1
  if (l >= r)
  {
    return;
  }
  
  // Choose pivot to be the last element in the subarray
  double pivot = x[r];
  
  // Index indicating the "split" between elements smaller than pivot and 
  // elements greater than pivot
  int cnt = l;
  
  // Traverse through array from l to r
  for (int i = l; i <= r; i++)
  {
    // If an element less than or equal to the pivot is found...
    if (x[i] <= pivot)
    {
      // Rcout << "ok" << std::endl;
      // Then swap arr[cnt] and arr[i] so that the smaller element arr[i] 
      // is to the left of all elements greater than pivot
      // swap(&arr[cnt], &arr[i]);
      int ordercnt = order[cnt];
      order[cnt] = order[i];
      order[i] = ordercnt; 
      
      double xcnt = x[cnt];
      x[cnt] = x[i];
      x[i] = xcnt; 
      
      // Make sure to increment cnt so we can keep track of what to swap
      // arr[i] with
      cnt++;
      
    }
  }
  
  // NOTE: cnt is currently at one plus the pivot's index 
  // (Hence, the cnt-2 when recursively sorting the left side of pivot)
  quicksort(order, x, l, cnt-2); // Recursively sort the left side of pivot
  quicksort(order, x, cnt, r);   // Recursively sort the right side of pivot
}

arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + arma::chol(sigma) * Y;
}

// [[Rcpp::export]]
bool checkPointIsInRegion(arma::vec x, double a1, double b1, 
                          double a2, double b2){
  
  if(x[0] > a1 & x[0] < b1 & x[1] > a2 & x[1] < b2){
    return(true);
  }
  
  return(false);
}

// [[Rcpp::export]]
arma::vec proposeNewPointPP(double a1, double b1, double a2, double b2){
  
  arma::vec proposedPoint = arma::zeros(2);
  
  proposedPoint[0] = R::runif(a1, b1);
  proposedPoint[1] = R::runif(a2, b2);
  
  return(proposedPoint);
}

// UNIVARIATE SOFTCORE

// [[Rcpp::export]]
double log_f_softcore_cpp(arma::mat x, int N, 
                          double beta, double theta){
  
  double loglikelihood = 0;
  
  loglikelihood += N * log(beta) - lfactorial_single(N);
  
  // within x1
  
  double loglikelihood1 = 0;
  
  for(int n = 1; n < N; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      double r = sqrt( pow(x(n,0) - x(n2,0), 2) + pow(x(n,1) - x(n2,1), 2) );
      
      loglikelihood1 += log(1 - exp(- r*r / (theta)));
      
    }
    
  }  
  
  return(loglikelihood + loglikelihood1);
}

// [[Rcpp::export]]
double log_f_softcore_cpp_quick_birth(arma::vec xnew, double oldLikelihood, 
                                      arma::mat x, int N, 
                                      double beta, double theta){
  
  double loglikelihood = oldLikelihood;
  
  loglikelihood += log(beta);
  loglikelihood -= log(N + 1);
  
  for(int j = 0; j < N; j++){
    
    double r_new = sqrt(pow(x(j,0) - xnew(0), 2) + pow(x(j,1) - xnew(1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_softcore_cpp_quick_death(int i, double oldLikelihood, 
                                      arma::mat x, int N, 
                                      double beta, double theta){
  
  double loglikelihood = oldLikelihood;
  
  loglikelihood -= log(beta);
  loglikelihood += log(N);
  
  // between
  
  for(int j = 0; j < N; j++){
    
    if(j != i){
      
      double r_old = sqrt(pow(x(j,0) - x(i,0), 2) + pow(x(j,1) - x(i,1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_softcore_cpp_quick_move(int i, double oldLikelihood, 
                                     arma::mat x, int N,
                                     arma::vec xold, double theta){
  
  double loglikelihood = oldLikelihood;
  
  // between
  
  for(int j = 0; j < N; j++){
    
    if(j != i){
      
      double r_old = sqrt(pow(x(j,0) - xold(0), 2) + pow(x(j,1) - xold(1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta)));
      
      double r_new = sqrt(pow(x(j,0) - x(i,0), 2) + pow(x(j,1) - x(i,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
List simulate_softcore_cpp(double theta, double beta,
                           int niter, arma::mat Sigma_prop, 
                           int Nmax, double lambda,
                           double a1, double b1,
                           double a2, double b2){
  
  double q = .3333333;
  double p = .5;
  
  arma::mat data;
  data.zeros(Nmax, 2);
  
  int N = R::rpois(lambda);
  for(int i = 0; i < N; i++){
    data.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPointPP(a1, b1, a2, b2));
  }
  
  double log_hastings_ratio_den = log_f_softcore_cpp(data, N, beta, theta);
  
  arma::vec likelihood = arma::zeros(niter);
  
  for(int iter = 0; iter < niter; iter++){
    
    likelihood[iter] = log_hastings_ratio_den;
    
    // x1
    
    if(R::runif(0,1) < q){
      
      // move
      
      // sample point
      int n = sample_int(N) - 1;
      
      arma::rowvec old_xi = data.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        data.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        double log_hastings_ratio_num = log_f_softcore_cpp_quick_move(n, log_hastings_ratio_den,
                                                                      data, N, oldxivec, theta);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data.row(n) = old_xi;
        }
        
      }
      
    } else {
      
      if(R::runif(0,1) < p){
        
        // birth
        
        arma::vec xi = proposeNewPointPP(a1, b1, a2, b2);
        
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          double log_hastings_ratio_num = log_f_softcore_cpp_quick_birth(xi, log_hastings_ratio_den,
                                                                         data, N, beta, theta);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data.row(N) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N += 1;
          }
          
        }
        
      } else {
        
        // death
        
        if (N > 0){
          
          // sample point
          int itemToRemove = sample_int(N) - 1;
          
          double log_hastings_ratio_num = log_f_softcore_cpp_quick_death(itemToRemove, log_hastings_ratio_den,
                                                                         data, N, beta, theta);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data.row(itemToRemove) = data.row(N - 1);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N -= 1;
          }
          
        }
        
      }
      
    }
    
  }
  
  return(List::create(_["data"] = data,
                      _["likelihood"] = likelihood,
                      _["N"] = N));
}

// [[Rcpp::export]]
List simulate_cond_softcore_cpp(double theta, double beta,
                           int niter, arma::mat Sigma_prop, 
                           int N, double lambda,
                           double a1, double b1,
                           double a2, double b2){
  
  arma::mat data;
  data.zeros(N, 2);
  
  for(int i = 0; i < N; i++){
    data.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPointPP(a1, b1, a2, b2));
  }
  
  double log_hastings_ratio_den = log_f_softcore_cpp(data, N, beta, theta);
  
  arma::vec likelihood = arma::zeros(niter);
  
  for(int iter = 0; iter < niter; iter++){
    
    // likelihood[iter] = log_hastings_ratio_den;
    
    // x1
    
    for(int n = 0; n < N; n++){
      
      arma::rowvec old_xi = data.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        data.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        double log_hastings_ratio_num = log_f_softcore_cpp_quick_move(n, log_hastings_ratio_den,
                                                                      data, N, oldxivec, theta);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data.row(n) = old_xi;
        }
        
      }
      
        
    }
    
  }
  
  return(List::create(_["data"] = data,
                      _["N"] = N));
}



//// NORMAL SOFTCORE PROCESS

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp(arma::mat x1, arma::mat x2, int N1, int N2,
                             double beta1, double beta2,
                             double theta1, double theta2, double theta3){
  
  double loglikelihood = 0;
  
  loglikelihood += N1 * log(beta1) - lfactorial_single(N1);
  loglikelihood += N2 * log(beta2) - lfactorial_single(N2);
  
  // within x1
  
  double loglikelihood1 = 0;
  
  for(int n = 1; n < N1; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      double r = sqrt( pow(x1(n,0) - x1(n2,0), 2) + pow(x1(n,1) - x1(n2,1), 2) );
      
      loglikelihood1 += log(1 - exp(- r*r / (theta1)));
      
    }
    
  }  
  
  // within x2
  
  double loglikelihood2 = 0;
  
  for(int n = 1; n < N2; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      double r = sqrt( pow(x2(n,0) - x2(n2,0), 2) + pow(x2(n,1) - x2(n2,1), 2) );
      
      loglikelihood2 += log(1 - exp(- r*r / (theta2)));
      
    }
    
  }
  
  // between
  
  double loglikelihood3 = 0;
  
  for(int n = 0; n < N1; n++){
    
    for(int n2 = 0; n2 < N2; n2++){
      
      double r = sqrt( pow(x1(n,0) - x2(n2,0), 2) + pow(x1(n,1) - x2(n2,1), 2) );
      
      loglikelihood3 += log(1 - exp(- r*r / (theta3)));
      
    }
    
  }
  
  return(loglikelihood + loglikelihood1 + loglikelihood2 + loglikelihood3);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_birth(arma::vec xnew,  
                                         arma::mat x, int N, 
                                         double beta, double theta){
  
  double loglikelihood = 0;
  
  loglikelihood += log(beta);
  loglikelihood -= log(N + 1);
  
  // between
  
  for(int j = 0; j < N; j++){
    
    double r_new = sqrt(pow(x(j,0) - xnew(0), 2) + pow(x(j,1) - xnew(1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_death(int i, arma::mat x, int N, 
                                         double beta, double theta){
  
  double loglikelihood = 0;
  
  loglikelihood -= log(beta);
  loglikelihood += log(N);
  
  // between
  
  for(int j = 0; j < N; j++){
    
    if(j != i){
      
      double r_old = sqrt(pow(x(j,0) - x(i,0), 2) + pow(x(j,1) - x(i,1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_move(int i, int dim, double oldLikelihood, 
                                        arma::mat x1, arma::mat x2, 
                                        int N1, int N2,
                                        arma::vec xold, 
                                        double theta1, double theta2, double theta3){
  
  double loglikelihood = oldLikelihood;
  
  if(dim == 1){
    
    // between
    
    for(int j = 0; j < N1; j++){
      
      if(j != i){
        
        double r_old = sqrt(pow(x1(j,0) - xold(0), 2) + pow(x1(j,1) - xold(1), 2));
        
        loglikelihood -= log(1 - exp(- r_old*r_old / (theta1)));
        
        double r_new = sqrt(pow(x1(j,0) - x1(i,0), 2) + pow(x1(j,1) - x1(i,1), 2));
        
        loglikelihood += log(1 - exp(- r_new*r_new / (theta1)));
        
      }
      
    }
    
    // within
    
    for(int j = 0; j < N2; j++){
      
      double r_old = sqrt(pow(x2(j,0) - xold(0), 2) + pow(x2(j,1) - xold(1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta3)));
      
      double r_new = sqrt(pow(x2(j,0) - x1(i,0), 2) + pow(x2(j,1) - x1(i,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta3)));
      
    }
    
  } else {
    
    // between
    
    for(int j = 0; j < N2; j++){
      
      if(j != i){
        
        double r_old = sqrt(pow(x2(j,0) - xold(0), 2) + pow(x2(j,1) - xold(1), 2));
        
        loglikelihood -= log(1 - exp(- r_old*r_old / (theta2)));
        
        double r_new = sqrt(pow(x2(j,0) - x2(i,0), 2) + pow(x2(j,1) - x2(i,1), 2));
        
        loglikelihood += log(1 - exp(- r_new*r_new / (theta2)));
        
      }
      
    }
    
    // within
    
    for(int j = 0; j < N1; j++){
      
      double r_old = sqrt(pow(x1(j,0) - xold(0), 2) + pow(x1(j,1) - xold(1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta3)));
      
      double r_new = sqrt(pow(x1(j,0) - x2(i,0), 2) + pow(x1(j,1) - x2(i,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta3)));
      
    }
    
  }
  
  
  return(loglikelihood);
}

double log_interfun_softcore(double r, double gamma){
  
  return(log(1 - exp(- r * r / (gamma))));
  
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_move_cond(int i, arma::mat x1, arma::mat x2, 
                                        int N1, int N2,
                                        arma::vec xold, 
                                        double theta1, double theta12){
  
  double loglikelihood = 0;
  
  // between
  
  for(int j = 0; j < N1; j++){
    
    if(j != i){
      
      double r_old = sqrt(pow(x1(j,0) - xold(0), 2) + pow(x1(j,1) - xold(1), 2));
      
      loglikelihood -= log_interfun_softcore(r_old, theta1);
      // loglikelihood -= log(1 - exp(- r_old * r_old / (theta1)));
      
      double r_new = sqrt(pow(x1(j,0) - x1(i,0), 2) + pow(x1(j,1) - x1(i,1), 2));
      
      loglikelihood += log_interfun_softcore(r_new, theta1);
      // loglikelihood += log(1 - exp(- r_new * r_new / (theta1)));
      
    }
    
  }
  
  // within
  
  for(int j = 0; j < N2; j++){
    
    double r_old = sqrt(pow(x2(j,0) - xold(0), 2) + pow(x2(j,1) - xold(1), 2));
    
    loglikelihood -= log_interfun_softcore(r_old, theta12);
    // loglikelihood -= log(1 - exp(- r_old * r_old / (theta3)));
    
    double r_new = sqrt(pow(x2(j,0) - x1(i,0), 2) + pow(x2(j,1) - x1(i,1), 2));
    
    loglikelihood += log_interfun_softcore(r_new, theta12);
    // loglikelihood += log(1 - exp(- r_new * r_new / (theta3)));
    
  }
  
  return(loglikelihood);
}

arma::mat proposePointProcess(arma::mat data1, arma::mat data2, int N1, int N2,
                              double theta1, double theta12,
                              double a1, double b1,
                              double a2, double b2,
                              arma::mat Sigma_prop){
  
  for(int n = 0; n < N1; n++){
    
    // propose point until accepted
    arma::rowvec old_xi = data1.row(n);
    arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
    arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
    
    bool pointOutside = !checkPointIsInRegion(xi, a1, b1, a2, b2);
    
    while(pointOutside){
      xi = mvrnormArma(oldxivec, Sigma_prop);
      pointOutside = !checkPointIsInRegion(xi, a1, b1, a2, b2);
    }
    
    if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
      
      data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
      
      double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move_cond(n, data1, data2,
                                                                            N1, N2,
                                                                            oldxivec,
                                                                            theta1, theta12);
      
      if(R::runif(0, 1) < exp(log_hastings_ratio_num)){
        data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
      } else {
        data1.row(n) = old_xi;
      }
      // Rcout << exp(log_hastings_ratio_num) << std::endl;
    }
    
  }
  
  return(data1);
}

// [[Rcpp::export]]
List simulate_conditional_bivsoftcore_cpp(int N1, int N2, 
                                          double theta1, double theta2, double theta12,
                                          int niter, arma::mat Sigma_prop,
                                          double a1, double b1,
                                          double a2, double b2){
  
  arma::mat data1;
  data1.zeros(N1, 2);
  for(int i = 0; i < N1; i++){
    data1.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPointPP(a1, b1, a2, b2));
  }
  
  arma::mat data2;
  data2.zeros(N2, 2);
  for(int i = 0; i < N2; i++){
    data2.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPointPP(a1, b1, a2, b2));
  }
  
  for(int iter = 0; iter < niter; iter++){
    
    // Rcout << iter << std::endl;
    
    // propose x1
    
    data1 = proposePointProcess(data1, data2, N1, N2, theta1, theta12, 
                                a1, b1, a2, b2, Sigma_prop);
    
    // for(int n = 0; n < N1; n++){
    //   
    //   // propose point until accepted
    //   arma::rowvec old_xi = data1.row(n);
    //   arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
    //   arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
    //   
    //   bool pointOutside = !checkPointIsInRegion(xi, a, b);
    //   
    //   while(pointOutside){
    //     xi = mvrnormArma(oldxivec, Sigma_prop);
    //     pointOutside = !checkPointIsInRegion(xi, a, b);
    //   }
    //   
    //   if(checkPointIsInRegion(xi, a, b)){
    //     
    //     data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
    //     
    //     double log_hastings_ratio_num = log_f_bivtan_cpp_quick_move(n, 1, 
    //                                                                      data1, data2,
    //                                                                      N1, N2,
    //                                                                      oldxivec,
    //                                                                      theta1, theta2, theta3);
    //     
    //     if(R::runif(0, 1) < exp(log_hastings_ratio_num)){
    //       data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
    //     } else {
    //       data1.row(n) = old_xi;
    //     }
    //     // Rcout << exp(log_hastings_ratio_num) << std::endl;
    //   }
    //   
    // }
    
    // propose x2
    
    // for(int n = 0; n < N2; n++){
    //   
    //   arma::rowvec old_xi = data2.row(n);
    //   arma::rowvec xi(2);
    //   xi[0] = R::rnorm(data2(n,0), sigma_prop);
    //   xi[1] =  R::rnorm(data2(n,1), sigma_prop);
    //   
    //   if(checkPointIsInRegion(xi)){
    //     
    //     data2.row(n) = xi;
    //     
    //     arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
    //     double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick(n, 2, log_hastings_ratio_den,
    //                                                                 data1, data2, oldxivec, 
    //                                                                 theta1, theta2, theta3);
    //     
    //     if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
    //       data2.row(n) = xi;
    //       log_hastings_ratio_den = log_hastings_ratio_num;
    //     } else {
    //       data2.row(n) = old_xi;
    //     }
    //     
    //   }
    //   
    // }
    
    data2 = proposePointProcess(data2, data1, N2, N1, theta2, theta12, 
                                a1, b1, a2, b2, Sigma_prop);
    
    // for(int n = 0; n < N2; n++){
    //   
    //   arma::rowvec old_xi = data2.row(n);
    //   arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
    //   arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
    //   
    //   bool pointOutside = !checkPointIsInRegion(xi, a, b);
    //   
    //   while(pointOutside){
    //     xi = mvrnormArma(oldxivec, Sigma_prop);
    //     pointOutside = !checkPointIsInRegion(xi, a, b);
    //   }
    //   
    //   if(checkPointIsInRegion(xi, a, b)){
    //     
    //     data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
    //     
    //     arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
    //     double log_hastings_ratio_num = log_f_bivtan_cpp_quick_move(n, 2, 
    //                                                                 data1, data2, N1, N2,
    //                                                                 oldxivec, 
    //                                                                 theta1, theta2, theta3);
    //     
    //     if(R::runif(0, 1) < exp(log_hastings_ratio_num)){
    //       data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
    //     } else {
    //       data2.row(n) = old_xi;
    //     }
    //     
    //   }
    //   
    // }
    
  }
  
  return(List::create(_["data1"] = data1,
                      _["data2"] = data2));
}

// // [[Rcpp::export]]
// arma::vec proposeNewPoint_old(arma::mat Sigma_prop){
//   
//   arma::vec zeroVec = arma::zeros(2);
//   
//   arma::vec proposedPoint = mvrnormArma(zeroVec, Sigma_prop);
//   
//   while(!checkPointIsInRegionNew(proposedPoint)){
//     proposedPoint = mvrnormArma(zeroVec, Sigma_prop);
//   }
//   
//   return(proposedPoint);
// }

// // [[Rcpp::export]]
// List simulate_homogeneous_PP(double beta1, int niter, arma::mat Sigma_prop, 
//                              int Nmax, double lambda,
//                              arma::mat Sigma_newpoint,
//                              double a1, double b1, double a2, double b2){
//   
//   arma::vec loglikelihoods = arma::zeros(niter);
//   
//   double q = .33333;
//   double p = .5;
//   
//   arma::mat data1;
//   data1.zeros(Nmax, 2);
//   
//   int N1 = R::rpois(lambda);
//   for(int i = 0; i < N1; i++){
//     data1.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPointPP(a1, b1, a2, b2));
//   }
//   
//   arma::mat data2;
//   data2.zeros(Nmax, 2);
//   
//   int N2 = R::rpois(lambda);
//   for(int i = 0; i < N2; i++){
//     data2.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPoint(a1, b1, a2, b2, traps, R));
//   }
//   
//   double log_hastings_ratio_den = log_f_bivsoftcore_cpp(data1, data2, N1, N2,
//                                                         beta1, beta2,
//                                                         theta1, theta2, theta3);
//   
//   for(int iter = 0; iter < niter; iter++){
//     
//     // x1
//     
//     if(R::runif(0,1) < q){
//       
//       // move
//       
//       // sample point
//       int n = sample_int(N1) - 1;
//       
//       arma::rowvec old_xi = data1.row(n);
//       arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
//       arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
//       
//       if(checkPointIsInRegionTraps(xi, traps, R)){
//         
//         data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 1, log_hastings_ratio_den,
//                                                                          data1, data2,
//                                                                          N1, N2,
//                                                                          oldxivec,
//                                                                          theta1, theta2, theta3);
//         
//         if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//           data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//         } else {
//           data1.row(n) = old_xi;
//         }
//         
//       }
//       
//     } else {
//       
//       if(R::runif(0,1) < p){
//         
//         // birth
//         
//         arma::vec xi = proposeNewPoint(a1, b1, a2, b2, traps, R);
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 1, log_hastings_ratio_den,
//                                                                           data1, data2,
//                                                                           N1, N2,
//                                                                           beta1, beta2,
//                                                                           theta1, theta2, theta3);
//         
//         if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//           data1.row(N1) = arma::conv_to<arma::rowvec>::from(xi);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//           N1 += 1;
//         }
//         
//       } else {
//         
//         // death
//         
//         if (N1 > 0){
//           
//           // sample point
//           int itemToRemove = sample_int(N1) - 1;
//           
//           double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
//                                                                             1, log_hastings_ratio_den,
//                                                                             data1, data2,
//                                                                             N1, N2,
//                                                                             beta1, beta2,
//                                                                             theta1, theta2, theta3);
//           
//           if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//             data1.row(itemToRemove) = data1.row(N1 - 1);
//             log_hastings_ratio_den = log_hastings_ratio_num;
//             N1 -= 1;
//           }
//           
//         }
//         
//       }
//       
//     }
//     
//     // x2
//     
//     if(R::runif(0,1) < q){ 
//       
//       // move
//       
//       // sample point
//       int n = sample_int(N2) - 1;
//       
//       arma::rowvec old_xi = data2.row(n);
//       arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
//       arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
//       
//       if(checkPointIsInRegionTraps(xi, traps, R)){
//         
//         data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
//         
//         arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 2, log_hastings_ratio_den,
//                                                                          data1, data2, N1, N2,
//                                                                          oldxivec, 
//                                                                          theta1, theta2, theta3);
//         
//         if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//           data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//         } else {
//           data2.row(n) = old_xi;
//         }
//         
//       }
//       
//     } else {
//       
//       if(R::runif(0,1) < p){
//         
//         // birth
//         
//         arma::vec xi = proposeNewPoint(a1, b1, a2, b2, traps, R);
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 2, log_hastings_ratio_den,
//                                                                           data1, data2, N1, N2,
//                                                                           beta1, beta2,
//                                                                           theta1, theta2, theta3);
//         
//         if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//           data2.row(N2) = arma::conv_to<arma::rowvec>::from(xi);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//           N2 += 1;
//         }
//         
//       } else {
//         
//         // death
//         
//         if(N2 > 0){
//           
//           // sample point
//           int itemToRemove = sample_int(N2) - 1;
//           
//           double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
//                                                                             2, log_hastings_ratio_den,
//                                                                             data1, data2,
//                                                                             N1, N2,
//                                                                             beta1, beta2,
//                                                                             theta1, theta2, theta3);
//           
//           if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
//             data2.row(itemToRemove) = data2.row(N2 - 1);
//             log_hastings_ratio_den = log_hastings_ratio_num;
//             N2 -= 1;
//           }
//           
//         }
//         
//       }
//       
//     }
//     
//     loglikelihoods[iter] = log_hastings_ratio_den;
//     
//   }  
//   
//   return(List::create(_["data1"] = data1,
//                       _["data2"] = data2,
//                       _["loglikelihoods"] = loglikelihoods,
//                       _["N1"] = N1,
//                       _["N2"] = N2));
// }

/////////////// MODEL FUNCTIONS

// [[Rcpp::export]]
double loglikelihood_xi_captured(arma::vec si, int i, int K, arma::vec S_usage,
                                 arma::vec trapsX, arma::vec trapsY,
                                 arma::mat CH, double p0, double sigma){
  
  double loglikelihood = 0;
  
  for (int j = 0; j < K; j++) {
    double p = p0 * exp(- (1 / (2 * sigma * sigma)) * (pow(si[0] - trapsX[j],2) + pow(si[1] - trapsY[j],2)));
    loglikelihood += R::dbinom(CH(i,j), S_usage[j], p, 1);
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double loglikelihood_xi_uncaptured(arma::vec si, int K, arma::vec S_usage,
                                   arma::vec trapsX, arma::vec trapsY,
                                   double p0, double sigma){
  
  double loglikelihood = 0;
  
  for (int j = 0; j < K; j++) {
    double p = p0 * exp(- (1 / (2 * sigma * sigma)) * (pow(si[0] - trapsX[j],2) + pow(si[1] - trapsY[j],2)));
    loglikelihood += R::dbinom(0, S_usage[j], p, 1);
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double loglikelihood_p_cpp(double p0, double sigma,
                           arma::mat CH, arma::mat s_n,
                           arma::vec trapsX, arma::vec trapsY,
                           arma::vec S_usage, int N){
  
  double loglikelihood = 0;
  
  int K = CH.n_cols;
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < K; j++){
      double p = p0 * exp(- (1 / (2 * sigma * sigma)) * (pow(s_n(i,0) - trapsX[j], 2) + pow(s_n(i,1) - trapsY[j],2)));
      
      loglikelihood += R::dbinom(CH(i,j), S_usage[j], p, 1);
      
    }
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double probAcceptingIndividual(arma::vec s_new, double p0, double sigma,
                               int D, int N, int K,
                               arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double p_notcaptured = exp(loglikelihood_xi_uncaptured(s_new, K, S_usage, trapsX, trapsY, p0, sigma));
  double prob = (p_notcaptured * (N + 1)) / (N + 1 - D);
  
  return(prob);
}

// [[Rcpp::export]]
double probProposingIndividual(arma::vec point, double probMixture,
                               double a1, double b1, double a2, double b2,
                               arma::mat mixtureMeans, arma::mat mixtureSd, arma::vec norm_const){
  
  double likelihood = 0;
  
  int numCenters = mixtureMeans.n_rows;
  
  likelihood += (1 - probMixture) * (R::dunif(point[0], a1, b1, 0)) *
    (R::dunif(point[1], a2, b2, 0)) / norm_const[0];
  Rcout << likelihood << std::endl;
  for(int l = 0; l < numCenters; l++){
    
    likelihood += (probMixture) * (1.0 / numCenters) * (R::dnorm(point[0], mixtureMeans(l,0), mixtureSd(l,0), 0) *
      R::dnorm(point[1], mixtureMeans(l,1), mixtureSd(l,1), 0)) / norm_const[l + 1];
    
  }
  Rcout << likelihood << std::endl;
  return(likelihood);
}

// [[Rcpp::export]]
List update_N_withacceptance(arma::mat s, double beta, double theta,
                             arma::mat CH, int N, int D,
                             double p0, double sigma,
                             arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
                             arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                             double a1, double b1, double a2, double b2){
  
  int K = CH.n_cols;
  
  bool isPointAccepted = false;
  arma::vec acceptedPoint = arma::zeros(2);
  
  if(R::runif(0,1) < .5){
    
    // birth
    
    arma::vec x1_new = proposeNewPointPP(a1, b1, a2, b2);
    
    
    double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(x1_new,
                                                                      s, N,
                                                                      beta, theta);
    
    double loglikelihood_xi_star = log(probAcceptingIndividual(x1_new, p0, sigma, D, N,
                                                               K, S_usage, trapsX, trapsY));
    
    double logNumerator = log_hastings_ratio_num + loglikelihood_xi_star;
    // double logposterior_xi = 0;
    
    // double logAddition = log(beta) - log(N + 1);
    
    // double loglikelihood_xi_star = log(probAcceptingIndividual(x1_new, p0, sigma, D, N,
    // K, S_usage, trapsX, trapsY));
    
    // double logNumerator = logAddition + loglikelihood_xi_star;
    double logDenominator = 0;
    
    if(R::runif(0,1) < exp(logNumerator - logDenominator)){
      
      s.row(N) = arma::conv_to<arma::rowvec>::from(x1_new);
      N += 1;
      isPointAccepted = true;
      acceptedPoint = x1_new;
    }
    
  } else {
    
    // death
    
    if(N > D){
      int itemToRemove = sample_int(N - D) + D - 1;
      
      arma::vec xi = arma::conv_to<arma::vec>::from(s.row(itemToRemove));
      
      double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove, s, N,
                                                                        beta, theta);
      
      double loglikelihood_xi_star = log(probAcceptingIndividual(xi, p0, sigma, D, N - 1, 
                                                                 K, S_usage, trapsX, trapsY));
      
      // double logdeletion = - log(beta) + log(N);
      // 
      // double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
      //                                                                   1, log_hastings_ratio_den,
      //                                                                   s1, s2,
      //                                                                   N1, N2,
      //                                                                   beta1, beta2,
      //                                                                   theta1, theta2, theta3);
      // 
      // double loglikelihood_xi_star = log(probAcceptingIndividual(xi, p0, sigma, D, N - 1,
      //                                                            K, S_usage, trapsX, trapsY));
      // 
      // double logNumerator = logdeletion;
      // double logDenominator = loglikelihood_xi_star;
      
      double logNumerator = log_hastings_ratio_num;
      double logDenominator = loglikelihood_xi_star;
      
      if(R::runif(0,1) < exp(logNumerator - logDenominator)){
        s.row(itemToRemove) = s.row(N - 1);
        N -= 1;
      }
    }
    
  }
  
  return(List::create(_["s"] = s,
                      _["N"] = N,
                      _["isPointAccepted"] = isPointAccepted,
                      _["acceptedPoint"] = acceptedPoint));
}

// [[Rcpp::export]]
double densProposal(arma::vec x, double normConst){
  
  double prob = 1 / normConst;
  
  return(prob);
}

// [[Rcpp::export]]
List simulate_CH_from_p_cpp(double p0_1, double sigma_1,
                            double p0_2, double sigma_2,
                            arma::mat s1, int N1, arma::mat s2, int N2,
                            arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  int K = trapsX.size();
  
  arma::mat CH_1 = arma::zeros(N1, K);
  
  for(int i = 0; i < N1; i++){
    for(int k = 0; k < K; k++){
      double p = p0_1 * exp(- (1 / (2 * sigma_1 * sigma_1)) * (pow(s1(i,0) - trapsX[k], 2) + pow(s1(i,1) - trapsY[k], 2)));
      CH_1(i, k) = R::rbinom(S_usage[k], p);
    }
  }
  
  arma::mat CH_2 = arma::zeros(N2, K);
  
  for(int i = 0; i < N2; i++){
    for(int k = 0; k < K; k++){
      double p = p0_2 * exp(- (1 / (2 * sigma_2 * sigma_2)) * (pow(s2(i,0) - trapsX[k], 2) + pow(s2(i,1) - trapsY[k], 2)));
      CH_2(i, k) = R::rbinom(S_usage[k], p);
    }
  }
  
  return(List::create(_["CH_1"] = CH_1,
                      _["CH_2"] = CH_2));
}

// NEW PROPOSAL FUNCTION (DEFINITIVE)

// [[Rcpp::export]]
double logDensityBirth(arma::vec point, double probMixture,
                       double a1, double b1, double a2, double b2,
                       arma::mat mixtureMeans, arma::mat mixtureSd, arma::vec norm_const){
  
  double likelihood = 0;
  
  int numCenters = mixtureMeans.n_rows;
  
  likelihood += (1 - probMixture) * (R::dunif(point[0], a1, b1, 0)) *
    (R::dunif(point[1], a2, b2, 0)) / norm_const[0];
  
  for(int l = 0; l < numCenters; l++){
    
    likelihood += (probMixture) * (1.0 / numCenters) * (R::dnorm(point[0], mixtureMeans(l,0), mixtureSd(l,0), 0) *
      R::dnorm(point[1], mixtureMeans(l,1), mixtureSd(l,1), 0)) / norm_const[l + 1];
    
  }
  
  return(log(likelihood));
}

// [[Rcpp::export]]
arma::vec proposeNewPointNewMixture(double probMixture, arma::mat mixtureMeans, arma::mat mixtureSd,
                                    double a1, double b1, double a2, double b2){
  arma::vec proposedPoint;
  
  if(R::runif(0,1) < (1 - probMixture)){
    
    proposedPoint = proposeNewPointPP(a1, b1, a2, b2);
    
  } else {
    
    int numCenters = mixtureMeans.n_rows;
    
    proposedPoint = arma::zeros(2);
    
    int indMixture = sample_int(numCenters) - 1;
    
    proposedPoint[0] = R::rnorm(mixtureMeans(indMixture, 0), mixtureSd(indMixture, 0));
    proposedPoint[1] = R::rnorm(mixtureMeans(indMixture, 1), mixtureSd(indMixture, 1));
    
    while(!checkPointIsInRegion(proposedPoint, a1, b1, a2, b2)){
      
      proposedPoint[0] = R::rnorm(mixtureMeans(indMixture, 0), mixtureSd(indMixture, 0));
      proposedPoint[1] = R::rnorm(mixtureMeans(indMixture, 1), mixtureSd(indMixture, 1));
    }
    
  }
  
  return proposedPoint;
  
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_birth_intonly(arma::vec xnew,
                                                 arma::mat x, int N,
                                                 double theta){
  
  double loglikelihood = 0;
  
  // between
  
  for(int j = 0; j < N; j++){
    
    double r_new = sqrt(pow(x(j,0) - xnew(0), 2) + pow(x(j,1) - xnew(1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double loglikelihood_addition(arma::mat x_new, int numProposedPoints, 
                              arma::mat x, int N, 
                              double beta, double theta){
  
  double loglikelihood = 0;
  
  loglikelihood += numProposedPoints * log(beta);
  
  // with existing points
  
  for(int l = 0; l < numProposedPoints; l++){
    arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
    loglikelihood += log_f_bivsoftcore_cpp_quick_birth_intonly(x_current, x, N,
                                                               theta);
  }
  
  // between new points
  
  for(int i = 0; i < numProposedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x_new(i,0) - x_new(j,0), 2) + pow(x_new(i,1) - x_new(j,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double loglikelihood_removal(IntegerVector idxPointsToRemove, int numRemovedPoints, 
                             arma::mat x, int N, double beta, double theta){
  
  idxPointsToRemove.sort(false);
  
  double loglikelihood = 0;
  
  loglikelihood -= numRemovedPoints * log(beta);
  
  // with points not removed
  
  for(int j = 0; j < numRemovedPoints; j++){
    
    int l = 0; // index of removed points
    
    for(int i = 0; i < N; i++){
      
      if(i != idxPointsToRemove[l]){
        
        double r_new = sqrt(pow(x(idxPointsToRemove[j],0) - x(i,0), 2) + 
                            pow(x(idxPointsToRemove[j],1) - x(i,1), 2));
        
        loglikelihood -= log(1 - exp(- r_new*r_new / (theta)));
        
      } else {
        
        l += 1;
        
      }
      
    }
    
  }
  
  // between them
  
  for(int i = 0; i < numRemovedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x(idxPointsToRemove[i],0) - x(idxPointsToRemove[j],0), 2) + 
                          pow(x(idxPointsToRemove[i],1) - x(idxPointsToRemove[j],1), 2));
      
      loglikelihood -= log(1 - exp(- r_new*r_new / (theta)));
      
    }
    
  }
  
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double probIndividualNotCaptured(arma::vec s_new, double p0, double sigma,
                                 int D, int N, int K,
                                 arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double p_notcaptured = exp(loglikelihood_xi_uncaptured(s_new, K, S_usage, trapsX, trapsY, p0, sigma));
  double prob = (p_notcaptured * (N + 1)) / (N + 1 - D);
  
  return(prob);
}

// [[Rcpp::export]]
double loglikelihood_move(int i, arma::mat x, int N, arma::vec xold, 
                          double theta){
  
  double loglikelihood = 0;
  
  // between
  
  for(int j = 0; j < N; j++){
    
    if(j != i){
      
      double r_old = sqrt(pow(x(j,0) - xold(0), 2) + pow(x(j,1) - xold(1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta)));
      
      double r_new = sqrt(pow(x(j,0) - x(i,0), 2) + pow(x(j,1) - x(i,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double computePointsMoveRatio(int i, arma::vec old_xi, arma::mat s, int D, 
                              int N, arma::vec xi, 
                              int K, arma::vec S_usage,
                              arma::mat trapsX, arma::mat trapsY, arma::mat CH,
                              double p0, double sigma, double theta){
  
  double log_move = loglikelihood_move(i, s, N, old_xi, theta);
  
  double loglikelihood_xi_star;
  double loglikelihood_xi_old;
  
  if(i < D){
    loglikelihood_xi_star = loglikelihood_xi_captured(xi, i, K, S_usage, trapsX, trapsY, CH, p0, sigma);
    loglikelihood_xi_old = loglikelihood_xi_captured(old_xi, i, K, S_usage, trapsX, trapsY, CH, p0, sigma);
  } else {
    loglikelihood_xi_star = loglikelihood_xi_uncaptured(xi, K, S_usage, trapsX, trapsY, p0, sigma);
    loglikelihood_xi_old = loglikelihood_xi_uncaptured(old_xi, K, S_usage, trapsX, trapsY, p0, sigma);
  }
  
  double logNumerator = log_move + loglikelihood_xi_star;
  double logDenominator = loglikelihood_xi_old;
  
  return(exp(logNumerator - logDenominator));
}

// [[Rcpp::export]]
double computePointsAdditionRatio(arma::mat x_new, int ptsToPropose, int N,
                                  arma::mat s, double beta, double theta,
                                  double probMixture, arma::mat mixtureMeans,
                                  arma::mat mixtureSd, arma::vec normConst,
                                  double a1, double b1, double a2, double b2, double S_area,
                                  int D, int K, double p0, double sigma,
                                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double logRatio_f = loglikelihood_addition(x_new, ptsToPropose, s, N, 
                                             beta / S_area, theta);
  
  double logbirthDensity = 0;
  for(int l = 0; l < ptsToPropose; l++){
    arma::vec x1_new = arma::conv_to<arma::vec>::from(x_new.row(l));
    logbirthDensity += logDensityBirth(x1_new, probMixture,
                                       a1, b1, a2, b2,
                                       mixtureMeans, mixtureSd, normConst);
  }
  
  double logDeathDensity = 0;
  for(int l = 0; l < ptsToPropose; l++){
    logDeathDensity -= log(N + 1 + l);
  }
  
  double loglikelihood_xi_star = 0;
  for(int l = 0; l < ptsToPropose; l++){
    arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
    loglikelihood_xi_star += log(probIndividualNotCaptured(x_current, p0, sigma, D, N + l,
                                                           K, S_usage, trapsX, trapsY));
  }
  
  double logNumerator = logRatio_f + loglikelihood_xi_star + logDeathDensity;
  double logDenominator = logbirthDensity;
  
  return(exp(logNumerator - logDenominator));
}

// // [[Rcpp::export]]
// double computePointsAdditionRatioSimple(arma::mat x_new, int ptsToPropose, int N,
//                                         arma::mat s, double beta,
//                                         double probMixture, arma::mat mixtureMeans,
//                                         arma::mat mixtureSd, arma::vec normConst,
//                                         double a1, double b1, double a2, double b2, double S_area,
//                                         int D, int K, double p0, double sigma,
//                                         arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
//   
//   double logRatio_f = loglikelihood_addition(ptsToPropose, beta);
//   
//   double logbirthDensity = 0;
//   // for(int l = 0; l < ptsToPropose; l++){
//   //   arma::vec x1_new = arma::conv_to<arma::vec>::from(x_new.row(l));
//   //   logbirthDensity += logDensityBirth(x1_new, probMixture,
//   //                                      a1, b1, a2, b2,
//   //                                      mixtureMeans, mixtureSd, normConst);
//   // }
//   
//   double logDeathDensity = 0;
//   for(int l = 0; l < ptsToPropose; l++){
//     logDeathDensity -= log(N + 1 + l);
//   }
//   
//   double loglikelihood_xi_star = 0;
//   for(int l = 0; l < ptsToPropose; l++){
//     arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
//     loglikelihood_xi_star += log(probIndividualNotCaptured(x_current, p0, sigma, D, N + l,
//                                                            K, S_usage, trapsX, trapsY));
//   }
//   
//   double logNumerator = logRatio_f + loglikelihood_xi_star + logDeathDensity;
//   double logDenominator = logbirthDensity;
//   
//   return(exp(logNumerator - logDenominator));
// }

// [[Rcpp::export]]
double computePointsDeletionRatio(IntegerVector idxPointsToRemove, int pointsToRemove, int N, arma::mat s,
                                  double beta, double theta, 
                                  double probMixture, arma::mat mixtureMeans,
                                  arma::mat mixtureSd, arma::vec normConst,
                                  double a1, double b1, double a2, double b2, double S_area,
                                  int D, int K, double p0, double sigma,
                                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double logRatio_f = loglikelihood_removal(idxPointsToRemove, pointsToRemove, s,
                                            N, beta / S_area, theta);
  
  double logbirthDensity = 0;
  for(int l = 0; l < pointsToRemove; l++){
    
    arma::vec xi = arma::conv_to<arma::vec>::from(s.row(idxPointsToRemove[l]));
    
    logbirthDensity += logDensityBirth(xi, probMixture,
                                       a1, b1, a2, b2,
                                       mixtureMeans, mixtureSd, normConst);
  }
  
  double logDeathDensity = 0;
  for(int l = 0; l < pointsToRemove; l++){
    logDeathDensity -= log(N - l);
  }
  
  double loglikelihood_xi_star = 0;
  for(int l = 0; l < pointsToRemove; l++){
    
    arma::vec xi = arma::conv_to<arma::vec>::from(s.row(idxPointsToRemove[l]));
    
    loglikelihood_xi_star += log(probIndividualNotCaptured(xi, p0, sigma, D, N - l - 1,
                                                           K, S_usage, trapsX, trapsY));
  }
  
  double logNumerator = logRatio_f + logbirthDensity;
  double logDenominator = logDeathDensity + loglikelihood_xi_star;
  
  return(exp(logNumerator - logDenominator));
}

// // [[Rcpp::export]]
// double computePointsDeletionRatioSimple(IntegerVector idxPointsToRemove, int pointsToRemove, int N, arma::mat s,
//                                         double beta, double probMixture, arma::mat mixtureMeans,
//                                         arma::mat mixtureSd, arma::vec normConst,
//                                         double a1, double b1, double a2, double b2, double S_area,
//                                         int D, int K, double p0, double sigma,
//                                         arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
//   
//   double logRatio_f = loglikelihood_removal(pointsToRemove, beta);
//   
//   double logbirthDensity = 0;
//   // for(int l = 0; l < pointsToRemove; l++){
//   //   
//   //   arma::vec xi = arma::conv_to<arma::vec>::from(s.row(idxPointsToRemove[l]));
//   //   
//   //   logbirthDensity += logDensityBirth(xi, probMixture,
//   //                                      a1, b1, a2, b2,
//   //                                      mixtureMeans, mixtureSd, normConst);
//   // }
//   
//   double logDeathDensity = 0;
//   for(int l = 0; l < pointsToRemove; l++){
//     logDeathDensity -= log(N - l);
//   }
//   
//   double loglikelihood_xi_star = 0;
//   for(int l = 0; l < pointsToRemove; l++){
//     
//     arma::vec xi = arma::conv_to<arma::vec>::from(s.row(idxPointsToRemove[l]));
//     
//     loglikelihood_xi_star += log(probIndividualNotCaptured(xi, p0, sigma, D, N - l - 1,
//                                                            K, S_usage, trapsX, trapsY));
//   }
//   
//   double logNumerator = logRatio_f + logbirthDensity;
//   double logDenominator = logDeathDensity + loglikelihood_xi_star;
//   
//   return(exp(logNumerator - logDenominator));
// }

// [[Rcpp::export]]
arma::mat update_s(arma::mat s, double beta, double theta,
                   arma::mat CH,
                   int N, int D,
                   double p0, double sigma, 
                   double lambda_newPoints,
                   arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                   arma::mat traps, //double R, double S_area,
                   arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
                   double a1, double b1, double a2, double b2){
  
  int K = CH.n_cols;
  
  int pointsToMove = R::rpois(lambda_newPoints);
  if(pointsToMove > N) pointsToMove = N;
  IntegerVector seqND = seq(1, N);
  
  IntegerVector idxPointsToMove = RcppArmadillo::sample(seqND, pointsToMove, 0) - 1;
  
  for(int j = 0; j < pointsToMove; j++) {
    // for(int i = 0; i < N1; i++) {
    
    int i = idxPointsToMove[j];
    
    arma::vec old_xi = arma::conv_to<arma::vec>::from(s.row(i));
    arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
    
    if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
      
      s.row(i) = arma::conv_to<arma::rowvec>::from(xi);
      
      double mh_ratio = computePointsMoveRatio(i, old_xi, s,
                                               D, N, xi, K, S_usage,
                                               trapsX, trapsY, CH,
                                               p0, sigma, theta);
      
      if(R::runif(0,1) > mh_ratio){
        s.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
      }
    }
    
  }
  
  return(s);
}

// [[Rcpp::export]]
arma::mat update_s_da(arma::mat s, double beta, double theta,
                      arma::vec z, arma::mat CH, int N, int D,
                      double p0, double sigma,
                      arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                      arma::mat traps, double R, double S_area,
                      arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
                      double a1, double b1, double a2, double b2){
  
  int K = CH.n_cols;
  
  int M = z.size();
  
  for(int i = 0; i < M; i++) {
    
    if(z[i] == 1){
      
      arma::vec old_xi = arma::conv_to<arma::vec>::from(s.row(i));
      arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
      
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        s.row(i) = arma::conv_to<arma::rowvec>::from(xi);
        
        double mh_ratio = computePointsMoveRatio(i, old_xi, s,
                                                 D, N, xi, K, S_usage,
                                                 trapsX, trapsY, CH,
                                                 p0, sigma, theta);
        
        if(R::runif(0,1) > mh_ratio){
          s.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
        }
      }
      
    } else {
      
      arma::vec proposedPoint = proposeNewPointPP(a1, b1, a2, b2);
      
      s.row(i) = arma::conv_to<arma::rowvec>::from(proposedPoint);
      
    }
    
    
    
  }
  
  return(s);
}

// [[Rcpp::export]]
List update_N(arma::mat s, double beta, double theta,
              arma::mat CH, int N, int D,
              double p0, double sigma,
              arma::mat Sigma_prop, arma::mat Sigma_newpoint,
              double lambda_newPoints, double probMixture,
              arma::mat mixtureMeans, arma::mat mixtureSd,
              arma::vec normConst, double S_area,
              arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
              double a1, double b1, double a2, double b2){
  
  int K = CH.n_cols;
  
  double mh_ratio_out;
  
  if(R::runif(0,1) < .5){
    
    // birth
    
    int pointsToPropose = 1 +  R::rpois(lambda_newPoints);
    
    arma::mat x_new = arma::zeros(pointsToPropose, 2);
    for(int l = 0; l < pointsToPropose; l++){
      arma::vec x1_new = proposeNewPointPP(a1, b1, a2, b2);
      // arma::vec x1_new = proposeNewPointNewMixture(probMixture, mixtureMeans, mixtureSd,
      //                                              a1, b1, a2, b2);
      
      x_new.row(l) = arma::conv_to<arma::rowvec>::from(x1_new);
    }
    
    double mh_ratio = computePointsAdditionRatio(x_new, pointsToPropose, N,
                                                 s, beta, theta, probMixture, mixtureMeans,
                                                 mixtureSd, normConst,
                                                 a1, b1, a2, b2, S_area,
                                                 D, K, p0, sigma,
                                                 S_usage, trapsX, trapsY);
    // Rcout << mh_ratio << std::endl;
    if(R::runif(0,1) < mh_ratio){
      
      for(int l = 0; l < pointsToPropose; l++){
        s.row(N + l) = x_new.row(l);
      }
      N += pointsToPropose;
      
    }
    
    mh_ratio_out = mh_ratio;
    
  } else {
    
    // death
    
    int pointsToRemove = 1 +  R::rpois(lambda_newPoints);
    
    if(N - pointsToRemove >= D){
      
      IntegerVector seqND = seq(1, N - D);
      
      IntegerVector idxPointsToRemove = D - 1 + RcppArmadillo::sample(seqND, pointsToRemove, 0);
      
      double mhRatio = computePointsDeletionRatio(idxPointsToRemove, pointsToRemove, N,
                                                  s, beta, theta, probMixture,
                                                  mixtureMeans, mixtureSd, normConst,
                                                  a1, b1, a2, b2, S_area,
                                                  D, K, p0, sigma,
                                                  S_usage, trapsX, trapsY);
      
      if(R::runif(0,1) < mhRatio){
        
        for(int l = 0; l < pointsToRemove; l++){
          s.row(idxPointsToRemove[l]) = s.row(N - 1);
          N -= 1;
        }
        
      }
      
      mh_ratio_out = mhRatio;
      
    }
    
  }
  
  return(List::create(_["s"] = s,
                      _["N"] = N,
                      _["mh_ratio"] = mh_ratio_out));
}

// [[Rcpp::export]]
List update_N_all(arma::mat s, double beta, double theta,
                  arma::mat CH, int N, int D,
                  double p0, double sigma, 
                  arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                  double lambda_newPoints, double probMixture,
                  arma::mat mixtureMeans, arma::mat mixtureSd,
                  arma::vec normConst, double S_area,
                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
                  double a1, double b1, double a2, double b2){
  
  int K = CH.n_cols;
  
  double mh_ratio_out;
  
  double p = .000000333333333;
  double q = .5;
  
  if(R::runif(0,1) < p){
    
    for(int i = 0; i < N; i++) {
      
      arma::vec old_xi = arma::conv_to<arma::vec>::from(s.row(i));
      arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
      
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        s.row(i) = arma::conv_to<arma::rowvec>::from(xi);
        
        double mh_ratio = computePointsMoveRatio(i, old_xi, s,
                                                 D, N, xi, K, S_usage,
                                                 trapsX, trapsY, CH,
                                                 p0, sigma, theta);
        
        if(R::runif(0,1) > mh_ratio){
          s.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
        }
      }
      
    }
    
    mh_ratio_out = -1;
    
  } else {
    
    if(R::runif(0,1) < q){
      
      // birth
      
      int pointsToPropose = 1 +  R::rpois(lambda_newPoints);
      
      arma::mat x_new = arma::zeros(pointsToPropose, 2);
      for(int l = 0; l < pointsToPropose; l++){
        arma::vec x1_new = proposeNewPointNewMixture(probMixture, mixtureMeans, mixtureSd,
                                                     a1, b1, a2, b2);
        
        x_new.row(l) = arma::conv_to<arma::rowvec>::from(x1_new);
      }
      
      double mh_ratio = computePointsAdditionRatio(x_new, pointsToPropose, N,
                                                   s, beta, theta, probMixture, mixtureMeans,
                                                   mixtureSd, normConst,
                                                   a1, b1, a2, b2, S_area,
                                                   D, K, p0, sigma,
                                                   S_usage, trapsX, trapsY);
      // Rcout << mh_ratio << std::endl;
      if(R::runif(0,1) < mh_ratio){
        
        for(int l = 0; l < pointsToPropose; l++){
          s.row(N + l) = x_new.row(l);
        }
        N += pointsToPropose;
        
      }
      
      mh_ratio_out = mh_ratio;
      
    } else {
      
      // death
      
      int pointsToRemove = 1 +  R::rpois(lambda_newPoints);
      
      if(N - pointsToRemove >= D){
        
        IntegerVector seqND = seq(1, N - D);
        
        IntegerVector idxPointsToRemove = D - 1 + RcppArmadillo::sample(seqND, pointsToRemove, 0);
        
        double mhRatio = computePointsDeletionRatio(idxPointsToRemove, pointsToRemove, N,
                                                    s, beta, theta, probMixture,
                                                    mixtureMeans, mixtureSd, normConst,
                                                    a1, b1, a2, b2, S_area,
                                                    D, K, p0, sigma,
                                                    S_usage, trapsX, trapsY);
        Rcout << mhRatio << std::endl;
        if(R::runif(0,1) < mhRatio){
          
          for(int l = 0; l < pointsToRemove; l++){
            s.row(idxPointsToRemove[l]) = s.row(N - 1);
            N -= 1;
          }
          
        }
        
        mh_ratio_out = mhRatio;
        
      }
      
    }
    
  }
  
  return(List::create(_["s"] = s,
                      _["N"] = N,
                      _["mh_ratio"] = mh_ratio_out));
}

// [[Rcpp::export]]
List update_z_cpp(arma::vec z, arma::mat s, int N, int D, double psi, double theta,
                  int K, arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
                  double p0, double sigma){
  
  int M = z.size();
  
  for(int j = D; j < M; j++){
    
    if(z[j] == 1){
      N -= 1;
      z[j] = 0;
    }
    
    arma::vec s_j = arma::conv_to<arma::vec>::from(s.row(j));
    
    double log_probAddition = 0;
    
    for(int i = 0; i < M; i++){
      
      if(z[i] == 1){
        
        if(i == j){
          
          double r_new = sqrt(pow(s(j,0) - s(i,0), 2) + pow(s(j,1) - s(i,1), 2));
          
          log_probAddition += log(1 - exp(- r_new*r_new / (theta)));
          
        }
        
      }
      
    }
    
    
    
    double prob_numerator = exp(loglikelihood_xi_uncaptured(s_j, K, S_usage, trapsX, trapsY, p0, sigma)) * 
      exp(log_probAddition);
    double prob_denominator = 1;//#(1 - exp(loglikelihood_xi_uncaptured(s[j,], K, S_usage, trapsX, trapsY, p0, sigma)))
    
    double R = prob_numerator / prob_denominator;
    
    double mh_ratio = psi * R / (1 - psi + psi * R);
    
    if(R::runif(0,1) < mh_ratio){
      z[j] = 1;
      N += 1;
    } 
    
  }
  
  return(List::create(_["z"] = z,
                      _["N"] = N));
  
} 
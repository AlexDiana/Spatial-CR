#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

int sample_int(int n) {
  Rcpp::IntegerVector pool = Rcpp::seq(1, n);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[0];
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

/// TRAPS FUNCTIONS

// [[Rcpp::export]]
bool checkPointIsInRegion(arma::vec x, 
                          double a1,
                          double b1,
                          double a2,
                          double b2){
  
  
  if(x[0] > a1 & x[0] < b1 & 
     x[1] > a2 & x[1] < b2){
    return(true);     
  }
  
  return(false);
}

// [[Rcpp::export]]
bool checkPointIsInRegionTraps(arma::vec x, arma::mat traps, double R){
  
  
  for(int k = 0; k < traps.n_rows; k++){
    
    if(sqrt(pow(x[0] - traps(k, 0), 2) + pow(x[1] - traps(k, 1), 2)) < R){
      return(true);
    }
    
  }
  
  return(false);
}


// [[Rcpp::export]]
bool insideSinglePolygon(arma::vec point, NumericMatrix vs, int start, int end) {
  
  double x = point[0], y = point[1];
  
  bool inside = false;
  
  int len = end - start;
  for (int i = 0, j = len - 1; i < len; j = i++) {
    double xi = vs(start + i,0), yi = vs(start + i,1);
    double xj = vs(start + j,0), yj = vs(start + j,1);
    bool intersect = ((yi > y) != (yj > y))
      && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  
  return inside;
  
}

// [[Rcpp::export]]
bool checkPointIsInRegionPolygons(arma::vec point, 
                                  NumericMatrix coord, 
                                  IntegerVector polygonHole,
                                  NumericVector polygonStarts){
  
  bool pointInsideHole;
  int holePolygons = sum(polygonHole);
  
  for (int i = 0; i < holePolygons; i++){
    int currentPolygonStart = polygonStarts[i];
    int currentPolygonEnd = polygonStarts[i + 1] - 1;
    
    bool inside = insideSinglePolygon(point, coord, currentPolygonStart, currentPolygonEnd);
    if(inside) return false;
  }
  
  for (int i = holePolygons; i < polygonStarts.length(); i++){
    
    int currentPolygonStart = polygonStarts[i];
    int currentPolygonEnd = (i < polygonStarts.length()) ? polygonStarts[i + 1] - 1 : (coord.nrow() - 1);
    
    bool inside = insideSinglePolygon(point, coord, currentPolygonStart, currentPolygonEnd);
    
    if(inside) return true;
  }
  
  return false;
}

// [[Rcpp::export]]
bool checkPointIsInRegionPolygonsAndTraps(arma::vec point, 
                                          NumericMatrix coord, 
                                          IntegerVector polygonHole,
                                          NumericVector polygonStarts,
                                          arma::mat traps, double R){
  
  if(checkPointIsInRegionPolygons(point, coord, polygonHole, polygonStarts) & 
     checkPointIsInRegionTraps(point, traps, R)){
    return(true);
  }
  
  return(false);
}

//// NORMAL SOFTCORE PROCESS

double lfactorial_single(int n){
  
  NumericVector xx(1); xx(0) = n;
  return(lfactorial(xx)[0]);
}

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
      
      double r2 = pow(x1(n,0) - x1(n2,0), 2) + pow(x1(n,1) - x1(n2,1), 2) ;
      
      loglikelihood1 += log(1 - exp(- r2 / (theta1)));
      
      // double r = sqrt( pow(x1(n,0) - x1(n2,0), 2) + pow(x1(n,1) - x1(n2,1), 2) );
      // 
      // loglikelihood1 += log(1 - exp(- r*r / (theta1)));
      
    }
    
  }  
  
  // within x2
  
  double loglikelihood2 = 0;
  
  for(int n = 1; n < N2; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      // double r = sqrt( pow(x2(n,0) - x2(n2,0), 2) + pow(x2(n,1) - x2(n2,1), 2) );
      // 
      // loglikelihood2 += log(1 - exp(- r*r / (theta2)));
      
      double r2 = pow(x2(n,0) - x2(n2,0), 2) + pow(x2(n,1) - x2(n2,1), 2) ;
      
      loglikelihood2 += log(1 - exp(- r2 / (theta2)));
      
    }
    
  }
  
  // between
  
  double loglikelihood3 = 0;
  
  for(int n = 0; n < N1; n++){
    
    for(int n2 = 0; n2 < N2; n2++){
      
      double r2 =  pow(x1(n,0) - x2(n2,0), 2) + pow(x1(n,1) - x2(n2,1), 2) ;
      
      loglikelihood3 += log(1 - exp(- r2 / (theta3)));
      
      // double r = sqrt( pow(x1(n,0) - x2(n2,0), 2) + pow(x1(n,1) - x2(n2,1), 2) );
      // 
      // loglikelihood3 += log(1 - exp(- r*r / (theta3)));
      
    }
    
  }
  
  return(loglikelihood + loglikelihood1 + loglikelihood2 + loglikelihood3);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_birth(arma::vec xnew, 
                                         int dim, double oldLikelihood, 
                                         arma::mat x1, arma::mat x2, 
                                         int N1, int N2, 
                                         double beta1, double beta2,
                                         double theta1, double theta2, double theta3){
  
  double loglikelihood = oldLikelihood;
  
  if(dim == 1){
    
    loglikelihood += log(beta1);
    loglikelihood -= log(N1 + 1);
    
    // between
    
    for(int j = 0; j < N1; j++){
      
      double r_new = sqrt(pow(x1(j,0) - xnew(0), 2) + pow(x1(j,1) - xnew(1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta1)));
      
    }
    
    // within
    
    for(int j = 0; j < N2; j++){
      
      double r_new = sqrt(pow(x2(j,0) - xnew(0), 2) + pow(x2(j,1) - xnew(1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta3)));
      
    }
    
  } else {
    
    loglikelihood += log(beta2);
    loglikelihood -= log(N2 + 1);
    
    // between
    
    for(int j = 0; j < N2; j++){
      
      double r_new = sqrt(pow(x2(j,0) - xnew(0), 2) + pow(x2(j,1) - xnew(1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta2)));
      
    }
    
    // within
    
    for(int j = 0; j < N1; j++){
      
      double r_new = sqrt(pow(x1(j,0) - xnew(0), 2) + pow(x1(j,1) - xnew(1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta3)));
      
    }
    
  }
  
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_death(int i,
                                         int dim, double oldLikelihood, 
                                         arma::mat x1, arma::mat x2, 
                                         int N1, int N2, 
                                         double beta1, double beta2,
                                         double theta1, double theta2, double theta3){
  
  double loglikelihood = oldLikelihood;
  
  if(dim == 1){
    
    loglikelihood -= log(beta1);
    loglikelihood += log(N1);
    
    // between
    
    for(int j = 0; j < N1; j++){
      
      if(j != i){
        
        double r_old = sqrt(pow(x1(j,0) - x1(i,0), 2) + pow(x1(j,1) - x1(i,1), 2));
        
        loglikelihood -= log(1 - exp(- r_old*r_old / (theta1)));
        
      }
      
    }
    
    // within
    
    for(int j = 0; j < N2; j++){
      
      // if(j != i){
      
      double r_old = sqrt(pow(x2(j,0) - x1(i,0), 2) + pow(x2(j,1) - x1(i,1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta3)));
      
      // }
      
    }
    
  } else {
    
    loglikelihood -= log(beta2);
    loglikelihood += log(N2);
    
    // between
    
    for(int j = 0; j < N2; j++){
      
      if(j != i){
        
        double r_old = sqrt(pow(x2(j,0) - x2(i,0), 2) + pow(x2(j,1) - x2(i,1), 2));
        
        loglikelihood -= log(1 - exp(- r_old*r_old / (theta2)));
        
      }
      
    }
    
    // within
    
    for(int j = 0; j < N1; j++){
      
      // if(j != i){
      
      double r_old = sqrt(pow(x1(j,0) - x2(i,0), 2) + pow(x1(j,1) - x2(i,1), 2));
      
      loglikelihood -= log(1 - exp(- r_old*r_old / (theta3)));
      
      // }
      
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

// [[Rcpp::export]]
arma::vec proposeNewPoint(double a1, 
                          double b1, 
                          double a2, 
                          double b2){
  
  arma::vec proposedPoint = arma::zeros(2);
  
  proposedPoint[0] = R::runif(a1, b1);
  proposedPoint[1] = R::runif(a2, b2);
  
  return(proposedPoint);
}

// [[Rcpp::export]]
List simulate_bivsoftcore_cpp(double theta1, double theta2, double theta3,
                              double beta1, double beta2,
                              int niter, arma::mat Sigma_prop, 
                              int Nmax, double lambda,
                              arma::mat Sigma_newpoint, 
                              double a1,
                              double b1,
                              double a2, 
                              double b2){
  
  arma::vec loglikelihoods = arma::zeros(niter);
  
  double q = .33333;
  double p = .5;
  
  arma::mat data1;
  data1.zeros(Nmax, 2);
  
  int N1 = R::rpois(lambda);
  for(int i = 0; i < N1; i++){
    data1.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPoint(a1, b1, a2, b2));
  }
  
  arma::mat data2;
  data2.zeros(Nmax, 2);
  
  int N2 = R::rpois(lambda);
  for(int i = 0; i < N2; i++){
    data2.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPoint(a1, b1, a2, b2));
  }
  
  double log_hastings_ratio_den = 0;//log_f_bivsoftcore_cpp(data1, data2, N1, N2,
  //                  beta1, beta2,
  //                theta1, theta2, theta3);
  
  for(int iter = 0; iter < niter; iter++){
    
    // x1
    
    if(R::runif(0,1) < q){
      
      // move
      
      // sample point
      int n = sample_int(N1) - 1;
      
      arma::rowvec old_xi = data1.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      // if(checkPointIsInRegionPolygons(xi, polycoord1, 
      // polygonHole1, polygonStarts1)){
      // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord1,
      // polygonHole1, polygonStarts1, traps, R)){
      if(checkPointIsInRegion(xi, 
                              a1, b1, a2, b2)){
        
        data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 1, log_hastings_ratio_den,
                                                                         data1, data2,
                                                                         N1, N2,
                                                                         oldxivec,
                                                                         theta1, theta2, theta3);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data1.row(n) = old_xi;
        }
        
      }
      
    } else {
      
      if(R::runif(0,1) < p){
        
        // birth
        
        arma::vec xi = proposeNewPoint(a1, b1, a2, b2);
        
        if(checkPointIsInRegion(xi, 
                                a1, b1, a2, b2)){
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 1, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data1.row(N1) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N1 += 1;
          }
          
        }
        
      } else {
        
        // death
        
        if (N1 > 0){
          
          // sample point
          int itemToRemove = sample_int(N1) - 1;
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
                                                                            1, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data1.row(itemToRemove) = data1.row(N1 - 1);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N1 -= 1;
          }
          
        }
        
      }
      
    }
    
    // x2
    
    if(R::runif(0,1) < q){ 
      
      // move
      
      // sample point
      int n = sample_int(N2) - 1;
      
      arma::rowvec old_xi = data2.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord2, 
      // polygonHole2, polygonStarts2, traps, R)){
      // if(checkPointIsInRegionPolygons(xi, polycoord2, 
      //                                 polygonHole2, polygonStarts2)){
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
        double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 2, log_hastings_ratio_den,
                                                                         data1, data2, N1, N2,
                                                                         oldxivec, 
                                                                         theta1, theta2, theta3);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data2.row(n) = old_xi;
        }
        
      }
      
    } else {
      
      if(R::runif(0,1) < p){
        
        // birth
        
        arma::vec xi = proposeNewPoint(a1, b1, a2, b2);
        
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 2, log_hastings_ratio_den,
                                                                            data1, data2, N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data2.row(N2) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N2 += 1;
          }
          
        }
        
      } else {
        
        // death
        
        if(N2 > 0){
          
          // sample point
          int itemToRemove = sample_int(N2) - 1;
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
                                                                            2, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data2.row(itemToRemove) = data2.row(N2 - 1);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N2 -= 1;
          }
          
        }
        
      }
      
    }
    
    loglikelihoods[iter] = log_hastings_ratio_den;
    
  }  
  
  return(List::create(_["data1"] = data1,
                      _["data2"] = data2,
                      _["loglikelihoods"] = loglikelihoods,
                      _["N1"] = N1,
                      _["N2"] = N2));
}

// [[Rcpp::export]]
List simulate_bivsoftcore_from_startingpoint(arma::mat data1, arma::mat data2,
                                             int N1, int N2, 
                                             double theta1, double theta2, double theta3,
                                             double beta1, double beta2,
                                             int niter, 
                                             arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                                             double a1,
                                             double b1,
                                             double a2, 
                                             double b2){
  
  double q = .33333;
  double p = .5;
  
  double log_hastings_ratio_den = 0;
  
  for(int iter = 0; iter < niter; iter++){
    
    // x1
    
    if(R::runif(0,1) < q){
      
      if(N1 > 0){
        
        // move
        
        // sample point
        int n = sample_int(N1) - 1;
        
        arma::rowvec old_xi = data1.row(n);
        arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
        arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
        
        // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord1, 
        // polygonHole1, polygonStarts1, traps, R)){
        // if(checkPointIsInRegionPolygons(xi, polycoord1, 
        //                                 polygonHole1, polygonStarts1)){
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 1, log_hastings_ratio_den,
                                                                           data1, data2,
                                                                           N1, N2,
                                                                           oldxivec,
                                                                           theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
          } else {
            data1.row(n) = old_xi;
          }
          
        }
        
      }
      
    } else {
      
      if(R::runif(0,1) < p){
        
        // birth
        
        arma::vec xi = proposeNewPoint(a1, b1, a2, b2);
        
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 1, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data1.row(N1) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N1 += 1;
          }
          
          
        }
        
      } else {
        
        // death
        
        if (N1 > 0){
          
          // sample point
          int itemToRemove = sample_int(N1) - 1;
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
                                                                            1, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data1.row(itemToRemove) = data1.row(N1 - 1);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N1 -= 1;
          }
          
        }
        
      }
      
      
    }
    
    // x2
    
    if(R::runif(0,1) < q){ 
      
      if(N2 > 0){
        
        // move
        
        // sample point
        int n = sample_int(N2) - 1;
        
        arma::rowvec old_xi = data2.row(n);
        arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
        arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
        
        // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord2, 
        //                                         polygonHole2, polygonStarts2, 
        //                                         traps, R)){
        // if(checkPointIsInRegionPolygons(xi, polycoord2, 
        //                                 polygonHole2, polygonStarts2)){
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          
          arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 2, log_hastings_ratio_den,
                                                                           data1, data2, N1, N2,
                                                                           oldxivec, 
                                                                           theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
          } else {
            data2.row(n) = old_xi;
          }
          
        }
        
      }
      
      
    } else {
      
      if(R::runif(0,1) < p){
        
        // birth
        
        arma::vec xi = proposeNewPoint(a1, b1, a2, b2);
        
        if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(xi, 2, log_hastings_ratio_den,
                                                                            data1, data2, N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data2.row(N2) = arma::conv_to<arma::rowvec>::from(xi);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N2 += 1;
          }
          
        }
        
      } else {
        
        // death
        
        if(N2 > 0){
          
          // sample point
          int itemToRemove = sample_int(N2) - 1;
          
          double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
                                                                            2, log_hastings_ratio_den,
                                                                            data1, data2,
                                                                            N1, N2,
                                                                            beta1, beta2,
                                                                            theta1, theta2, theta3);
          
          if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
            data2.row(itemToRemove) = data2.row(N2 - 1);
            log_hastings_ratio_den = log_hastings_ratio_num;
            N2 -= 1;
          }
          
        }
        
      }
      
    }
    
  }  
  
  return(List::create(_["data1"] = data1,
                      _["data2"] = data2,
                      _["N1"] = N1,
                      _["N2"] = N2));
}

// [[Rcpp::export]]
List update_x_all_cpp(arma::cube x_all1, arma::cube x_all2, arma::mat N_all, int niter,
                      double theta1, double theta2, double theta3, double beta1, double beta2,
                      arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                      double a1, double b1, double a2, double b2){
  
  int M = x_all1.n_rows;
  
  for(int m = 0; m < M; m++){
    
    int N1 = N_all(m,0);
    int N2 = N_all(m,1);
    
    arma::mat x1 = x_all1.subcube(arma::span(m), arma::span(), arma::span());
    arma::mat x2 = x_all2.subcube(arma::span(m), arma::span(), arma::span());
    
    List list_sims = simulate_bivsoftcore_from_startingpoint(x1, x2, N1, N2, 
                                                             theta1, theta2, theta3, beta1, beta2,
                                                             niter, Sigma_prop, Sigma_newpoint,
                                                             a1, b1, a2, b2);
    N1 = list_sims["N1"];
    N2 = list_sims["N2"];
    x1 = as<arma::mat>(list_sims["data1"]);
    x2 = as<arma::mat>(list_sims["data2"]);
    x_all1.subcube(arma::span(m),arma::span(),arma::span()) = x1;
    x_all2.subcube(arma::span(m),arma::span(),arma::span()) = x2;
    N_all(m,0) = N1;
    N_all(m,1) = N2;
    
  }
  
  return(List::create(_["x_all1"] = x_all1,
                      _["x_all2"] = x_all2,
                      _["N_all"] = N_all));
}


// [[Rcpp::export]]
List simulate_cond_bivsoftcore_cpp(double theta1, double theta2, double theta3,
                                   double N1, double N2,
                                   int niter, arma::mat Sigma_prop, 
                                   int Nmax, double lambda,
                                   arma::mat Sigma_newpoint, 
                                   double a1,
                                   double b1,
                                   double a2, 
                                   double b2){
  
  arma::vec loglikelihoods = arma::zeros(niter);
  
  double q = .33333;
  double p = .5;
  
  arma::mat data1;
  data1.zeros(Nmax, 2);
  
  for(int i = 0; i < N1; i++){
    data1.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPoint(a1, b1, a2, b2));
  }
  
  arma::mat data2;
  data2.zeros(Nmax, 2);
  
  for(int i = 0; i < N2; i++){
    data2.row(i) = arma::conv_to<arma::rowvec>::from(proposeNewPoint(a1, b1, a2, b2));
  }
  
  double log_hastings_ratio_den = 0;
  
  for(int iter = 0; iter < niter; iter++){
    
    // x1
    
    for(int n = 0; n < N1; n++){
      
      arma::rowvec old_xi = data1.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      // if(checkPointIsInRegionPolygons(xi, polycoord1, 
      // polygonHole1, polygonStarts1)){
      // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord1,
      // polygonHole1, polygonStarts1, traps, R)){
      if(checkPointIsInRegion(xi, 
                              a1, b1, a2, b2)){
        
        data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 1, log_hastings_ratio_den,
                                                                         data1, data2,
                                                                         N1, N2,
                                                                         oldxivec,
                                                                         theta1, theta2, theta3);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data1.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data1.row(n) = old_xi;
        }
        
      }
      
    }
    
    // x2
    
    for(int n = 0; n < N2; n++){
      
      arma::rowvec old_xi = data2.row(n);
      arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
      arma::vec xi = mvrnormArma(oldxivec, Sigma_prop);
      
      // if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord2, 
      // polygonHole2, polygonStarts2, traps, R)){
      // if(checkPointIsInRegionPolygons(xi, polycoord2, 
      //                                 polygonHole2, polygonStarts2)){
      if(checkPointIsInRegion(xi, a1, b1, a2, b2)){
        
        data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
        
        arma::vec oldxivec = arma::conv_to<arma::vec>::from(old_xi);
        double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(n, 2, log_hastings_ratio_den,
                                                                         data1, data2, N1, N2,
                                                                         oldxivec, 
                                                                         theta1, theta2, theta3);
        
        if(R::runif(0, 1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
          data2.row(n) = arma::conv_to<arma::rowvec>::from(xi);
          log_hastings_ratio_den = log_hastings_ratio_num;
        } else {
          data2.row(n) = old_xi;
        }
        
      }
      
    }
    
    loglikelihoods[iter] = log_hastings_ratio_den;
    
  }  
  
  return(List::create(_["data1"] = data1,
                      _["data2"] = data2));
}

/////////////// MODEL FUNCTIONS

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

// double loglikelihood_p_cpp_old(double p0, double sigma,
//                            arma::mat CH, arma::mat s_n, 
//                            arma::vec trapsX, arma::vec trapsY,
//                            arma::vec S_usage, int N){
//   
//   double loglikelihood = 0;
//   
//   int K = CH.n_cols;
//   
//   for(int i = 0; i < N; i++){
//     for(int j = 0; j < K; j++){
//       double p = p0 * exp(- (1 / (2 * sigma * sigma)) * (pow(s_n(i,0) - trapsX[j], 2) + pow(s_n(i,1) - trapsY[j],2)));
//       
//       loglikelihood += R::dbinom(CH(i,j), S_usage[j], p, 1);
//       
//     }
//   }
//   
//   return(loglikelihood);
// }

double loglikelihood_p_cpp(double p0, double sigma,
                           arma::mat CH, arma::mat s_n, 
                           arma::vec trapsX, arma::vec trapsY,
                           arma::vec S_usage, int N){
  
  double loglikelihood = 0;
  
  int K = CH.n_cols;
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < K; j++){
      
      double d = (pow(s_n(i,0) - trapsX[j], 2) + pow(s_n(i,1) - trapsY[j],2));
      
      double p = p0 * exp(- (1 / (2 * sigma * sigma)) * d);
      
      if(p > 0){
        // loglikelihood += R::dbinom(CH(i,j), S_usage[j], p, 1);
        loglikelihood += CH(i,j) * log(p) + (S_usage[j] - CH(i,j)) * log(1 - p);
        // Rcout << d  << " - " << loglikelihood << std::endl;
        
      } else if (p == 0 & CH(i,j) > 0){
        loglikelihood += R_NegInf;
      }
      
    }
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
arma::vec update_psigma_cpp(double p0, 
                            double sigma,
                            arma::mat CH,
                            arma::mat s,
                            int N, 
                            int K,
                            arma::vec S_usage,
                            arma::vec trapsX,
                            arma::vec trapsY,
                            double sigma_0,
                            double sd_sigma0,
                            double sd_sigma_prop,
                            double a_p,
                            double b_p){
  
  // update p0
  
  double a_fp = 0;
  double b_fp = 0;
  double c_fp = 0;
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < K; j++){
      double k_i = S_usage[j] - CH(i, j);
      double c_i =  exp(- (1 / (2 * sigma * sigma)) * ( pow(s(i, 0) - trapsX[j], 2) + pow(s(i, 1) - trapsY[j], 2)));
      double ci_sq = c_i * c_i;
      a_fp += CH(i, j);
      b_fp -= c_i * k_i; 
      c_fp += k_i * ci_sq;
    }
  }
  
  double p_star = (- b_fp - sqrt(b_fp * b_fp - 4 * a_fp * c_fp)) / (2 * c_fp);
  double h_star = 1 / (a_fp / (p_star * p_star));
  
  double p_new = R::rnorm(p_star, sqrt(h_star));
  
  if(p_new > 0 & p_new < 1){
    
    double loglikelihood_p0_star = loglikelihood_p_cpp(p_new, sigma, CH, s,
                                                       trapsX, trapsY, S_usage, N);
    
    double loglikelihood_p0 = loglikelihood_p_cpp(p0, sigma, CH, s,
                                                  trapsX, trapsY, S_usage, N);
    
    double logprior_p0_star = R::dbeta(p_new, a_p, b_p, 1);
    
    double logprior_p0 = R::dbeta(p0, a_p, b_p, 1);
    
    double logproposal_p0 = R::dnorm(p0, p_star, sqrt(h_star), 1);
    
    double logproposal_p_new = R::dnorm(p_new, p_star, sqrt(h_star), 1);
    
    double logposterior_p0_star = logprior_p0_star + loglikelihood_p0_star;
    
    double logposterior_p0 = logprior_p0 + loglikelihood_p0;
    
    double MH_ratio = exp(logposterior_p0_star + logproposal_p0 - 
                          logposterior_p0 - logproposal_p_new);
    
    if(R::runif(0, 1) < MH_ratio){
      p0 = p_new;
    }
    
  }
  
  // update sigma
  
  double sigma_star = R::rnorm(sigma, sd_sigma_prop);
  
  if(sigma_star > 0){
    
    double loglikelihood_sigma_star = loglikelihood_p_cpp(p0, sigma_star, CH, s,
                                                          trapsX, trapsY, S_usage, N);
    
    double loglikelihood_sigma = loglikelihood_p_cpp(p0, sigma, CH, s,
                                                     trapsX, trapsY, S_usage, N);
    
    double logprior_sigma_star = R::dnorm(sigma_star, sigma_0, sd_sigma0, 1);
    
    double logprior_sigma = R::dnorm(sigma, sigma_0, sd_sigma0, 1);
    
    double logposterior_sigma_star = loglikelihood_sigma_star + logprior_sigma_star;
    
    double logposterior_sigma = loglikelihood_sigma + logprior_sigma;
    
    if(R::runif(0, 1) < exp(loglikelihood_sigma_star - loglikelihood_sigma)){
      sigma = sigma_star;
    }
    
  }
  
  arma::vec out = arma::zeros(2);
  
  out[0] = p0;
  out[1] = sigma;
  
  return(out);
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


///////// 
/// S
/////////

double probAcceptingIndividual(arma::mat s, arma::vec s_new, double p0, double sigma,
                               int D, int N, int K,
                               arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  // // compute mean probability of capture
  // arma::vec p_all = arma::zeros(N);
  // for(int i = 0; i < N; i++){
  //   arma::vec si = arma::conv_to<arma::vec>::from(s.row(i));
  //   p_all[i] = 1 - exp(loglikelihood_xi_uncaptured(si, K, S, trapsX, trapsY, p0, sigma));
  // }  
  // 
  // double p_tilde = mean(p_all);
  // 
  // double ratio_probs = D * (1 - p_tilde) / ((N - D + 1) * p_tilde);
  
  double p_notcaptured = exp(loglikelihood_xi_uncaptured(s_new, K, S_usage, trapsX, trapsY, p0, sigma));
  double prob = (p_notcaptured * (N + 1)) / (N + 1 - D);
  
  return(prob);
}

double probProposingIndividual(arma::vec point, double probMixture,
                               double a1, double b1, double a2, double b2,
                               arma::mat mixtureMeans, arma::mat mixtureSd, arma::vec norm_const){
  
  double likelihood = 0;
  
  int numCenters = mixtureMeans.n_rows;
  
  likelihood += (1 - probMixture) * (R::dunif(point[0], a1, b1, 0)) * 
    (R::dunif(point[1], a2, b2, 0)) / norm_const[0];
  // Rcout << likelihood << std::endl;
  for(int l = 0; l < numCenters; l++){
    
    likelihood += (probMixture) * (1.0 / numCenters) * (R::dnorm(point[0], mixtureMeans(l,0), mixtureSd(l,0), 0) * 
      R::dnorm(point[1], mixtureMeans(l,1), mixtureSd(l,1), 0)) / norm_const[l + 1];
    
  }
  // Rcout << likelihood << std::endl;
  return(likelihood);
}

// // [[Rcpp::export]]
// List update_s_cpp_withacceptance(arma::mat s1, arma::mat s2, 
//                                  double theta1, double theta2, double theta3, 
//                                  double beta1, double beta2, 
//                                  arma::mat CH1, arma::mat CH2,
//                                  int N1, int N2, int D1, int D2,
//                                  double p0_1, double sigma_1,
//                                  double p0_2, double sigma_2,
//                                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY,
//                                  arma::mat Sigma_prop, arma::mat Sigma_newpoint,
//                                  NumericMatrix polycoord1, 
//                                  IntegerVector polygonHole1,
//                                  NumericVector polygonStarts1,
//                                  arma::vec polyBoundaries1,
//                                  NumericMatrix polycoord2, 
//                                  IntegerVector polygonHole2,
//                                  NumericVector polygonStarts2,
//                                  arma::vec polyBoundaries2,
//                                  arma::mat traps,
//                                  double R){
//   
//   int K = CH1.n_cols;
//   
//   // double q = .3333333333;
//   double q = 0;
//   double p = .5;
//   
//   double log_hastings_ratio_den = log_f_bivsoftcore_cpp(s1, s2, N1, N2,
//                                                         beta1, beta2,
//                                                         theta1, theta2, theta3);
//   
//   bool isPoint1Accepted = false;
//   bool isPoint2Accepted = false;
//   arma::vec acceptedPoint1 = arma::zeros(2);
//   arma::vec acceptedPoint2 = arma::zeros(2);
//   
//   // first process
//   
//   if(R::runif(0,1) < q){
//     
//     // move
//     
//     for(int i = 0; i < N1; i++) {
//       
//       arma::vec old_xi = arma::conv_to<arma::vec>::from(s1.row(i));
//       arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
//       
//       if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord1, 
//                                               polygonHole1, polygonStarts1, 
//                                               traps, R)){
//         // if(checkPointIsInRegionPolygons(xi, polycoord1, 
//         //                                 polygonHole1, polygonStarts1)){
//         
//         s1.row(i) = arma::conv_to<arma::rowvec>::from(xi);
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(i, 1, log_hastings_ratio_den,
//                                                                          s1, s2,
//                                                                          N1, N2,
//                                                                          old_xi,
//                                                                          theta1, theta2, theta3);
//         
//         double loglikelihood_xi_star;
//         double loglikelihood_xi_old;
//         
//         if(i < D1){
//           loglikelihood_xi_star = loglikelihood_xi_captured(xi, i, K, S_usage, trapsX, trapsY, CH1, p0_1, sigma_1);
//           loglikelihood_xi_old = loglikelihood_xi_captured(old_xi, i, K, S_usage, trapsX, trapsY, CH1, p0_1, sigma_1);
//         } else {
//           loglikelihood_xi_star = loglikelihood_xi_uncaptured(xi, K, S_usage, trapsX, trapsY, p0_1, sigma_1);
//           loglikelihood_xi_old = loglikelihood_xi_uncaptured(old_xi, K, S_usage, trapsX, trapsY, p0_1, sigma_1);
//         }
//         
//         double logposterior_xistar = log_hastings_ratio_num + loglikelihood_xi_star;
//         double logposterior_xi = log_hastings_ratio_den + loglikelihood_xi_old;
//         
//         if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//           log_hastings_ratio_den = log_hastings_ratio_num;
//         } else {
//           s1.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
//         }
//       }
//       
//     }
//     
//   } else {
//     
//     if(R::runif(0,1) < p){
//       
//       // birth
//       
//       arma::vec x1_new = proposeNewPoint(polyBoundaries1[0], 
//                                          polyBoundaries1[1], 
//                                                         polyBoundaries1[2], 
//                                                                        polyBoundaries1[3], 
//                                                                                       polycoord1, 
//                                                                                       polygonHole1,
//                                                                                       polygonStarts1,
//                                                                                       traps, R);
//       
//       double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(x1_new, 1, log_hastings_ratio_den,
//                                                                         s1, s2,
//                                                                         N1, N2,
//                                                                         beta1, beta2,
//                                                                         theta1, theta2, theta3);
//       
//       double loglikelihood_xi_star = log(probAcceptingIndividual(s1, x1_new, p0_1, sigma_1, D1, N1,
//                                                                  K, S_usage, trapsX, trapsY));
//       
//       double logposterior_xistar = log_hastings_ratio_num + loglikelihood_xi_star;
//       double logposterior_xi = log_hastings_ratio_den;
//       
//       // double birthDeathRatio = - log(N1);
//       
//       if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//         log_hastings_ratio_den = log_hastings_ratio_num;
//         s1.row(N1) = arma::conv_to<arma::rowvec>::from(x1_new);
//         N1 += 1;
//         isPoint1Accepted = true;
//         acceptedPoint1 = x1_new;
//       }
//       
//     } else {
//       
//       // death
//       
//       if(N1 > D1){
//         int itemToRemove = sample_int(N1 - D1) + D1 - 1;
//         
//         arma::vec xi = arma::conv_to<arma::vec>::from(s1.row(itemToRemove));
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
//                                                                           1, log_hastings_ratio_den,
//                                                                           s1, s2,
//                                                                           N1, N2,
//                                                                           beta1, beta2,
//                                                                           theta1, theta2, theta3);
//         
//         double loglikelihood_xi_star = log(probAcceptingIndividual(s1, xi, p0_1, sigma_1, D1, N1 - 1, 
//                                                                    K, S_usage, trapsX, trapsY));
//         
//         double logposterior_xistar = log_hastings_ratio_num;
//         double logposterior_xi = log_hastings_ratio_den + loglikelihood_xi_star;
//         
//         if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//           s1.row(itemToRemove) = s1.row(N1 - 1);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//           N1 -= 1;
//         }
//       }
//       
//     }
//     
//   }
//   
//   // second process
//   
//   if(R::runif(0,1) < q){
//     
//     // move
//     
//     for(int i = 0; i < N2; i++) {
//       
//       arma::vec old_xi = arma::conv_to<arma::vec>::from(s2.row(i));
//       arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
//       
//       if(checkPointIsInRegionPolygonsAndTraps(xi, polycoord2, 
//                                               polygonHole2, polygonStarts2, 
//                                               traps, R)){
//         // if(checkPointIsInRegionPolygons(xi, polycoord2, 
//         //                                 polygonHole2, polygonStarts2)){
//         
//         s2.row(i) = arma::conv_to<arma::rowvec>::from(xi);
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_move(i, 2, log_hastings_ratio_den,
//                                                                          s1, s2,
//                                                                          N1, N2,
//                                                                          old_xi,
//                                                                          theta1, theta2, theta3);
//         
//         double loglikelihood_xi_star;
//         double loglikelihood_xi_old;
//         
//         if(i < D2){
//           loglikelihood_xi_star = loglikelihood_xi_captured(xi, i, K, S_usage, trapsX, trapsY, CH2, p0_2, sigma_2);
//           loglikelihood_xi_old = loglikelihood_xi_captured(old_xi, i, K, S_usage, trapsX, trapsY, CH2, p0_2, sigma_2);
//         } else {
//           loglikelihood_xi_star = loglikelihood_xi_uncaptured(xi, K, S_usage, trapsX, trapsY, p0_2, sigma_2);
//           loglikelihood_xi_old = loglikelihood_xi_uncaptured(old_xi, K, S_usage, trapsX, trapsY, p0_2, sigma_2);
//         }
//         
//         double logposterior_xistar = log_hastings_ratio_num + loglikelihood_xi_star;
//         double logposterior_xi = log_hastings_ratio_den + loglikelihood_xi_old;
//         
//         if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//           log_hastings_ratio_den = log_hastings_ratio_num;
//         } else {
//           s2.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
//         }
//       }
//       
//     }
//     
//   } else {
//     
//     if(R::runif(0,1) < p){
//       
//       // birth
//       
//       arma::vec x2_new = proposeNewPoint(polyBoundaries2[0], 
//                                          polyBoundaries2[1], 
//                                                         polyBoundaries2[2], 
//                                                                        polyBoundaries2[3], 
//                                                                                       polycoord2, 
//                                                                                       polygonHole2,
//                                                                                       polygonStarts2,
//                                                                                       traps, R);
//       
//       double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_birth(x2_new, 2, log_hastings_ratio_den,
//                                                                         s1, s2,
//                                                                         N1, N2,
//                                                                         beta1, beta2,
//                                                                         theta1, theta2, theta3);
//       
//       double loglikelihood_xi_star = log(probAcceptingIndividual(s2, x2_new, p0_2, sigma_2, D2, N2, 
//                                                                  K, S_usage, trapsX, trapsY));
//       
//       double logposterior_xistar = log_hastings_ratio_num + loglikelihood_xi_star;
//       double logposterior_xi = log_hastings_ratio_den;
//       
//       if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//         log_hastings_ratio_den = log_hastings_ratio_num;
//         s2.row(N2) = arma::conv_to<arma::rowvec>::from(x2_new);
//         N2 += 1;
//         isPoint2Accepted = true;
//         acceptedPoint2 = x2_new;
//       }
//       
//     } else {
//       
//       // death
//       
//       if(N2 > D2){
//         int itemToRemove = sample_int(N2 - D2) + D2 - 1;
//         
//         arma::vec xi = arma::conv_to<arma::vec>::from(s2.row(itemToRemove));
//         
//         double log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove,
//                                                                           2, log_hastings_ratio_den,
//                                                                           s1, s2,
//                                                                           N1, N2,
//                                                                           beta1, beta2,
//                                                                           theta1, theta2, theta3);
//         
//         double loglikelihood_xi_star = log(probAcceptingIndividual(s2, xi, p0_2, sigma_2, D2, N2 - 1, 
//                                                                    K, S_usage, trapsX, trapsY));
//         
//         double logposterior_xistar = log_hastings_ratio_num;
//         double logposterior_xi = log_hastings_ratio_den + loglikelihood_xi_star;
//         
//         if(R::runif(0,1) < exp(logposterior_xistar - logposterior_xi)){
//           s2.row(itemToRemove) = s2.row(N2 - 1);
//           log_hastings_ratio_den = log_hastings_ratio_num;
//           N2 -= 1;
//         }
//       }
//       
//     }
//     
//   }
//   
//   
//   return(List::create(_["s1"] = s1,
//                       _["s2"] = s2,
//                       _["N1"] = N1,
//                       _["N2"] = N2,
//                       _["isPoint1Accepted"] = isPoint1Accepted,
//                       _["isPoint2Accepted"] = isPoint2Accepted,
//                       _["acceptedPoint1"] = acceptedPoint1,
//                       _["acceptedPoint2"] = acceptedPoint2));
// }

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_birth_intonly(arma::vec xnew,
                                                 arma::mat x1, arma::mat x2, 
                                                 int N1, int N2, 
                                                 double theta1, double theta12){
  
  double loglikelihood = 0;
  
  // between
  
  for(int j = 0; j < N1; j++){
    
    double r_new = sqrt(pow(x1(j,0) - xnew(0), 2) + pow(x1(j,1) - xnew(1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta1)));
    
  }
  
  // within
  
  for(int j = 0; j < N2; j++){
    
    double r_new = sqrt(pow(x2(j,0) - xnew(0), 2) + pow(x2(j,1) - xnew(1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta12)));
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_addition(arma::mat x_new, int numProposedPoints,
                                            arma::mat x1, arma::mat x2, 
                                            int N1, int N2, 
                                            double beta, double theta1, double theta12){
  
  double loglikelihood = 0;
  
  loglikelihood += numProposedPoints * log(beta);
  for(int l = 0; l < numProposedPoints; l++){
    
    loglikelihood -= log(N1 + 1 + l);
    
  }
  
  // with existing points
  
  for(int l = 0; l < numProposedPoints; l++){
    arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
    loglikelihood += log_f_bivsoftcore_cpp_quick_birth_intonly(x_current, x1, x2,
                                                               N1, N2,  
                                                               theta1, theta12);
  }
  
  // between new points
  
  for(int i = 0; i < numProposedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x_new(i,0) - x_new(j,0), 2) + pow(x_new(i,1) - x_new(j,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta1)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_quick_removal(IntegerVector idxPointsToRemove,
                                           int numRemovedPoints, 
                                           arma::mat x1, arma::mat x2, 
                                           int N1, int N2, double beta,
                                           double theta1, double theta12){
  
  idxPointsToRemove.sort(false);
  
  double loglikelihood = 0;
  
  loglikelihood -= numRemovedPoints * log(beta);
  for(int l = 0; l < numRemovedPoints; l++){
    loglikelihood += log(N1 - l);
  }
  
  // with points not removed
  
  for(int j = 0; j < numRemovedPoints; j++){
    
    int l = 0; // index of removed points
    
    for(int i = 0; i < N1; i++){
      
      if(i != idxPointsToRemove[l]){
        
        double r_new = sqrt(pow(x1(idxPointsToRemove[j],0) - x1(i,0), 2) + 
                            pow(x1(idxPointsToRemove[j],1) - x1(i,1), 2));
        
        loglikelihood -= log(1 - exp(- r_new*r_new / (theta1)));
        
      } else {
        
        l += 1;
        
      }
      
    }
    
  }
  
  // between them
  
  for(int i = 0; i < numRemovedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x1(idxPointsToRemove[i],0) - x1(idxPointsToRemove[j],0), 2) + 
                          pow(x1(idxPointsToRemove[i],1) - x1(idxPointsToRemove[j],1), 2));
      
      loglikelihood -= log(1 - exp(- r_new*r_new / (theta1)));
      
    }
    
  }
  
  // with point of other group
  
  for(int i = 0; i < numRemovedPoints; i++){
    
    for(int j = 0; j < N2; j++){
      
      double r_new = sqrt(pow(x1(idxPointsToRemove[i],0) - x2(j,0), 2) + 
                          pow(x1(idxPointsToRemove[i],1) - x2(j,1), 2));
      
      loglikelihood -= log(1 - exp(- r_new*r_new / (theta12)));
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double densProposal(arma::vec x, double normConst){
  
  double prob = 1 / normConst;
  
  return(prob);
}

// NEW PROPOSAL FUNCTION (DEFINITIVE)

// [[Rcpp::export]]
double logDensityBirth(arma::vec point, 
                       double a1, double b1, double a2, double b2){
  
  double likelihood = 0;
  
  likelihood += (R::dunif(point[0], a1, b1, 0)) * (R::dunif(point[1], a2, b2, 0));
  
  return(log(likelihood));
}

// [[Rcpp::export]]
double logDensityBirthMixture(arma::vec point, double probMixture,
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

// // [[Rcpp::export]]
// arma::vec proposeNewPointNewMixture(double probMixture, arma::mat mixtureMeans, arma::mat mixtureSd, 
//                                     NumericMatrix polycoord, 
//                                     IntegerVector polygonHole,
//                                     NumericVector polygonStarts,
//                                     arma::vec polyBoundaries){
//   arma::vec proposedPoint;
//   
//   if(R::runif(0,1) < (1 - probMixture)){
//     
//     proposedPoint = proposeNewPoint(a1, b1, a2, b2, traps, R);
//     
//     // proposedPoint[0] = R::runif(a1, b1);
//     // proposedPoint[1] = R::runif(a2, b2);
//     // 
//     // while(!checkPointIsInRegionTraps(proposedPoint, traps, R)){
//     //   proposedPoint[0] = R::runif(a1, b1);
//     //   proposedPoint[1] = R::runif(a2, b2);
//     // }
//     
//   } else {
//     
//     int numCenters = mixtureMeans.n_rows;
//     
//     proposedPoint = arma::zeros(2);
//     
//     int indMixture = sample_int(numCenters) - 1;
//     
//     proposedPoint[0] = R::rnorm(mixtureMeans(indMixture, 0), mixtureSd(indMixture, 0));
//     proposedPoint[1] = R::rnorm(mixtureMeans(indMixture, 1), mixtureSd(indMixture, 1));
//     
//     while(!checkPointIsInRegionTraps(proposedPoint, traps, R)){
//       
//       proposedPoint[0] = R::rnorm(mixtureMeans(indMixture, 0), mixtureSd(indMixture, 0));
//       proposedPoint[1] = R::rnorm(mixtureMeans(indMixture, 1), mixtureSd(indMixture, 1));
//     }
//     
//   }
//   
//   return proposedPoint;
//   
// }

double loglikelihood_move(int i, arma::mat x1, arma::mat x2, 
                          int N1, int N2, arma::vec xold, 
                          double theta1, double theta12){
  
  double loglikelihood = 0;
  
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
    
    loglikelihood -= log(1 - exp(- r_old*r_old / (theta12)));
    
    double r_new = sqrt(pow(x2(j,0) - x1(i,0), 2) + pow(x2(j,1) - x1(i,1), 2));
    
    loglikelihood += log(1 - exp(- r_new*r_new / (theta12)));
    
  }
  
  return(loglikelihood);
}

double loglikelihood_addition(arma::mat x_new, int numProposedPoints,
                              arma::mat x1, arma::mat x2, 
                              int N1, int N2, 
                              double beta, double theta1, double theta12){
  
  double loglikelihood = 0;
  
  loglikelihood += numProposedPoints * log(beta);
  
  // with existing points
  
  for(int l = 0; l < numProposedPoints; l++){
    arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
    loglikelihood += log_f_bivsoftcore_cpp_quick_birth_intonly(x_current, x1, x2,
                                                               N1, N2,  
                                                               theta1, theta12);
  }
  
  // between new points
  
  for(int i = 0; i < numProposedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x_new(i,0) - x_new(j,0), 2) + pow(x_new(i,1) - x_new(j,1), 2));
      
      loglikelihood += log(1 - exp(- r_new*r_new / (theta1)));
      
    }
    
  }
  
  return(loglikelihood);
}

double loglikelihood_removal(IntegerVector idxPointsToRemove,
                             int numRemovedPoints, 
                             arma::mat x1, arma::mat x2, 
                             int N1, int N2, double beta,
                             double theta1, double theta12){
  
  idxPointsToRemove.sort(false);
  
  double loglikelihood = 0;
  
  loglikelihood -= numRemovedPoints * log(beta);
  
  // with points not removed
  
  for(int j = 0; j < numRemovedPoints; j++){
    
    int l = 0; // index of removed points
    
    for(int i = 0; i < N1; i++){
      
      if(i != idxPointsToRemove[l]){
        
        double r_new = sqrt(pow(x1(idxPointsToRemove[j],0) - x1(i,0), 2) + 
                            pow(x1(idxPointsToRemove[j],1) - x1(i,1), 2));
        
        loglikelihood -= log(1 - exp(- r_new*r_new / (theta1)));
        
      } else {
        
        l += 1;
        
      }
      
    }
    
  }
  
  // between them
  
  for(int i = 0; i < numRemovedPoints; i++){
    
    for(int j = 0; j < i; j++){
      
      double r_new = sqrt(pow(x1(idxPointsToRemove[i],0) - x1(idxPointsToRemove[j],0), 2) + 
                          pow(x1(idxPointsToRemove[i],1) - x1(idxPointsToRemove[j],1), 2));
      
      loglikelihood -= log(1 - exp(- r_new*r_new / (theta1)));
      
    }
    
  }
  
  // with point of other group
  
  for(int i = 0; i < numRemovedPoints; i++){
    
    for(int j = 0; j < N2; j++){
      
      double r_new = sqrt(pow(x1(idxPointsToRemove[i],0) - x2(j,0), 2) + 
                          pow(x1(idxPointsToRemove[i],1) - x2(j,1), 2));
      
      loglikelihood -= log(1 - exp(- r_new*r_new / (theta12)));
      
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

double computePointsMoveRatio(int i, arma::mat s1, arma::mat s2, int N1, int N2,
                              arma::vec old_xi, double theta1, double theta12,
                              int D, arma::vec xi, int K, arma::vec S_usage,
                              arma::mat trapsX, arma::mat trapsY, arma::mat CH,
                              double p0, double sigma){
  
  double log_hastings_ratio_num = loglikelihood_move(i, s1, s2, N1, N2,
                                                     old_xi, theta1, theta12);
  
  double loglikelihood_xi_star;
  double loglikelihood_xi_old;
  
  if(i < D){
    loglikelihood_xi_star = loglikelihood_xi_captured(xi, i, K, S_usage, trapsX, trapsY, CH, p0, sigma);
    loglikelihood_xi_old = loglikelihood_xi_captured(old_xi, i, K, S_usage, trapsX, trapsY, CH, p0, sigma);
  } else {
    loglikelihood_xi_star = loglikelihood_xi_uncaptured(xi, K, S_usage, trapsX, trapsY, p0, sigma);
    loglikelihood_xi_old = loglikelihood_xi_uncaptured(old_xi, K, S_usage, trapsX, trapsY, p0, sigma);
  }
  
  // loglikelihood_xi_star = 0;
  // loglikelihood_xi_old = 0;
  
  double logNumerator = log_hastings_ratio_num + loglikelihood_xi_star;
  double logDenominator = loglikelihood_xi_old;  
  
  return(exp(logNumerator - logDenominator));
}

double computePointsAdditionRatio(arma::mat x_new, int ptsToPropose, int N, int N2,
                                  arma::mat s1, arma::mat s2, 
                                  double beta, double theta1, double theta12,
                                  double a1, double b1, double a2, double b2, double S_area,
                                  int D, int K, double p0, double sigma,
                                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double logRatio_f = loglikelihood_addition(x_new, ptsToPropose, 
                                             s1, s2, N, N2, 
                                             beta / S_area, theta1, theta12);
  
  double logbirthDensity = 0;
  for(int l = 0; l < ptsToPropose; l++){
    arma::vec x1_new = arma::conv_to<arma::vec>::from(x_new.row(l));
    logbirthDensity += logDensityBirth(x1_new, 
                                       a1, b1, a2, b2);
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
  // loglikelihood_xi_star = 0;
  
  double logNumerator = logRatio_f + loglikelihood_xi_star + logDeathDensity;
  double logDenominator = logbirthDensity;
  
  return(exp(logNumerator - logDenominator));
}

double computePointsDeletionRatio(IntegerVector idxPointsToRemove, int pointsToRemove,
                                  int N, int N2, arma::mat s1, arma::mat s2, 
                                  double beta, double theta1, double theta12,
                                  double a1, double b1, double a2, double b2, double S_area,
                                  int D, int K, double p0, double sigma,
                                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double logRatio_f = loglikelihood_removal(idxPointsToRemove, pointsToRemove, 
                                            s1, s2, N, N2, 
                                            beta / S_area, theta1, theta12);
  
  double logbirthDensity = 0;
  for(int l = 0; l < pointsToRemove; l++){
    
    arma::vec xi = arma::conv_to<arma::vec>::from(s1.row(idxPointsToRemove[l]));
    
    logbirthDensity += logDensityBirth(xi, 
                                       a1, b1, a2, b2);
  }
  
  double logDeathDensity = 0;
  for(int l = 0; l < pointsToRemove; l++){
    logDeathDensity -= log(N - l);
  }
  
  double loglikelihood_xi_star = 0;
  for(int l = 0; l < pointsToRemove; l++){
    
    arma::vec xi = arma::conv_to<arma::vec>::from(s1.row(idxPointsToRemove[l]));
    
    loglikelihood_xi_star += log(probIndividualNotCaptured(xi, p0, sigma, D, N - l - 1,
                                                           K, S_usage, trapsX, trapsY));
  }
  // loglikelihood_xi_star = 0;
  
  double logNumerator = logRatio_f + logbirthDensity;
  double logDenominator = logDeathDensity + loglikelihood_xi_star;
  
  return(exp(logNumerator - logDenominator));
}

// [[Rcpp::export]]
void ProposePointProcess(arma::mat& s1, arma::mat& s2,
                         double theta1, double theta2, double theta12,
                         double beta,
                         arma::mat CH, double p0, double sigma,
                         int& N1, int N2, int D1,
                         arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                         double lambda_movedPoints, double lambda_newPoints, 
                         double S_area, arma::vec S_usage, 
                         double a1, double b1, double a2, double b2,
                         arma::mat traps, double R,
                         arma::vec trapsX, arma::vec trapsY,
                         double q, double p, double& mhRatio){
  
  int K = CH.n_cols;
  
  // double q = .3333333333;
  // double p = .5;
  
  if(R::runif(0,1) < q){
    
    // move
    
    int pointsToMove = R::rpois(lambda_movedPoints);
    IntegerVector seqN1D1 = seq(1, N1);
    
    IntegerVector idxPointsToMove = RcppArmadillo::sample(seqN1D1, pointsToMove, 0) - 1;
    
    for(int j = 0; j < pointsToMove; j++) {
      // for(int i = 0; i < N1; i++) {
      
      int i = idxPointsToMove[j];
      
      arma::vec old_xi = arma::conv_to<arma::vec>::from(s1.row(i));
      arma::vec xi = mvrnormArma(old_xi, Sigma_prop);
      
      if(checkPointIsInRegion(xi, 
                              a1, b1, a2, b2)){
        // if(checkPointIsInRegionPolygons(xi, polycoord, 
        //                                 polygonHole, polygonStarts)){
        
        s1.row(i) = arma::conv_to<arma::rowvec>::from(xi);
        
        double mh_ratio = computePointsMoveRatio(i, s1, s2, N1, N2,
                                                 old_xi, theta1, theta12,
                                                 D1, xi, K, S_usage,
                                                 trapsX, trapsY, CH,
                                                 p0, sigma);
        // Rcout << mh_ratio << std::endl;
        if(R::runif(0,1) > mh_ratio){
          s1.row(i) = arma::conv_to<arma::rowvec>::from(old_xi);
        }
        
        mhRatio = mh_ratio;
        
      }
      
    }
    
  } else {
    
    if(R::runif(0,1) < p){
      
      // birth
      
      int pointsToPropose = 1 +  R::rpois(lambda_newPoints);
      
      arma::mat x_new = arma::zeros(pointsToPropose, 2);
      for(int l = 0; l < pointsToPropose; l++){
        arma::vec x1_new = proposeNewPoint(a1, b1, a2, b2);
        
        x_new.row(l) = arma::conv_to<arma::rowvec>::from(x1_new);
      }
      
      double mh_ratio = computePointsAdditionRatio(x_new, pointsToPropose, N1, N2,
                                                   s1, s2, beta, theta1, theta12,
                                                   a1, b1, a2, b2,                                                                                                         S_area,
                                                   D1, K, p0, sigma,
                                                   S_usage, trapsX, trapsY);
      
      
      if(R::runif(0,1) < mh_ratio){
        
        for(int l = 0; l < pointsToPropose; l++){
          s1.row(N1 + l) = x_new.row(l);
        }
        N1 += pointsToPropose;
      }
      
      mhRatio = mh_ratio;
      
    } else {
      
      // death
      
      int pointsToRemove = 1 +  R::rpois(lambda_newPoints);
      
      if(N1 - pointsToRemove >= D1){
        
        IntegerVector seqN1D1 = seq(1, N1 - D1);
        
        IntegerVector idxPointsToRemove = D1 - 1 + RcppArmadillo::sample(seqN1D1, pointsToRemove, 0);
        
        double mh_ratio = computePointsDeletionRatio(idxPointsToRemove, pointsToRemove, N1, N2,
                                                     s1, s2, beta, theta1, theta12,
                                                     a1, b1, a2, b2,
                                                     S_area,
                                                     D1, K, p0, sigma,
                                                     S_usage, trapsX, trapsY);
        
        
        if(R::runif(0,1) < mh_ratio){
          
          for(int l = 0; l < pointsToRemove; l++){
            s1.row(idxPointsToRemove[l]) = s1.row(N1 - 1);
            N1 -= 1;
          }
          
        }      
        
        mhRatio = mh_ratio;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
List update_s_cpp(arma::mat s1, arma::mat s2,
                  double theta1, double theta2, double theta12,
                  double beta1, double beta2,
                  arma::mat CH1, arma::mat CH2,
                  int N1, int N2, int D1, int D2,
                  double p0_1, double sigma_1,
                  double p0_2, double sigma_2,
                  arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                  double lambda_movedPoints1, double lambda_movedPoints2,
                  double lambda_newPoints1, double lambda_newPoints2,
                  double a1, double b1, double a2, double b2,
                  arma::mat traps, double R,
                  double S_area1, double S_area2,
                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double q = 1;
  double p = .5;
  double mh_ratio = 0;
  
  // FIRST PROCESS
  
  ProposePointProcess(s1, s2, theta1, theta2, theta12,
                      beta1, CH1, p0_1, sigma_1, N1, N2, D1,
                      Sigma_prop, Sigma_newpoint, 
                      lambda_movedPoints1, lambda_newPoints1, 
                      S_area1, S_usage, 
                      a1, b1, a2, b2,
                      traps, R,
                      trapsX, trapsY,
                      q, p, mh_ratio);
  
  // SECOND PROCESS
  
  ProposePointProcess(s2, s1, theta2, theta1, theta12,
                      beta2, CH2, p0_2, sigma_2, N2, N1, D2,
                      Sigma_prop, Sigma_newpoint, 
                      lambda_movedPoints2, lambda_newPoints2, 
                      S_area2, S_usage, 
                      a1, b1, a2, b2,
                      traps, R,
                      trapsX, trapsY,
                      q, p, mh_ratio);
  
  return(List::create(_["s1"] = s1,
                      _["s2"] = s2,
                      _["N1"] = N1,
                      _["N2"] = N2));
}

// [[Rcpp::export]]
List update_N_cpp(arma::mat s1, arma::mat s2,
                  double theta1, double theta2, double theta12,
                  double beta1, double beta2,
                  arma::mat CH1, arma::mat CH2,
                  int N1, int N2, int D1, int D2,
                  double p0_1, double sigma_1,
                  double p0_2, double sigma_2,
                  arma::mat Sigma_prop, arma::mat Sigma_newpoint,
                  double lambda_movedPoints1, double lambda_movedPoints2,
                  double lambda_newPoints1, double lambda_newPoints2,
                  double a1, double b1, double a2, double b2,
                  arma::mat traps, double R,
                  double S_area1, double S_area2,
                  arma::vec S_usage, arma::vec trapsX, arma::vec trapsY){
  
  double q = 0;
  double p = .5;
  
  double mh_ratio1 = 0;
  double mh_ratio2 = 0;
  
  // FIRST PROCESS
  
  ProposePointProcess(s1, s2, theta1, theta2, theta12,
                      beta1, CH1, p0_1, sigma_1, N1, N2, D1,
                      Sigma_prop, Sigma_newpoint, 
                      lambda_movedPoints1, lambda_newPoints1, 
                      S_area1, S_usage, 
                      a1, b1, a2, b2,
                      traps, R,
                      trapsX, trapsY, 
                      q, p, mh_ratio1);
  
  // SECOND PROCESS
  
  ProposePointProcess(s2, s1, theta2, theta1, theta12,
                      beta2, CH2, p0_2, sigma_2, N2, N1, D2,
                      Sigma_prop, Sigma_newpoint, 
                      lambda_movedPoints2, lambda_newPoints2, 
                      S_area2, S_usage, 
                      a1, b1, a2, b2,
                      traps, R,
                      trapsX, trapsY, 
                      q, p, mh_ratio2);
  
  return(List::create(_["s1"] = s1,
                      _["s2"] = s2,
                      _["N1"] = N1,
                      _["N2"] = N2,
                      _["mh_ratio1"] = mh_ratio1,
                      _["mh_ratio2"] = mh_ratio2));
}


/// OLD 


// [[Rcpp::export]]
arma::cube computeNewPointsDensity(arma::vec gridLength_x, arma::vec gridLength_y,
                                   arma::mat x1, arma::mat x2, int N1, int N2,
                                   double theta1, double theta2, double theta3,
                                   double beta1, double beta2) {
  
  int n_x = gridLength_x.size();
  int n_y = gridLength_y.size();
  arma::cube newPointsDensity;
  newPointsDensity.zeros(2, n_x, n_y);
  
  for(int i = 0; i < n_x; i++){
    for(int j = 0; j < n_y; j++){
      
      arma::vec xnew(2);
      xnew(0) = gridLength_x[i];
      xnew(1) = gridLength_y[j];
      newPointsDensity(0, i, j) = log_f_bivsoftcore_cpp_quick_birth(xnew, 0, 0, x1, 
                       x2, N1, N2, beta1, beta2, theta1, theta2, theta3);
      
      newPointsDensity(1, i, j) = log_f_bivsoftcore_cpp_quick_birth(xnew, 1, 0, x1, 
                       x2, N1, N2, beta1, beta2, theta1, theta2, theta3);
      
    }
  }
  
  return(newPointsDensity);
}

/////////////////////////////////
///////////  HASTINGS RATIO 
/////////////////////////////////


// [[Rcpp::export]]
double log_f_bivsoftcore_cpp_ratio(arma::mat x1, arma::mat x2, int N1, int N2,
                                   double beta1, double beta2,
                                   double theta1, double theta2, double theta3,
                                   double beta1_star, double beta2_star,
                                   double theta1_star, double theta2_star, double theta3_star){
  
  double loglikelihood = 0;
  
  loglikelihood += N1 * log(beta1_star);
  loglikelihood += N2 * log(beta2_star);
  
  loglikelihood -= N1 * log(beta1);
  loglikelihood -= N2 * log(beta2);
  
  // within x1
  
  double loglikelihood1 = 0;
  
  for(int n = 1; n < N1; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      double r2 = pow(x1(n,0) - x1(n2,0), 2) + pow(x1(n,1) - x1(n2,1), 2) ;
      
      loglikelihood1 += log(1 - exp(- r2 / (theta1_star))) - log(1 - exp(- r2 / (theta1)));
      
      // double r = sqrt( pow(x1(n,0) - x1(n2,0), 2) + pow(x1(n,1) - x1(n2,1), 2) );
      // 
      // loglikelihood1 += log(1 - exp(- r*r / (theta1)));
      
    }
    
  }  
  
  // within x2
  
  double loglikelihood2 = 0;
  
  for(int n = 1; n < N2; n++){
    
    for(int n2 = 0; n2 < n; n2++){
      
      // double r = sqrt( pow(x2(n,0) - x2(n2,0), 2) + pow(x2(n,1) - x2(n2,1), 2) );
      // 
      // loglikelihood2 += log(1 - exp(- r*r / (theta2)));
      
      double r2 = pow(x2(n,0) - x2(n2,0), 2) + pow(x2(n,1) - x2(n2,1), 2) ;
      
      loglikelihood2 += log(1 - exp(- r2 / (theta2_star))) - log(1 - exp(- r2 / (theta2)));
      
    }
    
  }
  
  // between
  
  double loglikelihood3 = 0;
  
  for(int n = 0; n < N1; n++){
    
    for(int n2 = 0; n2 < N2; n2++){
      
      double r2 =  pow(x1(n,0) - x2(n2,0), 2) + pow(x1(n,1) - x2(n2,1), 2) ;
      
      loglikelihood3 += log(1 - exp(- r2 / (theta3_star))) - log(1 - exp(- r2 / (theta3)));
      
      // double r = sqrt( pow(x1(n,0) - x2(n2,0), 2) + pow(x1(n,1) - x2(n2,1), 2) );
      // 
      // loglikelihood3 += log(1 - exp(- r*r / (theta3)));
      
    }
    
  }
  
  
  return(loglikelihood + loglikelihood1 + loglikelihood2 + loglikelihood3);
}

// [[Rcpp::export]]
double importanceSamplingEstimate(arma::cube &x_all1,
                                  arma::cube &x_all2,
                                  NumericVector params_current, 
                                  NumericVector params_star, 
                                  NumericMatrix N_all){
  
  int M = N_all.nrow();
  
  double expNumDen = 0;
  
  double beta1 = params_current[0];
  double beta2 = params_current[1];
  double theta1 = params_current[2];
  double theta2 = params_current[3];
  double theta12 = params_current[4];
  
  double beta1_star = params_star[0];
  double beta2_star = params_star[1];
  double theta1_star = params_star[2];
  double theta2_star = params_star[3];
  double theta12_star = params_star[4];
  
  for (int i = 0; i < M; i++) {
    
    int N1 = N_all(i, 0);
    int N2 = N_all(i, 1);
    
    arma::mat x1 = x_all1.subcube(arma::span(i), arma::span(), arma::span());
    arma::mat x2 = x_all2.subcube(arma::span(i), arma::span(), arma::span());
    
    double den = log_f_bivsoftcore_cpp(x1, x2, N1, N2,
                                       beta1, beta2,
                                       theta1, theta2, theta12);
    
    double num = log_f_bivsoftcore_cpp(x1, x2, N1, N2,
                                       beta1_star, beta2_star,
                                       theta1_star, theta2_star, theta12_star);
    
    // double num_minus_den = log_f_bivsoftcore_cpp_ratio(x1, x2, N1, N2, 
    //                                                    beta1, beta2,
    //                                                    theta1, theta2, theta12,
    //                                                    beta1_star, beta2_star,
    //                                                    theta1_star, theta2_star, theta12_star);
    
    // expNumDen += exp(num_minus_den) / M;
    expNumDen += exp(num - den) / M;
    
  }
  
  return expNumDen;
}


double log_prior_cpp(double theta1, double theta2, 
                     double theta12, double beta1, 
                     double beta2,
                     double mu_logtheta, double sd_logtheta, 
                     double a_beta, double b_beta){
  
  return (R::dnorm(log(theta1), mu_logtheta, sd_logtheta, 1) + 
          R::dnorm(log(theta2), mu_logtheta, sd_logtheta, 1) + 
          R::dnorm(log(theta12), mu_logtheta, sd_logtheta, 1) + 
          R::dgamma(beta1, a_beta, 1 / b_beta, 1) + 
          R::dgamma(beta2, a_beta, 1 / b_beta, 1));
  
}

struct ImpSampling : public Worker {
  
  // inputs
  const arma::cube X1;
  const arma::cube X2;
  const NumericVector params_current;
  const NumericVector params_star;
  const NumericMatrix N_all;
  
  // output matrix to write to
  RVector<double> expNumDen;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ImpSampling(arma::cube X1,
              arma::cube X2,
              NumericVector params_current,
              NumericVector params_star,
              NumericMatrix N_all,
              NumericVector expNumDen)
    : X1(X1), X2(X2),
      params_current(params_current),
      params_star(params_star),
      N_all(N_all),
      expNumDen(expNumDen) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      
      double beta1 = params_current[0];
      double beta2 = params_current[1];
      double theta1 = params_current[2];
      double theta2 = params_current[3];
      double theta12 = params_current[4];
      
      double beta1_star = params_star[0];
      double beta2_star = params_star[1];
      double theta1_star = params_star[2];
      double theta2_star = params_star[3];
      double theta12_star = params_star[4];
      
      int N1 = N_all(i, 0);
      int N2 = N_all(i, 1);
      
      arma::mat x1 = X1.subcube(arma::span(i), arma::span(), arma::span());
      arma::mat x2 = X2.subcube(arma::span(i), arma::span(), arma::span());
      
      // double den = log_f_bivsoftcore_cpp(x1, x2, N1, N2,
      //                                    beta1, beta2,
      //                                    theta1, theta2, theta12);
      //
      // double num = log_f_bivsoftcore_cpp(x1, x2, N1, N2,
      //                                    beta1_star, beta2_star,
      //                                    theta1_star, theta2_star, theta12_star);
      
      double numMinusDen = log_f_bivsoftcore_cpp_ratio(x1, x2, N1, N2,
                                                       beta1, beta2,
                                                       theta1, theta2, theta12,
                                                       beta1_star, beta2_star,
                                                       theta1_star, theta2_star, theta12_star);
      
      // expNumDen[i] = exp(num - den);
      expNumDen[i] = exp(numMinusDen);
      
    }
  }
};

// [[Rcpp::export]]
NumericVector importanceSamplingEstimateParallel(arma::cube &x_all1,
                                                 arma::cube &x_all2,
                                                 NumericVector params_current,
                                                 NumericVector params_star,
                                                 NumericMatrix N_all) {
  
  // allocate the matrix we will return
  NumericVector expNumDen(N_all.nrow());
  
  // create the worker
  ImpSampling impSampling(x_all1, x_all2, params_current, params_star, N_all, expNumDen);
  
  // call it with parallelFor
  parallelFor(0, expNumDen.length(), impSampling);
  
  return expNumDen;
}

// [[Rcpp::export]]
double hastings_ratio_cpp(arma::mat data1,
                          arma::mat data2,
                          int N1,
                          int N2,
                          arma::cube &x_all1,
                          arma::cube &x_all2,
                          NumericVector params_current,
                          NumericVector params_star,
                          NumericMatrix N_all,
                          double mu_logtheta,
                          double sd_logtheta,
                          double a_beta,
                          double b_beta) {
  
  double beta1 = params_current[0];
  double beta2 = params_current[1];
  double theta1 = params_current[2];
  double theta2 = params_current[3];
  double theta12 = params_current[4];
  
  double beta1_star = params_star[0];
  double beta2_star = params_star[1];
  double theta1_star = params_star[2];
  double theta2_star = params_star[3];
  double theta12_star = params_star[4];
  
  double ratio_constants = 1 /  mean(importanceSamplingEstimateParallel(x_all1, x_all2,
                                                                        params_current,
                                                                        params_star,
                                                                        N_all));
  
  double den = log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1, beta2, theta1, theta2, theta12) +
    log_prior_cpp(theta1, theta2, theta12, beta1, beta2, mu_logtheta, sd_logtheta, a_beta, b_beta);
  
  double num = log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1_star, beta2_star, theta1_star, theta2_star, theta12_star) +
    log_prior_cpp(theta1_star, theta2_star, theta12_star, beta1_star, beta2_star, mu_logtheta, sd_logtheta, a_beta, b_beta);
  
  return exp(num - den) * ratio_constants;
}

// [[Rcpp::export]]
double hastings_ratio_cpp_basic(arma::mat data1, 
                                arma::mat data2, 
                                int N1, 
                                int N2,
                                arma::cube &x_all1,
                                arma::cube &x_all2,
                                NumericVector params_current, 
                                NumericVector params_star, 
                                NumericMatrix N_all,
                                double mu_logtheta,
                                double sd_logtheta,
                                double a_beta,
                                double b_beta) {
  
  double beta1 = params_current[0];
  double beta2 = params_current[1];
  double theta1 = params_current[2];
  double theta2 = params_current[3];
  double theta12 = params_current[4];
  
  double beta1_star = params_star[0];
  double beta2_star = params_star[1];
  double theta1_star = params_star[2];
  double theta2_star = params_star[3];
  double theta12_star = params_star[4];
  
  double ratio_constants = 1 /  importanceSamplingEstimate(x_all1, x_all2, 
                                                           params_current,
                                                           params_star,
                                                           N_all);
  
  
  
  // double den = log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1, beta2, theta1, theta2, theta12) + 
  //   log_prior_cpp(theta1, theta2, theta12, beta1, beta2, mu_logtheta, sd_logtheta, a_beta, b_beta);
  // 
  // double num = log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1_star, beta2_star, theta1_star, theta2_star, theta12_star) + 
  //   log_prior_cpp(theta1_star, theta2_star, theta12_star, beta1_star, beta2_star, mu_logtheta, sd_logtheta, a_beta, b_beta);
  
  double num_minus_den = 
    log_f_bivsoftcore_cpp_ratio(data1, data2, N1, N2, 
                                beta1, beta2,
                                theta1, theta2, theta12,
                                beta1_star, beta2_star,
                                theta1_star, theta2_star, theta12_star) + 
                                  log_prior_cpp(theta1_star, theta2_star, theta12_star, beta1_star, beta2_star, mu_logtheta, sd_logtheta, a_beta, b_beta) - 
                                  log_prior_cpp(theta1, theta2, theta12, beta1, beta2, mu_logtheta, sd_logtheta, a_beta, b_beta);
  
  Rcout << num_minus_den << " - " << ratio_constants << std::endl;
  
  return exp(num_minus_den) * ratio_constants;
}




/////////////
// UPDATE X 
/////////////

// struct SimSoftCore : public Worker {
//   
//   // inputs
//   const arma::cube X1;
//   const arma::cube X2;
//   const NumericVector params_current;
//   const NumericVector params_star;
//   const NumericMatrix N_all;
//   
//   // output matrix to write to
//   RVector<double> expNumDen;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   SimSoftCore(arma::cube X1, 
//               arma::cube X2,
//               NumericMatrix N_all,
//               NumericVector params, 
//               NumericVector expNumDen)
//     : X1(X1), X2(X2), 
//       params_current(params_current), 
//       params_star(params_star),
//       N_all(N_all), 
//       expNumDen(expNumDen) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       int N1 = N_all(m,0);
//       int N2 = N_all(m,1);
//       
//       arma::mat x1 = x_all1.subcube(arma::span(m), arma::span(), arma::span());
//       arma::mat x2 = x_all2.subcube(arma::span(m), arma::span(), arma::span());
//       
//       List list_sims = simulate_bivsoftcore_from_startingpoint(x1, x2, N1, N2, 
//                                                                theta1, theta2, theta3, beta1, beta2,
//                                                                niter, Sigma_prop, Sigma_newpoint,
//                                                                traps, R,
//                                                                a1, b1, a2, b2);
//       outSims[i] = list_sims;
//       
//       // N1 = list_sims["N1"];
//       // N2 = list_sims["N2"];
//       // x1 = as<arma::mat>(list_sims["data1"]);
//       // x2 = as<arma::mat>(list_sims["data2"]);
//       // x_all1.subcube(arma::span(m),arma::span(),arma::span()) = x1;
//       // x_all2.subcube(arma::span(m),arma::span(),arma::span()) = x2;
//       // N_all(m,0) = N1;
//       // N_all(m,1) = N2;
//       
//     }
//   }
// }; 
// 
// // [[Rcpp::export]]
// List updateXparallel(arma::cube x_all1,
//                      arma::cube x_all2,
//                      NumericVector params_current,
//                      NumericVector params_star,
//                      NumericMatrix N_all) {
// 
//   // allocate the matrix we will return
//   List outSims(N_all.nrow());
// 
//   // create the worker
//   SimSoftCore simSoftCore(x_all1, x_all2, params_current, params_star, N_all, expNumDen);
// 
//   // call it with parallelFor
//   parallelFor(0, outSims.length(), simSoftCore);
// 
//   return outSims;
// }
// 
// // [[Rcpp::export]]
// List update_x_all_cpp(arma::cube x_all1, arma::cube x_all2, arma::mat N_all, int niter,
//                       double theta1, double theta2, double theta3, double beta1, double beta2,
//                       arma::mat Sigma_prop, arma::mat Sigma_newpoint,
//                       arma::mat traps, double R,
//                       double a1, double b1, double a2, double b2){
//   
//   int M = x_all1.n_rows;
//   
//   for(int m = 0; m < M; m++){
//     
//     int N1 = N_all(m,0);
//     int N2 = N_all(m,1);
//     
//     arma::mat x1 = x_all1.subcube(arma::span(m), arma::span(), arma::span());
//     arma::mat x2 = x_all2.subcube(arma::span(m), arma::span(), arma::span());
//     
//     List list_sims = simulate_bivsoftcore_from_startingpoint(x1, x2, N1, N2, 
//                                                              theta1, theta2, theta3, beta1, beta2,
//                                                              niter, Sigma_prop, Sigma_newpoint,
//                                                              traps, R,
//                                                              a1, b1, a2, b2);
//     N1 = list_sims["N1"];
//     N2 = list_sims["N2"];
//     x1 = as<arma::mat>(list_sims["data1"]);
//     x2 = as<arma::mat>(list_sims["data2"]);
//     x_all1.subcube(arma::span(m),arma::span(),arma::span()) = x1;
//     x_all2.subcube(arma::span(m),arma::span(),arma::span()) = x2;
//     N_all(m,0) = N1;
//     N_all(m,1) = N2;
//     
//   }
//   
//   return(List::create(_["x_all1"] = x_all1,
//                       _["x_all2"] = x_all2,
//                       _["N_all"] = N_all));
// }
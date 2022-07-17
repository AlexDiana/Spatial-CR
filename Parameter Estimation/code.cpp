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
    
    if(iter % 1000 == 0) Rcout << iter << std::endl;
    
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
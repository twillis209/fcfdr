#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec approxExtrap_rcpp(arma::vec x, arma::vec y, arma::vec xout) {
  if(!x.is_sorted()) {
    stop("Input vector x is not sorted in ascending order");
  }

  arma::vec yout(xout.size());

  for(int i = 0; i < xout.size(); ++i) {
    auto upper = std::lower_bound(x.begin(), x.end(), xout(i));

    // Handle boundaries
    if(upper == x.begin()) {
      yout(i) = (y(1) - y(0))/(x(1) - x(0))*(xout(i)-x(0))+y(0);
    } else if(upper == x.end()) {
      yout(i) = (y(y.size()-1) - y(y.size()-2))/(x(x.size()-1) - x(x.size()-2))*(xout(i)-x(x.size()-2))+y(y.size()-2);
    } else {
      uword x1_idx = std::distance(x.begin(), upper);
      double x1 = *upper;
      double x0 = x(x1_idx-1);
      double y1 = y(x1_idx);
      double y0 = y(x1_idx-1);

      yout(i) = y0+((y1-y0)/(x1-x0))*(xout(i)-x0);
    }
  }

  return yout;
}

//' @title Linear interpolation function
//'
//' @description This function carries out linear interpolation, reproducing a subset of the behaviour of stats::approx in R.
//'
//' @details Linear interpolation is performed with f(x)=f(x0)+[(f(x1)-f(x0))/(x1-x0)]*(x-x0). x is assumed to be sorted; if sorted in descending order, the order of elements in x and y is reversed. If a value of xout lies outside the range [min(x), max(x)], the value is interpolated using the nearest extremum (this corresponds to the behaviour of stats::approx with method=2).
//'
//' @param x x coordinates of points to be interpolated
//' @param y y coordinates of points to be interpolated
//' @param xout x coordiantes at which to interpolate
//'
//' @return interpolated values 
//' 
//' @author Tom Willis
// [[Rcpp::export]]
arma::vec approxfun_rcpp(arma::vec x, arma::vec y, arma::vec xout) {

  if(!x.is_sorted()) {
    stop("Input vector x is not sorted in ascending order");
  }

  arma::vec yout(xout.size());

  for(int i = 0; i < xout.size(); ++i) {
    auto upper = std::lower_bound(x.begin(), x.end(), xout(i));

    // Handle boundaries
    if(upper == x.begin()) {
      yout(i) = y(0);
      continue;
    } else if(upper == x.end()) {
      yout(i) = y(y.size()-1);
      continue;
    }

    uword x1_idx = std::distance(x.begin(), upper);
    double x1 = *upper;
    double x0 = x(x1_idx-1);
    double y1 = y(x1_idx);
    double y0 = y(x1_idx-1);

    yout(i) = y0+((y1-y0)/(x1-x0))*(xout(i)-x0);
  }

  return yout;

}

// Define the Rcpp function
// [[Rcpp::export]]
arma::vec per_group_binary_cfdr(arma::vec p_loo, arma::vec q_loo, arma::vec ps, arma::vec qs, arma::vec x) {
//    int n = p_loo.size();
//    int m = ps.size();
//
//    double q0 = (double) arma::intersect(find(q_loo == 1), find(p_loo > 0.5)).size() / (double) find(p_loo > 0.5).size();
//    double mult = (double) arma::intersect(find(q_loo == 0), find(p_loo > 0.5)).size() / (double) arma::intersect(find(q_loo == 1), find(p_loo > 0.5)).size();
//
//    arma::vec q0_sol(m);
//    arma::vec q1_sol(m);
//
//    for (int i = 0; i < m; ++i) {
//        double p = ps(i);
//        q0_sol(i) = arma::max(arma::intersect(find(p_loo <= p), find(q_loo == 0)).size(), 1);
//        q1_sol(i) = arma::max(arma::intersect(find(p_loo <= p), find(q_loo == 1)).size(), 1);
//    }
//
//    arma::vec sol = arma::zeros<arma::vec>(m);
//
//    for (int i = 0; i < m; ++i) {
//        if (qs(i) == 0) {
//            sol(i) = mult * ps(i) / q0_sol(i);
//        } else {
//            sol(i) = (1. / mult) * ps(i) / q1_sol(i);
//        }
//    }
//
//    arma::vec y(x.size());
//
//    for(int i = 0; i < x.size(); ++i) {
//      y(i) = (double) x(i) / (double) arma::max(arma::intersect(arma::find(p_loo <= x(i)), arma::find(q_loo == 0)).size(), 1);
//    }
//
//    arma::vec extr_x = sol.elem(arma::find_unique(sol));
//    arma::vec extr_y = approxExtrap_rcpp(y, x, extr_x);
//
//    arma::vec invg0_y = arma::clamp(arma::max(arma::min(extr_y, 1.0), 0.0), 0.0, 1.0);
//    arma::vec invg0_x = extr_x;
//
//    arma::vec y1 = x / arma::max(arma::mat(p_loo <= x, n, m) * arma::mat(q_loo == 1, n, 1), 1);
//    arma::vec extr1_x = unique(sol);
//    arma::vec extr1_y = approxExtrap_rcpp(y1, x, extr_x);
//    arma::vec invg1_y = arma::clamp(arma::max(arma::min(extr1_y, 1.0), 0.0), 0.0, 1.0);
//    arma::vec invg1_x = extr1_x;
//
//    arma::vec p1 = arma::zeros<arma::vec>(m);
//    arma::vec p0 = arma::zeros<arma::vec>(m);
//
//    for (int i = 0; i < m; ++i) {
//        if (qs(i) == 0) {
//          p1(i) = Rcpp::as<double>(approxfun_rcpp(invg1_x, invg1_y, sol(i)));
//        } else {
//            p1(i) = ps(i);
//        }
//
//        if (qs(i) == 1) {
//          p0(i) = Rcpp::as<double>(approxfun_rcpp(invg0_x, invg0_y, sol(i)));
//        } else {
//            p0(i) = ps(i);
//        }
//    }
//
//    return p0 * (1 - q0) + p1 * q0;
}


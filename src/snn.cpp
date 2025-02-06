#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <iomanip>
#include <numeric>  // for std::iota
#include <algorithm>  // for std::stable_sort

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0);
  return SNN;
}

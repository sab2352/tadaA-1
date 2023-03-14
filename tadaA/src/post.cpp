#include<Rcpp.h>
using namespace Rcpp;
using namespace std;
//' @export
// [[Rcpp::export]]
NumericVector post(NumericVector index,NumericVector Zg1,NumericVector Zg0,NumericVector prior) //aggregate 69891 into 17478 items
{
  int N=index.length();
  //cout<<"N:"<<N<<endl;
  int latter=0;int former=0;
  int i; int j;//long double res1=0;
  NumericVector prior_new=NumericVector(N-1);
  for(i=0;i<N-1;i++)
  {
    long double temp0=0;
    long double temp1=0;
    former=index[i]-1;        //这次的位置
    latter=index[i+1]-1;      //下次的位置
    for(j=former;j<latter;j++)
    {
      temp1=temp1+Zg1[j];
      temp0=temp0+Zg0[j];
    }
    prior_new[i]=exp(log(prior[i])+temp1)/(exp(log(prior[i])+temp1)+exp(log(1-prior[i])+temp0));
  }
  return prior_new;
}



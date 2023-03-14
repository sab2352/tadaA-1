#include<Rcpp.h>
using namespace Rcpp;
using namespace std;
//' @export
// [[Rcpp::export]]
long double ll_sum(NumericVector index,NumericVector Zg1,NumericVector Zg0,NumericVector prior)
{
  int N=index.length()-1;
  int latter=0;int former=0;
  int i; int j;long double res=0;//long double res1=0;
  for(i=0;i<N;i++)
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
    res=res+log(exp(log(prior[i])+temp1)+exp(log(1-prior[i])+temp0));
  }
//  cout<<"res:"<<res<<endl;
  return res;
}

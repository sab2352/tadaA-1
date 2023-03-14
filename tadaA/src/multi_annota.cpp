#include <Rcpp.h>
#include<chrono>
using namespace Rcpp;
using namespace std;
//' @export
// [[Rcpp::export]]
long double cal1(List x,NumericVector selected_annotations,NumericVector all_rr)  //计算一个annotation的logP_Zg1
{
  int i;long double temp=0;int index=0;
  NumericVector index1=as<NumericVector>(x["feature_vector"]);
  long double index2=as<long double>(x["sum_mut_rate_count"]);
  long double index3=as<long double>(x["sum_mut_rate"]);
  long double index4=as<long double>(x["sum_mut_count"]);
  long double index5=as<long double>(x["log_fcount"]);
  long double last;
  int N=selected_annotations.length();
  for(i=0;i<N;i++)
  {
    index=selected_annotations[i];
    temp=temp+index1[index-1]*all_rr[i];
  }
  last=index2+index4*temp-index3*exp(temp)-index5;
  return last;
}
//' @export
// [[Rcpp::export]]
long double cal2(List x)  //计算一个annotation的logP_Zg0
{
  long double index2=as<long double>(x["sum_mut_rate_count"]);
  long double index3=as<long double>(x["sum_mut_rate"]);
  long double index5=as<long double>(x["log_fcount"]);
  long double result=index2-index3-index5;
  return result;
}
//' @export
// [[Rcpp::export]]
NumericVector sumall1(List DATA,NumericVector selected_annotations,NumericVector all_rr) //合并所有annotations的logP_Zg1
{
  int N;
  long double res;
  List gene;
  List annotation;
  int i;
  int j;
  //auto start=chrono::system_clock::now();
  NumericVector Res=NumericVector(DATA.length());//69891
  for(i=0;i<DATA.length();i++)
  {
    gene=DATA[i];
    N=gene.length();
    res=0;
    for(j=0;j<N;j++)
    {
      annotation=gene[j];
      res=res+cal1(annotation,selected_annotations,all_rr);
    }
    Res[i]=res;
  }
  //auto end = std::chrono::system_clock::now();
  //auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  //cout <<0.001*elapsed.count() <<"s" << '\n';
  return Res;
}
//' @export
// [[Rcpp::export]]
NumericVector sumall0(List DATA) //合并所有annotations的logP_Zg
{
  int N;
  long double res;
  List gene;
  List annotation;
  int i;int j;
  //auto start=chrono::system_clock::now();
  NumericVector Res=NumericVector(DATA.length());//69891
  for(i=0;i<DATA.length();i++)
  {
    gene=DATA[i];
    N=gene.length();
    res=0;
    for(j=0;j<N;j++)
    {
      annotation=gene[j];
      res=res+cal2(annotation);
    }
    Res[i]=res;
  }
  //auto end = std::chrono::system_clock::now();
  //auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  //cout <<0.001*elapsed.count() <<"s" << '\n';
  return Res;
}

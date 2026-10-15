#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Exact trimmed-mean operation for nonnegative finite expression, evaluated
// across the SAME permutation indices supplied by CellChat. No RNG is used.
// [[Rcpp::export]]
NumericMatrix fast_bootstrap_means(NumericMatrix x, IntegerVector group,
                                   IntegerMatrix permutations, double trim) {
  int genes=x.nrow(), cells=x.ncol(), B=permutations.ncol();
  int G=*std::max_element(group.begin(),group.end());
  if(G>20 || permutations.nrow()!=cells || group.size()!=cells)stop("Invalid dimensions");
  std::vector<int> sizes(G,0), cut(G,0), keep(G,0);
  for(int c=0;c<cells;c++)sizes[group[c]-1]++;
  for(int g=0;g<G;g++){cut[g]=(int)(trim*sizes[g]);keep[g]=sizes[g]-2*cut[g];if(keep[g]<=0)stop("Empty trimmed mean");}
  typedef std::pair<double,int> Entry;
  std::vector<std::vector<Entry> > nz(genes);
  for(int c=0;c<cells;c++)for(int j=0;j<genes;j++){
    double v=x(j,c);if(!R_finite(v)||v<0)stop("Requires finite nonnegative expression");
    if(v>0)nz[j].push_back(Entry(v,c));
  }
  for(int j=0;j<genes;j++)std::sort(nz[j].begin(),nz[j].end());
  NumericMatrix out(genes,G*B);
  std::vector<int> boot(cells),npositive(G),position(G),lo(G),hi(G),zeros(G);
  std::vector<long double> sums(G),means(G),correction(G);
  for(int b=0;b<B;b++){
    if(b%100==0)checkUserInterrupt();
    for(int c=0;c<cells;c++)boot[c]=group[permutations(c,b)-1]-1;
    for(int j=0;j<genes;j++){
      std::fill(npositive.begin(),npositive.end(),0);
      for(const Entry &e:nz[j])npositive[boot[e.second]]++;
      for(int g=0;g<G;g++){
        int nzero=sizes[g]-npositive[g];
        lo[g]=std::max(0,cut[g]-nzero);
        hi[g]=std::max(0,sizes[g]-cut[g]-nzero);
        zeros[g]=keep[g]-(hi[g]-lo[g]);position[g]=0;sums[g]=0;
      }
      for(const Entry &e:nz[j]){int g=boot[e.second],k=position[g]++;if(k>=lo[g]&&k<hi[g])sums[g]+=(long double)e.first;}
      for(int g=0;g<G;g++){
        means[g]=sums[g]/keep[g];correction[g]=0;position[g]=0;
        // Match R mean's second accumulation over its sorted trimmed vector.
        for(int k=0;k<zeros[g];k++)correction[g]-=means[g];
      }
      for(const Entry &e:nz[j]){int g=boot[e.second],k=position[g]++;if(k>=lo[g]&&k<hi[g])correction[g]+=(long double)e.first-means[g];}
      for(int g=0;g<G;g++)out(j,G*b+g)=(double)(means[g]+correction[g]/keep[g]);
    }
  }
  rownames(out)=rownames(x);
  return out;
}

#include "TH1.h"
#include "TRandom.h"

class RejectionSampler {
  public:
  RejectionSampler(TH1F* dist);
  ~RejectionSampler(){}
  float getFloatValue();
  int getIntValue();
  private:
  float getValue();
  TH1F* d;
  TRandomMixMax r;
  bool int_flag;
};

RejectionSampler::RejectionSampler(TH1F* dist){
  d = dist;
}

float RejectionSampler::getFloatValue(){
  int_flag = false;
  return getValue();
}

int RejectionSampler::getIntValue(){
  int_flag = true;
  return round(getValue());
}

float RejectionSampler::getValue(){
  // form box
  float min_x = d->GetBinLowEdge(1); // ignore underflow bin
  int nbins = d->GetNbinsX();
  float max_x = d->GetBinLowEdge(nbins)+d->GetBinWidth(nbins); // ignore overflow bin
  if(d->GetMinimum()<0) std::cout << "RejectionSampler: !!WARNING!! Negative bin content detected. Sampler may not work correctly!" << std::endl;
  float min_y = 0.; // rejection sampling cannot deal with negative bin content
  float max_y = d->GetMaximum();
  
  
  float vx;
  float vy;
  bool found = false;
  while(!found){
    // generate random point
    if(int_flag){
      int min_int = round(min_x);
      int randint = r.Integer(round(max_x))+min_int;
      vx = (float)randint;
    }
    else{
      vx = r.Uniform(max_x)+min_x;
    }
    vy = r.Uniform(max_y);
  
    // check whether point is inside histogram
    float hval = d->GetBinContent(d->FindBin(vx));
    if(vy<hval) found = true;
  }
  return vx;
}

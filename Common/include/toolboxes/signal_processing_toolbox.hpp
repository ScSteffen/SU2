#pragma once

#include "../datatype_structure.hpp"
#include "../../include/option_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(std::vector<su2double> &data);
  
  class RunningAverage{
  protected:
    su2double val;
    unsigned long count;

  public:
    RunningAverage(){
      this->Reset();
    }
    virtual ~RunningAverage(){}
 
    virtual su2double Update(su2double valIn){ //su2double Signal_Processing::RunningAverage::Update(su2double valIn)
      su2double scaling = 1. / (su2double)(count + 1);
      val = valIn * scaling + val * (1. - scaling);
      count++;
      return val;
    }

    su2double Get(){
      return val;
    }
  
    unsigned long Count(){
      return count;
    }
  
    void Reset(){
      val = 0.;
      count = 0;
    }
  };

  class RunningWindowedAverage:public RunningAverage{
  private:
    std::vector<su2double> values;
    su2double timeStep;
    int functionIndex;

   public:
    RunningWindowedAverage(su2double timeStepIn, int functionIndexIn){
      timeStep = timeStepIn;
      functionIndex = functionIndexIn;
    }

    su2double Update(su2double valIn){

      values.push_back(valIn);

      switch (functionIndex){
        case 0: return HannWindowing();
        case 1: return HannSquaredWindowing();
        case 2: return BumpWindowing();
      }

        return 0;
      }

    su2double HannWindowing(){
      su2double wnd_timeAvg = 0.0;
      for(unsigned i=0; i<count; i++){
          wnd_timeAvg+=timeStep*values[i]*HannWindow( (su2double)i);
        }
      return wnd_timeAvg/(su2double) count;
    }

    su2double HannSquaredWindowing(){
      su2double wnd_timeAvg = 0.0;
      for(unsigned i=0; i<count; i++){
          wnd_timeAvg+=timeStep*values[i]*HannSquaredWindow( (su2double)i);
        }
      return wnd_timeAvg/(su2double) count;
    }

    su2double BumpWindowing(){
      su2double wnd_timeAvg = 0.0;
      for(unsigned i=0; i<count; i++){
          wnd_timeAvg+=timeStep*values[i]*BumpWindow( (su2double)i);
        }
      return wnd_timeAvg/(su2double) count;
    }

    su2double HannWindow(su2double i){
      return 1-cos(2*PI_NUMBER*i/(su2double)count);
    }


    su2double HannSquaredWindow(su2double i){
      return 2/3*(1-cos(2*PI_NUMBER*i/(su2double)count))*(1-cos(2*PI_NUMBER*i/(su2double)count));
    }

    su2double BumpWindow(su2double i){
      su2double tau = i/(su2double) count;
      return 1/0.00702986*(exp(-1/(tau-tau*tau)));
    }

  };
};

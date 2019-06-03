#pragma once

#include "../datatype_structure.hpp"
#include "../../include/option_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(std::vector<su2double> &data);
  
  class RunningAverage{
  private:
    su2double val;
    su2double wndVal;
    unsigned long count;
    std::vector<su2double> values;
    su2double timeStep;
    //int functionIndex;

  public:
    RunningAverage(){  this->Reset();} //timestep and functionIndex only mandatory for windowed averaging

    su2double GetVal();
    su2double GetWndVal();
    unsigned long Count();
    void      Reset();
    su2double Update(su2double valIn);         //Computes the arithmetic mean over the output values
    su2double WindowedUpdate(su2double valIn, su2double timeStep, int fctIdx); //Computes a (windowed) time average (integral)
    su2double NoWindowing();
    su2double HannWindowing();
    su2double HannSquaredWindowing();
    su2double BumpWindowing();
    su2double HannWindow(su2double i);
    su2double HannSquaredWindow(su2double i);
    su2double BumpWindow(su2double i);
  };

  /*
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
        case 1: return HannWindowing();
        case 2: return HannSquaredWindowing();
        case 3: return BumpWindowing();
        default: return NoWindowing();
      }

        return 0;
      }

    su2double NoWindowing(){
      su2double wnd_timeAvg = 0.0;
      for(unsigned i=0; i<count; i++){
          wnd_timeAvg+=timeStep*values[i];
        }
      return wnd_timeAvg/(su2double) count;
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

  }; */
};

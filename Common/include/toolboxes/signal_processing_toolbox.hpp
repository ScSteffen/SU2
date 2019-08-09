    #pragma once

#include "../datatype_structure.hpp"
#include "../../include/option_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(std::vector<su2double> &data);
  
  class RunningAverage{
  private:
    su2double val;
    su2double hannWndVal,hannSqWndVal,bumpWndVal,noWndAvgVal;
    unsigned long count;
    std::vector<su2double> values;

  public:
    RunningAverage(){  this->Reset();} //timestep and functionIndex only mandatory for windowed averaging

    su2double Update(su2double valIn){ //Computes the arithmetic mean of all values up to count
      su2double scaling = 1. / static_cast<su2double>(count +1);
      val = valIn * scaling + val * (1. - scaling);
      count++;
      return val;
    }

    su2double GetVal(){
      return val;
    }

    su2double GetWndVal(int fctIdx){
        switch (fctIdx){
          case 1: return hannWndVal;
          case 2: return hannSqWndVal;
          case 3: return bumpWndVal;
          default:return noWndAvgVal;
        }
    }

    su2double GetWndWeight(int fctIdx, su2double time, su2double endTime){ 	/*!< Window function evalutation at time with endTime */
        switch (fctIdx){
          case 1: return HannWindow(time, endTime);
          case 2: return HannSquaredWindow(time, endTime);
          case 3: return BumpWindow(time, endTime);
          default:return 1.0;
        }
    }

    unsigned long Count(){
      return count;
    }

    void Reset(){
      val         = 0.;
      count       = 0;
      hannWndVal  = 0.;
      hannSqWndVal= 0.;
      bumpWndVal  = 0.;
      noWndAvgVal = 0.;
    }

    void addValue(su2double valIn, unsigned long currentIter,unsigned long startIter = 0){
        if(currentIter >= startIter)values.push_back(valIn);
    }

    su2double WindowedUpdate(int fctIdx){ //Computes a windowed time average (integral)
      if(values.size()>1){
          switch (fctIdx){
            case 1: hannWndVal      = HannWindowing();        return hannWndVal;
            case 2: hannSqWndVal    = HannSquaredWindowing(); return hannSqWndVal;
            case 3: bumpWndVal      = BumpWindowing();        return bumpWndVal;
            default: noWndAvgVal    = NoWindowing();          return noWndAvgVal;
          }
      }
        return 0.0;
    }

  private:
    //Using Trapezoidal rule
    su2double NoWindowing(){
      su2double wnd_timeAvg = 0.0;
      for(unsigned i=0; i<values.size()-1; i++){
          wnd_timeAvg+=0.5*(values[i+1]+values[i]);
        }
      return wnd_timeAvg/static_cast<su2double>(values.size()-1);
    }

    su2double HannWindowing(){
      su2double wnd_timeAvg = 0.0;
      su2double endTime = static_cast<su2double>(values.size());
      for(unsigned i=0; i<values.size()-1; i++){
          wnd_timeAvg+=0.5*(values[i+1]*HannWindow(static_cast<su2double>(i+1),endTime)+values[i]*HannWindow(static_cast<su2double>(i), endTime));
      }
      return wnd_timeAvg/static_cast<su2double>(values.size()-1);
    }

    su2double HannSquaredWindowing(){
      su2double wnd_timeAvg = 0.0;
      su2double endTime = static_cast<su2double>(values.size());
      for(unsigned i=0; i<values.size()-1; i++){
          wnd_timeAvg+=0.5*(values[i+1]*HannSquaredWindow(static_cast<su2double>(i+1),endTime)+values[i]*HannSquaredWindow(static_cast<su2double>(i),endTime));
        }
      return wnd_timeAvg/static_cast<su2double>(values.size()-1);
    }

    su2double BumpWindowing(){
      su2double wnd_timeAvg = 0.0;
      su2double endTime = static_cast<su2double>(values.size());
      for(unsigned i=0; i<values.size()-1; i++){
          wnd_timeAvg+=0.5*(values[i+1]*BumpWindow(static_cast<su2double>(i+1),endTime)+values[i]*BumpWindow(static_cast<su2double>(i),endTime));
        }
      return wnd_timeAvg/static_cast<su2double>(values.size()-1);
    }

    su2double HannWindow(su2double i, su2double endTime){
      return 1.0-cos(2*PI_NUMBER*i/endTime);
    }

    su2double HannSquaredWindow(su2double i, su2double endTime){
     return 2.0/3.0*(1-cos(2*PI_NUMBER*i/endTime))*(1-cos(2*PI_NUMBER*i/endTime));
    }

    su2double BumpWindow(su2double i, su2double endTime){
      if(i==0.0) return 0;
      if(i==1.0) return 0;
      su2double tau = i/endTime;
      return 1.0/0.00702986*(exp(-1/(tau-tau*tau)));
    }
  };
};

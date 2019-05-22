#pragma once

#include "../datatype_structure.hpp"
#include "../../include/option_structure.hpp"
#include <vector>

namespace Signal_Processing {
  
  su2double Average(std::vector<su2double> &data);
  
  class RunningAverage{
  private:
    su2double val;
    unsigned long count;

  public:
    RunningAverage(){
      this->Reset();
    }
 
    su2double Update(su2double valIn){ //su2double Signal_Processing::RunningAverage::Update(su2double valIn)
      su2double scaling = 1. / (su2double)(count + 1);
      val = valIn * scaling + val * (1. - scaling);
      count++;
      return val;
    }

    su2double WindowedUpdate(su2double valIn, su2double endTime, int functionIndex){
        switch (functionIndex){
          case 0: valIn = valIn*HannWindow(endTime);
                break;
          case 1: valIn = valIn*HannSquaredWindow(endTime);
                break;
          case 2: valIn = valIn*BumpWindow(endTime);
            break;
          default: break;
        }

        su2double scaling = 1. / (su2double)(count + 1);
        val = valIn * scaling + val * (1. - scaling);
        count++;
        return val;
      }

    su2double HannWindow(su2double endTime){
      return 1-cos(2*PI_NUMBER*(count+1)/endTime);
    }


    su2double HannSquaredWindow(su2double endTime){
      return 2/3*(1-cos(2*PI_NUMBER*(count+1)/endTime))*(1-cos(2*PI_NUMBER*(count+1)/endTime));
    }

    su2double BumpWindow(su2double endTime){
      su2double tau = (count+1)/endTime;
      return 1/0.00702986*(exp(-1/(tau-tau*tau)));
    }

    /* TODO, separate into .cpp
    su2double WindowedUpdate(su2double valIn, su2double endTime, int functionIndex = 0);

    su2double HannWindow(su2double endTime);

    su2double HannSquaredWindow(su2double endTime);

    su2double BumpWindow(su2double endTime);
    */

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
};

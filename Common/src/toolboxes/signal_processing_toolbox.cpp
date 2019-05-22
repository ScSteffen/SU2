#include "../../include/toolboxes/signal_processing_toolbox.hpp"
#include "../../include/option_structure.hpp"

su2double Signal_Processing::Average(std::vector<su2double> &data){
  su2double avg = 0.0;
  for (std::vector<su2double>::iterator it = data.begin(); it != data.end(); it++){
    avg += (*it);
  }
  return avg/=data.size();
}

/*
su2double Signal_Processing::RunningAverage::WindowedUpdate(su2double valIn, su2double endTime, int functionIndex){
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

su2double Signal_Processing::RunningAverage::HannWindow(su2double endTime){
  return 1-cos(2*PI_NUMBER*(count+1)/endTime);
}


su2double Signal_Processing::RunningAverage::HannSquaredWindow(su2double endTime){
  return 2/3*(1-cos(2*PI_NUMBER*(count+1)/endTime))*(1-cos(2*PI_NUMBER*(count+1)/endTime));
}

su2double Signal_Processing::RunningAverage::BumpWindow(su2double endTime){
  su2double tau = (count+1)/endTime;
  return 1/0.00702986*(exp(-1/(tau-tau*tau)));
}
*/

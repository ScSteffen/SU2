#include "../../include/toolboxes/signal_processing_toolbox.hpp"
#include "../../include/option_structure.hpp"


using namespace Signal_Processing;

su2double Average(std::vector<su2double> &data){
  su2double avg = 0.0;
  for (std::vector<su2double>::iterator it = data.begin(); it != data.end(); it++){
    avg += (*it);
  }
  return avg/=data.size();
}


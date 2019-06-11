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

//RunningAverage Members
/*
su2double RunningAverage::Update(su2double valIn){ //Computes the arithmetic mean of all values up to count
  su2double scaling = 1. / static_cast<su2double>(count +1);
  val = valIn * scaling + val * (1. - scaling);
  count++;
  return val;
};

su2double RunningAverage::GetVal(){
  return val;
};

su2double RunningAverage::GetWndVal(){
  return wndVal;
}

unsigned long RunningAverage::Count(){
  return count;
}

void RunningAverage::Reset(){
  val    = 0.;
  count  = 0;
  wndVal = 0.;
}

su2double RunningAverage::WindowedUpdate(su2double valIn,su2double timeStepIn, int fctIdx){ //Computes a (windowed) time average (integral)

  values.push_back(valIn);
  this->timeStep = timeStepIn;

  switch (fctIdx){
    case 1: wndVal = HannWindowing();        return wndVal;
    case 2: wndVal = HannSquaredWindowing(); return wndVal;
    case 3: wndVal = BumpWindowing();        return wndVal;
    default: wndVal = NoWindowing();         return wndVal;
  }
    return 0;
}

su2double RunningAverage::NoWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned i=0; i<count; i++){
      wnd_timeAvg+=timeStep*values[i];
    }
  return wnd_timeAvg/static_cast<su2double>(count);
}

su2double RunningAverage::HannWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned i=0; i<count; i++){
      wnd_timeAvg+=timeStep*values[i]*HannWindow( static_cast<su2double>(i));
    }
  return wnd_timeAvg/static_cast<su2double>(count);
}

su2double RunningAverage::HannSquaredWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned i=0; i<count; i++){
      wnd_timeAvg+=timeStep*values[i]*HannSquaredWindow( static_cast<su2double>(i));
    }
  return wnd_timeAvg/static_cast<su2double>(count);
}

su2double RunningAverage::BumpWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned i=0; i<count; i++){
      wnd_timeAvg+=timeStep*values[i]*BumpWindow( static_cast<su2double>(i));
    }
  return wnd_timeAvg/static_cast<su2double>(count);
}

su2double RunningAverage::HannWindow(su2double i){
  return 1-cos(2*PI_NUMBER*i/static_cast<su2double>(count));
}

su2double RunningAverage::HannSquaredWindow(su2double i){
  return 2/3*(1-cos(2*PI_NUMBER*i/static_cast<su2double>(count)))*(1-cos(2*PI_NUMBER*i/static_cast<su2double>(count)));
}

su2double RunningAverage::BumpWindow(su2double i){
  su2double tau = i/static_cast<su2double>(count);
  return 1/0.00702986*(exp(-1/(tau-tau*tau)));
}
*/

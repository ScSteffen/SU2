/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/output/COutput.hpp"
#include "../../include/output/filewriter/CFVMDataSorter.hpp"
#include "../../include/output/filewriter/CFEMDataSorter.hpp"
#include "../../include/output/filewriter/CSurfaceFVMDataSorter.hpp"
#include "../../include/output/filewriter/CSurfaceFEMDataSorter.hpp"
#include "../../include/output/filewriter/CParaviewFileWriter.hpp"
#include "../../include/output/filewriter/CParaviewBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CCSVFileWriter.hpp"
#include "../../include/output/filewriter/CSU2FileWriter.hpp"
#include "../../include/output/filewriter/CSU2BinaryFileWriter.hpp"
#include "../../include/output/filewriter/CSU2MeshFileWriter.hpp"


#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

COutput::COutput(CConfig *config, unsigned short nDim, bool fem_output): femOutput(fem_output) {
  
  this->nDim = nDim;

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  fieldWidth = 12;
  
  convergenceTable = new PrintingToolbox::CTablePrinter(&std::cout);
  multiZoneHeaderTable = new PrintingToolbox::CTablePrinter(&std::cout);
  fileWritingTable = new PrintingToolbox::CTablePrinter(&std::cout);
  
  /*--- Set default filenames ---*/
  
  surfaceFilename = "surface";
  volumeFilename  = "volume";
  restartFilename = "restart";
  
  /*--- Retrieve the history filename ---*/
 
  historyFilename = config->GetConv_FileName();
  
  /*--- Add the correct file extension depending on the file format ---*/
  
  string hist_ext = ".csv";
  if (config->GetTabular_FileFormat() == TAB_TECPLOT) hist_ext = ".dat";
       
  /*--- Append the zone ID ---*/
  
  historyFilename = config->GetMultizone_HistoryFileName(historyFilename, config->GetiZone(), hist_ext);

  /*--- Append the restart iteration ---*/
  
  if (config->GetTime_Domain() && config->GetRestart()) {
    historyFilename = config->GetUnsteady_FileName(historyFilename, config->GetRestart_Iter(), hist_ext);
  }
  
  historySep = ","; 
  
  /*--- Initialize residual ---*/
  
  rhoResNew = EPS;
  rhoResOld = EPS;

  nRequestedHistoryFields = config->GetnHistoryOutput();
  for (unsigned short iField = 0; iField < nRequestedHistoryFields; iField++){
    requestedHistoryFields.push_back(config->GetHistoryOutput_Field(iField));
  }
  
  nRequestedScreenFields = config->GetnScreenOutput();
  for (unsigned short iField = 0; iField < nRequestedScreenFields; iField++){
    requestedScreenFields.push_back(config->GetScreenOutput_Field(iField));
  }
  
  nRequestedVolumeFields = config->GetnVolumeOutput();
  for (unsigned short iField = 0; iField < nRequestedVolumeFields; iField++){
    requestedVolumeFields.push_back(config->GetVolumeOutput_Field(iField));
  }
  
  gridMovement = config->GetGrid_Movement(); 
  
  multiZone     = config->GetMultizone_Problem();

  /*--- Default is to write history to file and screen --- */

  noWriting = false;

  cauchySerie = new su2double[config->GetCauchy_Elems()];

  convField = config->GetConv_Field();

  cauchyValue = 0.0;
  for (unsigned short iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
    cauchySerie[iCounter] = 0.0;

  convergence = false;
  TimeConvergence = false;

  WndCauchy_Serie = new su2double[config->GetWnd_Cauchy_Elems()];
  for (unsigned short iCounter = 0; iCounter < config->GetWnd_Cauchy_Elems(); iCounter++)
    WndCauchy_Serie[iCounter] = 0.0;

  WndCauchy_Value = 0.0;

  /*--- Initialize all convergence flags to false. ---*/
  
  convergence        = false;
  
  buildFieldIndexCache = false;
  
  curInnerIter = 0;
  curOuterIter = 0;
  curTimeIter  = 0;
  
  volumeDataSorter = nullptr;
  surfaceDataSorter = nullptr;
  
  headerNeeded = false;
  
}

COutput::~COutput(void) {
  
  delete convergenceTable;
  delete multiZoneHeaderTable;

  delete [] cauchySerie;
  delete [] WndCauchy_Serie;
  
    if (volumeDataSorter != nullptr)
      delete volumeDataSorter;
    
    volumeDataSorter = nullptr;
    
    if (surfaceDataSorter != nullptr)
      delete surfaceDataSorter;
    
    surfaceDataSorter = nullptr;
  
  
}



void COutput::SetHistory_Output(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CConfig *config,
                                  unsigned long TimeIter,
                                  unsigned long OuterIter,
                                  unsigned long InnerIter) {

  curTimeIter  = TimeIter;
  curAbsTimeIter = TimeIter - config->GetRestart_Iter();
  curOuterIter = OuterIter;
  curInnerIter = InnerIter;
  
  bool write_header, write_history, write_screen;
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadCommonHistoryData(config);
  
  LoadHistoryData(config, geometry, solver_container);
  
  Convergence_Monitoring(config, curInnerIter);
  
  Postprocess_HistoryData(config);

  MonitorTimeConvergence(config, curTimeIter);

  /*--- Output using only the master node ---*/
  
  if (rank == MASTER_NODE && !noWriting) {
    
    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(config);
    if (write_history) SetHistoryFile_Output(config);
    
    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(config);
    if (write_header) SetScreen_Header(config);
    
    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(config);
    if (write_screen) SetScreen_Output(config);
    
  }

}

void COutput::SetHistory_Output(CGeometry *geometry,
                                CSolver **solver_container,
                                CConfig *config) {
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadCommonHistoryData(config);
  
  LoadHistoryData(config, geometry, solver_container);
  
  Convergence_Monitoring(config, curInnerIter);  
  
  Postprocess_HistoryData(config);

}

void COutput::SetMultizoneHistory_Output(COutput **output, CConfig **config, CConfig *driver_config, unsigned long TimeIter, unsigned long OuterIter){
  
  curTimeIter  = TimeIter;
  curAbsTimeIter = TimeIter - driver_config->GetRestart_Iter();  
  curOuterIter = OuterIter;
  
  bool write_header, write_screen, write_history;
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadMultizoneHistoryData(output, config);
  
  Convergence_Monitoring(driver_config, curOuterIter);  
  
  /*--- Output using only the master node ---*/
  
  if (rank == MASTER_NODE && !noWriting) {
    
    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(driver_config);
    if (write_history) SetHistoryFile_Output(driver_config);
    
    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(driver_config);
    if (write_header) SetScreen_Header(driver_config);
    
    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(driver_config);
    if (write_screen) SetScreen_Output(driver_config);
    
  }
  
}
void COutput::SetCFL_Number(CSolver ****solver_container, CConfig *config) {
  
  su2double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
  unsigned short iMesh;
  
  unsigned short FinestMesh = config->GetFinestMesh();
  unsigned short nVar = 1;

  bool energy = config->GetEnergy_Equation();
  bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  switch( config->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES : case RANS:
    case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS:
      if (energy) {
        nVar = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetnVar();
        rhoResNew = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
      }
      else if (weakly_coupled_heat) {
        rhoResNew = solver_container[INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      }
      else {
        rhoResNew = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
      }
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
      rhoResNew = solver_container[INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
      break;
    case HEAT_EQUATION_FVM:
      rhoResNew = solver_container[INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      break;
  }
  
  if (rhoResNew < EPS) rhoResNew = EPS;
  if (rhoResOld < EPS) rhoResOld = rhoResNew;

  Div = rhoResOld/rhoResNew;
  Diff = rhoResNew-rhoResOld;

  /*--- Compute MG factor ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (iMesh == MESH_0) MGFactor[iMesh] = 1.0;
    else MGFactor[iMesh] = MGFactor[iMesh-1] * config->GetCFL(iMesh)/config->GetCFL(iMesh-1);
  }

  if (Div < 1.0) power = config->GetCFL_AdaptParam(0);
  else power = config->GetCFL_AdaptParam(1);

  /*--- Detect a stall in the residual ---*/

  if ((fabs(Diff) <= rhoResNew*1E-8) && (curInnerIter != 0)) { Div = 0.1; power = config->GetCFL_AdaptParam(1); }

  CFLMin = config->GetCFL_AdaptParam(2);
  CFLMax = config->GetCFL_AdaptParam(3);

  CFLFactor = pow(Div, power);

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    CFL = config->GetCFL(iMesh);
    CFL *= CFLFactor;

    if ((iMesh == MESH_0) && (CFL <= CFLMin)) {
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        config->SetCFL(iMesh, 1.001*CFLMin*MGFactor[iMesh]);
      }
      break;
    }
    if ((iMesh == MESH_0) && (CFL >= CFLMax)) {
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
        config->SetCFL(iMesh, 0.999*CFLMax*MGFactor[iMesh]);
      break;
    }

    config->SetCFL(iMesh, CFL);
  }

  switch( config->GetKind_Solver()) {
  case EULER : case NAVIER_STOKES : case RANS:
  case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS:      
    nVar = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetnVar();
    if (energy) rhoResOld = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
    else if (weakly_coupled_heat) rhoResOld = solver_container[INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    else rhoResOld = solver_container[INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
    break;
  case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
    rhoResOld = solver_container[INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
    break;
  case HEAT_EQUATION_FVM:
    rhoResOld = solver_container[INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    break;
  }
  
}

void COutput::Load_Data(CGeometry *geometry, CConfig *config, CSolver** solver_container){
  
  /*---- Construct a data sorter object to partition and distribute
   *  the local data into linear chunks across the processors ---*/
  
  
  if (femOutput){
    
    if (volumeDataSorter == nullptr)
      volumeDataSorter = new CFEMDataSorter(config, geometry, nVolumeFields);
    
    if (surfaceDataSorter == nullptr)
      surfaceDataSorter = new CSurfaceFEMDataSorter(config, geometry, nVolumeFields, 
                                                  dynamic_cast<CFEMDataSorter*>(volumeDataSorter));
    
  }  else {
   
    if (volumeDataSorter == nullptr)    
      volumeDataSorter = new CFVMDataSorter(config, geometry, nVolumeFields);
    
    if (surfaceDataSorter == nullptr)
      surfaceDataSorter = new CSurfaceFVMDataSorter(config, geometry, nVolumeFields,
                                                  dynamic_cast<CFVMDataSorter*>(volumeDataSorter));  
    
  }

  CollectVolumeData(config, geometry, solver_container);

  volumeDataSorter->SortOutputData();
  
}

void COutput::WriteToFile(CConfig *config, CGeometry *geometry, unsigned short format, string fileName){

  CFileWriter *fileWriter = NULL;
  
  unsigned short lastindex = fileName.find_last_of(".");
  fileName = fileName.substr(0, lastindex);
  
  /*--- Write files depending on the format --- */
  
  switch (format) {
    
    case SURFACE_CSV:
      
      if (fileName == "")
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);
       
      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();
      
      if (rank == MASTER_NODE) {
        (*fileWritingTable) << "CSV file" << fileName + CSU2FileWriter::fileExt;     
      }
      
      fileWriter = new CSU2FileWriter(volumeFieldNames, nDim, fileName, surfaceDataSorter);
      
      break;
      
    case RESTART_ASCII: case CSV:
     
      if (fileName == "")
        fileName = config->GetFilename(restartFilename, "", curTimeIter);
      
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 ASCII restart" << fileName + CSU2FileWriter::fileExt;     
      }
      
      fileWriter = new CSU2FileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter);
      
      break;
      
    case RESTART_BINARY:
      
      if (fileName == "")
        fileName = config->GetFilename(restartFilename, "", curTimeIter);
      
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 restart" << fileName + CSU2BinaryFileWriter::fileExt;   
      }
      
      fileWriter = new CSU2BinaryFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter);
      
      break;
      
    case MESH:
      
      if (fileName == "")
        fileName = volumeFilename;
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      volumeDataSorter->SortConnectivity(config, geometry, true);
      
      /*--- Set the mesh ASCII format ---*/
      
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 mesh" << fileName + CSU2MeshFileWriter::fileExt;
      }
      
      fileWriter = new CSU2MeshFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter,
                                          config->GetiZone(), config->GetnZone());
      
      
      break;    
      
    case TECPLOT_BINARY:
      
      if (fileName == "")
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      volumeDataSorter->SortConnectivity(config, geometry, true);
      
      /*--- Write tecplot binary ---*/
      
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot" << fileName + CTecplotBinaryFileWriter::fileExt;
      }
      
      fileWriter = new CTecplotBinaryFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter,
                                                curTimeIter, GetHistoryFieldValue("TIME_STEP"));
      
      break;
      
    case TECPLOT:
      
      if (fileName == "")
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      volumeDataSorter->SortConnectivity(config, geometry, true);
      
      /*--- Write tecplot binary ---*/
      
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot ASCII" << fileName + CTecplotFileWriter::fileExt;
      }
      
      fileWriter = new CTecplotFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter,
                                          curTimeIter, GetHistoryFieldValue("TIME_STEP"));
      
      break;
      
    case PARAVIEW_BINARY:
      
      if (fileName == "")
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      volumeDataSorter->SortConnectivity(config, geometry, true);
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview" << fileName + CParaviewBinaryFileWriter::fileExt;
      }
      
      fileWriter = new CParaviewBinaryFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter);
      
      break;
      
    case PARAVIEW:
      
      if (fileName == "")
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      volumeDataSorter->SortConnectivity(config, geometry, true);
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview ASCII" << fileName + CParaviewFileWriter::fileExt;      
      }
      
      fileWriter = new CParaviewFileWriter(volumeFieldNames, nDim, fileName, volumeDataSorter);
      
      break;
      
    case SURFACE_PARAVIEW:
      
      if (fileName == "")
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      surfaceDataSorter->SortConnectivity(config, geometry);      
      surfaceDataSorter->SortOutputData();
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview ASCII surface" << fileName + CParaviewFileWriter::fileExt;      
      }
      
      fileWriter = new CParaviewFileWriter(volumeFieldNames, nDim, fileName, surfaceDataSorter);
      
      break;
      
    case SURFACE_PARAVIEW_BINARY:
      
      if (fileName == "")
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      surfaceDataSorter->SortConnectivity(config, geometry);            
      surfaceDataSorter->SortOutputData();
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview surface" << fileName + CParaviewBinaryFileWriter::fileExt;      
      }
      
      fileWriter = new CParaviewBinaryFileWriter(volumeFieldNames, nDim, fileName, surfaceDataSorter);
      
      break;
      
    case SURFACE_TECPLOT:
      
      if (fileName == "")
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      surfaceDataSorter->SortConnectivity(config, geometry);      
      surfaceDataSorter->SortOutputData();
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot ASCII surface" << fileName + CTecplotFileWriter::fileExt;      
      }
      
      fileWriter = new CTecplotFileWriter(volumeFieldNames, nDim, fileName, surfaceDataSorter,
                                          curTimeIter, GetHistoryFieldValue("TIME_STEP"));
      
      break;
      
    case SURFACE_TECPLOT_BINARY:
      
      if (fileName == "")
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);
      
      /*--- Load and sort the output data and connectivity. ---*/
      
      surfaceDataSorter->SortConnectivity(config, geometry);      
      surfaceDataSorter->SortOutputData();
      
      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot surface" << fileName + CTecplotBinaryFileWriter::fileExt;      
      }
      
      fileWriter = new CTecplotBinaryFileWriter(volumeFieldNames, nDim, fileName, surfaceDataSorter,
                                                curTimeIter, GetHistoryFieldValue("TIME_STEP"));
      
      break;
      
    default:
      fileWriter = NULL;
      break;
  } 
  
  /*--- Write data to file ---*/
  
  if (fileWriter != NULL){
    
    fileWriter->Write_Data();
    
    delete fileWriter;
   
  }
}



bool COutput::SetResult_Files(CGeometry *geometry, CConfig *config, CSolver** solver_container, unsigned long Iter, bool force_writing){
  
  /*--- Load the volume data from the solver and sort it ---*/
  
  Load_Data(geometry, config, solver_container);
  
  if (WriteVolume_Output(config, Iter) || force_writing){
    
    if (rank == MASTER_NODE){
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);    
      fileWritingTable->PrintHeader();
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);    
    }
    unsigned short nVolumeFiles = config->GetnVolumeOutputFiles();
    unsigned short *VolumeFiles = config->GetVolumeOutputFiles();
        
    /*--- Loop through all requested output files ---*/
    
    for (unsigned short iFile = 0; iFile < nVolumeFiles; iFile++){
      
      WriteToFile(config, geometry, VolumeFiles[iFile]);

    }
    
    if (rank == MASTER_NODE){
      fileWritingTable->PrintFooter();
      headerNeeded = true;
    }
    
    /*--- Write any additonal files ----*/
    
    WriteAdditionalFiles(config, geometry, solver_container);
    
    return true;
  }
  
  return false;
}

bool COutput::Convergence_Monitoring(CConfig *config, unsigned long Iteration) {

  unsigned short iCounter;
    
  convergence = false;
  
  if( historyOutput_Map.count(convField) > 0 ){
    
    su2double monitor = historyOutput_Map[convField].value;
    
    /*--- Cauchy based convergence criteria ---*/
    
    if (historyOutput_Map[convField].fieldType == TYPE_COEFFICIENT) {
      
      if (Iteration == 0){
        for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++){
          cauchySerie[iCounter] = 0.0;
        }
        newFunc = monitor;
      }
      
      oldFunc = newFunc;
      newFunc = monitor;
      cauchyFunc = fabs(newFunc - oldFunc);
      
      cauchySerie[Iteration % config->GetCauchy_Elems()] = cauchyFunc;
      cauchyValue = 1.0;
      if (Iteration >= config->GetCauchy_Elems()){     
        cauchyValue = 0.0;
        for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
          cauchyValue += cauchySerie[iCounter];
      }
      
      if (cauchyValue >= config->GetCauchy_Eps()) { convergence = false;}
      else { convergence = true; newFunc = 0.0;}
      
      SetHistoryOutputValue("CAUCHY", cauchyValue);
      
    }
    
    /*--- Residual based convergence criteria ---*/
    
    if (historyOutput_Map[convField].fieldType == TYPE_RESIDUAL || historyOutput_Map[convField].fieldType == TYPE_AUTO_RESIDUAL) {
      
      /*--- Check the convergence ---*/
      
      if (Iteration != 0 && (monitor <= config->GetMinLogResidual())) { convergence = true;  }
      else { convergence = false; }
      
    }
    
    /*--- Do not apply any convergence criteria of the number
     of iterations is less than a particular value ---*/
    
    if (Iteration < config->GetStartConv_Iter()) {
      convergence = false;
    }
    
    /*--- Apply the same convergence criteria to all the processors ---*/
    
#ifdef HAVE_MPI
    
    unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
    sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
    rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
    
    /*--- Convergence criteria ---*/
    
    sbuf_conv[0] = convergence;
    SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
    
    /*-- Compute global convergence criteria in the master node --*/
    
    sbuf_conv[0] = 0;
    if (rank == MASTER_NODE) {
      if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
      else sbuf_conv[0] = 0;
    }
    
    SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (sbuf_conv[0] == 1) { convergence = true; }
    else { convergence = false;  }
    
    delete [] sbuf_conv;
    delete [] rbuf_conv;
    
#endif
    
    /*--- Stop the simulation in case a nan appears, do not save the solution ---*/
    
    if (monitor != monitor) {
      SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);
    }
    
  }
  return convergence;
}

bool COutput::MonitorTimeConvergence(CConfig *config, unsigned long TimeIteration) {
  TimeConvergence = false;

  bool Inner_IterConv = GetConvergence() || config->GetnInner_Iter()-1 <= curInnerIter; //Check, if Inner_Iter is converged

  if(Inner_IterConv && TimeIteration >= config->GetStartWindowIteration()){

     //wnd_values.push_back(GetHistoryFieldValue("BUMP_WND_AVG_"+ config->GetWndConv_Field()));

     unsigned short iCounter;

     string WndConv_Field = config->GetWndConv_Field();

     if(GetHistoryFields().count(WndConv_Field) > 0 ){

       su2double monitor = GetHistoryFields()["BUMP_WND_AVG_" + WndConv_Field].value; // Use Bump window

       /*--- Cauchy based convergence criteria ---*/

       if (GetHistoryFields()[WndConv_Field].fieldType == GetHistoryFieldType(2)/*TYPE_COEFFICIENT*/) {

         if (TimeIteration == config->GetStartWindowIteration()){
           WndNew_Func = monitor;
         }

         WndOld_Func = WndNew_Func;
         WndNew_Func = monitor;
         WndCauchy_Func = fabs(WndNew_Func - WndOld_Func);

         WndCauchy_Serie[TimeIteration % config->GetWnd_Cauchy_Elems()] = WndCauchy_Func;
         WndCauchy_Value = 1.0;
         if (TimeIteration >= config->GetWnd_Cauchy_Elems()+config->GetStartWindowIteration()){
           WndCauchy_Value = 0.0;
           for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
             WndCauchy_Value += WndCauchy_Serie[iCounter];
         }

         if (WndCauchy_Value >= config->GetWnd_Cauchy_Eps()) { TimeConvergence = false;}
         else { TimeConvergence = true;}

         SetHistoryOutputValue("TIME_WND_CAUCHY", WndCauchy_Value);
       }

       /*--- Residual based convergence criteria ---*/

       if (GetHistoryFields()[WndConv_Field].fieldType == GetHistoryFieldType(0) /*TYPE_RESIDUAL*/ || GetHistoryFields()[WndConv_Field].fieldType == GetHistoryFieldType(1)/*TYPE_AUTORESIDUAL*/) {

         /*--- Check the convergence ---*/
         if ((monitor <= config->GetWndMinLogResidual())) { TimeConvergence = true; }
         else { TimeConvergence = false; }

       }

       /*--- Do not apply any convergence criterion if the number
        of iterations is less than a particular value ---*/

       if (TimeIteration <= config->GetStartWindowIteration() + config->GetWnd_StartConv_Iter())TimeConvergence = false;

       /*--- Do not apply any convergence criterion if the option is disabled. */

       if(!config->GetWnd_Cauchy_Crit()) TimeConvergence = false;

       /*--- Stop the simulation in case a nan appears, do not save the solution ---*/

       if (monitor != monitor) {
         SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);
       }

     }
  }
  else if(TimeIteration == 0){
      SetHistoryOutputValue("TIME_WND_CAUCHY", 1.0);
    }
  return TimeConvergence;
}

void COutput::SetHistoryFile_Header(CConfig *config) { 

  unsigned short iField_Output = 0, 
      iReqField = 0,
      iMarker = 0;
  stringstream out;
  string RequestedField;
  int width = 20;
  
  for (iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    const HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField_Output]];
    for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = requestedHistoryFields[iReqField];   
      if (RequestedField == Field.outputGroup || (RequestedField == historyOutput_List[iField_Output])){
        if (Field.screenFormat == FORMAT_INTEGER) width = std::max((int)Field.fieldName.size()+2, 10);  
        else{ width = std::max((int)Field.fieldName.size()+2, 20);}
        historyFileTable->AddColumn("\"" + Field.fieldName + "\"", width);
      }
    }  
  }
  
  for (iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    for (iMarker = 0; iMarker < historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]].size(); iMarker++){  
      const HistoryOutputField &Field = historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]][iMarker];
      for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = requestedHistoryFields[iReqField];   
        if (RequestedField == Field.outputGroup || (RequestedField == historyOutputPerSurface_List[iField_Output])){
          if (Field.screenFormat == FORMAT_INTEGER) width = std::max((int)Field.fieldName.size()+2, 10);  
          else{ width = std::max((int)Field.fieldName.size()+2, 20);}
          historyFileTable->AddColumn("\"" + Field.fieldName + "\"", 20);          
        }
      }
    }
  }
  
  if (config->GetTabular_FileFormat() == TAB_TECPLOT) {
    histFile << "VARIABLES = \\" << endl;
  }
  historyFileTable->PrintHeader();
  histFile.flush();
}

void COutput::SetHistoryFile_Output(CConfig *config) { 
  
  unsigned short iField_Output = 0, 
      iReqField = 0,
      iMarker = 0;
  stringstream out;
  string RequestedField;
  
  for (iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    const HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField_Output]];
    for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = requestedHistoryFields[iReqField];   
      if (RequestedField == Field.outputGroup){
        (*historyFileTable) << Field.value; 
      }
    }
  }
  
  for (iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    for (iMarker = 0; iMarker < historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]].size(); iMarker++){
      const HistoryOutputField &Field = historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]][iMarker];
      for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = requestedHistoryFields[iReqField];   
        if (RequestedField == Field.outputGroup){
          (*historyFileTable) << Field.value;  
        }
      }
    }
  }
  
  /*--- Print the string to file and remove the last two characters (a separator and a space) ---*/

  histFile.flush();
}

void COutput::SetScreen_Header(CConfig *config) {
  if (config->GetMultizone_Problem()) 
    multiZoneHeaderTable->PrintHeader();
  convergenceTable->PrintHeader();
}

void COutput::SetScreen_Output(CConfig *config) {
  
  string RequestedField;
  
  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    stringstream out;
    RequestedField = requestedScreenFields[iReqField]; 
    if (historyOutput_Map.count(RequestedField) > 0){  
      switch (historyOutput_Map[RequestedField].screenFormat) {
        case FORMAT_INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(historyOutput_Map[RequestedField].value), fieldWidth);
          break;
        case FORMAT_FIXED:
          PrintingToolbox::PrintScreenFixed(out, historyOutput_Map[RequestedField].value, fieldWidth);
          break;
        case FORMAT_SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, historyOutput_Map[RequestedField].value, fieldWidth);
          break;      
      }
    }
    if (historyOutputPerSurface_Map.count(RequestedField) > 0){
      switch (historyOutputPerSurface_Map[RequestedField][0].screenFormat) {
        case FORMAT_INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(historyOutputPerSurface_Map[RequestedField][0].value), fieldWidth);
          break;
        case FORMAT_FIXED:
          PrintingToolbox::PrintScreenFixed(out, historyOutputPerSurface_Map[RequestedField][0].value, fieldWidth);
          break;
        case FORMAT_SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, historyOutputPerSurface_Map[RequestedField][0].value, fieldWidth);
          break;   
      }
    }      
    (*convergenceTable) << out.str();
  }
}

void COutput::PreprocessHistoryOutput(CConfig *config, bool wrt){
   
    noWriting = !wrt;

    /*--- Set the common output fields ---*/
    
    SetCommonHistoryFields(config);
    
    /*--- Set the History output fields using a virtual function call to the child implementation ---*/
    
    SetHistoryOutputFields(config);
    
    /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/
    
    Postprocess_HistoryFields(config);
    
    /*--- We use a fixed size of the file output summary table ---*/
    
    int total_width = 72;    
    fileWritingTable->AddColumn("File Writing Summary", (total_width-1)/2); 
    fileWritingTable->AddColumn("Filename", total_width/2);  
    fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    
    if (rank == MASTER_NODE && !noWriting){
      
      /*--- Check for consistency and remove fields that are requested but not available --- */
      
      CheckHistoryOutput();
      
      /*--- Open history file and print the header ---*/
      
      PrepareHistoryFile(config);
      
      total_width = nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1);
      
      /*--- Set the multizone screen header ---*/
      
      if (config->GetMultizone_Problem()){
        multiZoneHeaderTable->AddColumn(multiZoneHeaderString, total_width);      
        multiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
        multiZoneHeaderTable->SetPrintHeaderBottomLine(false);
      }

    }
    

    
}

void COutput::PreprocessMultizoneHistoryOutput(COutput **output, CConfig **config, bool wrt){
    
  noWriting = !wrt;

  /*--- Set the History output fields using a virtual function call to the child implementation ---*/
  
  SetMultizoneHistoryOutputFields(output, config);
  
  /*--- We use a fixed size of the file output summary table ---*/
  
  int total_width = 72;    
  fileWritingTable->AddColumn("File Writing Summary", (total_width-1)/2); 
  fileWritingTable->AddColumn("Filename", total_width/2);  
  fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);

  if (rank == MASTER_NODE && !noWriting){
    
    /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/
   
    //Postprocess_HistoryFields(config[ZONE_0]);
    
    /*--- Check for consistency and remove fields that are requested but not available --- */
    
    CheckHistoryOutput();
    
    /*--- Open history file and print the header ---*/
    
    PrepareHistoryFile(config[ZONE_0]);
    
    total_width = nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1);    
    
    /*--- Set the multizone screen header ---*/

    if (config[ZONE_0]->GetMultizone_Problem()){
      multiZoneHeaderTable->AddColumn(multiZoneHeaderString, nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1));      
      multiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
      multiZoneHeaderTable->SetPrintHeaderBottomLine(false);
    }
    
  }
  
}

void COutput::PrepareHistoryFile(CConfig *config){
  
  /*--- Open the history file ---*/
  
  histFile.open(historyFilename.c_str(), ios::out);
  
  /*--- Create and format the history file table ---*/
  
  historyFileTable = new PrintingToolbox::CTablePrinter(&histFile, "");
  historyFileTable->SetInnerSeparator(historySep);
  historyFileTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
  historyFileTable->SetPrintHeaderTopLine(false);
  historyFileTable->SetPrintHeaderBottomLine(false);
  historyFileTable->SetPrecision(14);
  
  /*--- Add the header to the history file. ---*/
  
  SetHistoryFile_Header(config);    
  
}

void COutput::CheckHistoryOutput(){
  
  
  /*--- Set screen convergence output header and remove unavailable fields ---*/
  
  string RequestedField;
  vector<string> FieldsToRemove;
  vector<bool> FoundField(nRequestedHistoryFields, false);
  
  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    RequestedField = requestedScreenFields[iReqField];  
    if (historyOutput_Map.count(RequestedField) > 0){ 
      convergenceTable->AddColumn(historyOutput_Map[RequestedField].fieldName, fieldWidth);
    }
    else if (historyOutputPerSurface_Map.count(RequestedField) > 0){
      convergenceTable->AddColumn(historyOutputPerSurface_Map[RequestedField][0].fieldName, fieldWidth);
    }else {
      FieldsToRemove.push_back(RequestedField);
    }
  }
  
  /*--- Remove fields which are not defined --- */
  
  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following screen output fields:" << endl;
        cout << "  ";
      }        cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedScreenFields.erase(std::find(requestedScreenFields.begin(), requestedScreenFields.end(), FieldsToRemove[iReqField]));
  }
  
  nRequestedScreenFields = requestedScreenFields.size();
  
  if (rank == MASTER_NODE){
    cout <<"Screen output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
      RequestedField = requestedScreenFields[iReqField];            
      cout << requestedScreenFields[iReqField];
      if (iReqField != nRequestedScreenFields - 1) cout << ", ";
    }
    cout << endl;
  }
  
  /*--- Remove unavailable fields from the history file output ---*/
  
  FieldsToRemove.clear();
  FoundField = vector<bool>(nRequestedHistoryFields, false);
  
  for (unsigned short iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = requestedHistoryFields[iReqField];   
      if (RequestedField == Field.outputGroup){
        FoundField[iReqField] = true;
      }
    }
  }
  
  for (unsigned short iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]].size(); iMarker++){
      HistoryOutputField &Field = historyOutputPerSurface_Map[historyOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = requestedHistoryFields[iReqField];   
        if (RequestedField == Field.outputGroup){
          FoundField[iReqField] = true;
        }
      }
    }
  }
  
  for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(requestedHistoryFields[iReqField]);
    }
  }
  
  /*--- Remove fields which are not defined --- */    
  
  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following history output groups:" << endl;
        cout << "  ";
      }        cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedHistoryFields.erase(std::find(requestedHistoryFields.begin(), requestedHistoryFields.end(), FieldsToRemove[iReqField]));
  }
  
  nRequestedHistoryFields = requestedHistoryFields.size();
  
  if (rank == MASTER_NODE){
    cout <<"History output groups: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = requestedHistoryFields[iReqField];            
      cout << requestedHistoryFields[iReqField];
      if (iReqField != nRequestedHistoryFields - 1) cout << ", ";
    }
    cout << endl;
  }
  
  /*--- Check that the requested convergence monitoring field is available ---*/

  if (historyOutput_Map.count(convField) == 0){
    SU2_MPI::Error(string("Convergence monitoring field ") + convField + string(" not available"), CURRENT_FUNCTION);
  }
}

void COutput::PreprocessVolumeOutput(CConfig *config){

  /*--- Set the volume output fields using a virtual function call to the child implementation ---*/  
  
  SetVolumeOutputFields(config);
   
  nVolumeFields = 0;
  
  string RequestedField;
  std::vector<bool> FoundField(nRequestedVolumeFields, false);
  vector<string> FieldsToRemove;
  
  
  /*--- Loop through all fields defined in the corresponding SetVolumeOutputFields(). 
 * If it is also defined in the config (either as part of a group or a single field), the field 
 * object gets an offset so that we know where to find the data in the Local_Data() array.
 *  Note that the default offset is -1. An index !=-1 defines this field as part of the output. ---*/

  for (unsigned short iField_Output = 0; iField_Output < volumeOutput_List.size(); iField_Output++){
    
    VolumeOutputField &Field = volumeOutput_Map[volumeOutput_List[iField_Output]];
    
    /*--- Loop through all fields specified in the config ---*/
    
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      
      RequestedField = requestedVolumeFields[iReqField];  
            
      if (((RequestedField == Field.outputGroup) || (RequestedField == volumeOutput_List[iField_Output])) && (Field.offset == -1)){
        Field.offset = nVolumeFields;
        volumeFieldNames.push_back(Field.fieldName);
        nVolumeFields++;
        
        FoundField[iReqField] = true;
      }
    }    
  }
  
  for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(requestedVolumeFields[iReqField]);
    }
  }
  
  /*--- Remove fields which are not defined --- */    
  
  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following volume output fields/groups:" << endl;
        cout << "  ";
      }
      cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedVolumeFields.erase(std::find(requestedVolumeFields.begin(), requestedVolumeFields.end(), FieldsToRemove[iReqField]));
  }
  
  if (rank == MASTER_NODE){
    cout <<"Volume output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      RequestedField = requestedVolumeFields[iReqField];            
      cout << requestedVolumeFields[iReqField];
      if (iReqField != nRequestedVolumeFields - 1) cout << ", ";
    }
    cout << endl;
  }
}

void COutput::CollectVolumeData(CConfig* config, CGeometry* geometry, CSolver** solver){
  
  unsigned short iMarker = 0;
  unsigned long iPoint = 0, jPoint = 0;
  unsigned long iVertex = 0;
  
  /*--- Reset the offset cache and index --- */
  curFieldIndex = 0;
  fieldIndexCache.clear();
  curGetFieldIndex = 0;
  fieldGetIndexCache.clear();
  
  if (femOutput){
    
    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/
  
    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
    
    /*--- Access the solution by looping over the owned volume elements. ---*/
  
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        
        buildFieldIndexCache = !fieldIndexCache.size() ? true : false;
        
        LoadVolumeDataFEM(config, geometry, solver, l, jPoint, j);
        
        jPoint++;
        
      }
    }
    
  } else {
    
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      /*--- Load the volume data into the Local_Data() array. --- */
      
      buildFieldIndexCache = !fieldIndexCache.size() ? true : false;

      LoadVolumeData(config, geometry, solver, iPoint);

    }
    
    /*--- Reset the offset cache and index --- */
    curFieldIndex = 0;
    fieldIndexCache.clear();
    curGetFieldIndex = 0;
    fieldGetIndexCache.clear();
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  
        if(geometry->node[iPoint]->GetDomain()){
          
          buildFieldIndexCache = !fieldIndexCache.size() ? true : false;
   
          LoadSurfaceData(config, geometry, solver, iPoint, iMarker, iVertex);
          
        }
      }   
    }
  }
}

void COutput::SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){
  
  if (buildFieldIndexCache){ 
    
    /*--- Build up the offset cache to speed up subsequent 
     * calls of this routine since the order of calls is 
     * the same for every value of iPoint --- */
    
    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map[name].offset;
      fieldIndexCache.push_back(Offset);        
      if (Offset != -1){
        volumeDataSorter->SetUnsorted_Data(iPoint, Offset, value);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
    }
  } else {
    
    /*--- Use the offset cache for the access ---*/
    
    const short Offset = fieldIndexCache[curFieldIndex++];
    if (Offset != -1){
      volumeDataSorter->SetUnsorted_Data(iPoint, Offset, value);
    }   
    if (curFieldIndex == fieldIndexCache.size()){
      curFieldIndex = 0;
    }
  }
  
}

su2double COutput::GetVolumeOutputValue(string name, unsigned long iPoint){
  
  if (buildFieldIndexCache){ 
    
    /*--- Build up the offset cache to speed up subsequent 
     * calls of this routine since the order of calls is 
     * the same for every value of iPoint --- */
    
    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map[name].offset;
      fieldGetIndexCache.push_back(Offset);        
      if (Offset != -1){
        return volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
    }
  } else {
    
    /*--- Use the offset cache for the access ---*/
    
    const short Offset = fieldGetIndexCache[curGetFieldIndex++];
  
    if (curGetFieldIndex == fieldGetIndexCache.size()){
      curGetFieldIndex = 0;
    }
    if (Offset != -1){
      return volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
    } 
  }
  
  return 0.0;
}

void COutput::SetAvgVolumeOutputValue(string name, unsigned long iPoint, su2double value){
  
  const su2double scaling = 1.0 / su2double(curAbsTimeIter + 1);
  
  if (buildFieldIndexCache){ 
    
    /*--- Build up the offset cache to speed up subsequent 
     * calls of this routine since the order of calls is 
     * the same for every value of iPoint --- */
    
    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map[name].offset;
      fieldIndexCache.push_back(Offset);        
      if (Offset != -1){
        
        const su2double old_value = volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
        const su2double new_value = value * scaling + old_value *( 1.0 - scaling);
        
        volumeDataSorter->SetUnsorted_Data(iPoint, Offset, new_value);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
    }
  } else {
    
    /*--- Use the offset cache for the access ---*/
    
    const short Offset = fieldIndexCache[curFieldIndex++];
    if (Offset != -1){
      
      const su2double old_value = volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
      const su2double new_value = value * scaling + old_value *( 1.0 - scaling);
      
      volumeDataSorter->SetUnsorted_Data(iPoint, Offset, new_value);
    }   
    if (curFieldIndex == fieldIndexCache.size()){
      curFieldIndex = 0;
    }
  }
  
}





void COutput::Postprocess_HistoryData(CConfig *config){
   
  map<string, su2double> Average;
  map<string, int> Count;
  
  for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    HistoryOutputField &currentField = historyOutput_Map[historyOutput_List[iField]];
    if (currentField.fieldType == TYPE_RESIDUAL){
      if ( SetInit_Residuals(config) || (currentField.value > initialResiduals[historyOutput_List[iField]]) ) {
        initialResiduals[historyOutput_List[iField]] = currentField.value;
      }
      SetHistoryOutputValue("REL_" + historyOutput_List[iField], 
                            currentField.value - initialResiduals[historyOutput_List[iField]]);
      
      Average[currentField.outputGroup] += currentField.value;
      Count[currentField.outputGroup]++;
           
    }
    if (currentField.fieldType == TYPE_COEFFICIENT){
      if (config->GetTime_Domain()){
      if(SetUpdate_Averages(config)){
        SetHistoryOutputValue("TAVG_" + historyOutput_List[iField], runningAverages[historyOutput_List[iField]].Update(currentField.value));
        runningAverages[historyOutput_List[iField]].addValue(currentField.value,config->GetTimeIter(), config->GetStartWindowIteration()); //Setting Start-Iteration for Windowing
        SetHistoryOutputValue("SQ_WND_AVG_" + historyOutput_List[iField],     runningAverages[historyOutput_List[iField]].WindowedUpdate( 0));
        SetHistoryOutputValue("HANN_WND_AVG_" + historyOutput_List[iField],   runningAverages[historyOutput_List[iField]].WindowedUpdate( 1));
        SetHistoryOutputValue("HANNSQ_WND_AVG_" + historyOutput_List[iField], runningAverages[historyOutput_List[iField]].WindowedUpdate( 2));
        SetHistoryOutputValue("BUMP_WND_AVG_" + historyOutput_List[iField],   runningAverages[historyOutput_List[iField]].WindowedUpdate( 3));
      }
      if (config->GetDirectDiff() != NO_DERIVATIVE){

        SetHistoryOutputValue("D_TAVG_" + historyOutput_List[iField]            , SU2_TYPE::GetDerivative(runningAverages[historyOutput_List[iField]].GetVal()));
        SetHistoryOutputValue("D_SQ_WND_AVG_" + historyOutput_List[iField]      , SU2_TYPE::GetDerivative(runningAverages[historyOutput_List[iField]].GetWndVal(0)));
        SetHistoryOutputValue("D_HANN_WND_AVG_" + historyOutput_List[iField]    , SU2_TYPE::GetDerivative(runningAverages[historyOutput_List[iField]].GetWndVal(1)));
        SetHistoryOutputValue("D_HANNSQ_WND_AVG_" + historyOutput_List[iField]  , SU2_TYPE::GetDerivative(runningAverages[historyOutput_List[iField]].GetWndVal(2)));
        SetHistoryOutputValue("D_BUMP_WND_AVG_" + historyOutput_List[iField]    , SU2_TYPE::GetDerivative(runningAverages[historyOutput_List[iField]].GetWndVal(3)));
      }
    }
    if (config->GetDirectDiff() != NO_DERIVATIVE){
        SetHistoryOutputValue("D_" + historyOutput_List[iField], SU2_TYPE::GetDerivative(currentField.value));
    }
}
  }
  
  map<string, su2double>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    SetHistoryOutputValue("AVG_" + it->first, it->second/Count[it->first]);
  }
  
}

void COutput::Postprocess_HistoryFields(CConfig *config){
  
  map<string, bool> Average;
  map<string, string> AverageGroupName =  CCreateMap<string, string>("BGS_RES", "bgs")("RMS_RES","rms")("MAX_RES", "max");
  
  for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    HistoryOutputField &currentField = historyOutput_Map[historyOutput_List[iField]];
    if (currentField.fieldType == TYPE_RESIDUAL){
      AddHistoryOutput("REL_" + historyOutput_List[iField], "rel" + currentField.fieldName, currentField.screenFormat,
                       "REL_" + currentField.outputGroup,  "Relative residual.", TYPE_AUTO_RESIDUAL);
      Average[currentField.outputGroup] = true;
    }
  }
  
  map<string, bool>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    AddHistoryOutput("AVG_" + it->first, "avg[" + AverageGroupName[it->first] + "]", FORMAT_FIXED,
                     "AVG_" + it->first , "Average residual over all solution variables.", TYPE_AUTO_RESIDUAL);
  }  
  if (config->GetTime_Domain()){
   for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    HistoryOutputField &currentField = historyOutput_Map[historyOutput_List[iField]];
    if (currentField.fieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("TAVG_"   + historyOutput_List[iField], "tavg["  + currentField.fieldName + "]", currentField.screenFormat, "TAVG_"   + currentField.outputGroup, "Time averaged values.", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("SQ_WND_AVG_" + historyOutput_List[iField]       , "sq_wnd_avg[" + currentField.fieldName + "]"      , currentField.screenFormat, "SQ_WND_AVG_" + currentField.outputGroup, "Time averaged square window weighted values.", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("HANN_WND_AVG_" + historyOutput_List[iField]     , "hann_wnd_avg[" + currentField.fieldName + "]"    , currentField.screenFormat, "HANN_WND_AVG_" + currentField.outputGroup, "Time averaged hann window weighted values.", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("HANNSQ_WND_AVG_" + historyOutput_List[iField]   , "hannSq_wnd_avg[" + currentField.fieldName + "]"  , currentField.screenFormat, "HANNSQ_WND_AVG_" + currentField.outputGroup, "Time averaged hann-square window weighted values.", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("BUMP_WND_AVG_" + historyOutput_List[iField]     , "bump_wnd_avg[" + currentField.fieldName + "]"    , currentField.screenFormat, "BUMP_WND_AVG_" + currentField.outputGroup, "Time averaged bump window weighted values.", TYPE_AUTO_COEFFICIENT);
    }
  }}
  
  if (config->GetDirectDiff()){
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      HistoryOutputField &currentField = historyOutput_Map[historyOutput_List[iField]];
      if (currentField.fieldType == TYPE_COEFFICIENT){
        AddHistoryOutput("D_"      + historyOutput_List[iField], "d["     + currentField.fieldName + "]",
                         currentField.screenFormat, "D_"      + currentField.outputGroup, 
                         "Derivative value (DIRECT_DIFF=YES)", TYPE_AUTO_COEFFICIENT);  
      }
    }
  }
  
if (config->GetTime_Domain() && config->GetDirectDiff()){
  for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    HistoryOutputField &currentField = historyOutput_Map[historyOutput_List[iField]];
    if (currentField.fieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("D_TAVG_" + historyOutput_List[iField], "dtavg[" + currentField.fieldName + "]", currentField.screenFormat, "D_TAVG_" + currentField.outputGroup, "Derivative of the time averaged value (DIRECT_DIFF=YES).", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("D_SQ_WND_AVG_" + historyOutput_List[iField]     , "d_sq_wnd_avg[" + currentField.fieldName + "]"    , currentField.screenFormat, "D_SQ_WND_AVG_" + currentField.outputGroup, "Time averaged square window weighted values (DIRECT_DIFF=YES).", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("D_HANN_WND_AVG_" + historyOutput_List[iField]   , "d_hann_wnd_avg[" + currentField.fieldName + "]"  , currentField.screenFormat, "D_HANN_WND_AVG_" + currentField.outputGroup, "Time averaged hann window weighted values (DIRECT_DIFF=YES).", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("D_HANNSQ_WND_AVG_" + historyOutput_List[iField] , "d_hannsq_wnd_avg[" + currentField.fieldName + "]", currentField.screenFormat, "D_HANNSQ_WND_AVG_" + currentField.outputGroup, "Time averaged hann-square window weighted values (DIRECT_DIFF=YES).", TYPE_AUTO_COEFFICIENT);
      AddHistoryOutput("D_BUMP_WND_AVG_" + historyOutput_List[iField]   , "d_bump_wnd_avg[" + currentField.fieldName + "]"  , currentField.screenFormat, "D_BUMP_WND_AVG_" + currentField.outputGroup, "Time averaged bump window weighted values (DIRECT_DIFF=YES).", TYPE_AUTO_COEFFICIENT);
    }
  }
 }
//if (historyOutput_Map[convField].fieldType == TYPE_COEFFICIENT){
  AddHistoryOutput("CAUCHY", "C["  + historyOutput_Map[convField].fieldName + "]", FORMAT_SCIENTIFIC, "CAUCHY","Cauchy residual value of field set with CONV_FIELD." ,TYPE_AUTO_COEFFICIENT);
  AddHistoryOutput("TIME_WND_CAUCHY", "t_wnd_C[" + historyOutput_Map[config->GetWndConv_Field()].fieldName+ "]", FORMAT_SCIENTIFIC, "TIME_WND_CAUCHY", "Cauchy residual value of field set with WND_CONV_FIELD.", TYPE_AUTO_COEFFICIENT);
//}

}

bool COutput::WriteScreen_Header(CConfig *config) {  
  
  unsigned long RestartIter = 0;
  
  if (config->GetRestart() && config->GetTime_Domain()){
    RestartIter = config->GetRestart_Iter();
  }
  
  unsigned long ScreenWrt_Freq_Inner = config->GetScreen_Wrt_Freq(2);
  unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
  unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);
  
  /*--- Always print header if it is forced ---*/
  
  if (headerNeeded){
    headerNeeded = false;
    return true;
  }
  
  /*--- Header is always disabled for multizone problems unless explicitely requested --- */
  
  if (config->GetMultizone_Problem() && !config->GetWrt_ZoneConv()){
    return false;
  }

  /* --- Always print header in the first iteration --- */
  
  if ((curInnerIter == 0) && 
      (curOuterIter == 0) && 
      (curTimeIter == RestartIter)){
    return true;
  }
  
  if (!PrintOutput(curTimeIter, ScreenWrt_Freq_Time)&& 
      !(curTimeIter == config->GetnTime_Iter() - 1)){
    return false;
  }
   
  /*--- If there is no inner or outer iteration, don't print header ---*/
  if (ScreenWrt_Freq_Outer == 0 && ScreenWrt_Freq_Inner == 0){
    return false;
  }
  
  /*--- Print header if we are at the first inner iteration ---*/
  
  if (curInnerIter == 0){
    return true;
  }
  
  return false;
}

bool COutput::WriteScreen_Output(CConfig *config) {
  
  unsigned long ScreenWrt_Freq_Inner = config->GetScreen_Wrt_Freq(2);
  unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
  unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);    
  
  if (config->GetMultizone_Problem() && !config->GetWrt_ZoneConv()){
    
    return false;
    
  }
  
  /*--- Check if screen output should be written --- */
  
  if (!PrintOutput(curTimeIter, ScreenWrt_Freq_Time)&& 
      !(curTimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (convergence) {return true;}
  
  if (!PrintOutput(curOuterIter, ScreenWrt_Freq_Outer) && 
      !(curOuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
  
  if (!PrintOutput(curInnerIter, ScreenWrt_Freq_Inner) &&
      !(curInnerIter == config->GetnInner_Iter() - 1)){
    
    return false;
    
  }
 
  return true;
  
}

bool COutput::WriteHistoryFile_Output(CConfig *config) { 

  unsigned long HistoryWrt_Freq_Inner = config->GetHistory_Wrt_Freq(2);
  unsigned long HistoryWrt_Freq_Outer = config->GetHistory_Wrt_Freq(1);
  unsigned long HistoryWrt_Freq_Time  = config->GetHistory_Wrt_Freq(0);    
    
  /*--- Check if screen output should be written --- */
  
  if (!PrintOutput(curTimeIter, HistoryWrt_Freq_Time)&& 
      !(curTimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (convergence) {return true;}
  
  if (!PrintOutput(curOuterIter,HistoryWrt_Freq_Outer) && 
      !(curOuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
  
  if (!PrintOutput(curInnerIter, HistoryWrt_Freq_Inner) &&
      !(curInnerIter == config->GetnInner_Iter() - 1)){
    
    return false;
    
  }
 
  return true;

}

bool COutput::WriteVolume_Output(CConfig *config, unsigned long Iter){
  return ((Iter > 0) && (Iter % config->GetVolume_Wrt_Freq() == 0));
}

void COutput::SetCommonHistoryFields(CConfig *config){
  
  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER", "Time iteration index"); 
  /// DESCRIPTION: The outer iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER", "Outer iteration index"); 
  /// DESCRIPTION: The inner iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", FORMAT_INTEGER,  "ITER", "Inner iteration index"); 
  /// END_GROUP
  
  /// BEGIN_GROUP: TIME_DOMAIN, DESCRIPTION: Time integration information
  /// Description: The current time
  AddHistoryOutput("CUR_TIME", "Cur_Time", FORMAT_SCIENTIFIC, "TIME_DOMAIN", "Current physical time (s)");
  /// Description: The current time step
  AddHistoryOutput("TIME_STEP", "Time_Step", FORMAT_SCIENTIFIC, "TIME_DOMAIN", "Current time step (s)");
 
  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("PHYS_TIME",   "Time(sec)", FORMAT_SCIENTIFIC, "PHYS_TIME", "Average wall-clock time"); 
  
}

void COutput::LoadCommonHistoryData(CConfig *config){
  
  SetHistoryOutputValue("TIME_ITER",  curTimeIter);  
  SetHistoryOutputValue("INNER_ITER", curInnerIter);
  SetHistoryOutputValue("OUTER_ITER", curOuterIter); 
  
  if (config->GetTime_Domain()){
    SetHistoryOutputValue("TIME_STEP", config->GetDelta_UnstTimeND()*config->GetTime_Ref());           
    if (curInnerIter == 0){
      SetHistoryOutputValue("CUR_TIME",  GetHistoryFieldValue("CUR_TIME") + GetHistoryFieldValue("TIME_STEP"));      
    }
  }
  
  su2double StopTime, UsedTime;
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  
  UsedTime = (StopTime - config->Get_StartTime())/((curOuterIter + 1) * (curInnerIter+1));
  
  SetHistoryOutputValue("PHYS_TIME", UsedTime);
  
}

void COutput::PrintHistoryFields(){ 
  
  if (rank == MASTER_NODE){
    
    PrintingToolbox::CTablePrinter HistoryFieldTable(&std::cout);
    
    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;
    
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      
      HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField]];
      
      if (Field.description != ""){
        if (historyOutput_List[iField].size() > NameSize){
          NameSize = historyOutput_List[iField].size();
        }
        if (Field.outputGroup.size() > GroupSize){
          GroupSize = Field.outputGroup.size();
        }
        if (Field.description.size() > DescrSize){
          DescrSize = Field.description.size();
        }
      }
    }
    
    cout << "Available output fields for the current configuration in " << multiZoneHeaderString << ":" << endl;
    
    HistoryFieldTable.AddColumn("Name", NameSize);
    HistoryFieldTable.AddColumn("Group Name", GroupSize);
    HistoryFieldTable.AddColumn("Type",5);
    HistoryFieldTable.AddColumn("Description", DescrSize);
    HistoryFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    
    HistoryFieldTable.PrintHeader();
    
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      
      HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField]];
      
      if (Field.fieldType == TYPE_DEFAULT || Field.fieldType == TYPE_COEFFICIENT || Field.fieldType == TYPE_RESIDUAL){
        string type;
        switch (Field.fieldType) {
          case TYPE_COEFFICIENT:
            type = "C";
            break;
          case TYPE_RESIDUAL:
            type = "R";
            break;
          default:
            type = "D";
            break;
        }
        
        if (Field.description != "")
          HistoryFieldTable << historyOutput_List[iField] << Field.outputGroup << type << Field.description;
        
      }
    }
    
    HistoryFieldTable.PrintFooter();
    
    cout << "Type legend: Default (D), Residual (R), Coefficient (C)" << endl;
    
    cout << "Generated output fields (only first field of every group is shown):" << endl;
    
    PrintingToolbox::CTablePrinter ModifierTable(&std::cout);
    
    ModifierTable.AddColumn("Name", NameSize);
    ModifierTable.AddColumn("Group Name", GroupSize);
    ModifierTable.AddColumn("Type",5);
    ModifierTable.AddColumn("Description", DescrSize);
    ModifierTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    ModifierTable.PrintHeader();
    
    std::map<string, bool> GroupVisited;
    
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      
      HistoryOutputField &Field = historyOutput_Map[historyOutput_List[iField]];
      
      if ((Field.fieldType == TYPE_AUTO_COEFFICIENT || Field.fieldType == TYPE_AUTO_RESIDUAL) && (GroupVisited.count(Field.outputGroup) == 0)){
        string type;
        switch (Field.fieldType) {
          case TYPE_AUTO_COEFFICIENT:
            type = "AC";
            break;
          case TYPE_AUTO_RESIDUAL:
            type = "AR";
            break;
          default:
            type = "AD";
            break;
        }
        
        if (Field.description != "")
          ModifierTable << historyOutput_List[iField] << Field.outputGroup << type << Field.description;
        
        GroupVisited[Field.outputGroup] = true;
      }
    }   
    ModifierTable.PrintFooter();

  }
}

void COutput::PrintVolumeFields(){
  
  if (rank == MASTER_NODE){
    
    PrintingToolbox::CTablePrinter VolumeFieldTable(&std::cout);
    
    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;
    
    for (unsigned short iField = 0; iField < volumeOutput_List.size(); iField++){
      
      VolumeOutputField &Field = volumeOutput_Map[volumeOutput_List[iField]];
      
      if (Field.description != ""){
        if (volumeOutput_List[iField].size() > NameSize){
          NameSize = volumeOutput_List[iField].size();
        }
        if (Field.outputGroup.size() > GroupSize){
          GroupSize = Field.outputGroup.size();
        }
        if (Field.description.size() > DescrSize){
          DescrSize = Field.description.size();
        }
      }
    }
    
    cout << "Available output fields for the current configuration in " << multiZoneHeaderString << ":" << endl;
    
    VolumeFieldTable.AddColumn("Name", NameSize);
    VolumeFieldTable.AddColumn("Group Name", GroupSize);
    VolumeFieldTable.AddColumn("Description", DescrSize);
    VolumeFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    
    VolumeFieldTable.PrintHeader();
    
    for (unsigned short iField = 0; iField < volumeOutput_List.size(); iField++){
      
      VolumeOutputField &Field = volumeOutput_Map[volumeOutput_List[iField]];

      if (Field.description != "")
        VolumeFieldTable << volumeOutput_List[iField] << Field.outputGroup << Field.description;
      
    }
    
    VolumeFieldTable.PrintFooter();
  }
}

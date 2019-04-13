/*!
 * \file output_heat.hpp
 * \brief Headers of the main subroutines for generating the file outputs.
 *        The subroutines and functions are in the <i>output_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, M. Colonno
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

#pragma once

#include "output.hpp"


/*! \class CHeatOutput
 *  \brief Output class for heat problems.
 *  \author R. Sanchez, T. Albring.
 *  \date June 5, 2018.
 */
class CHeatOutput : public COutput {
private:
  bool multizone;

  char char_histfile[200];

public:


  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CHeatOutput(CConfig *config, CGeometry *geometry, unsigned short iZone);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CHeatOutput(void);

  /*!
   * \brief Set the history file header
   * \param[in] config - Definition of the particular problem.
   */
  void LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst);


  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteHistoryFile_Output(CConfig *config, bool write_dualtime);

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreen_Header(CConfig *config);


  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreen_Output(CConfig *config, bool write_dualtime);
  
  /*!
   * \brief SetHistoryOutputFields
   * \param config
   */
  void SetHistoryOutputFields(CConfig *config);
   
  /*!
   * \brief SetVolumeOutputFields
   * \param config
   */
  void SetVolumeOutputFields(CConfig *config);
  
  /*!
   * \brief LoadVolumeData
   * \param config
   * \param geometry
   * \param solver
   * \param iPoint
   */
  void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint);
 
};
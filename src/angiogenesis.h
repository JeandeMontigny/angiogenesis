// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project and Jean de Montigny.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef ANGIOGENESIS_H_
#define ANGIOGENESIS_H_

#include "biodynamo.h"
#include "biology_modules.h"
#include "core/substance_initializers.h"
#include "neuroscience/neuroscience.h"
#include "extended_objects.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  auto set_param = [](Param* param) {
    // Create an artificial bound for the simulation space
    param->bound_space = true;
    param->min_bound = -100;
    param->max_bound = 250;
  };
  neuroscience::InitModule();
  Simulation simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();

  ModelInitializer::DefineSubstance(kVEGF, "VEGF", 0.5, 0, 50);

  auto* cell = new neuroscience::NeuronSoma({0,0,0});
  cell->SetDiameter(0.1);
  rm->AddAgent(cell);
  Vessel my_vessel;
  auto* vessel = bdm_static_cast<Vessel*>(
        cell->ExtendNewNeurite({1, 0, 0}, &my_vessel));
  vessel->SetDiameter(2);
  vessel->SetCanBranch(true);
  vessel->AddBehavior(new VascularGrowth_BM());
  // elongate vessel to set up initial condition
  for (int i = 0; i < 300; i++) {
    vessel->ElongateTerminalEnd(100, {1, i*0.01, 0});
    vessel->RunDiscretization();
  }
  std::cout << "initial condition set" << std::endl;

  // create cancerous cells
  auto* c_cell = new Cell({70, 150, 0});
  c_cell->SetDiameter(10);
  c_cell->AddBehavior(new VegfSecretion_BM());
  c_cell->AddBehavior(new TumourGrowth_BM());
  rm->AddAgent(c_cell);

  // Run simulation
  std::cout << "simulating.." << std::endl;
  simulation.GetScheduler()->Simulate(1000);
  std::cout << "done" << std::endl;

  return 0;
}

}  // namespace bdm

#endif  // ANGIOGENESIS_H_

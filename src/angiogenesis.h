// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
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

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  auto set_param = [](Param* param) {
    // Create an artificial bound for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = -100;
    param->max_bound_ = 250;
    param->run_mechanical_interactions_ = true;
  };
  neuroscience::InitModule();
  Simulation simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();

  ModelInitializer::DefineSubstance(kVEGF, "VEGF", 0.5, 0, 50);

  auto* cell = new neuroscience::NeuronSoma({0,0,0});
  cell->SetDiameter(2);
  rm->push_back(cell);
  auto vessel = cell->ExtendNewNeurite({1,0,0})->GetSoPtr<NeuriteElement>();
  vessel->SetDiameter(2);
  vessel->AddBiologyModule(new VascularGrowth_BM());
  // elongate vessel to set up initial condition
  for (int i = 0; i < 300; i++) {
    vessel->ElongateTerminalEnd(100, {1, i*0.01, 0});
    vessel->RunDiscretization();
  }

  // create cancerous cells
  auto* c_cell = new Cell({50, 150, 0});
  c_cell->SetDiameter(10);
  c_cell->AddBiologyModule(new VegfSecretion_BM());
  c_cell->AddBiologyModule(new TumourGrowth_BM());
  rm->push_back(c_cell);

  // Run simulation
  simulation.GetScheduler()->Simulate(1);

  return 0;
}

}  // namespace bdm

#endif  // ANGIOGENESIS_H_

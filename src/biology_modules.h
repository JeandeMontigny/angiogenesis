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
#ifndef BIOLOGY_MODULES_H_
#define BIOLOGY_MODULES_H_

#include "core/biology_module/biology_module.h"
#include "core/diffusion_grid.h"
#include "core/sim_object/cell.h"
#include "neuroscience/neurite_element.h"
#include "extended_objects.h"

namespace bdm {

enum Substances { kVEGF };

// -----------------------------------------------------------------------------
/// grow vascularity depending on chemical substance
struct VascularGrowth_BM : public BaseBiologyModule {
  VascularGrowth_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  VascularGrowth_BM(const Event& event, BaseBiologyModule* other,
                      uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new VascularGrowth_BM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new VascularGrowth_BM(*this);
  }

  void Run(SimObject* so) override {
    // if initial condition not set

    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    if (!init_) {
      dg_vegf_ = rm->GetDiffusionGrid(kVEGF);
      init_ = true;
    }

    auto* vessel = bdm_static_cast<Vessel*>(so);
    double concentration = dg_vegf_->GetConcentration(vessel->GetPosition());

    // remove biology module of main vessel terminal to avoid growth
    if (vessel->CanBranch() && vessel->IsTerminal()) {
      vessel->RemoveBiologyModule(this);
    }

    // if part of the main vessel and concentraiton >
    if (vessel->CanBranch() && concentration > 1e-6) {
      Double3 gradient;
      // just a single chance to branch, depending on concentration
      if (random->Uniform(0,1) < 1e5 * concentration){
        dg_vegf_->GetGradient(vessel->GetPosition(), &gradient);

        auto* branch = bdm_static_cast<Vessel*>(vessel->Branch(gradient));
        dg_vegf_->GetGradient(branch->GetPosition(), &gradient);

        branch->SetDiameter(1);
        branch->SetCanBranch(false);
        branch->ElongateTerminalEnd(50, gradient);
      }
      // if didn't branch, remove biology module, so won't do anything anymore
      vessel->RemoveBiologyModule(this);
    } // end if part of main vessel and concentration >

    // if is an extension of the main vessel
    if (!vessel->CanBranch()) {
      Double3 gradient;
      dg_vegf_->GetGradient(vessel->GetPosition(), &gradient);
      Double3 rand_dir = {random->Uniform(0,1), random->Uniform(0,1), random->Uniform(0,1)};
      // Double3 prev_dir = vessel->GetSpringAxis();
      Double3 next_dir = rand_dir + gradient;
      vessel->ElongateTerminalEnd(50, next_dir);
      // chance to Bifurcate
      // note: no branching allowed on non terminal point
      if (vessel->IsTerminal() && concentration > 1e-2 && random->Uniform(0,1) < 1e-2 * concentration) {
        // stop elongating once it reached tumour
        if (concentration > 0.8) {
          vessel->RemoveBiologyModule(this);
        }
        // get left daughter
        auto* branch_l = bdm_static_cast<Vessel*>(vessel->Bifurcate()[0]);
        branch_l->SetCanBranch(false);
      }
    } // end if not part of main vessel

  } // end run

  private:
    bool init_ = false;
    DiffusionGrid* dg_vegf_ = nullptr;
}; // end VascularGrowth_BM

// -----------------------------------------------------------------------------
/// Secrete substance at SimObject position
struct VegfSecretion_BM : public BaseBiologyModule {
  VegfSecretion_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  VegfSecretion_BM(const Event& event, BaseBiologyModule* other,
                      uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new VegfSecretion_BM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new VegfSecretion_BM(*this);
  }

  void Run(SimObject* so) override {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_vegf_ = rm->GetDiffusionGrid(kVEGF);
      init_ = true;
    }

    Cell* cell = dynamic_cast<Cell*>(so);

    Double3 cell_pos = cell->GetPosition();
    bool vessel_at_vicinity = false;
    // check if vessels at vicinity (dist < 4)
    rm->ApplyOnAllElements([&](SimObject* so) {
      if (auto* vessel = dynamic_cast<Vessel*>(so)) {
        Double3 vessel_pos = vessel->GetPosition();
        // note: square distance
        if (pow(cell_pos[0] - vessel_pos[0], 2) +
            pow(cell_pos[1] - vessel_pos[1], 2) +
            pow(cell_pos[2] - vessel_pos[2], 2) < 16 ) {
          vessel_at_vicinity = true;
        }
      }
    });

    // diffuse vegf only if no vessels at vicinity
    if (!vessel_at_vicinity) {
      dg_vegf_->IncreaseConcentrationBy(cell_pos, 1);
    }
  }

  private:
    bool init_ = false;
    DiffusionGrid* dg_vegf_ = nullptr;
}; // end VegfSecretion_BM

// -----------------------------------------------------------------------------
struct TumourGrowth_BM : public BaseBiologyModule {
  TumourGrowth_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  TumourGrowth_BM(const Event& event, BaseBiologyModule* other,
                      uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new TumourGrowth_BM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new TumourGrowth_BM(*this);
  }

  void Run(SimObject* so) override {
    if (Cell* cell = dynamic_cast<Cell*>(so)) {
      if (cell->GetDiameter() <= 12) {
        cell->ChangeVolume(200);
      } else {
        cell->Divide();
      }
    }
  } // end Run
}; // end TumourGrowth_BM


}  // namespace bdm

#endif  // BIOLOGY_MODULES_H_

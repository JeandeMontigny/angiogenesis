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

    if (!init_) {
      dg_vegf_ = rm->GetDiffusionGrid(kVEGF);
      init_ = true;
    }

    auto* vessel = bdm_static_cast<NeuriteElement*>(so);

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

    auto* vessel = bdm_static_cast<NeuriteElement*>(so);

    auto& secretion_position = vessel->GetPosition();
    dg_vegf_->IncreaseConcentrationBy(secretion_position, 1);
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

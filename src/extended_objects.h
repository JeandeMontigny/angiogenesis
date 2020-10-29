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
#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

#include "core/sim_object/sim_object.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

// Define custom neurite Vessel extending NeuriteElement
  class Vessel : public neuroscience::NeuriteElement {
    BDM_SIM_OBJECT_HEADER(Vessel, neuroscience::NeuriteElement, 1);

   public:
    Vessel() : Base() {}
    virtual ~Vessel() {}

    // Default event constructor
    // copy can_branch_ member for all events
    Vessel(const Event& event, SimObject* other, uint64_t new_oid = 0)
      : Base(event, other, new_oid) {
      can_branch_ = static_cast<Vessel*>(other)->can_branch_;
  }

    void SetCanBranch(bool b) {can_branch_ = b; }
    bool CanBranch() { return can_branch_; }

  private:
     bool can_branch_;
  }; // end Vessel definition


}

#endif

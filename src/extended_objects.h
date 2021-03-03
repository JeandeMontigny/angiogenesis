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

#include "core/agent/agent.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

// Define custom neurite Vessel extending NeuriteElement
  class Vessel : public neuroscience::NeuriteElement {
    BDM_AGENT_HEADER(Vessel, neuroscience::NeuriteElement, 1);

   public:
    Vessel() {}
    virtual ~Vessel() {}

    void Initialize(const NewAgentEvent& event) override {
      Base::Initialize(event);
    }

    void SetCanBranch(bool b) {can_branch_ = b; }
    bool CanBranch() { return can_branch_; }

  private:
     bool can_branch_;
  }; // end Vessel definition


}

#endif

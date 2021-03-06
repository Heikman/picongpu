/**
 * Copyright 2013 Heiko Burau, René Widera, Richard Pausch
 *
 * This file is part of PIConGPU. 
 * 
 * PIConGPU is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * 
 * PIConGPU is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License 
 * along with PIConGPU.  
 * If not, see <http://www.gnu.org/licenses/>. 
 */ 
 


#pragma once

#include "particles/memory/boxes/TileDataBox.hpp"
#include "types.h"
#include "particles/memory/frames/NullFrame.hpp"

namespace picongpu
{
    using namespace PMacc;
    

    template<class Base = NullFrame>
    class Momentum_mt1 : public Base
    {
    private:

        typedef float3_X MyFloat;

        PMACC_ALIGN(momentum_mt1[Base::tileSize], MyFloat);

    public:

        template<class OTHERFRAME>
        HDINLINE void copy(lcellId_t myId, OTHERFRAME &other, lcellId_t otherId)
        {
            Base::copy(myId, other, otherId);
            this->getMomentum_mt1()[myId] = other.getMomentum_mt1()[otherId];
        }

        HDINLINE VectorDataBox<MyFloat> getMomentum_mt1()
        {
            return VectorDataBox<MyFloat > (momentum_mt1);
        }
    };


}//namespace picongpu



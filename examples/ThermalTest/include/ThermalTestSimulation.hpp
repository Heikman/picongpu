/**
 * Copyright 2013 Axel Huebl, Heiko Burau, René Widera
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
 


#ifndef TERMALTESTSIMULATION_HPP
#define	TERMALTESTSIMULATION_HPP

#include "simulation_defines.hpp"

#include "simulationControl/MySimulation.hpp"

#include "simulationControl/SimulationHelper.hpp"
#include "simulation_classTypes.hpp"

#include "fields/FieldE.hpp"

#include "dimensions/GridLayout.hpp"
#include "simulation_types.hpp"
#include "eventSystem/EventSystem.hpp"

#include "nvidia/memory/MemoryInfo.hpp"
#include "mappings/kernel/MappingDescription.hpp"
#include "ArgsParser.hpp"

#include <assert.h>

#include "particles/Species.hpp"

#include "plugins/PluginController.hpp"

#include "cuSTL/container/DeviceBuffer.hpp"
#include "cuSTL/container/HostBuffer.hpp"
#include "cuSTL/algorithm/kernel/Foreach.hpp"
#include "cuSTL/algorithm/kernel/Reduce.hpp"
#include "cuSTL/algorithm/mpi/Gather.hpp"
#include "lambda/Expression.hpp"
#include "math/vector/compile-time/Int.hpp"
#include "math/vector/compile-time/Size_t.hpp"
#include "cuSTL/cursor/tools/slice.hpp"
#include "cuSTL/cursor/FunctorCursor.hpp"

#include "cuSTL/container/allocator/DeviceMemEvenPitchAllocator.hpp"
#include "cuSTL/algorithm/host/Foreach.hpp"
#include "math/vector/math_functor/min.hpp"
#include "math/vector/math_functor/max.hpp"
#include "math/vector/math_functor/sqrtf.hpp"
#include "math/vector/math_functor/cosf.hpp"


namespace picongpu
{

using namespace PMacc;

class ThermalTestSimulation : public MySimulation
{
public:

    ThermalTestSimulation()
    : MySimulation()
    {
    }

    uint32_t init()
    {
        MySimulation::init();

        using namespace ::PMacc::math;

        BOOST_AUTO(fieldE_coreBorder,
             this->fieldE->getGridBuffer().getDeviceBuffer().cartBuffer().view(
                   GuardDim().vec(), -GuardDim().vec()));
                   
        this->eField_zt[0] = new container::DeviceBuffer<float, 2 > (Size_t < 2 > (fieldE_coreBorder.size().z(), this->collectTimesteps));
        this->eField_zt[1] = new container::DeviceBuffer<float, 2 >(this->eField_zt[0]->size());

        return 0;
    }

    void moduleRegisterHelp(po::options_description& desc)
    {
        MySimulation::moduleRegisterHelp(desc);
    }

    void moduleLoad()
    {
        MySimulation::moduleLoad();
    }

    virtual ~ThermalTestSimulation()
    {
        __delete(eField_zt[0]);
        __delete(eField_zt[1]);
    }

    void writeOutput()
    {
        using namespace ::PMacc::math;

        PMacc::GridController<SIMDIM>& con = PMacc::GridController<SIMDIM>::getInstance();
        Size_t<SIMDIM> gpuDim = (Size_t<SIMDIM>)con.getGpuNodes();
        zone::SphericZone<SIMDIM> gpuGatheringZone(Size_t<SIMDIM > (1, 1, gpuDim.z()));
        algorithm::mpi::Gather<SIMDIM> gather(gpuGatheringZone);
        if (!gather.participate()) return;

        container::HostBuffer<float, 2 > eField_zt_host(eField_zt[0]->size());
        for (int i = 0; i < 2; i++)
        {
            eField_zt_host = *(eField_zt[i]);
            container::HostBuffer<float, 2 > global_eField_zt
                    (Size_t < 2 > (gpuDim.z() * eField_zt[i]->size().x(), eField_zt[i]->size().y()));
            gather(global_eField_zt, eField_zt_host, 1);
            if (gather.root())
            {
                std::string filename;
                if (i == 0)
                    filename = "eField_zt_trans.dat";
                else
                    filename = "eField_zt_long.dat";
                std::ofstream eField_zt_dat(filename.data());
                eField_zt_dat << global_eField_zt;
                eField_zt_dat.close();
            }
        }

    }

    /**
     * Run one simulation step.
     *
     * @param currentStep iteration number of the current step
     */
    void runOneStep(uint32_t currentStep)
    {
        MySimulation::runOneStep(currentStep);

        if (currentStep > this->collectTimesteps + firstTimestep)
            return;
        if (currentStep < firstTimestep)
            return;

        using namespace math;

        BOOST_AUTO(fieldE_coreBorder,
           this->fieldE->getGridBuffer().getDeviceBuffer().cartBuffer().view(
                GuardDim().vec(), -GuardDim().vec()));

        for (size_t z = 0; z < eField_zt[0]->size().x(); z++)
        {
            zone::SphericZone < 2 > reduceZone(fieldE_coreBorder.size().shrink<2>());
            using namespace lambda;
            for (int i = 0; i < 2; i++)
            {
                algorithm::kernel::Reduce<BlockDim > ()
                        (eField_zt[i]->origin()(z, currentStep - firstTimestep), reduceZone,
                        cursor::make_FunctorCursor(
                                                   cursor::tools::slice(fieldE_coreBorder.origin()(0, 0, z)),
                                                   _1[i == 0 ? 0 : 2]),
                        _1 + _2);
            }
        }
        
        if (currentStep == this->collectTimesteps + firstTimestep)
            writeOutput();
    }

private:
    // number of timesteps which collect the data
    static const uint32_t collectTimesteps = 768;
    // first timestep which collects data
    //   you may like to let the plasma develope/thermalize a little bit
    static const uint32_t firstTimestep = 1024;

    container::DeviceBuffer<float, 2 >* eField_zt[2];

    typedef PMacc::math::CT::Size_t < 16, 16, 1 > BlockDim;
    typedef PMacc::math::CT::Size_t<TILE_WIDTH, TILE_HEIGHT, TILE_DEPTH> GuardDim;
};

} // namespace picongpu

#endif	/* TERMALTESTSIMULATION_HPP */


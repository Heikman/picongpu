/**
 * Copyright 2013 René Widera
 *
 * This file is part of libPMacc. 
 * 
 * libPMacc is free software: you can redistribute it and/or modify 
 * it under the terms of of either the GNU General Public License or 
 * the GNU Lesser General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * libPMacc is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License and the GNU Lesser General Public License 
 * for more details. 
 * 
 * You should have received a copy of the GNU General Public License 
 * and the GNU Lesser General Public License along with libPMacc. 
 * If not, see <http://www.gnu.org/licenses/>. 
 */ 
 
/* 
 * File:   Factory.hpp
 * Author: widera
 *
 * Created on 9. April 2010, 14:26
 */


#include "particles/tasks/ParticleFactory.hpp"
#include "particles/tasks/TaskSendParticlesExchange.hpp"
#include "particles/tasks/TaskReceiveParticlesExchange.hpp"
#include "particles/tasks/TaskParticlesReceive.hpp"
#include "particles/tasks/TaskParticlesSend.hpp"

#include "eventSystem/tasks/Factory.hpp"
#include "eventSystem/tasks/ITask.hpp"

namespace PMacc
{

    template<class ParBase>
    inline EventTask ParticleFactory::createTaskParticlesReceive(ParBase &parBase,
    ITask *registeringTask)
    {
        TaskParticlesReceive<ParBase>* task = new TaskParticlesReceive<ParBase > (parBase);

        return Factory::getInstance().startTask(*task, registeringTask);
    }

    template<class ParBase>
    inline EventTask ParticleFactory::createTaskReceiveParticlesExchange(ParBase &parBase, uint32_t exchange,
    ITask *registeringTask)
    {
        TaskReceiveParticlesExchange<ParBase>* task = new TaskReceiveParticlesExchange<ParBase > (parBase, exchange);

        return Factory::getInstance().startTask(*task, registeringTask);
    }

    template<class ParBase>
    inline EventTask ParticleFactory::createTaskParticlesSend(ParBase &parBase,
    ITask *registeringTask)
    {
        TaskParticlesSend<ParBase>* task = new TaskParticlesSend<ParBase > (parBase);

        return Factory::getInstance().startTask(*task, registeringTask);
    }

    template<class ParBase>
    inline EventTask ParticleFactory::createTaskSendParticlesExchange(ParBase &parBase, uint32_t exchange,
    ITask *registeringTask)
    {
        TaskSendParticlesExchange<ParBase>* task = new TaskSendParticlesExchange<ParBase > (parBase, exchange);

        return Factory::getInstance().startTask(*task, registeringTask);
    }



} //namespace PMacc





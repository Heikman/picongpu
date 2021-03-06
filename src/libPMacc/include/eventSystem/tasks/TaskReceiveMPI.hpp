/**
 * Copyright 2013 Felix Schmitt, René Widera, Wolfgang Hoenig
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
 * File:   TaskReceiveMPI.hpp
 * Author: wolfgang
 *
 * Created on 13. August 2010, 20:33
 */

#ifndef TASKRECEIVEMPI_HPP
#define	TASKRECEIVEMPI_HPP

#include <mpi.h>

#include "memory/buffers/Exchange.hpp"

#include "eventSystem/tasks/MPITask.hpp"

#include "communication/manager_common.h"
#include "communication/ICommunicator.hpp"


namespace PMacc
{

template <class TYPE, unsigned DIM>
class TaskReceiveMPI : public MPITask
{
public:

    TaskReceiveMPI(Exchange<TYPE, DIM> *exchange) :
    MPITask(),
    exchange(exchange)
    {

    }

    virtual void init()
    {
        __startAtomicTransaction();
        this->request = EnvironmentController::getInstance()
                .getCommunicator().startReceive(
                                                exchange->getExchangeType(),
                                                (char*) exchange->getHostBuffer().getBasePointer(),
                                                exchange->getHostBuffer().getDataSpace().getElementCount() * sizeof (TYPE),
                                                exchange->getCommunicationTag());
        __endTransaction();
    }

    bool executeIntern()
    {
        if (this->isFinished())
            return true;

        if (this->request == NULL)
            throw std::runtime_error("request was NULL (call executeIntern after freed");
        
        int flag=0;
        MPI_CHECK(MPI_Test(this->request, &flag, &(this->status)));

        if (flag) //finished
        {
            delete this->request;
            this->request = NULL;
            setFinished();
            return true;
        }
        return false;
    }

    virtual ~TaskReceiveMPI()
    {
        //\\todo: this make problems because we send bytes and not combined types
        int recv_data_count;
        MPI_CHECK(MPI_Get_count(&(this->status), MPI_CHAR, &recv_data_count));
        

        IEventData *edata = new EventDataReceive(NULL, recv_data_count);

        notify(this->myId, RECVFINISHED, edata); /*add notify her*/
        delete edata;

    }

    void event(id_t eventId, EventType type, IEventData* data)
    {


    }

    std::string toString()
    {
        return "TaskReceiveMPI";
    }

private:
    Exchange<TYPE, DIM> *exchange;
    MPI_Request *request;
    MPI_Status status;
};

} //namespace PMacc

#endif	/* TASKRECEIVEMPI_HPP */


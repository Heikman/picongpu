/**
 * Copyright 2013 Heiko Burau, René Widera
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
 * File:   AllReduce.hpp
 * Author: widera
 *
 * Created on 27. November 2012, 14:58
 */


#pragma once

#include "types.h"
#include <mpi.h>
#include "communication/manager_common.h"

namespace PMacc
{
namespace mpi
{

namespace reduceMethods
{

struct AllReduce
{

    HINLINE bool hasResult(int mpiRank) const
    {
        return mpiRank != -1;
    }

    template<class Functor, typename Type >
    HINLINE void operator()(Functor func, Type* dest, Type* src, const size_t count, MPI_Datatype type, MPI_Op op, MPI_Comm comm) const
    {
        MPI_CHECK(MPI_Allreduce((void*) src,
                                (void*) dest,
                                count,
                                type,
                                op, comm));
    }
};

} /*namespace reduceMethods*/

} /*namepsace mpi*/

} /*namespace PMacc*/



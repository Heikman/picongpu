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
 
#ifndef CURSOR_CTBUFFERNAVIGATOR_HPP
#define CURSOR_CTBUFFERNAVIGATOR_HPP

#include <boost/type_traits/remove_pointer.hpp>
#include "math/vector/Int.hpp"

namespace PMacc
{
namespace cursor
{
namespace CT
{
    
template<typename Pitch, int dim = Pitch::dim + 1>
struct BufferNavigator;
    
template<typename Pitch>
struct BufferNavigator<Pitch, 1>
{
    static const int dim = 1;

    template<typename Data>
    HDINLINE
    Data operator()(const Data& data, const math::Int<dim>& jump) const
    {
        char* result = (char*)data;
        result += jump.x() * sizeof(typename boost::remove_pointer<Data>::type);
        return (Data)result;
    }
};

template<typename Pitch>
struct BufferNavigator<Pitch, 2>
{
    static const int dim = 2;

    template<typename Data>
    HDINLINE
    Data operator()(const Data& data, const math::Int<dim>& jump) const
    {
        char* result = (char*)data;
        result += jump.x() * sizeof(typename boost::remove_pointer<Data>::type)
               + jump.y() * Pitch::x::value;
        return (Data)result;
    }
};

template<typename Pitch>
struct BufferNavigator<Pitch, 3>
{
    static const int dim = 3;

    template<typename Data>
    HDINLINE
    Data operator()(const Data& data, const math::Int<dim>& jump) const
    {
        char* result = (char*)data;
        result += jump.x() * sizeof(typename boost::remove_pointer<Data>::type)
               + jump.y() * Pitch::x::value
               + jump.z() * Pitch::y::value;
        return (Data)result;
    }
};
    
} // CT
} // cursor
} // PMacc

#endif // CURSOR_CTBUFFERNAVIGATOR_HPP

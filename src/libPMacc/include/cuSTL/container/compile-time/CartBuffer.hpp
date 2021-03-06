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
 
#ifndef CONTAINER_CT_CTCARTBUFFER_HPP
#define CONTAINER_CT_CTCARTBUFFER_HPP

#include "types.h"
#include "math/vector/compile-time/UInt.hpp"
#include "cuSTL/cursor/compile-time/BufferCursor.hpp"
#include "cuSTL/cursor/navigator/CartNavigator.hpp"
#include "cuSTL/cursor/accessor/PointerAccessor.hpp"
#include <math/vector/compile-time/Int.hpp>
#include <cuSTL/zone/compile-time/SphericZone.hpp>

namespace PMacc
{
namespace container
{
namespace CT
{
    
template<typename Type, typename _Size, typename Allocator, typename Copier, typename Assigner>
class CartBuffer
{
public:
    typedef Type type;
    typedef _Size Size;
    typedef typename Allocator::Pitch Pitch;
    typedef cursor::CT::BufferCursor<Type, Pitch> Cursor;
    static const int dim = Size::dim;
    typedef zone::CT::SphericZone<_Size, typename math::CT::make_Int<dim, 0>::type> Zone;
private:
    Type* dataPointer;
    //HDINLINE void init();
public:
    DINLINE CartBuffer();
    DINLINE CartBuffer(const CT::CartBuffer<Type, Size, Allocator, Copier, Assigner>& other);
    
    DINLINE CT::CartBuffer<Type, Size, Allocator, Copier, Assigner>& 
    operator=(const CT::CartBuffer<Type, Size, Allocator, Copier, Assigner>& rhs);
    
    DINLINE void assign(const Type& value);
    DINLINE Type* getDataPointer() const {return dataPointer;}
    
    DINLINE cursor::CT::BufferCursor<Type, Pitch> origin() const;
    /*
    HDINLINE Cursor<PointerAccessor<Type>, CartNavigator<dim>, char*>
    originCustomAxes(const math::UInt<dim>& axes) const;
    */
    DINLINE math::Size_t<dim> size() const {return math::Size_t<dim>(Size());}
};
 
} // CT
} // container
} // PMacc

#include "CartBuffer.tpp"

#endif // CONTAINER_CT_CTCARTBUFFER_HPP

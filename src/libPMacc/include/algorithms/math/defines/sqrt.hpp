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
 * File:   math.hpp
 * Author: widera
 *
 * Created on 30. Januar 2013, 09:55
 */

#pragma once

namespace PMacc
{
namespace algorithms
{
namespace math
{

namespace detail
{

template<typename Type>
struct Sqrt;

template<typename Type>
struct RSqrt;


} //namespace detail

template<typename T1>
HDINLINE static typename detail::Sqrt< T1 >::result sqrt(const T1& value)
{
    return detail::Sqrt< T1 > ()(value);
}

template<typename T1>
HDINLINE static typename detail::RSqrt< T1 >::result rsqrt(const T1& value)
{
    return detail::RSqrt< T1 > ()(value);
}

} //namespace math
} //namespace algorithms
}//namespace PMacc


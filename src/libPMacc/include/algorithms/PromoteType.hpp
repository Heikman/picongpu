/**
 * Copyright 2013 Axel Huebl, René Widera
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
 * File:   PromoteType.hpp
 * Author: huebl
 *
 * Created on Juli 12th, 2013 - 12:00pm
 */

#pragma once

namespace PMacc
{
namespace algorithms
{
namespace promoteType
{

    // general: use first type
    template<class T1, class T2>
    struct promoteType {
        typedef T1 ValueType;
    };
    
    // special: promote float to double
    template< >
    struct promoteType<float, double> {
        typedef double ValueType;
    };
    

} //namespace promoteType
} //namespace algorithms
} //namespace PMacc
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

#pragma once

#include "types.h"

#include <boost/mpl/pair.hpp>

namespace PMacc
{
namespace bmpl = boost::mpl;

template<typename T_Key>
struct KeyToAliasPair
{
    typedef
    bmpl::pair< T_Key,
            T_Key >
            type;
};

template<template<typename,typename> class T_Alias,typename T_Key>
struct KeyToAliasPair< T_Alias<T_Key,pmacc_isAlias> >
{
    typedef
    bmpl::pair< T_Alias<pmacc_void,pmacc_isAlias> ,
            T_Alias<T_Key,pmacc_isAlias> >
            type;
};


}//namespace PMacc

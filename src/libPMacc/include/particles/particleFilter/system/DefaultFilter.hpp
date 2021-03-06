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
 * File:   DefaultFilter.hpp
 * Author: widera
 *
 * Created on 18. Oktober 2011, 15:10
 */

#ifndef DEFAULTFILTER_HPP
#define	DEFAULTFILTER_HPP

#include "types.h"
#include "particles/frame_types.hpp"
#include "particles/memory/frames/NullFrame.hpp"

namespace PMacc
{
    
    
template<class Base = NullFrame>
class DefaultFilter : public Base
{
    private:
        bool filterActive;
    public:
        
        HDINLINE DefaultFilter() : filterActive(false)
        {}
        
        HDINLINE virtual ~DefaultFilter()
        {}
        
        template<class FRAME>
        HDINLINE bool operator()(FRAME & frame,lcellId_t id)
        {
            return (!filterActive)||Base::operator() (frame,id);
        }

        /*disable or enable filter
         * true == enable filters
         * false == disable filters
         */
        HDINLINE void setStatus(bool active)
        {
            filterActive=active;
        }
        
        HDINLINE bool getStatus()
        {
            return filterActive;
        }
};

template<>
class DefaultFilter<NullFrame> 
{
    private:
        bool alwaysTrue;
    public:
        
        HDINLINE DefaultFilter() : alwaysTrue(true)
        {}
        
        HDINLINE virtual ~DefaultFilter()
        {}
        
        template<class FRAME>
        HDINLINE bool operator()(FRAME & frame,lcellId_t id)
        {
            return alwaysTrue;
        }

        HDINLINE void setDefault(bool value)
        {
            alwaysTrue=value;
        }
        
        HDINLINE bool getDefault()
        {
            return alwaysTrue;
        }
};


} //namepsace Frame

#endif	/* DEFAULTFILTER_HPP */


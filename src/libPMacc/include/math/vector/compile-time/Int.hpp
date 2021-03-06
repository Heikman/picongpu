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
 
#ifndef STLPICCTINT_HPP
#define STLPICCTINT_HPP

#include <stdint.h>
#include "Vector.hpp"
#include <boost/mpl/integral_c.hpp>

namespace PMacc
{
namespace math
{
namespace CT
{
    
template<int x = 123456789, int y = 123456789, int z = 123456789,
         int dummy = 123456789>
struct Int;

template<>
struct Int<> : public CT::Vector<>
{};

template<int x>
struct Int<x> : public CT::Vector<mpl::integral_c<int, x> >
{
    typedef CT::Vector<mpl::integral_c<int, x> > vector_type;
};

template<int x, int y>
struct Int<x, y> : public CT::Vector<mpl::integral_c<int, x>,
                                     mpl::integral_c<int, y> >
{
    typedef CT::Vector<mpl::integral_c<int, x>,
                       mpl::integral_c<int, y> > vector_type;
};

template<int x, int y, int z>
struct Int<x, y, z> : public CT::Vector<mpl::integral_c<int, x>,
                                                   mpl::integral_c<int, y>,
                                                   mpl::integral_c<int, z> >
{
    typedef CT::Vector<mpl::integral_c<int, x>,
                       mpl::integral_c<int, y>,
                       mpl::integral_c<int, z> > vector_type;
};

template<int dim, int val>
struct make_Int;

template<int val>
struct make_Int<1, val>
{
    typedef Int<val> type;
};

template<int val>
struct make_Int<2, val>
{
    typedef Int<val, val> type;
};

template<int val>
struct make_Int<3, val>
{
    typedef Int<val, val, val> type;
};

} // CT
} // math
} // PMacc

#endif //STLPICCTINT_HPP

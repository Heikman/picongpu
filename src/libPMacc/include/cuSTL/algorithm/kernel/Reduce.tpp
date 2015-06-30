/**
 * Copyright 2013 Heiko Burau, Rene Widera
 *
 * This file is part of libPMacc.
 *
 * libPMacc is free software: you can redistribute it and/or modify
 * it under the terms of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
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

#include "math/vector/Size_t.hpp"
#include "cuSTL/algorithm/kernel/Foreach.hpp"
#include "lambda/Expression.hpp"
#include "cuSTL/container/compile-time/SharedBuffer.hpp"
#include "math/Vector.hpp"
#include "cuSTL/cursor/NestedCursor.hpp"
#include "cuSTL/zone/SphericZone.hpp"
#include <boost/type_traits/remove_reference.hpp>
#include <cuSTL/cursor/accessor/CursorAccessor.hpp>
#include <nvidia/reduce/Reduce.hpp>

namespace PMacc
{
namespace algorithm
{
namespace kernel
{
    
namespace detail
{
    
/**
 * @class MappedNavigator
 * @author Heiko Burau 
 * @date 30/06/15
 * @file Reduce.tpp
 * @brief Wraps a ndim-cursor into a 1D-cursor
 */
template<int T_dim>
class MappedNavigator
{
public:
    static const int dim = T_dim;
private:
    math::Size_t<dim> shape;
    int pos;
public:
    HDINLINE
    MappedNavigator(math::Size_t<dim> shape) 
     : shape(shape), pos(0) {}
     
    
    HDINLINE
    math::Int<dim> toNdim(int idx) const
    {
        math::Int<dim> result;
        int volume = 1;
        for(int i = 0; i < dim; i++)
        {
            result[i] = (idx / volume) % this->shape[i];
            volume *= this->shape[i];
        }
        return result;
    }

    template<typename Cursor>
    HDINLINE
    Cursor operator()(const Cursor& cursor, math::Int<1> jump)
    {
        math::Int<dim> ndstart = toNdim(this->pos);
        this->pos += jump.x();
        math::Int<dim> ndend = toNdim(this->pos);
        
        math::Int<dim> ndjump = ndend - ndstart;
        
        return cursor(ndjump);
    }

};
    
} // detail

template<typename SrcCursor, typename Zone, typename NVidiaFunctor>
typename SrcCursor::ValueType Reduce::operator()(const SrcCursor& srcCursor, const Zone& p_zone, const NVidiaFunctor& functor)
{
    SrcCursor srcCursor_shifted = srcCursor(p_zone.offset);
    
    detail::MappedNavigator<Zone::dim> myNavi(p_zone.size);
    
    BOOST_AUTO(_srcCursor, cursor::make_Cursor(cursor::CursorAccessor<SrcCursor>(),
                                               myNavi,
                                               srcCursor_shifted));
    
    PMacc::nvidia::reduce::Reduce reduce(1024);
    return reduce(functor, _srcCursor, p_zone.size.productOfComponents());
}

} // kernel
} // algorithm
} // PMacc

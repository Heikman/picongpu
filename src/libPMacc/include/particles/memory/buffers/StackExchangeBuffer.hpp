/**
 * Copyright 2013 Felix Schmitt, René Widera
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
 * File:   StackExchangeBuffer.hpp
 * Author: widera
 *
 * Created on 26. November 2010, 14:26
 */

#ifndef STACKEXCHANGEBUFFER_HPP
#define	STACKEXCHANGEBUFFER_HPP

#include "particles/memory/boxes/ExchangePopDataBox.hpp"
#include "particles/memory/boxes/ExchangePushDataBox.hpp"
#include "memory/buffers/Exchange.hpp"
#include <cassert>

namespace PMacc
{

    

    /**
     * Can be used for creating several DataBox types from an Exchange.
     *
     * @tparam FRAME frame datatype
     */
    template<class FRAME, class FRAMEINDEX, unsigned DIM>
    class StackExchangeBuffer
    {
    public:

        /**
         * Create a stack from any ExchangeBuffer<FRAME,DIM>.
         *
         * If the stack's internal GridBuffer has no sizeOnDevice, no device querys are allowed.
         *
         * @param stack Exchange
         */
        StackExchangeBuffer(Exchange<FRAME, DIM1> &stack, Exchange<FRAMEINDEX, DIM1> &stackIndexer) :
        stack(stack), stackIndexer(stackIndexer)
        {

        }

        /**
         * Returns a PushDataBox for the internal HostBuffer.
         *
         * @return PushDataBox for host buffer
         */
        ExchangePushDataBox<vint_t, FRAME, DIM> getHostExchangePushDataBox()
        {
            return ExchangePushDataBox<vint_t, FRAME, DIM > (
                                                             stack.getHostBuffer().getBasePointer(),
                                                             stack.getHostBuffer().getCurrentSizePointer(),
                                                             stack.getHostBuffer().getDataSpace().getElementCount(),
                                                             PushDataBox<vint_t, FRAMEINDEX > (
                                                                                               stackIndexer.getHostBuffer().getBasePointer(),
                                                                                               stackIndexer.getHostBuffer().getCurrentSizePointer()));
        }

        /**
         * Returns a PopDataBox for the internal HostBuffer.
         *
         * @return PopDataBox for host buffer
         */
        ExchangePopDataBox<vint_t, FRAME, DIM> getHostExchangePopDataBox()
        {
            return ExchangePopDataBox<vint_t, FRAME, DIM > (
                                                            stack.getHostBuffer().getBasePointer(),
                                                            PopDataBox<vint_t, FRAMEINDEX > (
                                                                                             stackIndexer.getHostBuffer().getBasePointer(),
                                                                                             (vint_t*) stackIndexer.getHostBuffer().getCurrentSizePointer(),
                                                                                             (vint_t) stackIndexer.getHostBuffer().getCurrentSize()));
        }

        /**
         * Returns a PushDataBox for the internal DeviceBuffer.
         *
         * @return PushDataBox for device buffer
         */
        ExchangePushDataBox<vint_t, FRAME, DIM> getDeviceExchangePushDataBox()
        {
            assert(stack.getDeviceBuffer().hasCurrentSizeOnDevice() == true);
            assert(stackIndexer.getDeviceBuffer().hasCurrentSizeOnDevice() == true);
            return ExchangePushDataBox<vint_t, FRAME, DIM > (
                                                             stack.getDeviceBuffer().getBasePointer(),
                                                             (vint_t*) stack.getDeviceBuffer().getCurrentSizeOnDevicePointer(),
                                                             stack.getDeviceBuffer().getDataSpace().getElementCount(),
                                                             PushDataBox<vint_t, FRAMEINDEX > (
                                                                                               stackIndexer.getDeviceBuffer().getBasePointer(),
                                                                                               (vint_t*) stackIndexer.getDeviceBuffer().getCurrentSizeOnDevicePointer()));
        }

        /**
         * Returns a PopDataBox for the internal DeviceBuffer.
         *
         * @return PopDataBox for device buffer
         */
        ExchangePopDataBox<vint_t, FRAME, DIM> getDeviceExchangePopDataBox()
        {
            return ExchangePopDataBox<vint_t, FRAME, DIM > (
                                                            stack.getDeviceBuffer().getBasePointer(),
                                                            PopDataBox<vint_t, FRAMEINDEX > (
                                                                                             stackIndexer.getDeviceBuffer().getBasePointer(),
                                                                                             (vint_t*) stackIndexer.getDeviceBuffer().getCurrentSizeOnDevicePointer(),
                                                                                             (vint_t) stackIndexer.getDeviceBuffer().getCurrentSize()));
        }

        void setCurrentSize(const size_t size)
        {
            // do host and device setCurrentSize parallel
            EventTask split = __getTransactionEvent();
            __startTransaction(split);
            stackIndexer.getHostBuffer().setCurrentSize(size);
            stack.getHostBuffer().setCurrentSize(size);
            EventTask e1 = __endTransaction();

            __startTransaction(split);
            stackIndexer.getDeviceBuffer().setCurrentSize(size);
            EventTask e2 = __endTransaction();
            __startTransaction(split);
            stack.getDeviceBuffer().setCurrentSize(size);
            EventTask e3 = __endTransaction();

            __setTransactionEvent(e1 + e2 + e3);
        }

        size_t getHostCurrentSize()
        {
            return stackIndexer.getHostBuffer().getCurrentSize();
        }

        size_t getDeviceCurrentSize()
        {
            return stackIndexer.getDeviceBuffer().getCurrentSize();
        }

        size_t getDeviceParticlesCurrentSize()
        {
            return stack.getDeviceBuffer().getCurrentSize();
        }

        size_t getHostParticlesCurrentSize()
        {
            return stack.getHostBuffer().getCurrentSize();
        }

        size_t getMaxParticlesCount()
        {
            return stack.getHostBuffer().getDataSpace().getElementCount();
        }

    private:

        Exchange<FRAME, DIM1> &getExchangeBuffer()
        {
            return stack;
        }

        Exchange<FRAME, DIM1>& stack;
        Exchange<FRAMEINDEX, DIM1>& stackIndexer;
    };
}

#endif	/* STACKEXCHANGEBUFFER_HPP */


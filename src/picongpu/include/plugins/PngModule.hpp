/**
 * Copyright 2013 Axel Huebl, René Widera
 *
 * This file is part of PIConGPU. 
 * 
 * PIConGPU is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * 
 * PIConGPU is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License 
 * along with PIConGPU.  
 * If not, see <http://www.gnu.org/licenses/>. 
 */ 
 


#ifndef PNGMODULE_HPP
#define	PNGMODULE_HPP

#include "types.h"
#include "simulation_defines.hpp"
#include "simulation_types.hpp"
#include "dimensions/DataSpace.hpp"

#include "simulation_classTypes.hpp"
#include "plugins/IPluginModule.hpp"
#include <vector>
#include <list>


#include <cassert>

#include <stdexcept>


namespace picongpu
{
    using namespace PMacc;


    namespace po = boost::program_options;

    template<class VisClass>
    class PngModule : public IPluginModule
    {
    public:

        typedef VisClass VisType;
        typedef std::list<VisType*> VisPointerList;

        PngModule(std::string name, std::string prefix) :
        analyzerName(name),
        analyzerPrefix(prefix),
        cellDescription(NULL)
        {
            ModuleConnector::getInstance().registerModule(this);
        }

        virtual ~PngModule()
        {

        }

        std::string moduleGetName() const
        {
            return analyzerName;
        }

        void moduleRegisterHelp(po::options_description& desc)
        {
            desc.add_options()
                    ((analyzerPrefix + ".period").c_str(), po::value<std::vector<uint32_t> > (&notifyFrequencys)->multitoken(), "enable data output [for each n-th step]")
                    ((analyzerPrefix + ".axis").c_str(), po::value<std::vector<std::string > > (&axis)->multitoken(), "axis which are shown [valid values x,y,z] example: yz")
                    ((analyzerPrefix + ".slicePoint").c_str(), po::value<std::vector<float> > (&slicePoints)->multitoken(), "value range: 0 <= x <= 1 , point of the slice")
                    ((analyzerPrefix + ".folder").c_str(), po::value<std::vector<std::string> > (&folders)->multitoken(), "folder for output files");
        }

        void setMappingDescription(MappingDesc *cellDescription)
        {
            this->cellDescription = cellDescription;
        }


    private:

        void moduleLoad()
        {

            if (0 != notifyFrequencys.size())
            {
                if (0 != slicePoints.size() &&
                    0 != axis.size())
                {
                    for (int i = 0; i < (int) slicePoints.size(); ++i) /*!\todo: use vactor with max elements*/
                    {
                        uint32_t frequ = getValue(notifyFrequencys, i);
                        if (frequ != 0)
                        {

                            if (getValue(axis, i).length() == 2u)
                            {
                                std::stringstream o_slicePoint;
                                o_slicePoint << getValue(slicePoints, i);
                                /*add default value for folder*/
                                if (folders.empty())
                                {
                                    folders.push_back(std::string("."));
                                }
                                std::string filename(analyzerName + "_" + getValue(axis, i) + "_" + o_slicePoint.str());
                                typename VisType::CreatorType pngCreator(filename, getValue(folders, i));
                                DataSpace<DIM2 > transpose(
                                                           charToAxisNumber(getValue(axis, i)[0]),
                                                           charToAxisNumber(getValue(axis, i)[1])
                                                           );
                                VisType* tmp = new VisType(analyzerName, pngCreator, frequ, transpose, getValue(slicePoints, i));
                                visIO.push_back(tmp);
                                tmp->setMappingDescription(cellDescription);
                                tmp->init();
                            }
                            else
                                throw std::runtime_error((std::string("[Png Module] wrong charecter count in axis: ") + getValue(axis, i)).c_str());
                        }
                    }
                }
                else
                {
                    throw std::runtime_error("[Png Module] One parameter is missing");
                }
            }
        }

        void moduleUnload()
        {
            for (typename VisPointerList::iterator iter = visIO.begin();
                 iter != visIO.end();
                 ++iter)
            {
                __delete(*iter);
            }
            visIO.clear();
        }

        /*! Get value of the postition in a vector
         * @return value at id postition, if id >= size of vector last valid value is given back
         */
        template<class Vec>
        typename Vec::value_type getValue(Vec vec, size_t id)
        {
            if (vec.size() == 0)
                throw std::runtime_error("[Livew View] getValue is used with a parameter set with no parameters (count is 0)");
            if (id >= vec.size())
            {
                return vec[vec.size() - 1];
            }
            return vec[id];
        }

        int charToAxisNumber(char c)
        {
            if (c == 'x')
                return 0;
            if (c == 'y')
                return 1;
            return 2;
        }


        std::string analyzerName;
        std::string analyzerPrefix;

        std::vector<uint32_t> notifyFrequencys;
        std::vector<float> slicePoints;
        std::vector<std::string> folders;
        std::vector<std::string> axis;
        VisPointerList visIO;

        MappingDesc* cellDescription;

    };

}//namespace

#endif	/* PNGMODULE_HPP */

